"""Tests for the dataset loaders, statistics, and feature encoding.

The folding code is checked against brute force in ``test_mccaskill.py``. The
risk in this half of the project is different: a silent bug in ranking, fold
assignment, or feature indexing would not crash anything, it would just quietly
produce a wrong correlation. These tests pin down the pieces where that could
happen, especially the ones that protect against over-optimistic results
(grouped folds, honest out-of-fold prediction).
"""

from __future__ import annotations

import random
import tempfile
import unittest
from pathlib import Path

from datasets import Guide, assign_percentiles, load_crisprscan, percentile_ranks
from model_features import FeatureSpec, build_design_rows
from stats import (
    SparseDesign,
    mean_within_group_spearman,
    within_group_spearman_ci,
    assign_group_folds,
    benjamini_hochberg,
    cholesky_solve,
    cluster_bootstrap_ci,
    fit_full,
    nested_cv_predictions,
    paired_spearman_delta_ci,
    pearson,
    rank_average,
    spearman,
    spearman_permutation_p,
)


class TestRanking(unittest.TestCase):
    def test_percentile_ranks_span_zero_to_one(self):
        self.assertEqual(percentile_ranks([10.0, 20.0, 30.0]), [0.0, 0.5, 1.0])

    def test_percentile_ranks_average_ties(self):
        self.assertEqual(percentile_ranks([5.0, 5.0, 9.0]), [0.25, 0.25, 1.0])

    def test_percentile_ranks_edge_cases(self):
        self.assertEqual(percentile_ranks([]), [])
        self.assertEqual(percentile_ranks([7.0]), [0.5])

    def test_percentiles_are_computed_within_each_screen(self):
        """Two screens on different scales must not be ranked against each other."""
        guides = [
            Guide("a", "g1", "A" * 20, 1.0, "screenA", "d"),
            Guide("b", "g1", "C" * 20, 2.0, "screenA", "d"),
            Guide("c", "g2", "G" * 20, 500.0, "screenB", "d"),
            Guide("d", "g2", "T" * 20, 900.0, "screenB", "d"),
        ]
        ranked = assign_percentiles(guides)
        self.assertEqual([g.percentile for g in ranked], [0.0, 1.0, 0.0, 1.0])

    def test_rank_average_handles_ties(self):
        self.assertEqual(rank_average([3.0, 1.0, 1.0, 7.0]), [3.0, 1.5, 1.5, 4.0])


class TestCorrelation(unittest.TestCase):
    def test_perfect_monotone_relationships(self):
        self.assertAlmostEqual(spearman([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
        self.assertAlmostEqual(spearman([1, 2, 3, 4], [40, 30, 20, 10]), -1.0)

    def test_spearman_is_invariant_to_monotone_transforms(self):
        x = [0.1, 0.9, 0.4, 0.7, 0.2]
        y = [3.0, 1.0, 2.0, 8.0, 5.0]
        transformed = [value ** 3 for value in y]
        self.assertAlmostEqual(spearman(x, y), spearman(x, transformed))

    def test_constant_input_gives_zero(self):
        self.assertEqual(spearman([1, 1, 1, 1], [1, 2, 3, 4]), 0.0)
        self.assertEqual(pearson([2.0, 2.0], [1.0, 5.0]), 0.0)

    def test_length_mismatch_is_rejected(self):
        with self.assertRaises(ValueError):
            pearson([1.0, 2.0], [1.0])

    def test_benjamini_hochberg_is_monotone_and_bounded(self):
        raw = [0.001, 0.008, 0.02, 0.04, 0.3, 0.9]
        adjusted = benjamini_hochberg(raw)
        self.assertEqual(len(adjusted), len(raw))
        for value in adjusted:
            self.assertGreaterEqual(value, 0.0)
            self.assertLessEqual(value, 1.0)
        for earlier, later in zip(adjusted, adjusted[1:]):
            self.assertLessEqual(earlier, later + 1e-12)
        for raw_value, adjusted_value in zip(raw, adjusted):
            self.assertGreaterEqual(adjusted_value, raw_value - 1e-12)

    def test_permutation_p_is_small_for_a_real_association(self):
        rng = random.Random(3)
        x = [rng.random() for _ in range(120)]
        y = [value * 2 + rng.gauss(0, 0.05) for value in x]
        self.assertLess(spearman_permutation_p(x, y, n_perm=500, seed=1), 0.01)

    def test_permutation_p_is_calibrated_under_the_null(self):
        """Checking one random pair would be flaky: a single draw of pure noise
        has a genuine 5% chance of looking significant. The meaningful property
        is calibration across many draws - null p-values should be roughly
        uniform, so only about 5% land below 0.05."""
        pvalues = []
        for trial in range(40):
            rng = random.Random(100 + trial)
            x = [rng.random() for _ in range(120)]
            y = [rng.random() for _ in range(120)]
            pvalues.append(spearman_permutation_p(x, y, n_perm=300, seed=trial))
        significant = sum(1 for p in pvalues if p < 0.05) / len(pvalues)
        self.assertLess(significant, 0.20)
        pvalues.sort()
        median = pvalues[len(pvalues) // 2]
        self.assertGreater(median, 0.2)
        self.assertLess(median, 0.8)

    def test_cluster_bootstrap_brackets_the_estimate(self):
        rng = random.Random(5)
        n = 300
        x = [rng.random() for _ in range(n)]
        y = [value + rng.gauss(0, 0.4) for value in x]
        groups = [f"g{i % 12}" for i in range(n)]
        observed = spearman(x, y)
        low, high = cluster_bootstrap_ci(
            lambda idx: spearman([x[i] for i in idx], [y[i] for i in idx]),
            groups, n_boot=300, seed=7,
        )
        self.assertLessEqual(low, observed)
        self.assertGreaterEqual(high, observed)


class TestLinearAlgebra(unittest.TestCase):
    def test_cholesky_solves_a_known_system(self):
        matrix = [[4.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]]
        expected = [1.0, -2.0, 3.0]
        vector = [
            sum(matrix[i][j] * expected[j] for j in range(3)) for i in range(3)
        ]
        solution = cholesky_solve(matrix, vector)
        for got, want in zip(solution, expected):
            self.assertAlmostEqual(got, want, places=10)

    def test_cholesky_matches_a_random_spd_system(self):
        rng = random.Random(11)
        size = 25
        factor = [[rng.gauss(0, 1) for _ in range(size)] for _ in range(size)]
        matrix = [
            [sum(factor[k][i] * factor[k][j] for k in range(size)) for j in range(size)]
            for i in range(size)
        ]
        for i in range(size):
            matrix[i][i] += size
        expected = [rng.gauss(0, 1) for _ in range(size)]
        vector = [
            sum(matrix[i][j] * expected[j] for j in range(size)) for i in range(size)
        ]
        solution = cholesky_solve(matrix, vector)
        for got, want in zip(solution, expected):
            self.assertAlmostEqual(got, want, places=7)


class TestFolds(unittest.TestCase):
    def test_groups_are_never_split_across_folds(self):
        groups = [f"gene{i % 9}" for i in range(200)]
        folds = assign_group_folds(groups, 5, seed=2)
        seen: dict[str, int] = {}
        for group, fold in zip(groups, folds):
            if group in seen:
                self.assertEqual(seen[group], fold, f"{group} was split")
            seen[group] = fold
        self.assertEqual(len(set(folds)), 5)

    def test_folds_are_roughly_balanced(self):
        groups = [f"g{i}" for i in range(40) for _ in range(i % 7 + 1)]
        folds = assign_group_folds(groups, 4, seed=1)
        sizes = [folds.count(f) for f in range(4)]
        self.assertLess(max(sizes) - min(sizes), len(groups) * 0.25)

    def test_assignment_is_deterministic(self):
        groups = [f"g{i % 11}" for i in range(90)]
        self.assertEqual(
            assign_group_folds(groups, 3, seed=42), assign_group_folds(groups, 3, seed=42)
        )


class TestRidge(unittest.TestCase):
    def _linear_problem(self, n=400, seed=17):
        """A dataset where y really is a linear function of two features."""
        rng = random.Random(seed)
        rows, targets, groups = [], [], []
        for index in range(n):
            a = rng.random()
            b = rng.random()
            # Column 0 is the intercept.
            rows.append([(0, 1.0), (1, a), (2, b)])
            targets.append(0.3 + 2.0 * a - 1.0 * b + rng.gauss(0, 0.02))
            groups.append(f"g{index % 20}")
        return SparseDesign(rows, 3), targets, groups

    def test_fit_recovers_known_coefficients(self):
        design, targets, _ = self._linear_problem()
        weights = fit_full(design, targets, alpha=1e-6, intercept_index=0)
        self.assertAlmostEqual(weights[0], 0.3, places=1)
        self.assertAlmostEqual(weights[1], 2.0, places=1)
        self.assertAlmostEqual(weights[2], -1.0, places=1)

    def test_ridge_shrinks_coefficients(self):
        design, targets, _ = self._linear_problem()
        weak = fit_full(design, targets, alpha=1e-6, intercept_index=0)
        strong = fit_full(design, targets, alpha=1e5, intercept_index=0)
        self.assertLess(abs(strong[1]), abs(weak[1]))

    def test_intercept_is_not_penalised(self):
        """With a huge penalty the slopes vanish but the intercept still fits."""
        design, targets, _ = self._linear_problem()
        weights = fit_full(design, targets, alpha=1e8, intercept_index=0)
        self.assertLess(abs(weights[1]), 0.01)
        mean = sum(targets) / len(targets)
        self.assertAlmostEqual(weights[0], mean, places=2)

    def test_nested_cv_recovers_a_real_signal(self):
        design, targets, groups = self._linear_problem()
        predictions, alphas = nested_cv_predictions(
            design, targets, groups, [0.01, 1.0, 100.0],
            n_outer=4, n_inner=3, seed=3, intercept_index=0,
        )
        self.assertEqual(len(alphas), 4)
        self.assertGreater(spearman(predictions, targets), 0.9)

    def test_nested_cv_finds_nothing_in_pure_noise(self):
        """The critical guard: out-of-fold predictions must not learn noise.

        If cross-validation leaked, a model given random features would still
        appear to predict the labels. It must not.
        """
        rng = random.Random(23)
        rows, targets, groups = [], [], []
        for index in range(300):
            row = [(0, 1.0)] + [(f, rng.random()) for f in range(1, 25)]
            rows.append(row)
            targets.append(rng.random())
            groups.append(f"g{index % 15}")
        design = SparseDesign(rows, 25)
        predictions, _ = nested_cv_predictions(
            design, targets, groups, [0.1, 1.0, 10.0, 100.0],
            n_outer=4, n_inner=3, seed=5, intercept_index=0,
        )
        self.assertLess(abs(spearman(predictions, targets)), 0.25)

    def test_group_offsets_destroy_the_pooled_metric_but_not_the_within_group_one(self):
        """The bug this metric exists to avoid.

        Each group has its own baseline level that a model holding out whole
        groups cannot know. Here the predictions rank perfectly inside every
        group but carry no information about the group offsets. Pooling the
        predictions therefore looks like nothing; scoring within groups
        correctly reports a perfect ranking.
        """
        rng = random.Random(77)
        predicted, truth, groups = [], [], []
        for group in range(10):
            offset = rng.uniform(0, 10)  # this group's baseline, unknowable to the model
            for _ in range(30):
                signal = rng.random()
                predicted.append(signal)          # no offset: the model cannot see it
                truth.append(offset + signal)     # offset dominates the pooled ranking
                groups.append(f"g{group}")

        pooled = spearman(predicted, truth)
        within = mean_within_group_spearman(predicted, truth, groups)
        self.assertLess(abs(pooled), 0.35)
        self.assertGreater(within, 0.99)

    def test_within_group_spearman_skips_small_groups(self):
        predicted = [1.0, 2.0, 3.0] + [float(i) for i in range(10)]
        truth = [3.0, 2.0, 1.0] + [float(i) for i in range(10)]
        groups = ["tiny"] * 3 + ["big"] * 10
        # The tiny group correlates -1 and the big one +1; only the big one counts.
        self.assertAlmostEqual(
            mean_within_group_spearman(predicted, truth, groups, min_size=8), 1.0
        )

    def test_within_group_spearman_ci_brackets_the_mean(self):
        rng = random.Random(101)
        predicted, truth, groups = [], [], []
        for group in range(14):
            for _ in range(20):
                value = rng.random()
                predicted.append(value)
                truth.append(value + rng.gauss(0, 0.3))
                groups.append(f"g{group}")
        mean, low, high, n_groups = within_group_spearman_ci(
            predicted, truth, groups, n_boot=400, seed=2
        )
        self.assertEqual(n_groups, 14)
        self.assertLessEqual(low, mean)
        self.assertGreaterEqual(high, mean)
        self.assertGreater(mean, 0.5)

    def test_paired_delta_ci_contains_the_observed_difference(self):
        rng = random.Random(31)
        n = 300
        truth = [rng.random() for _ in range(n)]
        baseline = [value + rng.gauss(0, 0.5) for value in truth]
        better = [value + rng.gauss(0, 0.2) for value in truth]
        groups = [f"g{i % 15}" for i in range(n)]
        delta, low, high = paired_spearman_delta_ci(
            baseline, better, truth, groups, n_boot=300, seed=9
        )
        self.assertGreater(delta, 0.0)
        self.assertLessEqual(low, delta)
        self.assertGreaterEqual(high, delta)


class TestFeatureEncoding(unittest.TestCase):
    def _guide(self, spacer: str) -> Guide:
        return Guide("id", "gene", spacer, 0.5, "screen", "dataset")

    def test_position_and_dinucleotide_columns(self):
        spec = FeatureSpec(position=True, dinucleotide=True, gc=True)
        guide = self._guide("ACGT" * 5)
        rows = build_design_rows([guide], [{}], spec)
        entries = dict(rows[0])
        # 1 intercept + 20 positions + 19 dinucleotides + 1 G/C
        self.assertEqual(len(rows[0]), 41)
        self.assertEqual(entries[spec.intercept_index], 1.0)
        # Position 0 is 'A' -> offset + 0 * 4 + 0
        self.assertIn(spec.position_offset, entries)
        # G/C content of ACGTACGT... is exactly 0.5
        self.assertAlmostEqual(entries[spec.gc_offset], 0.5)

    def test_feature_count_matches_the_spec(self):
        spec = FeatureSpec(position=True, dinucleotide=True, gc=True)
        self.assertEqual(spec.n_features, 1 + 20 * 4 + 19 * 16 + 1)

    def test_structure_columns_are_appended_and_scaled(self):
        spec = FeatureSpec(
            position=False, dinucleotide=False, gc=False,
            structure=["seed_ensemble_full", "ensemble_energy_full"],
        )
        cached = {"seed_ensemble_full": "0.400", "ensemble_energy_full": "-20.0"}
        rows = build_design_rows([self._guide("A" * 20)], [cached], spec)
        entries = dict(rows[0])
        self.assertAlmostEqual(entries[spec.structure_offset], 0.4)
        # Energies are divided by 10 to bring them to order 1.
        self.assertAlmostEqual(entries[spec.structure_offset + 1], -2.0)

    def test_different_spacers_give_different_encodings(self):
        spec = FeatureSpec()
        rows = build_design_rows(
            [self._guide("A" * 20), self._guide("C" * 20)], [{}, {}], spec
        )
        self.assertNotEqual(set(rows[0]), set(rows[1]))

    def test_wrong_length_spacer_is_rejected(self):
        with self.assertRaises(ValueError):
            build_design_rows([self._guide("ACGT")], [{}], FeatureSpec())

    def test_misaligned_cache_is_rejected(self):
        with self.assertRaises(ValueError):
            build_design_rows([self._guide("A" * 20)], [], FeatureSpec())


class TestDatasetLoading(unittest.TestCase):
    def test_crisprscan_loads_and_validates(self):
        guides = load_crisprscan()
        self.assertGreater(len(guides), 900)
        self.assertTrue(all(len(g.spacer) == 20 for g in guides))
        self.assertTrue(all(g.screen == "MorenoMateos2015" for g in guides))
        self.assertGreater(len({g.group for g in guides}), 50)

    def test_bad_spacer_length_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.csv"
            path.write_text(
                "guide_id,gene,spacer_dna,pam,mod_freq\nx,g,ACGT,TGG,0.5\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "20 nt spacer"):
                load_crisprscan(path)

    def test_missing_file_is_reported_clearly(self):
        with self.assertRaises(FileNotFoundError):
            load_crisprscan(Path("does/not/exist.csv"))


if __name__ == "__main__":
    unittest.main()
