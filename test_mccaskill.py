"""Correctness tests for the energy model, the MFE folder, and the partition function.

The important tests here are the brute-force ones. For short sequences every
pseudoknot-free structure can be enumerated, its energy evaluated directly by
loop decomposition, and the exact Boltzmann statistics computed by hand. The
dynamic programming must reproduce those numbers to floating-point tolerance.
That is a much stronger claim than "the output looks like a hairpin", and it is
what makes the ensemble accessibility numbers in the study trustworthy.
"""

from __future__ import annotations

import math
import random
import unittest

import energy_model
import mccaskill
from energy_model import (
    RT37,
    BoltzmannTables,
    EnergyModel,
    enumerate_structures,
    pairs_from_dot_bracket,
    structure_energy,
)
from mccaskill import partition_fold, sample_structures, unpaired_from_samples
from zuker import zuker_fold


# Short enough to enumerate exhaustively, varied enough to exercise hairpins,
# bulges, internal loops, wobble pairs, and three- and four-way junctions.
BRUTE_FORCE_SEQUENCES = [
    "GGGAAACCC",
    "GGGGAAAACCCC",
    "GCGCAAAAGCGC",
    "AUGCUAGCUAGC",
    "ACGUACGUACGUAC",
    "GGGAAACCCUUGGG",
    "AAAAAAAAAAAA",
    "GCGCGCAAAGCGCGC",
    "GGAAACGAAACC",
    "GGGAAACGAAACCC",
    "CGGAAACCGAAACCG",
    "GGAAACGAAACGAAACC",
    "GGGCAAAGCCGAAACCC",
    "GUGAAACGAAACAC",
    "GGCAUAUGCCAAAGGC",
]


def brute_force_ensemble(sequence: str, rt: float = RT37):
    """Exact partition function and pair probabilities by enumeration."""
    total = 0.0
    pair_weight: dict[tuple[int, int], float] = {}
    for pairs in enumerate_structures(sequence):
        energy = structure_energy(sequence, pairs)
        weight = 0.0 if energy == math.inf else math.exp(-energy / rt)
        total += weight
        for pair in pairs:
            pair_weight[pair] = pair_weight.get(pair, 0.0) + weight
    return total, {pair: weight / total for pair, weight in pair_weight.items()}


class TestEnergyModel(unittest.TestCase):
    def test_dot_bracket_round_trip(self):
        structure = "((((....))))"
        pairs = pairs_from_dot_bracket(structure)
        self.assertEqual(pairs, [(0, 11), (1, 10), (2, 9), (3, 8)])

    def test_unbalanced_dot_bracket_is_rejected(self):
        with self.assertRaises(ValueError):
            pairs_from_dot_bracket("((..)")
        with self.assertRaises(ValueError):
            pairs_from_dot_bracket("(..))")

    def test_structure_energy_matches_hand_calculation(self):
        # GGGG/CCCC stem, four-base loop: three G-C stacks plus hairpin initiation.
        model = EnergyModel("GGGGAAAACCCC")
        expected = (
            model.stack(0, 11, 1, 10)
            + model.stack(1, 10, 2, 9)
            + model.stack(2, 9, 3, 8)
            + model.hairpin(3, 8)
        )
        self.assertAlmostEqual(structure_energy("GGGGAAAACCCC", "((((....))))"), expected, places=10)

    def test_overlapping_pairs_are_rejected(self):
        with self.assertRaises(ValueError):
            structure_energy("GGGGAAAACCCC", [(0, 11), (0, 10)])

    def test_boltzmann_tables_match_the_energy_model(self):
        """The tabulated weights must equal exp(-E/RT) computed the slow way."""
        for sequence in ["GGGGAAAACCCCAGGCAUAUGCCUU", "GCGCAAAAGCGCUUAGGCAAAGCC"]:
            model = EnergyModel(sequence)
            tables = BoltzmannTables(model, RT37)
            n = model.n
            checked = 0
            for i in range(n):
                for j in range(i + 5, n):
                    if not model.can_pair(i, j):
                        continue
                    for k in range(i + 1, j):
                        for l in range(k + 4, j):
                            if not model.can_pair(k, l):
                                continue
                            if (k - i - 1) + (j - l - 1) > energy_model.MAX_INTERNAL:
                                continue
                            expected = math.exp(-model.interior(i, j, k, l) / RT37)
                            self.assertAlmostEqual(
                                tables.interior(i, j, k, l), expected, places=12,
                                msg=f"{sequence} ({i},{j})<-({k},{l})",
                            )
                            checked += 1
            self.assertGreater(checked, 50)


class TestZuker(unittest.TestCase):
    def test_mfe_equals_brute_force_minimum(self):
        for sequence in BRUTE_FORCE_SEQUENCES:
            with self.subTest(sequence=sequence):
                best = min(structure_energy(sequence, p) for p in enumerate_structures(sequence))
                result = zuker_fold(sequence)
                self.assertAlmostEqual(result.energy, best, places=7)

    def test_traceback_structure_has_the_reported_energy(self):
        """A traceback bug shows up as a structure whose energy is not the MFE."""
        for sequence in BRUTE_FORCE_SEQUENCES:
            with self.subTest(sequence=sequence):
                result = zuker_fold(sequence)
                self.assertAlmostEqual(
                    structure_energy(sequence, result.pairs), result.energy, places=7
                )

    def test_trna_phe_recovers_the_accepted_cloverleaf(self):
        sequence = (
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
        )
        accepted = "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))...."
        result = zuker_fold(sequence)
        self.assertEqual(result.structure, accepted)

    def test_unpaired_flags_are_zero_one(self):
        result = zuker_fold("GGGGAAAACCCC")
        self.assertEqual(result.unpaired_flags(), [0.0] * 4 + [1.0] * 4 + [0.0] * 4)


class TestPartitionFunction(unittest.TestCase):
    def test_partition_function_matches_enumeration(self):
        for sequence in BRUTE_FORCE_SEQUENCES:
            with self.subTest(sequence=sequence):
                expected, _ = brute_force_ensemble(sequence)
                result = partition_fold(sequence)
                self.assertAlmostEqual(
                    result.partition_function / expected, 1.0, places=9
                )

    def test_pair_probabilities_match_enumeration(self):
        for sequence in BRUTE_FORCE_SEQUENCES:
            with self.subTest(sequence=sequence):
                _, expected = brute_force_ensemble(sequence)
                result = partition_fold(sequence)
                n = len(sequence)
                for i in range(n):
                    for j in range(i + 1, n):
                        self.assertAlmostEqual(
                            result.pair_prob[i][j], expected.get((i, j), 0.0), places=9,
                            msg=f"{sequence} pair ({i},{j})",
                        )

    def test_unpaired_probabilities_match_enumeration(self):
        for sequence in BRUTE_FORCE_SEQUENCES:
            with self.subTest(sequence=sequence):
                _, expected = brute_force_ensemble(sequence)
                result = partition_fold(sequence)
                for i in range(len(sequence)):
                    paired = sum(
                        weight for (a, b), weight in expected.items() if a == i or b == i
                    )
                    self.assertAlmostEqual(result.unpaired[i], 1.0 - paired, places=9)

    def test_multiloop_recursions_under_stressed_parameters(self):
        """Shipped parameters make multiloops rare, so force them to matter.

        With the defaults, multiloop structures carry well under 1% of the
        ensemble weight, and ML_UNPAIRED is 0 so every z^k factor equals 1. Both
        modules read these constants at call time, so patching them keeps brute
        force and dynamic programming describing the same model while pushing
        multiloop weight above 50%.
        """
        cases = [
            (-1.0, -0.9, 0.0, RT37),
            (0.5, -0.5, 0.4, RT37),
            (-1.5, -1.2, 0.3, 1.2),
            (1.0, 0.3, -0.2, RT37),
        ]
        saved = (energy_model.ML_CLOSE, energy_model.ML_BRANCH, energy_model.ML_UNPAIRED)
        try:
            for ml_close, ml_branch, ml_unpaired, rt in cases:
                for module in (energy_model, mccaskill):
                    module.ML_CLOSE = ml_close
                    module.ML_BRANCH = ml_branch
                    module.ML_UNPAIRED = ml_unpaired
                for sequence in BRUTE_FORCE_SEQUENCES:
                    with self.subTest(sequence=sequence, ml_close=ml_close, rt=rt):
                        expected_q, expected_p = brute_force_ensemble(sequence, rt=rt)
                        result = partition_fold(sequence, temperature_rt=rt)
                        self.assertAlmostEqual(
                            result.partition_function / expected_q, 1.0, places=9
                        )
                        n = len(sequence)
                        for i in range(n):
                            for j in range(i + 1, n):
                                self.assertAlmostEqual(
                                    result.pair_prob[i][j],
                                    expected_p.get((i, j), 0.0),
                                    places=9,
                                )
        finally:
            for module in (energy_model, mccaskill):
                module.ML_CLOSE, module.ML_BRANCH, module.ML_UNPAIRED = saved

    def test_ensemble_free_energy_is_at_most_the_mfe(self):
        """Summing over structures can only lower the free energy below the best one."""
        for sequence in BRUTE_FORCE_SEQUENCES:
            with self.subTest(sequence=sequence):
                ensemble = partition_fold(sequence).ensemble_free_energy
                self.assertLessEqual(ensemble, zuker_fold(sequence).energy + 1e-9)

    def test_probabilities_are_bounded(self):
        result = partition_fold("GGCAUAUGCCAAAGGCAUAUGCCUUAGCAUAUGCU")
        for value in result.unpaired:
            self.assertGreaterEqual(value, 0.0)
            self.assertLessEqual(value, 1.0)
        n = result.n
        for i in range(n):
            row_total = sum(result.pair_prob[i][j] for j in range(i + 1, n))
            row_total += sum(result.pair_prob[k][i] for k in range(0, i))
            self.assertLessEqual(row_total, 1.0 + 1e-9)

    def test_sampling_converges_to_exact_probabilities(self):
        """Stochastic traceback uses only the inside matrices; the outside pass is
        an independent computation. Agreement between them checks both."""
        sequence = "GGCAUAUGCCAAAGGCAUAUGCCUUAGGCAUAUGCCAAAUU"
        result = partition_fold(sequence)
        samples = sample_structures(result, 6000, rng=random.Random(20260724))
        approximate = unpaired_from_samples(samples, result.n)
        for exact, sampled in zip(result.unpaired, approximate):
            # 6000 draws gives a standard error of at most ~0.0065 per base.
            self.assertLess(abs(exact - sampled), 0.03)

    def test_sampled_structures_are_valid(self):
        result = partition_fold("GGCAUAUGCCAAAGGCAUAUGCCUU")
        model = EnergyModel(result.sequence)
        for pairs in sample_structures(result, 40, rng=random.Random(5)):
            used = set()
            for i, j in pairs:
                self.assertTrue(model.can_pair(i, j))
                self.assertNotIn(i, used)
                self.assertNotIn(j, used)
                used.update((i, j))
            for a, b in pairs:
                for c, d in pairs:
                    # No crossing pairs: intervals nest or stay disjoint.
                    self.assertFalse(a < c < b < d)

    def test_max_span_forbids_long_range_pairs(self):
        sequence = "GGGGGGAAAAAAAAAAAAAAAAAAAAAAAACCCCCC"
        loose = partition_fold(sequence)
        tight = partition_fold(sequence, max_span=10)
        self.assertGreater(max(loose.pair_prob[i][j] for i in range(len(sequence))
                               for j in range(i + 1, len(sequence))), 0.5)
        for i in range(len(sequence)):
            for j in range(i + 1, len(sequence)):
                if j - i > 10:
                    self.assertEqual(tight.pair_prob[i][j], 0.0)

    def test_empty_and_unfoldable_sequences(self):
        empty = partition_fold("")
        self.assertEqual(empty.partition_function, 1.0)
        self.assertEqual(empty.unpaired, [])

        flat = partition_fold("AAAAAAAAAA")
        self.assertAlmostEqual(flat.partition_function, 1.0, places=12)
        self.assertAlmostEqual(flat.ensemble_free_energy, 0.0, places=12)
        self.assertEqual(flat.unpaired, [1.0] * 10)

    def test_mean_unpaired_window(self):
        result = partition_fold("GGGGAAAACCCC")
        loop = result.mean_unpaired(4, 8)
        stem = result.mean_unpaired(0, 4)
        self.assertGreater(loop, stem)
        self.assertAlmostEqual(loop, sum(result.unpaired[4:8]) / 4, places=12)
        self.assertEqual(result.mean_unpaired(5, 5), 0.0)

    def test_dna_input_is_accepted(self):
        from_dna = partition_fold("GGGGTTTTCCCC")
        from_rna = partition_fold("GGGGUUUUCCCC")
        self.assertEqual(from_dna.sequence, "GGGGUUUUCCCC")
        self.assertAlmostEqual(from_dna.partition_function, from_rna.partition_function)


if __name__ == "__main__":
    unittest.main()
