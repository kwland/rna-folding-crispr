"""Statistics and ridge regression, using only the Python standard library.

Three things here are worth reading before trusting a number the study reports.

**Held-out evaluation is grouped by gene.** Guides that target the same gene
share chromatin state, expression, and local sequence composition. A random
split puts near-duplicates on both sides and inflates every score. Every fold in
this module holds out whole genes.

**Uncertainty is estimated by resampling genes, not guides.** For the same
reason, a confidence interval that treats 4,685 guides as 4,685 independent
observations is too narrow. ``cluster_bootstrap_ci`` resamples whole genes.

**Ridge fits are built from per-cell Gram matrices.** The normal-equation pieces
X'X and X'y are sums over samples, so they can be accumulated once per
cross-validation cell and then added or subtracted to get any training subset.
This makes nested cross-validation cheap enough to run in pure Python, and it is
exact rather than an approximation.
"""

from __future__ import annotations

import math
import random
from operator import mul


__all__ = [
    "rank_average",
    "spearman",
    "pearson",
    "mean_within_group_spearman",
    "within_group_scores",
    "within_group_spearman_ci",
    "benjamini_hochberg",
    "cluster_bootstrap_ci",
    "paired_spearman_delta_ci",
    "spearman_permutation_p",
    "assign_group_folds",
    "SparseDesign",
    "nested_cv_predictions",
    "select_alpha",
    "fit_full",
]


# --------------------------------------------------------------- correlation


def rank_average(values: list[float]) -> list[float]:
    """Ranks with ties replaced by their average rank."""
    count = len(values)
    order = sorted(range(count), key=lambda i: values[i])
    ranks = [0.0] * count
    position = 0
    while position < count:
        end = position
        while end + 1 < count and values[order[end + 1]] == values[order[position]]:
            end += 1
        shared = (position + end) / 2.0 + 1.0
        for index in range(position, end + 1):
            ranks[order[index]] = shared
        position = end + 1
    return ranks


def pearson(x: list[float], y: list[float]) -> float:
    count = len(x)
    if count != len(y):
        raise ValueError("pearson needs two equally long sequences")
    if count < 2:
        return 0.0
    mean_x = sum(x) / count
    mean_y = sum(y) / count
    dx = [value - mean_x for value in x]
    dy = [value - mean_y for value in y]
    denominator = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    if denominator == 0.0:
        return 0.0
    return sum(map(mul, dx, dy)) / denominator


def spearman(x: list[float], y: list[float]) -> float:
    """Spearman rank correlation, tie-corrected."""
    if len(x) < 2:
        return 0.0
    return pearson(rank_average(list(x)), rank_average(list(y)))


def mean_within_group_spearman(
    predicted: list[float],
    truth: list[float],
    groups: list[str],
    min_size: int = 8,
) -> float:
    """Average Spearman computed separately inside each group.

    This is the headline metric for every cross-validated model here, and the
    reason is worth stating plainly.

    Genes differ enormously in how editable they are. In this dataset the mean
    activity percentile of a held-out gene ranges from about 0.32 to 0.77. A
    model that has never seen a gene cannot know where that gene sits, so it
    predicts near the overall mean for all of its guides. Pooling out-of-fold
    predictions across genes then measures mostly one thing: whether the model
    guessed each gene's baseline level. It did not, so the pooled correlation
    collapses to roughly zero even when the model ranks guides *within* every
    gene at rho = 0.10 to 0.30.

    Ranking within a target is also the question a guide designer actually
    asks. They have one gene and need the best guide for it; they are never
    choosing between a guide for NF2 and a guide for HPRT1.

    Groups smaller than ``min_size`` are skipped, because a rank correlation on
    a handful of points is mostly noise.
    """
    members: dict[str, list[int]] = {}
    for index, group in enumerate(groups):
        members.setdefault(group, []).append(index)

    scores = []
    for indices in members.values():
        if len(indices) < min_size:
            continue
        scores.append(
            spearman([predicted[i] for i in indices], [truth[i] for i in indices])
        )
    if not scores:
        return float("nan")
    return sum(scores) / len(scores)


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values, in the input order.

    The position-resolved analysis runs one test per spacer position, so some
    position will look interesting by chance. This is the correction for that.
    """
    count = len(pvalues)
    if count == 0:
        return []
    order = sorted(range(count), key=lambda i: pvalues[i])
    adjusted = [0.0] * count
    previous = 1.0
    for rank in range(count - 1, -1, -1):
        index = order[rank]
        value = pvalues[index] * count / (rank + 1)
        previous = min(previous, value)
        adjusted[index] = min(1.0, previous)
    return adjusted


# ----------------------------------------------------------------- resampling


def cluster_bootstrap_ci(
    statistic,
    groups: list[str],
    n_boot: int = 2000,
    seed: int = 0,
    level: float = 0.95,
) -> tuple[float, float]:
    """Percentile confidence interval, resampling whole groups with replacement.

    ``statistic`` is called with a list of sample indices and returns a number.
    """
    index_by_group: dict[str, list[int]] = {}
    for index, group in enumerate(groups):
        index_by_group.setdefault(group, []).append(index)
    names = list(index_by_group)
    if len(names) < 2:
        return (float("nan"), float("nan"))

    rng = random.Random(seed)
    draws = []
    for _ in range(n_boot):
        indices: list[int] = []
        for _ in range(len(names)):
            indices.extend(index_by_group[names[rng.randrange(len(names))]])
        try:
            value = statistic(indices)
        except (ValueError, ZeroDivisionError):
            continue
        if value == value:  # skip NaN
            draws.append(value)
    if not draws:
        return (float("nan"), float("nan"))
    draws.sort()
    tail = (1.0 - level) / 2.0
    low = draws[max(0, int(math.floor(tail * len(draws))))]
    high = draws[min(len(draws) - 1, int(math.ceil((1.0 - tail) * len(draws))) - 1)]
    return (low, high)


def within_group_scores(
    predicted: list[float],
    truth: list[float],
    groups: list[str],
    min_size: int = 8,
) -> dict[str, float]:
    """Spearman inside each group, keyed by group name."""
    members: dict[str, list[int]] = {}
    for index, group in enumerate(groups):
        members.setdefault(group, []).append(index)
    return {
        group: spearman([predicted[i] for i in indices], [truth[i] for i in indices])
        for group, indices in members.items()
        if len(indices) >= min_size
    }


def _bootstrap_mean_ci(
    values: list[float], n_boot: int, seed: int, level: float
) -> tuple[float, float]:
    """Percentile interval for the mean of a list of independent values."""
    if len(values) < 2:
        return (float("nan"), float("nan"))
    rng = random.Random(seed)
    count = len(values)
    draws = []
    for _ in range(n_boot):
        total = 0.0
        for _ in range(count):
            total += values[rng.randrange(count)]
        draws.append(total / count)
    draws.sort()
    tail = (1.0 - level) / 2.0
    low = draws[max(0, int(math.floor(tail * len(draws))))]
    high = draws[min(len(draws) - 1, int(math.ceil((1.0 - tail) * len(draws))) - 1)]
    return (low, high)


def within_group_spearman_ci(
    predicted: list[float],
    truth: list[float],
    groups: list[str],
    n_boot: int = 2000,
    seed: int = 0,
    level: float = 0.95,
    min_size: int = 8,
) -> tuple[float, float, float, int]:
    """Mean within-group Spearman with an interval from resampling groups."""
    scores = list(within_group_scores(predicted, truth, groups, min_size).values())
    if not scores:
        return (float("nan"), float("nan"), float("nan"), 0)
    mean = sum(scores) / len(scores)
    low, high = _bootstrap_mean_ci(scores, n_boot, seed, level)
    return (mean, low, high, len(scores))


def paired_spearman_delta_ci(
    baseline: list[float],
    augmented: list[float],
    truth: list[float],
    groups: list[str],
    n_boot: int = 2000,
    seed: int = 0,
    level: float = 0.95,
    min_size: int = 8,
) -> tuple[float, float, float]:
    """Confidence interval for the *change* in accuracy from adding features.

    Both models are scored inside the same genes, so the shared difficulty of a
    gene cancels and what remains is the paired per-gene difference. Resampling
    those differences is a bootstrap of a mean over independent units, which is
    both cleaner and cheaper than resampling guides and recomputing.

    This is the quantity the study actually cares about: not whether structure
    correlates with activity, but whether it adds anything a sequence-only model
    did not already have.
    """
    base = within_group_scores(baseline, truth, groups, min_size)
    aug = within_group_scores(augmented, truth, groups, min_size)
    shared = sorted(set(base) & set(aug))
    deltas = [aug[group] - base[group] for group in shared]
    if not deltas:
        return (float("nan"), float("nan"), float("nan"))
    observed = sum(deltas) / len(deltas)
    low, high = _bootstrap_mean_ci(deltas, n_boot, seed, level)
    return (observed, low, high)


def spearman_permutation_p(
    x: list[float],
    y: list[float],
    n_perm: int = 5000,
    seed: int = 0,
) -> float:
    """Two-sided permutation p-value for a Spearman correlation.

    CAUTION: this shuffles individual guides, which destroys the gene structure
    in the data. Guides targeting one gene are correlated, so the effective
    sample size is closer to the number of genes than the number of guides, and
    a free shuffle therefore produces a null distribution that is too narrow.
    The p-values it returns are anti-conservative and should be read as a rough
    screen only. The cluster bootstrap intervals, which resample whole genes,
    are the trustworthy uncertainty statement, and they are what the study's
    conclusions rest on.
    """
    observed = abs(spearman(x, y))
    shuffled = list(y)
    rng = random.Random(seed)
    hits = 0
    for _ in range(n_perm):
        rng.shuffle(shuffled)
        if abs(spearman(x, shuffled)) >= observed:
            hits += 1
    return (hits + 1) / (n_perm + 1)


# ------------------------------------------------------------ linear algebra


def cholesky_solve(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """Solve a symmetric positive-definite system by Cholesky decomposition.

    The inner products run through ``sum(map(mul, ...))`` so the O(p^3) work
    happens in C rather than in interpreted loops.
    """
    size = len(vector)
    lower = [[0.0] * size for _ in range(size)]
    for i in range(size):
        row_i = lower[i]
        source = matrix[i]
        for j in range(i + 1):
            row_j = lower[j]
            total = source[j] - sum(map(mul, row_i[:j], row_j[:j]))
            if i == j:
                if total <= 0.0:
                    # Only reachable if the ridge penalty is too small to make
                    # the system positive definite; nudge it and continue.
                    total = 1e-12
                row_i[j] = math.sqrt(total)
            else:
                row_i[j] = total / row_j[j]

    # Forward substitution, then back substitution.
    forward = [0.0] * size
    for i in range(size):
        row_i = lower[i]
        forward[i] = (vector[i] - sum(map(mul, row_i[:i], forward[:i]))) / row_i[i]

    solution = [0.0] * size
    for i in range(size - 1, -1, -1):
        total = forward[i]
        for k in range(i + 1, size):
            total -= lower[k][i] * solution[k]
        solution[i] = total / lower[i][i]
    return solution


# ------------------------------------------------------- cross-validated ridge


def assign_group_folds(groups: list[str], n_folds: int, seed: int = 0) -> list[int]:
    """Assign each sample a fold index, keeping every group intact.

    Groups are shuffled and then handed to whichever fold currently holds the
    fewest samples, which keeps folds close to equal size even when one gene
    contributes far more guides than another.
    """
    counts: dict[str, int] = {}
    for group in groups:
        counts[group] = counts.get(group, 0) + 1
    names = sorted(counts)
    random.Random(seed).shuffle(names)
    names.sort(key=lambda name: -counts[name])

    load = [0] * n_folds
    fold_of_group: dict[str, int] = {}
    for name in names:
        target = min(range(n_folds), key=lambda f: load[f])
        fold_of_group[name] = target
        load[target] += counts[name]
    return [fold_of_group[group] for group in groups]


class SparseDesign:
    """A design matrix stored as sparse rows of (feature index, value).

    Guide features are mostly one-hot: 20 position indicators and 19
    dinucleotide indicators out of several hundred columns. Storing and
    multiplying only the nonzero entries turns the Gram accumulation from
    hundreds of millions of operations into a few million.
    """

    def __init__(self, rows: list[list[tuple[int, float]]], n_features: int) -> None:
        self.rows = rows
        self.n_features = n_features

    def __len__(self) -> int:
        return len(self.rows)

    def predict(self, weights: list[float], indices: list[int] | None = None) -> list[float]:
        chosen = range(len(self.rows)) if indices is None else indices
        out = []
        for index in chosen:
            out.append(sum(weights[f] * v for f, v in self.rows[index]))
        return out


def _cell_gram(
    design: SparseDesign,
    y: list[float],
    indices: list[int],
) -> tuple[list[list[float]], list[float]]:
    """Accumulate X'X and X'y over a subset of rows."""
    p = design.n_features
    gram = [[0.0] * p for _ in range(p)]
    rhs = [0.0] * p
    rows = design.rows
    for index in indices:
        row = rows[index]
        target = y[index]
        for a, (fa, va) in enumerate(row):
            gram_fa = gram[fa]
            rhs[fa] += va * target
            for fb, vb in row[: a + 1]:
                gram_fa[fb] += va * vb
    # Mirror the lower triangle into the upper triangle.
    for i in range(p):
        gram_i = gram[i]
        for j in range(i):
            gram[j][i] = gram_i[j]
    return gram, rhs


def _add_cells(cells, keys, p: int):
    gram = [[0.0] * p for _ in range(p)]
    rhs = [0.0] * p
    for key in keys:
        cell_gram, cell_rhs = cells[key]
        for i in range(p):
            row = gram[i]
            source = cell_gram[i]
            for j in range(p):
                row[j] += source[j]
            rhs[i] += cell_rhs[i]
    return gram, rhs


def _subtract(left, right, p: int):
    """(X'X, X'y) for a training set, as the whole minus the held-out part."""
    gram_l, rhs_l = left
    gram_r, rhs_r = right
    gram = [
        [gram_l[i][j] - gram_r[i][j] for j in range(p)]
        for i in range(p)
    ]
    rhs = [rhs_l[i] - rhs_r[i] for i in range(p)]
    return gram, rhs


def _fit(gram, rhs, alpha: float, p: int, intercept_index: int) -> list[float]:
    regularised = [row[:] for row in gram]
    for i in range(p):
        if i != intercept_index:
            regularised[i][i] += alpha
    return cholesky_solve(regularised, rhs)


def nested_cv_predictions(
    design: SparseDesign,
    y: list[float],
    groups: list[str],
    alphas: list[float],
    n_outer: int = 5,
    n_inner: int = 4,
    seed: int = 0,
    intercept_index: int = 0,
) -> tuple[list[float], list[float]]:
    """Out-of-fold predictions from ridge regression with nested model selection.

    The ridge penalty is chosen inside each training fold, never on the fold
    being predicted, so the returned predictions are honestly out-of-sample.
    Returns the predictions and the penalty chosen for each outer fold.
    """
    p = design.n_features
    outer = assign_group_folds(groups, n_outer, seed=seed)
    inner = assign_group_folds(groups, n_inner, seed=seed + 977)

    # One Gram per (outer, inner) cell; every training subset below is a sum of cells.
    cells: dict[tuple[int, int], tuple[list[list[float]], list[float]]] = {}
    members: dict[tuple[int, int], list[int]] = {}
    for index in range(len(design)):
        members.setdefault((outer[index], inner[index]), []).append(index)
    for key, indices in members.items():
        cells[key] = _cell_gram(design, y, indices)

    predictions = [0.0] * len(design)
    chosen_alphas = []
    for fold in range(n_outer):
        train_keys = [key for key in cells if key[0] != fold]
        train_system = _add_cells(cells, train_keys, p)

        # The Gram matrix does not depend on the ridge penalty, so build each
        # inner training system once and reuse it across the whole alpha grid.
        inner_systems = {}
        for inner_fold in range(n_inner):
            held_keys = [key for key in train_keys if key[1] == inner_fold]
            held_indices = [i for key in held_keys for i in members[key]]
            if not held_indices or len(held_keys) == len(train_keys):
                continue
            held_system = _add_cells(cells, held_keys, p)
            inner_systems[inner_fold] = (
                _subtract(train_system, held_system, p),
                held_indices,
            )

        best_alpha, best_score = alphas[0], -2.0
        for alpha in alphas:
            scores = []
            for (gram, rhs), held_indices in inner_systems.values():
                weights = _fit(gram, rhs, alpha, p, intercept_index)
                predicted = design.predict(weights, held_indices)
                truth = [y[i] for i in held_indices]
                scores.append(spearman(predicted, truth))
            if scores:
                mean_score = sum(scores) / len(scores)
                if mean_score > best_score:
                    best_alpha, best_score = alpha, mean_score
        chosen_alphas.append(best_alpha)

        weights = _fit(train_system[0], train_system[1], best_alpha, p, intercept_index)
        test_indices = [i for i in range(len(design)) if outer[i] == fold]
        for index, value in zip(test_indices, design.predict(weights, test_indices)):
            predictions[index] = value

    return predictions, chosen_alphas


def select_alpha(
    design: SparseDesign,
    y: list[float],
    groups: list[str],
    alphas: list[float],
    n_folds: int = 4,
    seed: int = 0,
    intercept_index: int = 0,
) -> float:
    """Choose a ridge penalty by grouped cross-validation within one screen.

    Used before transferring a model to a different screen: the penalty has to
    come from the training screen alone, never from the screen being predicted.
    """
    p = design.n_features
    folds = assign_group_folds(groups, n_folds, seed=seed)
    members: dict[int, list[int]] = {}
    for index, fold in enumerate(folds):
        members.setdefault(fold, []).append(index)
    cells = {fold: _cell_gram(design, y, indices) for fold, indices in members.items()}
    total = _add_cells(cells, list(cells), p)

    best_alpha, best_score = alphas[0], -2.0
    systems = {
        fold: (_subtract(total, cells[fold], p), members[fold])
        for fold in members
        if len(members) > 1
    }
    for alpha in alphas:
        scores = []
        for (gram, rhs), held in systems.values():
            weights = _fit(gram, rhs, alpha, p, intercept_index)
            scores.append(spearman(design.predict(weights, held), [y[i] for i in held]))
        if scores:
            mean_score = sum(scores) / len(scores)
            if mean_score > best_score:
                best_alpha, best_score = alpha, mean_score
    return best_alpha


def fit_full(
    design: SparseDesign,
    y: list[float],
    alpha: float,
    intercept_index: int = 0,
) -> list[float]:
    """Fit ridge on every row; used only for transferring a model between screens."""
    gram, rhs = _cell_gram(design, y, list(range(len(design))))
    return _fit(gram, rhs, alpha, design.n_features, intercept_index)
