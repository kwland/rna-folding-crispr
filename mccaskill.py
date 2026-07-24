"""McCaskill partition function: base-pair probabilities over the whole ensemble.

Why this module exists
----------------------
Nussinov and Zuker each return **one** structure. Asking "is the seed unpaired?"
of a single structure gives a 0/1 answer that can flip on an energy difference of
a tenth of a kcal/mol. An RNA molecule does not sit in one structure; it moves
through a Boltzmann-weighted ensemble of them.

McCaskill (1990) computes, exactly and in the same O(n^3) dynamic-programming
family as Zuker, the probability that each pair (i, j) is formed:

    p_ij = (inside weight of i..j given i-j paired)
           * (outside weight of everything else)
           / Q

The probability that base i is unpaired is then ``1 - sum_j p_ij``. Seed
accessibility becomes a continuous number - the mean unpaired probability over
the seed - rather than a count of positions that happened to be unpaired in one
arbitrarily chosen structure.

Structure of the algorithm
--------------------------
Inside (partition functions, computed over increasing span):

    QB[i][j]  = sum over structures on i..j **given (i, j) pair**
    QM1[i][j] = one multiloop branch that starts exactly at i, then unpaired to j
    QM[i][j]  = at least one multiloop branch somewhere in i..j
    Qf[t]     = exterior partition function of the prefix 0..t-1
    Qr[t]     = exterior partition function of the suffix t..n-1

Outside (pair probabilities, computed over decreasing span). A pair (i, j) is
formed in exactly one of three contexts, so p_ij is a sum of three terms:

    1. exterior      - nothing encloses it
    2. interior loop - a closer pair (k, l) encloses it with no other branch
    3. multiloop     - a closer pair (k, l) encloses it alongside other branches

Term 3 is the awkward one: written naively it is a sum over every enclosing pair,
which is O(n^4). The standard fix, used here, is to notice that the enclosing
sum factorises into two accumulator tables (``acc_unpaired`` and ``acc_branch``)
that are updated once per finalised pair, which brings the whole outside pass
back to O(n^3).

Correctness
-----------
``test_mccaskill.py`` checks the recursions against brute-force enumeration of
every pseudoknot-free structure for short sequences: the partition function, the
ensemble free energy, and every single pair probability must match to numerical
tolerance. A separate test checks that stochastic sampling from the same inside
matrices converges to the same probabilities on a longer sequence.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Sequence

from energy_model import (
    MAX_INTERNAL,
    MIN_LOOP,
    ML_BRANCH,
    ML_CLOSE,
    ML_UNPAIRED,
    RT37,
    BoltzmannTables,
    EnergyModel,
)


__all__ = ["EnsembleResult", "partition_fold", "sample_structures"]


@dataclass
class EnsembleResult:
    """Everything the partition function knows about one sequence."""

    sequence: str
    partition_function: float
    ensemble_free_energy: float
    pair_prob: list[list[float]]
    unpaired: list[float]
    #: Inside matrices, retained so a caller can draw Boltzmann samples.
    _qb: list[list[float]] = field(repr=False, default_factory=list)
    _qm: list[list[float]] = field(repr=False, default_factory=list)
    _qm1: list[list[float]] = field(repr=False, default_factory=list)
    _qr: list[float] = field(repr=False, default_factory=list)
    _model: EnergyModel | None = field(repr=False, default=None)

    @property
    def n(self) -> int:
        return len(self.sequence)

    def pair_probability(self, i: int, j: int) -> float:
        if i > j:
            i, j = j, i
        return self.pair_prob[i][j]

    def mean_unpaired(self, start: int, end: int) -> float:
        """Mean unpaired probability over the half-open window [start, end).

        This is the ensemble replacement for "fraction of seed positions that
        were unpaired in the one structure we looked at".
        """
        window = self.unpaired[start:end]
        if not window:
            return 0.0
        return sum(window) / len(window)

    def expected_paired_bases(self) -> float:
        """Expected number of paired bases across the ensemble."""
        return sum(1.0 - p for p in self.unpaired)


def _boltzmann(energy: float, rt: float) -> float:
    if energy == math.inf:
        return 0.0
    return math.exp(-energy / rt)


def partition_fold(
    sequence: str,
    allow_wobble: bool = True,
    temperature_rt: float = RT37,
    max_span: int | None = None,
) -> EnsembleResult:
    """Run the McCaskill inside and outside recursions.

    ``max_span`` optionally forbids pairs spanning more than that many bases,
    which is the local-folding restriction RNAplfold uses to model a molecule
    that folds while it is still being transcribed. ``None`` means global folding.
    """
    model = EnergyModel(sequence, allow_wobble=allow_wobble)
    seq = model.seq
    n = model.n
    rt = temperature_rt

    if n == 0:
        return EnsembleResult("", 1.0, 0.0, [], [])

    span_limit = n - 1 if max_span is None else min(max_span, n - 1)

    def pairable(i: int, j: int) -> bool:
        return j - i <= span_limit and model.can_pair(i, j)

    # Boltzmann weights that get reused constantly, computed once.
    z_unpaired = math.exp(-ML_UNPAIRED / rt)          # per unpaired base in a multiloop
    z_pow = [z_unpaired ** k for k in range(n + 2)]   # z_unpaired ** k, tabulated
    ml_open = math.exp(-(ML_CLOSE + ML_BRANCH) / rt)  # closing a multiloop
    ml_branch = math.exp(-ML_BRANCH / rt)             # one more branch in a multiloop

    # Loop weights, tabulated once. The recursions below inline the lookups from
    # BoltzmannTables.interior rather than calling it, because that call sits in
    # the innermost loop and the call overhead alone dominated the runtime.
    tables = BoltzmannTables(model, rt)
    term = tables.term
    int_term = tables.int_term
    w_hairpin = tables.hairpin
    w_stack = tables.stack
    w_bulge1_l = tables.bulge1_left
    w_bulge1_r = tables.bulge1_right
    w_bulge = tables.bulge
    w_internal = tables.internal

    # ------------------------------------------------------------ inside pass

    qb = [[0.0] * n for _ in range(n)]
    qm = [[0.0] * n for _ in range(n)]
    qm1 = [[0.0] * n for _ in range(n)]

    for span in range(MIN_LOOP + 1, n):
        for i in range(0, n - span):
            j = i + span

            # ---- QB: structures on i..j given that i pairs with j
            if pairable(i, j):
                total = w_hairpin[i][j]

                # stack / bulge / internal loop with one inner pair (k, l)
                k_max = min(i + MAX_INTERNAL + 1, j - MIN_LOOP - 2)
                term_ij = term[i][j]
                int_term_ij = int_term[i][j]
                for k in range(i + 1, k_max + 1):
                    left = k - i - 1
                    qb_k = qb[k]
                    term_k = term[k]
                    l_min = max(k + MIN_LOOP + 1, j - 1 - (MAX_INTERNAL - left))
                    if left == 0:
                        # Nothing unpaired on the left: stack, right bulge, or bulge.
                        for l in range(j - 1, l_min - 1, -1):
                            inner = qb_k[l]
                            if inner == 0.0:
                                continue
                            right = j - l - 1
                            if right == 0:
                                total += inner * w_stack[i][j]
                            elif right == 1:
                                total += inner * w_bulge1_r[i][j]
                            else:
                                total += inner * w_bulge[right] * term_ij * term_k[l]
                    else:
                        w_int_row = w_internal[left]
                        int_term_k = int_term[k]
                        for l in range(j - 1, l_min - 1, -1):
                            inner = qb_k[l]
                            if inner == 0.0:
                                continue
                            right = j - l - 1
                            if right == 0:
                                if left == 1:
                                    total += inner * w_bulge1_l[i][j]
                                else:
                                    total += inner * w_bulge[left] * term_ij * term_k[l]
                            else:
                                total += inner * w_int_row[right] * int_term_ij * int_term_k[l]

                # multiloop: at least two branches inside, last one starting at h
                if span >= 2 * (MIN_LOOP + 2) + 1:
                    ml = 0.0
                    qm_row = qm[i + 1]
                    for h in range(i + 2, j - MIN_LOOP - 1):
                        left = qm_row[h - 1]
                        if left == 0.0:
                            continue
                        right = qm1[h][j - 1]
                        if right == 0.0:
                            continue
                        ml += left * right
                    if ml:
                        total += ml * ml_open * term[i][j]

                qb[i][j] = total

            # ---- QM1: exactly one branch, opening at i, trailing bases unpaired
            value = qm1[i][j - 1] * z_unpaired if j - 1 >= i else 0.0
            if qb[i][j]:
                value += qb[i][j] * ml_branch * term[i][j]
            qm1[i][j] = value

            # ---- QM: at least one branch anywhere in i..j
            acc = 0.0
            qm_i = qm[i]
            for h in range(i, j - MIN_LOOP):
                branch = qm1[h][j]
                if branch == 0.0:
                    continue
                prefix = z_pow[h - i] + (qm_i[h - 1] if h - 1 >= i else 0.0)
                acc += prefix * branch
            qm_i[j] = acc

    # ---- exterior: prefix and suffix partition functions
    qf = [1.0] * (n + 1)  # qf[t] covers bases 0..t-1
    for t in range(1, n + 1):
        j = t - 1
        value = qf[t - 1]
        for h in range(0, j - MIN_LOOP):
            if qb[h][j]:
                value += qf[h] * qb[h][j] * term[h][j]
        qf[t] = value

    qr = [1.0] * (n + 2)  # qr[t] covers bases t..n-1
    for t in range(n - 1, -1, -1):
        value = qr[t + 1]
        qb_t = qb[t]
        term_t = term[t]
        for l in range(t + MIN_LOOP + 1, n):
            if qb_t[l]:
                value += qb_t[l] * term_t[l] * qr[l + 1]
        qr[t] = value

    q_total = qf[n]
    if not math.isfinite(q_total) or q_total <= 0.0:
        raise OverflowError(
            f"Partition function is not usable ({q_total}). The sequence is likely too long "
            "for unscaled double precision; fold a shorter window or set max_span."
        )
    # Two independent routes to the same number; disagreement means a broken recursion.
    if abs(qr[0] - q_total) > 1e-6 * max(1.0, q_total):
        raise AssertionError("Prefix and suffix partition functions disagree")

    # ----------------------------------------------------------- outside pass

    prob = [[0.0] * n for _ in range(n)]

    # Accumulators for the multiloop term. For an enclosing pair (k, l) with
    # weight w:
    #   acc_unpaired[k][q] += w * z^(l-1-q)      -> everything right of q is unpaired
    #   acc_branch[k][q]   += w * QM[q+1][l-1]   -> at least one more branch right of q
    # Folding each finalised pair into these tables costs O(n), so the whole
    # multiloop term costs O(n^3) instead of the naive O(n^4).
    acc_unpaired = [[0.0] * n for _ in range(n)]
    acc_branch = [[0.0] * n for _ in range(n)]
    pending: list[tuple[int, int, float]] = []

    for span in range(n - 1, MIN_LOOP, -1):
        # Pairs finalised at the previous (wider) span can now enclose things.
        for k, l, weight in pending:
            au_k = acc_unpaired[k]
            ab_k = acc_branch[k]
            for q in range(k + 1, l):
                au_k[q] += weight * z_pow[l - 1 - q]
                if q + 1 <= l - 1:
                    inner = qm[q + 1][l - 1]
                    if inner:
                        ab_k[q] += weight * inner
        pending = []

        for i in range(0, n - span):
            j = i + span
            qb_ij = qb[i][j]
            if qb_ij == 0.0:
                continue

            # 1. exterior context
            total = qf[i] * qb_ij * term[i][j] * qr[j + 1] / q_total

            # 2. enclosed by (k, l) as the only branch: stack, bulge or internal loop
            k_lo = max(0, i - MAX_INTERNAL - 1)
            term_ij = term[i][j]
            int_term_ij = int_term[i][j]
            for k in range(k_lo, i):
                left = i - k - 1
                l_hi = min(n - 1, j + 1 + (MAX_INTERNAL - left))
                prob_k = prob[k]
                qb_k = qb[k]
                term_k = term[k]
                int_term_k = int_term[k]
                w_int_row = w_internal[left] if left else None
                for l in range(j + 1, l_hi + 1):
                    outer = prob_k[l]
                    if outer == 0.0:
                        continue
                    right = l - j - 1
                    if left == 0:
                        if right == 0:
                            weight = w_stack[k][l]
                        elif right == 1:
                            weight = w_bulge1_r[k][l]
                        else:
                            weight = w_bulge[right] * term_k[l] * term_ij
                    elif right == 0:
                        if left == 1:
                            weight = w_bulge1_l[k][l]
                        else:
                            weight = w_bulge[left] * term_k[l] * term_ij
                    else:
                        weight = w_int_row[right] * int_term_k[l] * int_term_ij
                    total += (outer / qb_k[l]) * weight * qb_ij

            # 3. one branch of a multiloop closed by some enclosing pair
            multi = 0.0
            for k in range(0, i):
                branch = acc_branch[k][j]
                unpaired_right = acc_unpaired[k][j]
                if branch == 0.0 and unpaired_right == 0.0:
                    continue
                left_branches = qm[k + 1][i - 1] if k + 1 <= i - 1 else 0.0
                if left_branches:
                    multi += left_branches * (unpaired_right + branch)
                if branch:
                    multi += z_pow[i - 1 - k] * branch
            if multi:
                total += qb_ij * ml_branch * term[i][j] * multi

            if total <= 0.0:
                continue
            # Numerical guard: probabilities may drift a hair above 1.0.
            prob[i][j] = min(total, 1.0)
            pending.append((i, j, prob[i][j] / qb_ij * ml_open * term[i][j]))

    unpaired = [1.0] * n
    for i in range(n):
        used = 0.0
        for j in range(i + MIN_LOOP + 1, n):
            used += prob[i][j]
        for k in range(0, i - MIN_LOOP):
            used += prob[k][i]
        unpaired[i] = max(0.0, min(1.0, 1.0 - used))

    return EnsembleResult(
        sequence=seq,
        partition_function=q_total,
        ensemble_free_energy=-rt * math.log(q_total),
        pair_prob=prob,
        unpaired=unpaired,
        _qb=qb,
        _qm=qm,
        _qm1=qm1,
        _qr=qr,
        _model=model,
    )


# --------------------------------------------------------------- sampling


def sample_structures(
    result: EnsembleResult,
    count: int,
    rng: random.Random | None = None,
    temperature_rt: float = RT37,
) -> list[list[tuple[int, int]]]:
    """Draw structures from the Boltzmann ensemble by stochastic traceback.

    This is the Ding & Lawrence sampling idea: walk the same inside matrices the
    partition function built, but at every branch point choose an option with
    probability proportional to its weight instead of taking the maximum. The
    samples are exact draws from the ensemble, so their pair frequencies converge
    to the pair probabilities computed by the outside pass - which is how
    ``test_mccaskill.py`` cross-checks the two halves of the algorithm.
    """
    if result._model is None:
        raise ValueError("This result was not produced by partition_fold")
    model = result._model
    rt = temperature_rt
    n = result.n
    qb, qm, qm1, qr = result._qb, result._qm, result._qm1, result._qr
    rng = rng or random.Random(0)

    z_unpaired = math.exp(-ML_UNPAIRED / rt)
    z_pow = [z_unpaired ** k for k in range(n + 2)]
    ml_open = math.exp(-(ML_CLOSE + ML_BRANCH) / rt)
    ml_branch = math.exp(-ML_BRANCH / rt)

    def term(i: int, j: int) -> float:
        return math.exp(-model.terminal_penalty(i, j) / rt)

    def pick(options: list[tuple[float, object]], total: float) -> object:
        """Choose one option with probability proportional to its weight."""
        target = rng.random() * total
        running = 0.0
        for weight, payload in options:
            running += weight
            if running >= target:
                return payload
        return options[-1][1]

    def trace_pair(i: int, j: int, pairs: list[tuple[int, int]]) -> None:
        pairs.append((i, j))
        options: list[tuple[float, object]] = [
            (_boltzmann(model.hairpin(i, j), rt), ("hairpin",))
        ]
        k_max = min(i + MAX_INTERNAL + 1, j - MIN_LOOP - 2)
        for k in range(i + 1, k_max + 1):
            left = k - i - 1
            l_min = max(k + MIN_LOOP + 1, j - 1 - (MAX_INTERNAL - left))
            for l in range(j - 1, l_min - 1, -1):
                if qb[k][l] == 0.0:
                    continue
                weight = qb[k][l] * _boltzmann(model.interior(i, j, k, l), rt)
                if weight:
                    options.append((weight, ("interior", k, l)))
        for h in range(i + 2, j - MIN_LOOP - 1):
            left = qm[i + 1][h - 1]
            right = qm1[h][j - 1]
            if left and right:
                options.append((left * right * ml_open * term(i, j), ("multi", h)))

        total = sum(weight for weight, _ in options)
        choice = pick(options, total)
        if choice[0] == "hairpin":
            return
        if choice[0] == "interior":
            trace_pair(choice[1], choice[2], pairs)
            return
        h = choice[1]
        trace_qm(i + 1, h - 1, pairs)
        trace_qm1(h, j - 1, pairs)

    def trace_qm1(i: int, j: int, pairs: list[tuple[int, int]]) -> None:
        options = []
        for l in range(i + MIN_LOOP + 1, j + 1):
            if qb[i][l]:
                options.append((qb[i][l] * ml_branch * term(i, l) * z_pow[j - l], l))
        total = sum(weight for weight, _ in options)
        l = pick(options, total)
        trace_pair(i, l, pairs)

    def trace_qm(i: int, j: int, pairs: list[tuple[int, int]]) -> None:
        options = []
        for h in range(i, j - MIN_LOOP):
            branch = qm1[h][j]
            if branch == 0.0:
                continue
            left = qm[i][h - 1] if h - 1 >= i else 0.0
            options.append((z_pow[h - i] * branch, ("empty", h)))
            if left:
                options.append((left * branch, ("branch", h)))
        total = sum(weight for weight, _ in options)
        kind, h = pick(options, total)
        if kind == "branch":
            trace_qm(i, h - 1, pairs)
        trace_qm1(h, j, pairs)

    def trace_exterior(pairs: list[tuple[int, int]]) -> None:
        t = 0
        while t < n:
            total = qr[t]
            options: list[tuple[float, object]] = [(qr[t + 1], ("unpaired",))]
            for l in range(t + MIN_LOOP + 1, n):
                if qb[t][l]:
                    options.append((qb[t][l] * term(t, l) * qr[l + 1], ("pair", l)))
            choice = pick(options, total)
            if choice[0] == "unpaired":
                t += 1
            else:
                l = choice[1]
                trace_pair(t, l, pairs)
                t = l + 1

    samples = []
    for _ in range(count):
        pairs: list[tuple[int, int]] = []
        trace_exterior(pairs)
        pairs.sort()
        samples.append(pairs)
    return samples


def unpaired_from_samples(samples: Sequence[Sequence[tuple[int, int]]], n: int) -> list[float]:
    """Per-base unpaired frequency across sampled structures."""
    if not samples:
        return [1.0] * n
    counts = [0] * n
    for pairs in samples:
        for i, j in pairs:
            counts[i] += 1
            counts[j] += 1
    return [1.0 - counts[i] / len(samples) for i in range(n)]
