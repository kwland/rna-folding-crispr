"""Zuker minimum-free-energy folding, sharing the energy model with McCaskill.

This exists so that "single structure" and "whole ensemble" can be compared
without changing anything else. Nussinov answers a different question (most
pairs, not lowest energy), so comparing Nussinov to McCaskill would confound two
changes at once: the objective *and* the single-structure/ensemble distinction.
Running Zuker and McCaskill over the identical parameters in ``energy_model``
isolates the second one, which is the variable this project is actually testing.

    V[i][j]  = MFE of i..j given that i pairs with j
    WM[i][j] = MFE of i..j inside a multiloop, at least one branch
    W[i][j]  = MFE of i..j with no constraint (exterior context)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from energy_model import (
    MAX_INTERNAL,
    MIN_LOOP,
    ML_BRANCH,
    ML_CLOSE,
    ML_UNPAIRED,
    EnergyModel,
)
from nussinov import dot_bracket


__all__ = ["MfeResult", "zuker_fold"]

_EPS = 1e-7  # float tolerance for traceback equality tests


@dataclass(frozen=True)
class MfeResult:
    sequence: str
    energy: float
    structure: str
    pairs: list[tuple[int, int]]

    @property
    def n(self) -> int:
        return len(self.sequence)

    def unpaired_flags(self) -> list[float]:
        """1.0 where a base is unpaired in this one structure, else 0.0.

        Deliberately the same shape as ``EnsembleResult.unpaired`` so the two
        accessibility definitions are interchangeable downstream. The difference
        is that these are hard 0/1 values, which is exactly the coarseness the
        ensemble measure is meant to replace.
        """
        flags = [1.0] * self.n
        for i, j in self.pairs:
            flags[i] = 0.0
            flags[j] = 0.0
        return flags


def zuker_fold(sequence: str, allow_wobble: bool = True) -> MfeResult:
    """Fold one sequence to its minimum free energy structure."""
    model = EnergyModel(sequence, allow_wobble=allow_wobble)
    n = model.n
    if n == 0:
        return MfeResult("", 0.0, "", [])

    inf = math.inf
    v = [[inf] * n for _ in range(n)]
    wm = [[inf] * n for _ in range(n)]
    w = [[0.0] * n for _ in range(n)]

    for span in range(1, n):
        for i in range(0, n - span):
            j = i + span

            # ---- V
            best = inf
            if model.can_pair(i, j):
                best = model.hairpin(i, j)

                k_max = min(i + MAX_INTERNAL + 1, j - MIN_LOOP - 2)
                for k in range(i + 1, k_max + 1):
                    left = k - i - 1
                    l_min = max(k + MIN_LOOP + 1, j - 1 - (MAX_INTERNAL - left))
                    for l in range(j - 1, l_min - 1, -1):
                        if v[k][l] == inf:
                            continue
                        candidate = model.interior(i, j, k, l) + v[k][l]
                        if candidate < best:
                            best = candidate

                closing = ML_CLOSE + ML_BRANCH + model.terminal_penalty(i, j)
                for m in range(i + 2, j - 1):
                    left = wm[i + 1][m]
                    right = wm[m + 1][j - 1]
                    if left == inf or right == inf:
                        continue
                    candidate = left + right + closing
                    if candidate < best:
                        best = candidate
            v[i][j] = best

            # ---- WM
            best = inf
            if v[i][j] != inf:
                best = v[i][j] + ML_BRANCH + model.terminal_penalty(i, j)
            if wm[i + 1][j] != inf and wm[i + 1][j] + ML_UNPAIRED < best:
                best = wm[i + 1][j] + ML_UNPAIRED
            if wm[i][j - 1] != inf and wm[i][j - 1] + ML_UNPAIRED < best:
                best = wm[i][j - 1] + ML_UNPAIRED
            for q in range(i, j):
                if wm[i][q] == inf or wm[q + 1][j] == inf:
                    continue
                candidate = wm[i][q] + wm[q + 1][j]
                if candidate < best:
                    best = candidate
            wm[i][j] = best

            # ---- W (exterior loop: unpaired bases are free)
            best = 0.0
            if w[i + 1][j] < best:
                best = w[i + 1][j]
            if w[i][j - 1] < best:
                best = w[i][j - 1]
            if v[i][j] != inf:
                closed = v[i][j] + model.terminal_penalty(i, j)
                if closed < best:
                    best = closed
            for r in range(i, j):
                candidate = w[i][r] + w[r + 1][j]
                if candidate < best:
                    best = candidate
            w[i][j] = best

    pairs: list[tuple[int, int]] = []

    def eq(a: float, b: float) -> bool:
        if a == inf or b == inf:
            return False
        return abs(a - b) < _EPS

    def trace_w(i: int, j: int) -> None:
        if j - i < MIN_LOOP + 1:
            return
        if eq(w[i][j], w[i + 1][j]):
            return trace_w(i + 1, j)
        if eq(w[i][j], w[i][j - 1]):
            return trace_w(i, j - 1)
        if v[i][j] != inf and eq(w[i][j], v[i][j] + model.terminal_penalty(i, j)):
            pairs.append((i, j))
            return trace_v(i, j)
        for k in range(i, j):
            if eq(w[i][j], w[i][k] + w[k + 1][j]):
                trace_w(i, k)
                trace_w(k + 1, j)
                return

    def trace_v(i: int, j: int) -> None:
        value = v[i][j]
        if value == inf or eq(value, model.hairpin(i, j)):
            return

        k_max = min(i + MAX_INTERNAL + 1, j - MIN_LOOP - 2)
        for k in range(i + 1, k_max + 1):
            left = k - i - 1
            l_min = max(k + MIN_LOOP + 1, j - 1 - (MAX_INTERNAL - left))
            for l in range(j - 1, l_min - 1, -1):
                if v[k][l] == inf:
                    continue
                if eq(value, model.interior(i, j, k, l) + v[k][l]):
                    pairs.append((k, l))
                    return trace_v(k, l)

        closing = ML_CLOSE + ML_BRANCH + model.terminal_penalty(i, j)
        for m in range(i + 2, j - 1):
            if wm[i + 1][m] == inf or wm[m + 1][j - 1] == inf:
                continue
            if eq(value, wm[i + 1][m] + wm[m + 1][j - 1] + closing):
                trace_wm(i + 1, m)
                trace_wm(m + 1, j - 1)
                return

    def trace_wm(i: int, j: int) -> None:
        value = wm[i][j]
        if value == inf:
            return
        if v[i][j] != inf and eq(value, v[i][j] + ML_BRANCH + model.terminal_penalty(i, j)):
            pairs.append((i, j))
            return trace_v(i, j)
        if wm[i + 1][j] != inf and eq(value, wm[i + 1][j] + ML_UNPAIRED):
            return trace_wm(i + 1, j)
        if wm[i][j - 1] != inf and eq(value, wm[i][j - 1] + ML_UNPAIRED):
            return trace_wm(i, j - 1)
        for q in range(i, j):
            if wm[i][q] == inf or wm[q + 1][j] == inf:
                continue
            if eq(value, wm[i][q] + wm[q + 1][j]):
                trace_wm(i, q)
                trace_wm(q + 1, j)
                return

    trace_w(0, n - 1)
    pairs.sort()
    return MfeResult(model.seq, w[0][n - 1], dot_bracket(n, pairs), pairs)
