"""Nearest-neighbour free-energy model shared by the MFE and partition-function folders.

This is the Python counterpart of the energy model in ``docs/js/rna-fold.js``.
Keeping one set of parameters in both places means the browser demonstration and
the batch analysis describe the same physical model, so a disagreement between
them is a bug rather than a difference of opinion.

PARAMETER PROVENANCE. Read this before quoting numbers in a report.

* The ten Watson-Crick nearest-neighbour stacking energies are the published
  values (Xia et al. 1998; Turner & Mathews 2010), kcal/mol at 37 C. Exact.
* G-U wobble stacking, the loop initiation tables, and the multiloop model are
  SIMPLIFIED approximations of the Turner rules. They preserve the qualitative
  ordering but are not the exact published tables. Every one of them carries an
  APPROXIMATION comment below.

Consequence: predicted structures are usually sensible, but absolute kcal/mol
values will not match ViennaRNA exactly. ViennaRNA remains the reference for any
quantitative energy claim. What this module is for is a *self-consistent* model
in which the only thing that changes between two experiments is the thing being
tested.
"""

from __future__ import annotations

import math
from typing import Iterable

from nussinov import CANONICAL_PAIRS, WOBBLE_PAIRS, normalize_rna


__all__ = [
    "RT37",
    "MIN_LOOP",
    "MAX_INTERNAL",
    "ML_CLOSE",
    "ML_BRANCH",
    "ML_UNPAIRED",
    "EnergyModel",
    "BoltzmannTables",
    "normalize_rna",
    "structure_energy",
    "enumerate_structures",
    "pairs_from_dot_bracket",
]


# Gas constant times temperature at 37 C, kcal/mol. Every Boltzmann weight in
# the partition function is exp(-E / RT37).
RT37 = 0.6163
# Jacobson-Stockmayer coefficient for extrapolating loop initiation past the table.
LOOP_SCALE = 1.75 * RT37

MIN_LOOP = 3       # a hairpin must enclose at least this many unpaired bases
MAX_INTERNAL = 30  # Turner convention: ignore internal loops larger than this

# Linear multiloop model: a + b * branches + c * unpaired bases.
#
# b is NEGATIVE, so an extra branch is rewarded. That looks wrong but it is the
# Turner 2004 convention and it is load-bearing: with a positive b the model
# nests tRNA's arms into a chain instead of opening the four-way junction.
# APPROXIMATION: the real Turner multiloop term is not strictly linear.
ML_CLOSE = 3.4
ML_BRANCH = -0.9
ML_UNPAIRED = 0.0

# Terminal A-U / G-U penalty at the end of a helix (Turner: +0.45).
TERMINAL_AU = 0.45
# APPROXIMATION: Turner uses a distinct A-U/G-U closure cost inside internal loops.
INTERNAL_TERMINAL_AU = 0.7

# Loop initiation tables (Turner-style, kcal/mol); larger sizes are extrapolated.
HAIRPIN_INIT = {3: 5.4, 4: 5.6, 5: 5.7, 6: 5.4, 7: 6.0, 8: 5.5, 9: 6.4}
BULGE_INIT = {1: 3.8, 2: 2.8, 3: 3.2, 4: 3.6, 5: 4.0, 6: 4.4}
# APPROXIMATION: sizes 2-3 (1x1 and 1x2 loops) use special tables in Turner; flattened here.
INTERNAL_INIT = {2: 1.5, 3: 1.6, 4: 1.7, 5: 1.8, 6: 2.0}

# "WX/ZY" is the duplex 5'-WX-3' over 3'-ZY-5', i.e. outer pair (W,Z) carrying
# inner pair (X,Y) stacked directly on it.
_WC_STACKS = [
    ("AA/UU", -0.93),
    ("AU/UA", -1.10),
    ("UA/AU", -1.33),
    ("CU/GA", -2.08),
    ("CA/GU", -2.11),
    ("GU/CA", -2.24),
    ("GA/CU", -2.35),
    ("CG/GC", -2.36),
    ("GG/CC", -3.26),
    ("GC/CG", -3.42),
]

# APPROXIMATION. Real Turner wobble parameters are context dependent (tandem G-U
# has its own tables). Collapsed to two flat values that keep the ordering:
# a wobble stack is weaker than any Watson-Crick stack, tandem wobbles weaker still.
WOBBLE_SINGLE_STACK = -1.3
WOBBLE_TANDEM_STACK = -0.5


def _build_stack_table() -> dict[str, float]:
    table: dict[str, float] = {}
    for key, value in _WC_STACKS:
        outer, inner = key.split("/")
        w, x = outer[0], outer[1]
        z, y = inner[0], inner[1]
        # Key layout: outer pair (w, z) then inner pair (x, y).
        table[w + x + z + y] = value
        # Reading the same duplex from the other strand: WX/ZY == YZ/XW.
        table[y + z + x + w] = value
    return table


_STACK = _build_stack_table()
_CANONICAL = {a + b for a, b in CANONICAL_PAIRS}
_WOBBLE = {a + b for a, b in WOBBLE_PAIRS}


def _extrapolate(table: dict[int, float], size: int, anchor: int) -> float:
    if size in table:
        return table[size]
    return table[anchor] + LOOP_SCALE * math.log(size / anchor)


class EnergyModel:
    """Loop energies for one sequence, in kcal/mol at 37 C.

    The class exists so the folding routines can ask ``model.hairpin(i, j)``
    without re-deriving sequence lookups, and so a test can swap in a different
    parameter set without touching the recursions.
    """

    def __init__(self, sequence: str, allow_wobble: bool = True) -> None:
        self.seq = normalize_rna(sequence)
        self.n = len(self.seq)
        self.allow_wobble = allow_wobble

    # -------------------------------------------------------------- pairing

    def can_pair(self, i: int, j: int) -> bool:
        """True when bases i and j may form a pair with a legal enclosed loop."""
        if j - i - 1 < MIN_LOOP:
            return False
        key = self.seq[i] + self.seq[j]
        if key in _CANONICAL:
            return True
        return self.allow_wobble and key in _WOBBLE

    def _is_wobble(self, i: int, j: int) -> bool:
        return (self.seq[i] + self.seq[j]) in _WOBBLE

    # ------------------------------------------------------------- energies

    def terminal_penalty(self, i: int, j: int) -> float:
        """Helix-end penalty: free for G-C, charged for A-U and G-U."""
        a, b = self.seq[i], self.seq[j]
        if (a + b) in _CANONICAL and (a == "G" or a == "C"):
            return 0.0
        return TERMINAL_AU

    def _internal_terminal_penalty(self, i: int, j: int) -> float:
        a, b = self.seq[i], self.seq[j]
        if (a + b) in _CANONICAL and (a == "G" or a == "C"):
            return 0.0
        return INTERNAL_TERMINAL_AU

    def stack(self, i: int, j: int, k: int, l: int) -> float:
        """Stacking energy of inner pair (k, l) on outer pair (i, j)."""
        seq = self.seq
        wc = _STACK.get(seq[i] + seq[k] + seq[j] + seq[l])
        if wc is not None:
            return wc
        wobbles = (1 if self._is_wobble(i, j) else 0) + (1 if self._is_wobble(k, l) else 0)
        if wobbles >= 2:
            return WOBBLE_TANDEM_STACK
        if wobbles == 1:
            return WOBBLE_SINGLE_STACK
        return 0.0

    def hairpin(self, i: int, j: int) -> float:
        """Energy of the hairpin loop closed by (i, j)."""
        size = j - i - 1
        if size < MIN_LOOP:
            return math.inf
        return _extrapolate(HAIRPIN_INIT, size, 9) + self.terminal_penalty(i, j)

    def interior(self, i: int, j: int, k: int, l: int) -> float:
        """Energy of the loop closed by (i, j) with a single inner pair (k, l).

        Covers the three cases that share this shape: a stack (no unpaired
        bases), a bulge (unpaired on one side), and an internal loop (both).
        """
        left = k - i - 1
        right = j - l - 1
        total = left + right

        if total == 0:
            return self.stack(i, j, k, l)

        if left == 0 or right == 0:
            energy = _extrapolate(BULGE_INIT, total, 6)
            if total == 1:
                # A single bulged base does not break the helix, so the stack survives.
                energy += self.stack(i, j, k, l)
            else:
                energy += self.terminal_penalty(i, j) + self.terminal_penalty(k, l)
            return energy

        energy = _extrapolate(INTERNAL_INIT, total, 6) if total >= 2 else math.inf
        # APPROXIMATION: Turner's asymmetry term uses a fitted coefficient and cap.
        energy += min(3.0, 0.6 * abs(left - right))
        energy += self._internal_terminal_penalty(i, j) + self._internal_terminal_penalty(k, l)
        return energy


class BoltzmannTables:
    """Pre-exponentiated loop weights for one sequence at one temperature.

    The partition function evaluates the same handful of loop energies millions
    of times, and ``math.log``/``math.exp`` in that inner loop dominated the
    runtime. Every weight here is a pure function of a small key, so it can be
    tabulated once and then read as an array lookup.

    The factorisations below are exact, not approximations. They just exploit
    that the loop terms separate: a bulge of size >= 2 costs an initiation that
    depends only on its size, plus one helix-end penalty for each of the two
    pairs; an internal loop likewise splits into a size/asymmetry part and two
    per-pair parts.
    """

    def __init__(self, model: EnergyModel, rt: float = RT37) -> None:
        self.model = model
        self.rt = rt
        n = model.n
        seq = model.seq

        def w(energy: float) -> float:
            return 0.0 if energy == math.inf else math.exp(-energy / rt)

        # Per-pair helix-end weights.
        self.term = [[1.0] * n for _ in range(n)]
        self.int_term = [[1.0] * n for _ in range(n)]
        # Loop weights keyed by the closing pair.
        self.hairpin = [[0.0] * n for _ in range(n)]
        self.stack = [[0.0] * n for _ in range(n)]
        self.bulge1_left = [[0.0] * n for _ in range(n)]
        self.bulge1_right = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(i + MIN_LOOP + 1, n):
                if not model.can_pair(i, j):
                    continue
                self.term[i][j] = w(model.terminal_penalty(i, j))
                self.int_term[i][j] = w(model._internal_terminal_penalty(i, j))
                self.hairpin[i][j] = w(model.hairpin(i, j))
                if model.can_pair(i + 1, j - 1):
                    self.stack[i][j] = w(model.stack(i, j, i + 1, j - 1))
                # A single bulged base on the left: inner pair is (i+2, j-1).
                if i + 2 < j - 1 and model.can_pair(i + 2, j - 1):
                    self.bulge1_left[i][j] = w(
                        _extrapolate(BULGE_INIT, 1, 6) + model.stack(i, j, i + 2, j - 1)
                    )
                # A single bulged base on the right: inner pair is (i+1, j-2).
                if i + 1 < j - 2 and model.can_pair(i + 1, j - 2):
                    self.bulge1_right[i][j] = w(
                        _extrapolate(BULGE_INIT, 1, 6) + model.stack(i, j, i + 1, j - 2)
                    )

        # Bulges of size >= 2: initiation only; the two end penalties are per-pair.
        self.bulge = [0.0] * (MAX_INTERNAL + 2)
        for size in range(2, MAX_INTERNAL + 2):
            self.bulge[size] = w(_extrapolate(BULGE_INIT, size, 6))

        # Internal loops: initiation plus asymmetry, indexed by the two side sizes.
        limit = MAX_INTERNAL + 1
        self.internal = [[0.0] * limit for _ in range(limit)]
        for left in range(1, limit):
            for right in range(1, limit):
                total = left + right
                if total > MAX_INTERNAL or total < 2:
                    continue
                energy = _extrapolate(INTERNAL_INIT, total, 6)
                energy += min(3.0, 0.6 * abs(left - right))
                self.internal[left][right] = w(energy)

    def interior(self, i: int, j: int, k: int, l: int) -> float:
        """Boltzmann weight of the loop closed by (i, j) around inner pair (k, l).

        Kept for readability and testing. The folding routines inline this same
        arithmetic in their innermost loops, and ``test_mccaskill.py`` checks the
        two agree with ``EnergyModel.interior``.
        """
        left = k - i - 1
        right = j - l - 1
        total = left + right
        if total == 0:
            return self.stack[i][j]
        if left == 0 or right == 0:
            if total == 1:
                return self.bulge1_right[i][j] if right == 1 else self.bulge1_left[i][j]
            return self.bulge[total] * self.term[i][j] * self.term[k][l]
        return self.internal[left][right] * self.int_term[i][j] * self.int_term[k][l]


# ------------------------------------------------------------------ helpers


def pairs_from_dot_bracket(structure: str) -> list[tuple[int, int]]:
    """Convert dot-bracket notation to a sorted list of zero-based pairs."""
    stack: list[int] = []
    pairs: list[tuple[int, int]] = []
    for index, char in enumerate(structure):
        if char == "(":
            stack.append(index)
        elif char == ")":
            if not stack:
                raise ValueError(f"Unbalanced ')' at position {index + 1}")
            pairs.append((stack.pop(), index))
        elif char != ".":
            raise ValueError(f"Unexpected dot-bracket character {char!r}")
    if stack:
        raise ValueError(f"Unbalanced '(' at position {stack[-1] + 1}")
    pairs.sort()
    return pairs


def structure_energy(
    sequence: str,
    pairs: Iterable[tuple[int, int]] | str,
    allow_wobble: bool = True,
) -> float:
    """Free energy of one explicit structure, evaluated loop by loop.

    This does not use dynamic programming. It decomposes the structure into its
    loops and adds them up, which makes it an independent check on the folding
    recursions: the DP and this function must agree on any structure the DP
    returns. ``test_mccaskill.py`` relies on that.
    """
    model = EnergyModel(sequence, allow_wobble=allow_wobble)
    n = model.n
    if isinstance(pairs, str):
        pair_list = pairs_from_dot_bracket(pairs)
    else:
        pair_list = sorted(pairs)

    partner = [-1] * n
    for i, j in pair_list:
        if not (0 <= i < j < n):
            raise ValueError(f"Pair ({i}, {j}) is outside the sequence")
        if partner[i] != -1 or partner[j] != -1:
            raise ValueError(f"Base in pair ({i}, {j}) is already paired")
        partner[i] = j
        partner[j] = i

    def children_of(i: int, j: int) -> list[tuple[int, int]]:
        """Pairs directly enclosed by (i, j), skipping anything deeper."""
        out: list[tuple[int, int]] = []
        k = i + 1
        while k < j:
            if partner[k] > k:
                out.append((k, partner[k]))
                k = partner[k] + 1
            else:
                k += 1
        return out

    def loop_energy(i: int, j: int) -> float:
        kids = children_of(i, j)
        if not kids:
            return model.hairpin(i, j)
        if len(kids) == 1:
            k, l = kids[0]
            return model.interior(i, j, k, l) + loop_energy(k, l)
        unpaired = (j - i - 1) - sum(l - k + 1 for k, l in kids)
        total = ML_CLOSE + ML_BRANCH + model.terminal_penalty(i, j)
        total += ML_UNPAIRED * unpaired
        for k, l in kids:
            total += ML_BRANCH + model.terminal_penalty(k, l) + loop_energy(k, l)
        return total

    # Exterior loop: unpaired bases are free, each top-level helix pays its end penalty.
    energy = 0.0
    index = 0
    while index < n:
        j = partner[index]
        if j > index:
            energy += model.terminal_penalty(index, j) + loop_energy(index, j)
            index = j + 1
        else:
            index += 1
    return energy


def enumerate_structures(
    sequence: str,
    allow_wobble: bool = True,
) -> list[list[tuple[int, int]]]:
    """Every pseudoknot-free structure of a short sequence, for exact checks.

    The count grows quickly, so this is a validation tool for sequences of
    roughly 16 nt or less, not an analysis routine.
    """
    model = EnergyModel(sequence, allow_wobble=allow_wobble)
    n = model.n
    cache: dict[tuple[int, int], list[list[tuple[int, int]]]] = {}

    def solve(i: int, j: int) -> list[list[tuple[int, int]]]:
        if j - i < MIN_LOOP + 1:
            return [[]]
        key = (i, j)
        if key in cache:
            return cache[key]
        # Unambiguous split: either j is unpaired, or j pairs with exactly one k.
        results = [list(rest) for rest in solve(i, j - 1)]
        for k in range(i, j - MIN_LOOP):
            if not model.can_pair(k, j):
                continue
            for left in solve(i, k - 1) if k > i else [[]]:
                for inner in solve(k + 1, j - 1):
                    results.append([*left, (k, j), *inner])
        cache[key] = results
        return results

    if n == 0:
        return [[]]
    return [sorted(structure) for structure in solve(0, n - 1)]
