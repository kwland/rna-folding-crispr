"""Feature blocks for the predictive models.

The study's central question is *incremental*: does folding tell you anything a
sequence-only model did not already know? Answering that requires a sequence
baseline that is genuinely competitive, otherwise "structure helps" can just
mean "the baseline was weak". The baseline here is the standard one for guide
activity - position-specific nucleotides and dinucleotides plus G/C content -
which is the feature family behind published rule sets.

Feature scaling uses fixed constants (percentages divided by 100, energies by
10), never statistics estimated from the data. That keeps the design matrix
completely independent of the labels and of which fold a guide lands in, so
there is no route by which held-out information can leak into training.
"""

from __future__ import annotations

from datasets import Guide


__all__ = ["FeatureSpec", "BLOCKS", "build_design_rows", "structure_feature_columns"]

BASES = "ACGT"
BASE_INDEX = {base: index for index, base in enumerate(BASES)}
DINUCS = [a + b for a in BASES for b in BASES]
DINUC_INDEX = {pair: index for index, pair in enumerate(DINUCS)}

SPACER_LENGTH = 20

#: Scalar folding features, with the fixed divisor used to bring them to O(1).
STRUCTURE_FEATURES: list[tuple[str, float]] = [
    ("seed_nussinov", 1.0),
    ("seed_mfe_spacer", 1.0),
    ("seed_ensemble_spacer", 1.0),
    ("seed_mfe_full", 1.0),
    ("seed_ensemble_full", 1.0),
    ("mean_unpaired_spacer", 1.0),
    ("mean_unpaired_full", 1.0),
    ("mfe_energy_spacer", 10.0),
    ("mfe_energy_full", 10.0),
    ("ensemble_energy_spacer", 10.0),
    ("ensemble_energy_full", 10.0),
]

#: Subsets used to separate "does structure help at all" from "which structure view".
STRUCTURE_SUBSETS: dict[str, list[str]] = {
    "structure_all": [name for name, _ in STRUCTURE_FEATURES],
    "structure_nussinov": ["seed_nussinov"],
    "structure_mfe_spacer": ["seed_mfe_spacer", "mfe_energy_spacer"],
    "structure_ensemble_spacer": [
        "seed_ensemble_spacer", "mean_unpaired_spacer", "ensemble_energy_spacer",
    ],
    "structure_ensemble_full": [
        "seed_ensemble_full", "mean_unpaired_full", "ensemble_energy_full",
    ],
    "structure_profile_full": [f"unpaired_full_{i:02d}" for i in range(1, SPACER_LENGTH + 1)],
    "structure_profile_spacer": [
        f"unpaired_spacer_{i:02d}" for i in range(1, SPACER_LENGTH + 1)
    ],
}

BLOCKS = ["position", "dinucleotide", "gc"]


class FeatureSpec:
    """Which feature blocks to include, and the resulting column layout."""

    def __init__(
        self,
        position: bool = True,
        dinucleotide: bool = True,
        gc: bool = True,
        structure: list[str] | None = None,
    ) -> None:
        self.position = position
        self.dinucleotide = dinucleotide
        self.gc = gc
        self.structure = list(structure or [])

        # Index 0 is the intercept and is never penalised by the ridge.
        self.intercept_index = 0
        offset = 1
        self.position_offset = offset
        if position:
            offset += SPACER_LENGTH * len(BASES)
        self.dinuc_offset = offset
        if dinucleotide:
            offset += (SPACER_LENGTH - 1) * len(DINUCS)
        self.gc_offset = offset
        if gc:
            offset += 1
        self.structure_offset = offset
        offset += len(self.structure)
        self.n_features = offset

    def describe(self) -> str:
        parts = []
        if self.position:
            parts.append("position")
        if self.dinucleotide:
            parts.append("dinucleotide")
        if self.gc:
            parts.append("gc")
        if self.structure:
            parts.append(f"structure[{len(self.structure)}]")
        return "+".join(parts) or "intercept only"


def structure_feature_columns() -> list[str]:
    return [name for name, _ in STRUCTURE_FEATURES]


def _scale_for(name: str) -> float:
    for feature, divisor in STRUCTURE_FEATURES:
        if feature == name:
            return divisor
    return 1.0  # per-position unpaired probabilities are already in [0, 1]


def build_design_rows(
    guides: list[Guide],
    features: list[dict[str, str]],
    spec: FeatureSpec,
) -> list[list[tuple[int, float]]]:
    """Build sparse rows of (column index, value) for the given guides.

    ``features`` holds one cached row per guide, in the same order, as written
    by ``compute_features.py``.
    """
    if len(guides) != len(features):
        raise ValueError("guides and cached features are not aligned")

    rows: list[list[tuple[int, float]]] = []
    for guide, cached in zip(guides, features):
        spacer = guide.spacer.upper().replace("U", "T")
        if len(spacer) != SPACER_LENGTH:
            raise ValueError(f"{guide.guide_id}: expected {SPACER_LENGTH} nt spacer")
        row: list[tuple[int, float]] = [(spec.intercept_index, 1.0)]

        if spec.position:
            for position, base in enumerate(spacer):
                column = spec.position_offset + position * len(BASES) + BASE_INDEX[base]
                row.append((column, 1.0))

        if spec.dinucleotide:
            for position in range(SPACER_LENGTH - 1):
                pair = spacer[position : position + 2]
                column = (
                    spec.dinuc_offset
                    + position * len(DINUCS)
                    + DINUC_INDEX[pair]
                )
                row.append((column, 1.0))

        if spec.gc:
            gc = (spacer.count("G") + spacer.count("C")) / SPACER_LENGTH
            row.append((spec.gc_offset, gc))

        for offset, name in enumerate(spec.structure):
            value = float(cached[name]) / _scale_for(name)
            if value:
                row.append((spec.structure_offset + offset, value))

        rows.append(row)
    return rows
