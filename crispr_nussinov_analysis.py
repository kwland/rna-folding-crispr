"""Analyze CRISPR guide spacers with three folding models side by side.

The starting idea: a guide that folds back on itself, especially in the
PAM-proximal seed region, may have less spacer available to bind its DNA target.

Three different answers to "is the seed accessible?" are reported, because they
disagree and the disagreement is the interesting part:

* **Nussinov** - one structure with the most base pairs. Every seed position is
  scored 0 or 1. Transparent, but not a physical model.
* **Zuker MFE** - one structure, the lowest free energy under the nearest-
  neighbour model. Still 0 or 1 per position.
* **McCaskill ensemble** - the probability each position is unpaired, averaged
  over the whole Boltzmann ensemble rather than a single structure. This is what
  RNAplfold-style accessibility actually means.

``--context sgrna`` additionally folds each spacer joined to the 76 nt sgRNA
scaffold, because a spacer in a real guide can pair with the scaffold and a
spacer folded alone cannot. That is slower (the molecule is 96 nt rather than
20) but it is the biologically honest question.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from guide_features import SGRNA_SCAFFOLD, gc_fraction
from mccaskill import partition_fold
from nussinov import nussinov_fold
from zuker import zuker_fold


DEFAULT_SEED_LENGTH = 8
FEATURE_COLUMNS = [
    "spacer_rna",
    "length",
    "gc_percent",
    "nussinov_pairs",
    "self_pair_fraction",
    "dot_bracket",
    "seed_region_1_based",
    "seed_paired_bases",
    "seed_accessibility",
    "mfe_dot_bracket",
    "mfe_energy_kcal",
    "seed_accessibility_mfe",
    "ensemble_free_energy_kcal",
    "seed_accessibility_ensemble",
    "mean_unpaired_ensemble",
    "design_warning",
    "analysis_status",
]

#: Extra columns produced by --context sgrna.
SCAFFOLD_COLUMNS = [
    "sgrna_length",
    "sgrna_mfe_energy_kcal",
    "seed_accessibility_mfe_sgrna",
    "sgrna_ensemble_free_energy_kcal",
    "seed_accessibility_ensemble_sgrna",
    "mean_unpaired_ensemble_sgrna",
]


def paired_positions(pairs: list[tuple[int, int]]) -> set[int]:
    positions: set[int] = set()
    for i, j in pairs:
        positions.add(i)
        positions.add(j)
    return positions


def design_warnings(
    sequence: str,
    gc_percent: float,
    seed_accessibility: float,
    seed_ensemble: float | None = None,
) -> str:
    warnings = []
    dna = sequence.upper().replace("U", "T")
    if len(dna) != 20:
        warnings.append("not_20_nt_spacer")
    if gc_percent < 30:
        warnings.append("low_gc")
    elif gc_percent > 75:
        warnings.append("high_gc")
    if "TTTT" in dna:
        warnings.append("poly_t_u6_termination_risk")
    if seed_accessibility < 0.5:
        warnings.append("seed_mostly_paired")
    if seed_ensemble is not None and seed_ensemble < 0.5:
        # The ensemble version can fire when the single structure looks fine,
        # which is the whole reason for computing it.
        warnings.append("seed_low_ensemble_accessibility")
    return ";".join(warnings) or "none"


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def score_guide(
    spacer: str,
    seed_length: int = DEFAULT_SEED_LENGTH,
    min_loop_length: int = 3,
    context: str = "spacer",
    scaffold: str = SGRNA_SCAFFOLD,
) -> dict[str, str]:
    """Fold one spacer and summarise it under all three models.

    ``context="sgrna"`` also folds the spacer joined to the scaffold and reports
    the seed accessibility the spacer has inside the real molecule.
    """
    if seed_length <= 0:
        raise ValueError("seed_length must be a positive integer")
    if context not in ("spacer", "sgrna"):
        raise ValueError("context must be 'spacer' or 'sgrna'")

    result = nussinov_fold(spacer, min_loop_length=min_loop_length)
    paired = paired_positions(result.pairs)
    n = len(result.sequence)
    seed_start = max(0, n - seed_length)
    seed_positions = set(range(seed_start, n))
    seed_paired = len(seed_positions & paired)
    seed_unpaired = len(seed_positions) - seed_paired
    seed_accessibility = seed_unpaired / len(seed_positions) if seed_positions else 0.0
    gc_percent = gc_fraction(result.sequence) * 100
    self_pair_fraction = (2 * result.score / n) if n else 0.0

    mfe = zuker_fold(result.sequence)
    mfe_flags = mfe.unpaired_flags()
    ensemble = partition_fold(result.sequence)
    seed_ensemble = ensemble.mean_unpaired(seed_start, n)

    features = {
        "spacer_rna": result.sequence,
        "length": str(n),
        "gc_percent": f"{gc_percent:.1f}",
        "nussinov_pairs": str(result.score),
        "self_pair_fraction": f"{self_pair_fraction:.3f}",
        "dot_bracket": result.structure,
        "seed_region_1_based": f"{seed_start + 1}-{n}",
        "seed_paired_bases": str(seed_paired),
        "seed_accessibility": f"{seed_accessibility:.3f}",
        "mfe_dot_bracket": mfe.structure,
        "mfe_energy_kcal": f"{mfe.energy:.2f}",
        "seed_accessibility_mfe": f"{_mean(mfe_flags[seed_start:n]):.3f}",
        "ensemble_free_energy_kcal": f"{ensemble.ensemble_free_energy:.2f}",
        "seed_accessibility_ensemble": f"{seed_ensemble:.3f}",
        "mean_unpaired_ensemble": f"{ensemble.mean_unpaired(0, n):.3f}",
        "design_warning": design_warnings(
            result.sequence, gc_percent, seed_accessibility, seed_ensemble
        ),
        "analysis_status": "ok",
    }

    if context == "sgrna":
        full = result.sequence + scaffold
        full_mfe = zuker_fold(full)
        full_ensemble = partition_fold(full)
        features.update({
            "sgrna_length": str(len(full)),
            "sgrna_mfe_energy_kcal": f"{full_mfe.energy:.2f}",
            "seed_accessibility_mfe_sgrna":
                f"{_mean(full_mfe.unpaired_flags()[seed_start:n]):.3f}",
            "sgrna_ensemble_free_energy_kcal": f"{full_ensemble.ensemble_free_energy:.2f}",
            "seed_accessibility_ensemble_sgrna":
                f"{full_ensemble.mean_unpaired(seed_start, n):.3f}",
            "mean_unpaired_ensemble_sgrna": f"{full_ensemble.mean_unpaired(0, n):.3f}",
        })
    return features


def placeholder_features(status: str, context: str = "spacer") -> dict[str, str]:
    columns = list(FEATURE_COLUMNS)
    if context == "sgrna":
        columns += SCAFFOLD_COLUMNS
    return {column: (status if column == "analysis_status" else "") for column in columns}


def analyze_guides(
    input_csv: Path,
    output_csv: Path,
    seed_length: int,
    context: str = "spacer",
) -> None:
    if seed_length <= 0:
        raise ValueError("seed_length must be a positive integer")

    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        input_columns = reader.fieldnames or []
        if "spacer_dna" not in input_columns and "spacer_rna" not in input_columns:
            raise ValueError("Input CSV needs a spacer_dna or spacer_rna column.")

        rows = []
        for row in reader:
            spacer = row.get("spacer_rna") or row.get("spacer_dna") or ""
            if not spacer or "REPLACE_WITH" in spacer.upper():
                features = placeholder_features("needs_sequence", context)
            else:
                features = score_guide(spacer, seed_length=seed_length, context=context)
            rows.append({**row, **features})

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    default_columns = [*input_columns, *FEATURE_COLUMNS]
    if context == "sgrna":
        default_columns += SCAFFOLD_COLUMNS
    fieldnames = list(rows[0].keys()) if rows else default_columns
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fold CRISPR guide spacers with Nussinov, Zuker, and McCaskill."
    )
    parser.add_argument("--input", default="data/crispr_guide_examples.csv")
    parser.add_argument("--output", default="analysis_outputs/crispr_guide_nussinov_features.csv")
    parser.add_argument("--seed-length", type=int, default=DEFAULT_SEED_LENGTH)
    parser.add_argument(
        "--context", choices=("spacer", "sgrna"), default="spacer",
        help="Fold the spacer alone (fast) or joined to the 76 nt sgRNA scaffold (slower).",
    )
    args = parser.parse_args()

    analyze_guides(Path(args.input), Path(args.output), args.seed_length, args.context)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
