"""Check every headline result against ViennaRNA.

The study's folding code uses published Watson-Crick stacking values but
simplified wobble, loop, and multiloop terms. This module measures what that
costs, by recomputing the same quantities with ViennaRNA under standard Turner
parameters and asking four questions:

1. **Agreement.** How closely does the custom ensemble accessibility track
   ViennaRNA's, guide by guide and position by position?
2. **Do the conclusions survive?** Recompute the accessibility-to-activity
   correlations and the incremental model value using ViennaRNA features. A
   result that only holds under the simplified parameters is not a result.
3. **Temperature sensitivity.** Folding is temperature dependent, and screens
   are not all run at 37 C. If the conclusion flips at 25 C or 42 C, that needs
   saying.
4. **Folding context.** Global folding assumes the whole molecule equilibrates
   together. RNAplfold-style local folding, which forbids long-range pairs, is
   closer to a molecule folding as it is transcribed.

Run ``vienna_reference.py`` first to produce the inputs. This module needs only
the standard library: it reads the CSVs ViennaRNA wrote.

    python validate_vienna.py
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from model_features import FeatureSpec, build_design_rows
from run_study import ALPHA_GRID, LOGO_MAX_GENES, load_cached
from stats import (
    SparseDesign,
    mean_within_group_spearman,
    nested_cv_predictions,
    paired_spearman_delta_ci,
    partial_spearman,
    pearson,
    spearman,
    within_group_spearman_ci,
)


OUTPUTS = Path("analysis_outputs")

#: (custom column, ViennaRNA column, label) for the measures the study reports.
COMPARISONS = [
    ("seed_mfe_spacer", "vienna_seed_mfe_spacer", "MFE seed accessibility (spacer alone)"),
    ("seed_ensemble_spacer", "vienna_seed_ensemble_spacer", "Ensemble seed accessibility (spacer alone)"),
    ("seed_mfe_full", "vienna_seed_mfe_full", "MFE seed accessibility (spacer + scaffold)"),
    ("seed_ensemble_full", "vienna_seed_ensemble_full", "Ensemble seed accessibility (spacer + scaffold)"),
    ("mean_unpaired_spacer", "vienna_mean_unpaired_spacer", "Mean unpaired probability (spacer alone)"),
    ("mean_unpaired_full", "vienna_mean_unpaired_full", "Mean unpaired probability (in sgRNA)"),
    ("mfe_energy_full", "vienna_mfe_energy_full", "MFE energy (spacer + scaffold)"),
    ("ensemble_energy_full", "vienna_ensemble_energy_full", "Ensemble free energy (spacer + scaffold)"),
]

#: ViennaRNA structure features offered to the incremental model.
VIENNA_MODEL_FEATURES = [
    "vienna_seed_ensemble_full",
    "vienna_mean_unpaired_full",
    "vienna_ensemble_energy_full",
    "vienna_seed_mfe_full",
    "vienna_mfe_energy_full",
]


def vienna_path(dataset: str, temperature: float = 37.0, max_span: int | None = None) -> Path:
    stem = f"vienna_{dataset}"
    if abs(temperature - 37.0) > 1e-9:
        stem += f"_t{temperature:g}"
    if max_span:
        stem += f"_span{max_span}"
    return OUTPUTS / f"{stem}.csv"


def load_vienna(path: Path, guide_ids: list[str]) -> list[dict[str, str]]:
    """Load a ViennaRNA feature file and align it to the study's guide order."""
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. Run: python vienna_reference.py "
            "(with an interpreter that has ViennaRNA installed)"
        )
    with path.open(newline="", encoding="utf-8") as handle:
        rows = {row["guide_id"]: row for row in csv.DictReader(handle)}
    missing = [g for g in guide_ids if g not in rows]
    if missing:
        raise ValueError(f"{path} is missing {len(missing)} guides, e.g. {missing[:3]}")
    return [rows[g] for g in guide_ids]


def _column(rows: list[dict[str, str]], name: str) -> list[float]:
    return [float(row[name]) for row in rows]


# ------------------------------------------------------------- 1. agreement


def agreement(custom: list[dict], vienna: list[dict]) -> list[dict]:
    """How closely the two implementations track each other, measure by measure."""
    results = []
    for mine, theirs, label in COMPARISONS:
        a = _column(custom, mine)
        b = _column(vienna, theirs)
        differences = [x - y for x, y in zip(a, b)]
        n = len(differences)
        results.append({
            "feature": mine,
            "vienna_feature": theirs,
            "label": label,
            "pearson": pearson(a, b),
            "spearman": spearman(a, b),
            "mean_difference": sum(differences) / n,
            "mean_absolute_difference": sum(abs(d) for d in differences) / n,
            "max_absolute_difference": max(abs(d) for d in differences),
        })
    return results


def position_agreement(custom: list[dict], vienna: list[dict]) -> list[dict]:
    """Per-position agreement of the unpaired-probability profiles."""
    out = []
    for context, mine_prefix, their_prefix in (
        ("spacer", "unpaired_spacer", "vienna_unpaired_spacer"),
        ("full", "unpaired_full", "vienna_unpaired_full"),
    ):
        for position in range(1, 21):
            a = _column(custom, f"{mine_prefix}_{position:02d}")
            b = _column(vienna, f"{their_prefix}_{position:02d}")
            out.append({
                "context": context,
                "position": position,
                "pearson": pearson(a, b),
                "spearman": spearman(a, b),
                "mean_absolute_difference": sum(abs(x - y) for x, y in zip(a, b)) / len(a),
                "mean_custom": sum(a) / len(a),
                "mean_vienna": sum(b) / len(b),
            })
    return out


# ------------------------------------------- 2. do the conclusions survive?


def correlations_under_vienna(guides, custom, vienna, n_boot: int, seed: int) -> list[dict]:
    """Accessibility-to-activity correlations computed from ViennaRNA features."""
    activity = [g.percentile for g in guides]
    genes = [g.group for g in guides]
    gc = _column(custom, "gc_percent")

    results = []
    for mine, theirs, label in COMPARISONS:
        values = _column(vienna, theirs)
        within, low, high, n_groups = within_group_spearman_ci(
            values, activity, genes, n_boot=n_boot, seed=seed
        )
        custom_values = _column(custom, mine)
        custom_within, _, _, _ = within_group_spearman_ci(
            custom_values, activity, genes, n_boot=n_boot, seed=seed
        )
        results.append({
            "feature": theirs,
            "label": label,
            "vienna_within_gene": within,
            "vienna_ci_low": low,
            "vienna_ci_high": high,
            "custom_within_gene": custom_within,
            "difference": within - custom_within,
            "vienna_pooled": spearman(values, activity),
            "vienna_partial_gc": partial_spearman(values, activity, gc),
            "n_groups": n_groups,
        })
    return results


def incremental_under_vienna(guides, custom, vienna, n_boot: int, seed: int) -> dict:
    """Does adding ViennaRNA folding features help a sequence-only model?"""
    targets = [g.percentile for g in guides]
    genes = [g.group for g in guides]

    merged = [{**c, **v} for c, v in zip(custom, vienna)]

    # Same rule as run_study: leave-one-gene-out where the gene count allows it,
    # grouped five-fold otherwise.
    n_outer = 0 if len({g for g in genes}) <= LOGO_MAX_GENES else 5

    baseline_spec = FeatureSpec()
    baseline_design = SparseDesign(
        build_design_rows(guides, merged, baseline_spec), baseline_spec.n_features
    )
    baseline_predictions, _ = nested_cv_predictions(
        baseline_design, targets, genes, ALPHA_GRID, n_outer=n_outer, n_inner=3,
        seed=seed, intercept_index=baseline_spec.intercept_index,
    )
    baseline_score = mean_within_group_spearman(baseline_predictions, targets, genes)

    models = []
    for name, columns in (
        ("vienna_structure_all", VIENNA_MODEL_FEATURES),
        ("vienna_ensemble_only", ["vienna_seed_ensemble_full", "vienna_mean_unpaired_full",
                                  "vienna_ensemble_energy_full"]),
    ):
        spec = FeatureSpec(structure=columns)
        design = SparseDesign(build_design_rows(guides, merged, spec), spec.n_features)
        predictions, _ = nested_cv_predictions(
            design, targets, genes, ALPHA_GRID, n_outer=n_outer, n_inner=3,
            seed=seed, intercept_index=spec.intercept_index,
        )
        delta, low, high = paired_spearman_delta_ci(
            baseline_predictions, predictions, targets, genes, n_boot=n_boot, seed=seed
        )
        models.append({
            "name": name,
            "spearman": mean_within_group_spearman(predictions, targets, genes),
            "delta": delta,
            "delta_ci_low": low,
            "delta_ci_high": high,
        })

    return {"baseline_spearman": baseline_score, "models": models}


# --------------------------------------- 3 and 4. temperature and context


def sensitivity(guides, custom, dataset: str, n_boot: int, seed: int) -> list[dict]:
    """Repeat the headline correlation under other temperatures and folding contexts."""
    activity = [g.percentile for g in guides]
    genes = [g.group for g in guides]
    guide_ids = [g.guide_id for g in guides]

    settings = [
        ("37 C, global folding", 37.0, None),
        ("25 C, global folding", 25.0, None),
        ("42 C, global folding", 42.0, None),
        ("37 C, local folding (max span 40)", 37.0, 40),
    ]

    out = []
    for label, temperature, max_span in settings:
        path = vienna_path(dataset, temperature, max_span)
        if not path.exists():
            out.append({"setting": label, "available": False})
            continue
        rows = load_vienna(path, guide_ids)
        entry = {"setting": label, "available": True,
                 "temperature": temperature, "max_span": max_span}
        for feature in ("vienna_seed_ensemble_full", "vienna_mean_unpaired_full",
                        "vienna_ensemble_energy_full"):
            values = _column(rows, feature)
            within, low, high, _ = within_group_spearman_ci(
                values, activity, genes, n_boot=n_boot, seed=seed
            )
            entry[feature] = {
                "within_gene": within, "ci_low": low, "ci_high": high,
                "mean": sum(values) / len(values),
            }
        out.append(entry)
    return out


# ------------------------------------------------------------------ report


def _fmt(value, places: int = 3) -> str:
    if value is None or value != value:
        return "n/a"
    return f"{value:.{places}f}"


def report(results: dict) -> str:
    lines: list[str] = []
    add = lines.append

    for dataset, block in results["datasets"].items():
        add(f"\n### {block['label']}\n")
        add("**Agreement with ViennaRNA**\n")
        add("| Measure | Pearson r | Spearman | mean diff | mean abs diff |")
        add("|---|---:|---:|---:|---:|")
        for row in block["agreement"]:
            add(f"| {row['label']} | {_fmt(row['pearson'])} | {_fmt(row['spearman'])} | "
                f"{_fmt(row['mean_difference'])} | {_fmt(row['mean_absolute_difference'])} |")

        add("\n**Correlation with activity, custom against ViennaRNA** "
            "(mean within-gene Spearman)\n")
        add("| Measure | custom | ViennaRNA | 95% interval | difference |")
        add("|---|---:|---:|---|---:|")
        for row in block["correlations"]:
            add(f"| {row['label']} | {_fmt(row['custom_within_gene'])} | "
                f"{_fmt(row['vienna_within_gene'])} | "
                f"[{_fmt(row['vienna_ci_low'])}, {_fmt(row['vienna_ci_high'])}] | "
                f"{_fmt(row['difference'])} |")

        inc = block["incremental"]
        add(f"\n**Incremental value of ViennaRNA folding features.** Sequence-only "
            f"baseline {_fmt(inc['baseline_spearman'])}.\n")
        add("| Added features | held-out rho | change | 95% interval |")
        add("|---|---:|---:|---|")
        for model in inc["models"]:
            add(f"| {model['name']} | {_fmt(model['spearman'])} | {_fmt(model['delta'])} | "
                f"[{_fmt(model['delta_ci_low'])}, {_fmt(model['delta_ci_high'])}] |")

        add("\n**Temperature and folding context**\n")
        add("| Setting | seed accessibility rho | mean unpaired rho | ensemble energy rho |")
        add("|---|---:|---:|---:|")
        for row in block["sensitivity"]:
            if not row.get("available"):
                add(f"| {row['setting']} | not generated | | |")
                continue
            add(f"| {row['setting']} | "
                f"{_fmt(row['vienna_seed_ensemble_full']['within_gene'])} | "
                f"{_fmt(row['vienna_mean_unpaired_full']['within_gene'])} | "
                f"{_fmt(row['vienna_ensemble_energy_full']['within_gene'])} |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate the study against ViennaRNA.")
    parser.add_argument("--output", type=Path, default=OUTPUTS / "vienna_validation.json")
    parser.add_argument("--summary", type=Path, default=OUTPUTS / "vienna_validation.md")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260724)
    args = parser.parse_args()

    labels = {"doench": "Doench 2014 + 2016 pooled (human cells)",
              "crisprscan": "CRISPRscan / Moreno-Mateos 2015 (zebrafish, in vivo)"}

    results: dict = {"datasets": {}, "settings": {
        "bootstrap": args.bootstrap, "seed": args.seed,
        "note": "ViennaRNA uses standard Turner parameters; the custom model "
                "uses published stacking values with simplified loop terms.",
    }}

    for dataset in ("doench", "crisprscan"):
        guides, custom = load_cached(dataset)
        guide_ids = [g.guide_id for g in guides]
        vienna = load_vienna(vienna_path(dataset), guide_ids)
        print(f"[{dataset}] {len(guides)} guides", flush=True)

        results["datasets"][dataset] = {
            "label": labels[dataset],
            "n": len(guides),
            "agreement": agreement(custom, vienna),
            "position_agreement": position_agreement(custom, vienna),
            "correlations": correlations_under_vienna(
                guides, custom, vienna, args.bootstrap, args.seed
            ),
            "incremental": incremental_under_vienna(
                guides, custom, vienna, args.bootstrap, args.seed
            ),
            "sensitivity": sensitivity(guides, custom, dataset, args.bootstrap, args.seed),
        }
        print(f"[{dataset}] done", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    summary = report(results)
    args.summary.write_text(
        "# ViennaRNA validation\n\nGenerated by `validate_vienna.py`.\n" + summary + "\n",
        encoding="utf-8",
    )
    print(summary)
    print(f"\nWrote {args.output} and {args.summary}")


if __name__ == "__main__":
    main()
