"""Reduce the study results to the compact file the website reads.

``analysis_outputs/study_results.json`` is the full record, including bootstrap
settings and fields only the tables use. The site needs a smaller, rounded
subset, so this script derives ``docs/data/study.json`` from it rather than
having the browser parse the whole thing.

Keeping this separate from ``run_study.py`` means the site data can be rebuilt
in a second, without re-running the twenty-five minute analysis.

    python export_site_data.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def _round(value: float, places: int = 4) -> float | None:
    """Round, mapping NaN to null so JSON stays valid and JavaScript sees null."""
    if value != value:
        return None
    return round(value, places)


def _correlation(entry: dict) -> dict:
    return {
        "feature": entry["feature"],
        "label": entry["label"],
        "spearman": _round(entry["spearman"]),
        "ciLow": _round(entry["ci_low"]),
        "ciHigh": _round(entry["ci_high"]),
        "withinGene": _round(entry["within_gene_spearman"]),
        "withinLow": _round(entry["within_gene_ci_low"]),
        "withinHigh": _round(entry["within_gene_ci_high"]),
    }


def _position(entry: dict) -> dict:
    return {
        "position": entry["position"],
        "spearman": _round(entry["spearman"]),
        "ciLow": _round(entry["ci_low"]),
        "ciHigh": _round(entry["ci_high"]),
        "pAdjusted": _round(entry["p_adjusted"]),
        "meanUnpaired": _round(entry["mean_unpaired"]),
    }


def _model(entry: dict) -> dict:
    record = {
        "name": entry["name"],
        "spearman": _round(entry["spearman"]),
        "ciLow": _round(entry["ci_low"]),
        "ciHigh": _round(entry["ci_high"]),
        "nFeatures": entry["n_features"],
    }
    if "delta" in entry:
        record["delta"] = _round(entry["delta"])
        record["deltaLow"] = _round(entry["delta_ci_low"])
        record["deltaHigh"] = _round(entry["delta_ci_high"])
    return record


def build(results: dict) -> dict:
    site: dict = {"generatedBy": "run_study.py", "datasets": {}, "transfer": []}

    for name, block in results["datasets"].items():
        incremental = block["incremental"]
        site["datasets"][name] = {
            "label": block["label"],
            "n": block["n"],
            "nGenes": block["n_genes"],
            "correlations": [_correlation(c) for c in block["correlations"]],
            "positions": {
                context: [_position(p) for p in entries]
                for context, entries in block["positions"].items()
            },
            "baseline": {
                "spearman": _round(incremental["baseline_spearman"]),
                "ciLow": _round(incremental["baseline_ci"][0]),
                "ciHigh": _round(incremental["baseline_ci"][1]),
                "pooled": _round(incremental["baseline_pooled_spearman"]),
                "nGenesScored": incremental["n_scored_genes"],
            },
            "models": [
                _model(m) for m in incremental["models"] if m["name"] != "baseline"
            ],
        }

    for entry in results.get("transfer", []):
        structure = entry["structure_all"]
        site["transfer"].append({
            "transfer": entry["transfer"],
            "nTrain": entry["n_train"],
            "nTest": entry["n_test"],
            "baseline": _round(entry["baseline"]["spearman"]),
            "structure": _round(structure["spearman"]),
            "delta": _round(structure["delta"]),
            "deltaLow": _round(structure["delta_ci_low"]),
            "deltaHigh": _round(structure["delta_ci_high"]),
        })
    return site


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the website's study data file.")
    parser.add_argument("--input", type=Path, default=Path("analysis_outputs/study_results.json"))
    parser.add_argument("--output", type=Path, default=Path("docs/data/study.json"))
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"{args.input} is missing. Run run_study.py first.")
    results = json.loads(args.input.read_text(encoding="utf-8"))
    payload = build(results)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    size = args.output.stat().st_size / 1024
    print(f"Wrote {args.output} ({size:.1f} KB)")


if __name__ == "__main__":
    main()
