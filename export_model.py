"""Train the browser's efficiency model with this repository's own code.

The site needs a small model it can evaluate in JavaScript. This script fits one
using ``stats.py`` - ridge regression on position, dinucleotide, G/C, and two
folding features - and writes it to ``docs/data/model.json``.

Two design choices matter:

**Only features the browser can actually compute.** The study's best structure
features come from the partition function, and there is no McCaskill
implementation in the browser. So the exported model uses the two features the
site's own Zuker folder produces for the full sgRNA: the fraction of the seed
left unpaired in the MFE structure, and the MFE itself. A model whose inputs the
site cannot reproduce would be useless there.

**The reported accuracy is the honest one.** The weights are fitted on every
guide, but the quoted correlation comes from nested cross-validation that holds
out whole genes. Those two numbers are different, and quoting the fitted-on-
everything number would overstate the model badly. Guides targeting the same
gene are not independent, so a random split leaks.

    python export_model.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from model_features import FeatureSpec, build_design_rows
from run_study import ALPHA_GRID, load_cached
from stats import (
    SparseDesign,
    fit_full,
    mean_within_group_spearman,
    nested_cv_predictions,
    select_alpha,
    spearman,
)


#: Folding features the browser can recompute with js/rna-fold.js.
BROWSER_STRUCTURE = ["seed_mfe_full", "mfe_energy_full"]

BASES = "ACGT"


def build(seed: int = 20260724) -> dict:
    guides, cached = load_cached("doench")
    targets = [g.percentile for g in guides]
    genes = [g.group for g in guides]

    spec = FeatureSpec(structure=BROWSER_STRUCTURE)
    design = SparseDesign(build_design_rows(guides, cached, spec), spec.n_features)

    baseline_spec = FeatureSpec()
    baseline = SparseDesign(
        build_design_rows(guides, cached, baseline_spec), baseline_spec.n_features
    )

    print("Cross-validating (grouped by gene, nested penalty selection)...")
    predictions, _ = nested_cv_predictions(
        design, targets, genes, ALPHA_GRID, n_outer=5, n_inner=3,
        seed=seed, intercept_index=spec.intercept_index,
    )
    held_out = mean_within_group_spearman(predictions, targets, genes)

    baseline_predictions, _ = nested_cv_predictions(
        baseline, targets, genes, ALPHA_GRID, n_outer=5, n_inner=3,
        seed=seed, intercept_index=baseline_spec.intercept_index,
    )
    held_out_sequence_only = mean_within_group_spearman(baseline_predictions, targets, genes)

    gc_values = [float(row["gc_percent"]) for row in cached]
    gc_only = mean_within_group_spearman(gc_values, targets, genes)
    gc_pooled = spearman(gc_values, targets)

    print("Fitting the exported weights on every guide...")
    alpha = select_alpha(
        design, targets, genes, ALPHA_GRID, seed=seed,
        intercept_index=spec.intercept_index,
    )
    weights = fit_full(design, targets, alpha, intercept_index=spec.intercept_index)

    position = weights[spec.position_offset : spec.position_offset + 80]
    dinucleotide = weights[spec.dinuc_offset : spec.dinuc_offset + 304]

    return {
        "kind": "ridge-linear-v1",
        "intercept": weights[spec.intercept_index],
        "position": position,
        "dinucleotide": dinucleotide,
        "gc": weights[spec.gc_offset],
        "structure": {
            "seedOpenness": weights[spec.structure_offset],
            "mfeEnergyScaled": weights[spec.structure_offset + 1],
        },
        "meta": {
            "features": (
                "position one-hot (20x4) + dinucleotide one-hot (19x16) + G/C fraction "
                "+ seed openness and MFE of the full sgRNA"
            ),
            "nfeat": spec.n_features,
            "params": spec.n_features,
            "alpha": alpha,
            "trainedOn": "Doench 2014 + 2016 pooled (4,685 guides), activity as within-screen percentile",
            "trainN": len(guides),
            "nGenes": len({g.group for g in guides}),
            "validation": (
                "5-fold cross-validation holding out whole genes, penalty chosen inside "
                "each training fold, scored as mean Spearman within each gene"
            ),
            "heldOutSpearman": round(held_out, 4),
            "heldOutSpearmanSequenceOnly": round(held_out_sequence_only, 4),
            "baselineGcSpearman": round(gc_only, 4),
            "baselineGcSpearmanPooled": round(gc_pooled, 4),
            "energyScale": 10.0,
            "bases": BASES,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Train and export the browser model.")
    parser.add_argument("--output", type=Path, default=Path("docs/data/model.json"))
    parser.add_argument("--seed", type=int, default=20260724)
    args = parser.parse_args()

    payload = build(seed=args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    meta = payload["meta"]
    print(f"\nWrote {args.output} ({meta['params']} weights)")
    print(f"  held-out Spearman, with structure : {meta['heldOutSpearman']:.3f}")
    print(f"  held-out Spearman, sequence only  : {meta['heldOutSpearmanSequenceOnly']:.3f}")
    print(f"  G/C content alone                 : {meta['baselineGcSpearman']:.3f}")
    print(f"  ridge penalty                     : {meta['alpha']}")


if __name__ == "__main__":
    main()
