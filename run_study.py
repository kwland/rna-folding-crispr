"""The study: does RNA secondary structure predict CRISPR guide activity?

Reads the cached folding features written by ``compute_features.py`` and runs
four experiments, each aimed at one way the original single-number analysis
could have been wrong.

1. **Measurement.** Seed accessibility under five definitions - Nussinov pair
   counting, Zuker MFE, and the Boltzmann ensemble, each for the spacer alone
   and (for the energy models) for the spacer inside the full sgRNA. If the
   original near-zero correlation was an artefact of collapsing the ensemble to
   one structure, a better measure should recover a signal.

2. **Position.** "The last 8 nt" is an assumption, not a finding. With per-base
   unpaired probabilities the correlation can be computed at every spacer
   position and the assumption tested directly.

3. **Incremental value.** The publishable question is not whether structure
   correlates with activity, but whether it adds anything to a sequence-only
   baseline. Measured as the change in held-out Spearman under grouped nested
   cross-validation, with a paired confidence interval.

4. **Replication.** Any of the above could be a property of one screen. Every
   headline number is recomputed on an independent screen (CRISPRscan, zebrafish
   in vivo) and across screen-to-screen transfers.

    python run_study.py                      # everything
    python run_study.py --skip-models        # experiments 1 and 2 only (fast)
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from datasets import Guide, load_crisprscan, load_doench_pooled
from model_features import (
    STRUCTURE_SUBSETS,
    FeatureSpec,
    build_design_rows,
    structure_feature_columns,
)
from stats import (
    SparseDesign,
    benjamini_hochberg,
    cluster_bootstrap_ci,
    fit_full,
    mean_within_group_spearman,
    nested_cv_predictions,
    paired_spearman_delta_ci,
    select_alpha,
    spearman,
    spearman_permutation_p,
    within_group_spearman_ci,
)


CACHE = {
    "doench": Path("analysis_outputs/features_doench_pooled.csv"),
    "crisprscan": Path("analysis_outputs/features_crisprscan.csv"),
}
LOADERS = {"doench": load_doench_pooled, "crisprscan": load_crisprscan}

#: Ridge penalties searched inside each training fold.
ALPHA_GRID = [1.0, 10.0, 100.0, 1000.0]

ACCESSIBILITY_MEASURES = [
    ("seed_nussinov", "Nussinov seed accessibility (spacer alone)"),
    ("seed_mfe_spacer", "MFE seed accessibility (spacer alone)"),
    ("seed_ensemble_spacer", "Ensemble seed accessibility (spacer alone)"),
    ("seed_mfe_full", "MFE seed accessibility (spacer + scaffold)"),
    ("seed_ensemble_full", "Ensemble seed accessibility (spacer + scaffold)"),
    ("mean_unpaired_spacer", "Mean unpaired probability, whole spacer (alone)"),
    ("mean_unpaired_full", "Mean unpaired probability, whole spacer (in sgRNA)"),
    ("ensemble_energy_full", "Ensemble free energy (spacer + scaffold)"),
    ("gc_percent", "G/C percent"),
]


def load_cached(dataset: str) -> tuple[list[Guide], list[dict[str, str]]]:
    """Load a screen together with its cached folding features, keeping order."""
    path = CACHE[dataset]
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. Run: python compute_features.py --dataset {dataset}"
        )
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    guides = LOADERS[dataset]()
    if len(rows) != len(guides):
        raise ValueError(
            f"{path} has {len(rows)} rows but the screen has {len(guides)} guides; "
            "re-run compute_features.py"
        )
    for guide, row in zip(guides, rows):
        if guide.guide_id != row["guide_id"]:
            raise ValueError(f"Cached features are out of order at {guide.guide_id}")
    return guides, rows


# ----------------------------------------------------- 1. accessibility measures


def experiment_correlations(
    guides: list[Guide], rows: list[dict[str, str]], n_boot: int, seed: int, n_perm: int
) -> list[dict]:
    """Correlate each accessibility measure with activity, pooled and within gene.

    Both are reported because they answer different questions and can disagree.
    The pooled number mixes guides from genes with very different baseline
    editability; the within-gene number asks the question a designer asks.
    """
    activity = [g.percentile for g in guides]
    genes = [g.group for g in guides]
    results = []
    for column, label in ACCESSIBILITY_MEASURES:
        values = [float(row[column]) for row in rows]
        rho = spearman(values, activity)
        low, high = cluster_bootstrap_ci(
            lambda idx: spearman([values[i] for i in idx], [activity[i] for i in idx]),
            genes, n_boot=n_boot, seed=seed,
        )
        within, within_low, within_high, n_groups = within_group_spearman_ci(
            values, activity, genes, n_boot=n_boot, seed=seed
        )
        results.append({
            "feature": column,
            "label": label,
            "spearman": rho,
            "ci_low": low,
            "ci_high": high,
            "within_gene_spearman": within,
            "within_gene_ci_low": within_low,
            "within_gene_ci_high": within_high,
            "n_groups": n_groups,
            "permutation_p": spearman_permutation_p(values, activity, n_perm=n_perm, seed=seed),
            "n": len(values),
        })
    return results


# --------------------------------------------------------- 2. position-resolved


def experiment_positions(
    guides: list[Guide], rows: list[dict[str, str]], n_boot: int, seed: int, n_perm: int
) -> dict[str, list[dict]]:
    activity = [g.percentile for g in guides]
    genes = [g.group for g in guides]
    output: dict[str, list[dict]] = {}

    for context, prefix in (("spacer", "unpaired_spacer"), ("full", "unpaired_full")):
        entries = []
        for position in range(1, 21):
            values = [float(row[f"{prefix}_{position:02d}"]) for row in rows]
            spread = max(values) - min(values)
            if spread == 0.0:
                # A position that is unpaired in every guide carries no information.
                entries.append({
                    "position": position, "spearman": 0.0, "ci_low": float("nan"),
                    "ci_high": float("nan"), "p": 1.0, "mean_unpaired": values[0],
                    "within_gene_spearman": 0.0, "constant": True,
                })
                continue
            rho = spearman(values, activity)
            low, high = cluster_bootstrap_ci(
                lambda idx: spearman([values[i] for i in idx], [activity[i] for i in idx]),
                genes, n_boot=n_boot, seed=seed + position,
            )
            within, _, _, _ = within_group_spearman_ci(
                values, activity, genes, n_boot=1, seed=seed
            )
            entries.append({
                "position": position,
                "spearman": rho,
                "ci_low": low,
                "ci_high": high,
                "within_gene_spearman": within,
                "p": spearman_permutation_p(values, activity, n_perm=n_perm, seed=seed + position),
                "mean_unpaired": sum(values) / len(values),
                "constant": False,
            })
        adjusted = benjamini_hochberg([entry["p"] for entry in entries])
        for entry, value in zip(entries, adjusted):
            entry["p_adjusted"] = value
        output[context] = entries
    return output


# ------------------------------------------------------- 3. incremental value


def _evaluate_spec(payload):
    """Run nested grouped CV for one feature set. Top-level, for worker processes."""
    name, spec_kwargs, rows, targets, genes, seed = payload
    spec = FeatureSpec(**spec_kwargs)
    design = SparseDesign(rows, spec.n_features)
    started = time.perf_counter()
    predictions, alphas = nested_cv_predictions(
        design, targets, genes, ALPHA_GRID, n_outer=5, n_inner=3,
        seed=seed, intercept_index=spec.intercept_index,
    )
    return {
        "name": name,
        "description": spec.describe(),
        "n_features": spec.n_features,
        "predictions": predictions,
        "alphas": alphas,
        # Headline metric: ranking within a gene. See stats.mean_within_group_spearman.
        "within_gene_spearman": mean_within_group_spearman(predictions, targets, genes),
        # Reported alongside so the gap between the two is visible rather than hidden.
        "pooled_spearman": spearman(predictions, targets),
        "seconds": time.perf_counter() - started,
    }


def experiment_incremental(
    guides: list[Guide],
    rows: list[dict[str, str]],
    n_boot: int,
    seed: int,
    workers: int,
) -> dict:
    targets = [g.percentile for g in guides]
    genes = [g.group for g in guides]

    specs: dict[str, dict] = {"baseline": {}}
    for name, columns in STRUCTURE_SUBSETS.items():
        specs[name] = {"structure": columns}

    jobs = []
    for name, kwargs in specs.items():
        spec = FeatureSpec(**kwargs)
        design_rows = build_design_rows(guides, rows, spec)
        jobs.append((name, kwargs, design_rows, targets, genes, seed))

    if workers <= 1:
        evaluated = [_evaluate_spec(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            evaluated = list(pool.map(_evaluate_spec, jobs))

    by_name = {entry["name"]: entry for entry in evaluated}
    baseline = by_name["baseline"]

    summary = []
    for entry in evaluated:
        mean, low, high, n_groups = within_group_spearman_ci(
            entry["predictions"], targets, genes, n_boot=n_boot, seed=seed
        )
        record = {
            "name": entry["name"],
            "description": entry["description"],
            "n_features": entry["n_features"],
            "spearman": mean,
            "ci_low": low,
            "ci_high": high,
            "n_groups": n_groups,
            "pooled_spearman": entry["pooled_spearman"],
            "alphas": entry["alphas"],
            "seconds": entry["seconds"],
        }
        if entry["name"] != "baseline":
            delta, delta_low, delta_high = paired_spearman_delta_ci(
                baseline["predictions"], entry["predictions"], targets, genes,
                n_boot=n_boot, seed=seed,
            )
            record.update({
                "delta": delta, "delta_ci_low": delta_low, "delta_ci_high": delta_high
            })
        summary.append(record)

    mean, low, high, n_groups = within_group_spearman_ci(
        baseline["predictions"], targets, genes, n_boot=n_boot, seed=seed
    )
    return {
        "baseline_spearman": mean,
        "baseline_ci": [low, high],
        "baseline_pooled_spearman": baseline["pooled_spearman"],
        "n_scored_genes": n_groups,
        "models": summary,
    }


# ----------------------------------------------------------- 4. cross-screen


def _subset(guides, rows, predicate):
    keep = [i for i, g in enumerate(guides) if predicate(g)]
    return [guides[i] for i in keep], [rows[i] for i in keep]


def experiment_transfer(
    sources: dict[str, tuple[list[Guide], list[dict[str, str]]]],
    n_boot: int,
    seed: int,
) -> list[dict]:
    """Train on one screen, predict a different one, with no retuning."""
    doench_guides, doench_rows = sources["doench"]
    scan_guides, scan_rows = sources["crisprscan"]

    d2016 = _subset(doench_guides, doench_rows, lambda g: g.screen == "Doench2016")
    d2014 = _subset(doench_guides, doench_rows, lambda g: g.screen == "Doench2014")

    pairs = [
        ("Doench2016 -> Doench2014", d2016, d2014),
        ("Doench2014 -> Doench2016", d2014, d2016),
        ("Doench pooled -> CRISPRscan", (doench_guides, doench_rows), (scan_guides, scan_rows)),
        ("CRISPRscan -> Doench pooled", (scan_guides, scan_rows), (doench_guides, doench_rows)),
    ]
    variants = {
        "baseline": {},
        "structure_all": {"structure": STRUCTURE_SUBSETS["structure_all"]},
        "structure_ensemble_full": {
            "structure": STRUCTURE_SUBSETS["structure_ensemble_full"]
        },
    }

    results = []
    for label, (train_guides, train_rows), (test_guides, test_rows) in pairs:
        entry = {"transfer": label, "n_train": len(train_guides), "n_test": len(test_guides)}
        test_truth = [g.percentile for g in test_guides]
        test_genes = [g.group for g in test_guides]
        predictions_by_variant = {}
        for name, kwargs in variants.items():
            spec = FeatureSpec(**kwargs)
            train_design = SparseDesign(
                build_design_rows(train_guides, train_rows, spec), spec.n_features
            )
            test_design = SparseDesign(
                build_design_rows(test_guides, test_rows, spec), spec.n_features
            )
            alpha = select_alpha(
                train_design, [g.percentile for g in train_guides],
                [g.group for g in train_guides], ALPHA_GRID, seed=seed,
                intercept_index=spec.intercept_index,
            )
            weights = fit_full(
                train_design, [g.percentile for g in train_guides], alpha,
                intercept_index=spec.intercept_index,
            )
            predicted = test_design.predict(weights)
            predictions_by_variant[name] = predicted
            mean, low, high, n_groups = within_group_spearman_ci(
                predicted, test_truth, test_genes, n_boot=n_boot, seed=seed
            )
            entry[name] = {
                "alpha": alpha,
                "spearman": mean,
                "ci_low": low,
                "ci_high": high,
                "n_groups": n_groups,
                "pooled_spearman": spearman(predicted, test_truth),
            }

        for name in variants:
            if name == "baseline":
                continue
            delta, low, high = paired_spearman_delta_ci(
                predictions_by_variant["baseline"], predictions_by_variant[name],
                test_truth, test_genes, n_boot=n_boot, seed=seed,
            )
            entry[name].update({"delta": delta, "delta_ci_low": low, "delta_ci_high": high})
        results.append(entry)
    return results


# ------------------------------------------------------------------ reporting


def _fmt(value: float, places: int = 3) -> str:
    if value != value:
        return "n/a"
    return f"{value:.{places}f}"


def report(results: dict) -> str:
    lines: list[str] = []
    add = lines.append

    for dataset, block in results["datasets"].items():
        add(f"\n### {block['label']}  (n = {block['n']} guides, {block['n_genes']} genes)\n")
        add("| Accessibility measure | pooled rho | 95% CI | within-gene rho | 95% CI | perm. p |")
        add("|---|---:|---|---:|---|---:|")
        for row in block["correlations"]:
            add(
                f"| {row['label']} | {_fmt(row['spearman'])} | "
                f"[{_fmt(row['ci_low'])}, {_fmt(row['ci_high'])}] | "
                f"{_fmt(row['within_gene_spearman'])} | "
                f"[{_fmt(row['within_gene_ci_low'])}, {_fmt(row['within_gene_ci_high'])}] | "
                f"{_fmt(row['permutation_p'], 4)} |"
            )

        for context, entries in block["positions"].items():
            name = "spacer alone" if context == "spacer" else "spacer + scaffold"
            significant = [e for e in entries if e["p_adjusted"] < 0.05]
            strongest = max(entries, key=lambda e: abs(e["spearman"]))
            add(
                f"\n**Position-resolved ({name}):** strongest position "
                f"{strongest['position']} (rho = {_fmt(strongest['spearman'])}); "
                f"{len(significant)} of 20 positions significant after "
                f"Benjamini-Hochberg correction."
            )

        if "incremental" in block:
            inc = block["incremental"]
            add(
                f"\n**Held-out sequence baseline**, scored as mean within-gene Spearman over "
                f"{inc['n_scored_genes']} genes: **{_fmt(inc['baseline_spearman'])}** "
                f"[{_fmt(inc['baseline_ci'][0])}, {_fmt(inc['baseline_ci'][1])}]. "
                f"Pooling the same predictions across genes instead gives "
                f"{_fmt(inc['baseline_pooled_spearman'])}, because a model cannot guess the "
                f"baseline editability of a gene it has never seen.\n"
            )
            add("| Added structure features | held-out rho | change vs baseline | 95% CI of change |")
            add("|---|---:|---:|---|")
            for model in inc["models"]:
                if model["name"] == "baseline":
                    continue
                add(
                    f"| {model['name']} | {_fmt(model['spearman'])} | "
                    f"{_fmt(model.get('delta', float('nan')))} | "
                    f"[{_fmt(model.get('delta_ci_low', float('nan')))}, "
                    f"{_fmt(model.get('delta_ci_high', float('nan')))}] |"
                )

    if results.get("transfer"):
        add("\n### Cross-screen transfer\n")
        add("Trained on one screen, applied to another with no retuning. "
            "Scored as mean within-gene Spearman.\n")
        add("| Train -> test | baseline rho | +structure rho | change | 95% CI of change |")
        add("|---|---:|---:|---:|---|")
        for row in results["transfer"]:
            structure = row["structure_all"]
            add(
                f"| {row['transfer']} | {_fmt(row['baseline']['spearman'])} | "
                f"{_fmt(structure['spearman'])} | {_fmt(structure.get('delta', float('nan')))} | "
                f"[{_fmt(structure.get('delta_ci_low', float('nan')))}, "
                f"{_fmt(structure.get('delta_ci_high', float('nan')))}] |"
            )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the CRISPR structure study.")
    parser.add_argument("--output", type=Path, default=Path("analysis_outputs/study_results.json"))
    parser.add_argument("--summary", type=Path, default=Path("analysis_outputs/study_summary.md"))
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260724)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--skip-models", action="store_true",
                        help="Skip the cross-validated models (experiments 3 and 4).")
    args = parser.parse_args()

    labels = {"doench": "Doench 2014 + 2016 pooled (human cells)",
              "crisprscan": "CRISPRscan / Moreno-Mateos 2015 (zebrafish, in vivo)"}

    sources = {name: load_cached(name) for name in CACHE}
    results: dict = {"datasets": {}, "settings": {
        "bootstrap": args.bootstrap, "permutations": args.permutations,
        "seed": args.seed, "alpha_grid": ALPHA_GRID,
        "structure_features": structure_feature_columns(),
        "primary_metric": "mean within-gene Spearman (see stats.mean_within_group_spearman)",
    }}

    for name, (guides, rows) in sources.items():
        print(f"[{name}] {len(guides)} guides", flush=True)
        block = {
            "label": labels[name],
            "n": len(guides),
            "n_genes": len({g.group for g in guides}),
            "correlations": experiment_correlations(
                guides, rows, args.bootstrap, args.seed, args.permutations
            ),
            "positions": experiment_positions(
                guides, rows, args.bootstrap, args.seed, args.permutations
            ),
        }
        if not args.skip_models:
            print(f"[{name}] cross-validated models...", flush=True)
            started = time.perf_counter()
            block["incremental"] = experiment_incremental(
                guides, rows, args.bootstrap, args.seed, args.workers
            )
            print(f"[{name}] models done in {(time.perf_counter() - started) / 60:.1f} min", flush=True)
        results["datasets"][name] = block

    if not args.skip_models:
        print("cross-screen transfer...", flush=True)
        results["transfer"] = experiment_transfer(sources, args.bootstrap, args.seed)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    summary = report(results)
    args.summary.write_text(
        "# Study results\n\nGenerated by `run_study.py`.\n" + summary + "\n",
        encoding="utf-8",
    )
    print(summary)
    print(f"\nWrote {args.output} and {args.summary}")


if __name__ == "__main__":
    main()
