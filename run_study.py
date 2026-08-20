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
    partial_spearman,
    pearson,
    repeated_cv_predictions,
    select_alpha,
    spearman,
    spearman_permutation_p,
    within_group_permutation_p,
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
    gc = [float(row["gc_percent"]) for row in rows]
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
        # G/C content is the obvious confound: composition drives both pairing
        # and activity, so a raw structure correlation may just be reporting G/C.
        partial = (
            float("nan") if column == "gc_percent"
            else partial_spearman(values, activity, gc)
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
            "partial_spearman_gc": partial,
            "spearman_with_gc": spearman(values, gc),
            # Gene-blocked: activity is shuffled inside each gene, never across genes.
            "permutation_p": within_group_permutation_p(
                values, activity, genes, n_perm=n_perm, seed=seed
            ),
            # Reported only to contrast with the gene-blocked test above, which is
            # the one the conclusions use, so a coarser estimate is enough.
            "pooled_permutation_p": spearman_permutation_p(
                values, activity, n_perm=min(n_perm, 200), seed=seed
            ),
            "n": len(values),
        })
    return results


def experiment_context(rows: list[dict[str, str]]) -> dict:
    """How much folding the spacer alone differs from folding the real sgRNA.

    Cheap to compute (means and correlations, no resampling), so it is also
    recalculated in ``--report-only`` mode. It exists here rather than as a
    side calculation so the numbers quoted in ANALYSIS.md come from the pipeline.
    """
    def column(name: str) -> list[float]:
        return [float(row[name]) for row in rows]

    pairs = [
        ("seed_mfe", "seed_mfe_spacer", "seed_mfe_full",
         "MFE seed accessibility"),
        ("seed_ensemble", "seed_ensemble_spacer", "seed_ensemble_full",
         "Ensemble seed accessibility"),
        ("mean_unpaired", "mean_unpaired_spacer", "mean_unpaired_full",
         "Mean unpaired probability"),
        ("ensemble_energy", "ensemble_energy_spacer", "ensemble_energy_full",
         "Ensemble free energy (kcal/mol)"),
    ]
    measures = []
    for key, alone_col, full_col, label in pairs:
        alone = column(alone_col)
        full = column(full_col)
        measures.append({
            "key": key,
            "label": label,
            "mean_spacer_alone": sum(alone) / len(alone),
            "mean_in_sgrna": sum(full) / len(full),
            "pearson_between_contexts": pearson(alone, full),
            "spearman_between_contexts": spearman(alone, full),
        })

    unstructured = sum(1 for v in column("mean_unpaired_spacer") if v > 0.99)
    return {
        "measures": measures,
        "n": len(rows),
        "essentially_unstructured_alone": unstructured,
        "essentially_unstructured_fraction": unstructured / len(rows),
    }


# --------------------------------------------------------- 2. position-resolved


def experiment_positions(
    guides: list[Guide], rows: list[dict[str, str]], n_boot: int, seed: int, n_perm: int
) -> dict[str, list[dict]]:
    """Correlate unpaired probability with activity at each of the 20 positions.

    Everything here is gene-aware. The reported correlation is the mean within
    gene, the interval comes from resampling genes, and the p-value comes from a
    permutation that shuffles activity **inside** each gene. Benjamini-Hochberg
    is then applied to those gene-blocked p-values across the 20 positions.

    An earlier version reported pooled correlations with p-values from a free
    shuffle. That test ignores the clustering of guides within genes, so it
    rejects far too readily; the counts it produced were not trustworthy and
    have been replaced rather than merely annotated.
    """
    activity = [g.percentile for g in guides]
    genes = [g.group for g in guides]
    output: dict[str, list[dict]] = {}

    for context, prefix in (("spacer", "unpaired_spacer"), ("full", "unpaired_full")):
        entries = []
        for position in range(1, 21):
            values = [float(row[f"{prefix}_{position:02d}"]) for row in rows]
            spread = max(values) - min(values)
            if spread == 0.0:
                # A position unpaired in every guide carries no information.
                entries.append({
                    "position": position, "spearman": 0.0, "pooled_spearman": 0.0,
                    "ci_low": float("nan"), "ci_high": float("nan"), "p": 1.0,
                    "mean_unpaired": values[0], "constant": True,
                })
                continue

            within, within_low, within_high, _ = within_group_spearman_ci(
                values, activity, genes, n_boot=n_boot, seed=seed + position
            )
            entries.append({
                "position": position,
                # Headline: within gene, with a gene-resampled interval.
                "spearman": within,
                "ci_low": within_low,
                "ci_high": within_high,
                # Kept for comparison with the pooled view.
                "pooled_spearman": spearman(values, activity),
                "p": within_group_permutation_p(
                    values, activity, genes, n_perm=n_perm, seed=seed + position
                ),
                "mean_unpaired": sum(values) / len(values),
                "constant": False,
            })
        adjusted = benjamini_hochberg([entry["p"] for entry in entries])
        for entry, value in zip(entries, adjusted):
            entry["p_adjusted"] = value
        output[context] = entries
    return output


# ------------------------------------------------------- 3. incremental value


#: Above this many genes, leave-one-gene-out costs one full nested fit per gene
#: and the memory for one Gram matrix per gene, which is not worth it. Below it,
#: leave-one-gene-out is preferred because it removes fold assignment entirely.
LOGO_MAX_GENES = 25


def _evaluate_spec(payload):
    """Run grouped cross-validation for one feature set, in a worker process.

    Two schemes, because they answer different worries.

    *Leave-one-gene-out* removes any dependence on how genes happened to be
    assigned to folds. With the 18-gene human screen that is affordable and is
    used as the headline. With 111 zebrafish genes it would mean 111 nested fits
    and 111 resident Gram matrices per feature set, so five-fold stands in.

    *Repeated five-fold* re-runs over several gene-to-fold assignments. The
    spread across those assignments is reported, because with few genes that
    spread can be larger than the effect being tested, and a single split can
    make a null look like a finding.
    """
    name, spec_kwargs, rows, targets, genes, seeds = payload
    spec = FeatureSpec(**spec_kwargs)
    design = SparseDesign(rows, spec.n_features)
    started = time.perf_counter()

    n_genes = len({g for g in genes})
    use_logo = n_genes <= LOGO_MAX_GENES

    if use_logo:
        # n_outer <= 0 means one fold per gene. There is no fold assignment to
        # vary, so repeating over seeds would measure nothing: every repeat would
        # return the identical partition. The spread is reported only for the
        # screen that genuinely has to choose an assignment.
        predictions, alphas = nested_cv_predictions(
            design, targets, genes, ALPHA_GRID, n_outer=0, n_inner=3,
            seed=seeds[0], intercept_index=spec.intercept_index,
        )
        scheme = "leave-one-gene-out"
        repeated_scores = []
    else:
        runs = repeated_cv_predictions(
            design, targets, genes, ALPHA_GRID, n_outer=5, n_inner=3,
            seeds=seeds, intercept_index=spec.intercept_index,
        )
        repeated_scores = [
            mean_within_group_spearman(p, targets, genes) for p, _ in runs
        ]
        predictions, alphas = runs[0]
        scheme = "grouped 5-fold, repeated over 3 gene-to-fold assignments"

    return {
        "name": name,
        "description": spec.describe(),
        "n_features": spec.n_features,
        "scheme": scheme,
        "predictions": predictions,
        "alphas": alphas,
        "within_gene_spearman": mean_within_group_spearman(predictions, targets, genes),
        "pooled_spearman": spearman(predictions, targets),
        "repeated_scores": repeated_scores,
        "repeated_mean": (sum(repeated_scores) / len(repeated_scores)) if repeated_scores else None,
        "repeated_min": min(repeated_scores) if repeated_scores else None,
        "repeated_max": max(repeated_scores) if repeated_scores else None,
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

    fold_seeds = [seed, seed + 101, seed + 202]
    jobs = []
    for name, kwargs in specs.items():
        spec = FeatureSpec(**kwargs)
        design_rows = build_design_rows(guides, rows, spec)
        jobs.append((name, kwargs, design_rows, targets, genes, fold_seeds))

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
            "repeated_scores": entry["repeated_scores"],
            "repeated_mean": entry["repeated_mean"],
            "repeated_min": entry["repeated_min"],
            "repeated_max": entry["repeated_max"],
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
        "baseline_repeated": baseline["repeated_scores"],
        "n_scored_genes": n_groups,
        "scheme": baseline["scheme"],
        "validation": (
            f"{baseline['scheme']}; the ridge penalty is chosen by within-gene "
            "Spearman inside each training fold, never on the fold being predicted"
        ),
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
        add("Correlations with measured activity. The within-gene column is the headline: "
            "it compares guides only against others aimed at the same gene. The p-value "
            "comes from shuffling activity inside each gene, never across genes. The "
            "partial column removes G/C content, the obvious confound.\n")
        add("| Accessibility measure | pooled rho | within-gene rho | 95% CI | partial (G/C removed) | rho with G/C | gene-blocked p |")
        add("|---|---:|---:|---|---:|---:|---:|")
        for row in block["correlations"]:
            add(
                f"| {row['label']} | {_fmt(row['spearman'])} | "
                f"{_fmt(row['within_gene_spearman'])} | "
                f"[{_fmt(row['within_gene_ci_low'])}, {_fmt(row['within_gene_ci_high'])}] | "
                f"{_fmt(row.get('partial_spearman_gc', float('nan')))} | "
                f"{_fmt(row.get('spearman_with_gc', float('nan')))} | "
                f"{_fmt(row['permutation_p'], 4)} |"
            )

        if "context" in block:
            ctx = block["context"]
            add("\n**Folding context.** Each spacer folded alone, and folded as part "
                "of the full sgRNA.\n")
            add("| Measure | spacer alone | in the sgRNA | correlation between the two |")
            add("|---|---:|---:|---:|")
            for row in ctx["measures"]:
                add(f"| {row['label']} | {row['mean_spacer_alone']:.3f} | "
                    f"{row['mean_in_sgrna']:.3f} | {row['pearson_between_contexts']:.3f} |")
            add(f"\n{ctx['essentially_unstructured_alone']} of {ctx['n']} spacers "
                f"({ctx['essentially_unstructured_fraction']:.1%}) are essentially "
                f"unstructured when folded alone, meaning a mean unpaired probability "
                f"above 0.99.\n")

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
            add(f"Validation scheme: {inc['validation']}.\n")
            spread = inc.get("baseline_repeated") or []
            if spread:
                scores = ", ".join(_fmt(v) for v in spread)
                add(
                    f"Across {len(spread)} gene-to-fold assignments the baseline scored "
                    f"{scores}, so the spread from the split alone is about "
                    f"{_fmt(max(spread) - min(spread))}.\n"
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
    parser.add_argument("--report-only", action="store_true",
                        help="Rebuild the summary from an existing results JSON.")
    args = parser.parse_args()

    if args.report_only:
        if not args.output.exists():
            raise SystemExit(f"{args.output} is missing; run the study first.")
        saved = json.loads(args.output.read_text(encoding="utf-8"))
        # The context block needs no resampling, so it can be filled in here for
        # results produced before it existed.
        for name, block in saved["datasets"].items():
            if "context" not in block:
                _, cached_rows = load_cached(name)
                block["context"] = experiment_context(cached_rows)
        args.output.write_text(json.dumps(saved, indent=2), encoding="utf-8")
        text = report(saved)
        args.summary.write_text(
            "# Study results\n\nGenerated by `run_study.py`.\n" + text + "\n",
            encoding="utf-8",
        )
        print(text)
        print(f"\nRebuilt {args.summary} from {args.output}")
        return

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
            "context": experiment_context(rows),
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
