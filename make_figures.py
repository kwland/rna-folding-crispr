"""Create the repository's SVG figures with the Python standard library."""

from __future__ import annotations

import argparse
import csv
import json
import math
from html import escape
from pathlib import Path


INK = "#16211d"
MUTED = "#61706a"
PAPER = "#f7f4ec"
GREEN = "#2f7d62"
GOLD = "#d6a84a"
LINE = "#d8ded8"
RED = "#b4543a"
FONT = "Arial, sans-serif"


def _svg_text(x: float, y: float, text: str, **attrs: object) -> str:
    options = " ".join(f'{key.replace("_", "-")}="{escape(str(value))}"' for key, value in attrs.items())
    return f'<text x="{x}" y="{y}" {options}>{escape(text)}</text>'


def write_workflow_figure(path: Path) -> None:
    boxes = [
        (45, "1", "Sequence", "RNA or DNA input"),
        (300, "2", "Dynamic programming", "Maximize allowed pairs"),
        (555, "3", "Traceback", "Recover one optimum"),
        (810, "4", "Interpret", "Structure + guide features"),
    ]
    elements = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1080" height="330" viewBox="0 0 1080 330" role="img" aria-labelledby="title desc">',
        '<title id="title">Nussinov analysis workflow</title>',
        '<desc id="desc">Four steps: sequence input, dynamic programming, traceback, and interpretation.</desc>',
        f'<rect width="1080" height="330" fill="{PAPER}" rx="24"/>',
        _svg_text(45, 54, "How a sequence becomes an explainable fold", fill=INK, font_size=28, font_weight=700, font_family="Arial, sans-serif"),
        _svg_text(45, 82, "The Python path is deliberately small: each result can be inspected from input to output.", fill=MUTED, font_size=16, font_family="Arial, sans-serif"),
    ]
    for index, (x, number, title, subtitle) in enumerate(boxes):
        elements.extend(
            [
                f'<rect x="{x}" y="125" width="215" height="140" rx="16" fill="#ffffff" stroke="{LINE}" stroke-width="2"/>',
                f'<circle cx="{x + 28}" cy="156" r="15" fill="{GREEN}"/>',
                _svg_text(x + 28, 162, number, fill="#ffffff", font_size=15, font_weight=700, text_anchor="middle", font_family="Arial, sans-serif"),
                _svg_text(x + 22, 205, title, fill=INK, font_size=18, font_weight=700, font_family="Arial, sans-serif"),
                _svg_text(x + 22, 234, subtitle, fill=MUTED, font_size=14, font_family="Arial, sans-serif"),
            ]
        )
        if index < len(boxes) - 1:
            elements.append(f'<path d="M {x + 220} 195 H {x + 246}" stroke="{GOLD}" stroke-width="4" stroke-linecap="round"/>')
            elements.append(f'<path d="M {x + 240} 189 L {x + 248} 195 L {x + 240} 201" fill="none" stroke="{GOLD}" stroke-width="4" stroke-linecap="round" stroke-linejoin="round"/>')
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def read_features(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row.get("analysis_status") == "ok"]
    if not rows:
        raise ValueError(f"No analyzed rows found in {path}")
    return rows


def write_feature_figure(path: Path, rows: list[dict[str, str]]) -> None:
    width = 1080
    height = 190 + 74 * len(rows)
    plot_left = 300
    plot_width = 680
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Guide feature comparison</title>',
        '<desc id="desc">Horizontal bars compare GC fraction and Nussinov seed accessibility for illustrative guide sequences.</desc>',
        f'<rect width="{width}" height="{height}" fill="{PAPER}" rx="24"/>',
        _svg_text(45, 54, "Illustrative guide features", fill=INK, font_size=28, font_weight=700, font_family="Arial, sans-serif"),
        _svg_text(45, 82, "These synthetic examples demonstrate the pipeline; they are not measurements of editing activity.", fill=MUTED, font_size=16, font_family="Arial, sans-serif"),
        f'<rect x="45" y="108" width="14" height="14" rx="3" fill="{GREEN}"/>',
        _svg_text(68, 120, "G/C fraction", fill=INK, font_size=14, font_family="Arial, sans-serif"),
        f'<rect x="180" y="108" width="14" height="14" rx="3" fill="{GOLD}"/>',
        _svg_text(203, 120, "Seed accessibility", fill=INK, font_size=14, font_family="Arial, sans-serif"),
    ]
    for tick in range(0, 101, 25):
        x = plot_left + plot_width * tick / 100
        elements.append(f'<line x1="{x}" y1="148" x2="{x}" y2="{height - 34}" stroke="{LINE}" stroke-width="1"/>')
        elements.append(_svg_text(x, 143, f"{tick}%", fill=MUTED, font_size=12, text_anchor="middle", font_family="Arial, sans-serif"))

    for index, row in enumerate(rows):
        y = 175 + index * 74
        name = row.get("name") or f"guide_{index + 1}"
        gc_value = float(row["gc_percent"]) / 100
        seed_value = float(row["seed_accessibility"])
        elements.append(_svg_text(45, y + 18, name.replace("_", " "), fill=INK, font_size=15, font_weight=700, font_family="Arial, sans-serif"))
        elements.append(f'<rect x="{plot_left}" y="{y}" width="{plot_width * gc_value:.1f}" height="18" rx="5" fill="{GREEN}"/>')
        elements.append(f'<rect x="{plot_left}" y="{y + 25}" width="{plot_width * seed_value:.1f}" height="18" rx="5" fill="{GOLD}"/>')
        elements.append(_svg_text(plot_left + plot_width + 14, y + 15, f"{gc_value:.0%}", fill=GREEN, font_size=13, font_weight=700, font_family="Arial, sans-serif"))
        elements.append(_svg_text(plot_left + plot_width + 14, y + 40, f"{seed_value:.0%}", fill="#9a6e15", font_size=13, font_weight=700, font_family="Arial, sans-serif"))
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


# ------------------------------------------------------------- study figures


def _header(width: int, height: int, title: str, subtitle: str, description: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        f'<title id="title">{escape(title)}</title>',
        f'<desc id="desc">{escape(description)}</desc>',
        f'<rect width="{width}" height="{height}" fill="{PAPER}" rx="24"/>',
        _svg_text(45, 54, title, fill=INK, font_size=28, font_weight=700, font_family=FONT),
        _svg_text(45, 82, subtitle, fill=MUTED, font_size=15, font_family=FONT),
    ]


def _nice_bounds(values: list[float]) -> tuple[float, float]:
    """A symmetric, rounded axis range that comfortably contains the data."""
    finite = [v for v in values if v == v]
    limit = max(0.05, max(abs(v) for v in finite) if finite else 0.05)
    step = 0.05
    while limit / step > 6:
        step *= 2
    return -math.ceil(limit / step) * step, math.ceil(limit / step) * step


def write_position_figure(path: Path, results: dict, dataset: str) -> None:
    """Spearman correlation between unpaired probability and activity, per position."""
    block = results["datasets"][dataset]
    contexts = [
        ("spacer", "spacer folded alone", GREEN),
        ("full", "spacer folded with the 76 nt scaffold", GOLD),
    ]
    width, height = 1080, 470
    left, right = 90, 1035
    top, bottom = 150, 380
    plot_width = right - left

    values = [
        entry[key]
        for _, entries in block["positions"].items()
        for entry in entries
        for key in ("spearman", "ci_low", "ci_high")
        if entry[key] == entry[key]
    ]
    low_bound, high_bound = _nice_bounds(values)

    def y_of(value: float) -> float:
        return bottom - (value - low_bound) / (high_bound - low_bound) * (bottom - top)

    def x_of(position: int) -> float:
        return left + (position - 0.5) / 20 * plot_width

    elements = _header(
        width, height,
        "Does any spacer position matter?",
        f"{block['label']} — mean within-gene Spearman between unpaired probability and "
        f"editing efficiency, computed separately at every position.",
        "Line chart of per-position Spearman correlation with 95% confidence intervals, "
        "for the spacer folded alone and folded with the sgRNA scaffold.",
    )

    # Canonical seed region, shaded, so the reader can see it was an assumption.
    seed_x = x_of(13) - plot_width / 40
    elements.append(
        f'<rect x="{seed_x:.1f}" y="{top}" width="{right - seed_x:.1f}" '
        f'height="{bottom - top}" fill="{GOLD}" opacity="0.10"/>'
    )
    elements.append(
        _svg_text((seed_x + right) / 2, top - 12, "canonical seed (13-20)", fill=MUTED,
                  font_size=13, text_anchor="middle", font_family=FONT)
    )

    # Axes.
    tick = 0.05 if high_bound <= 0.2 else 0.1
    steps = int(round(high_bound / tick))
    for index in range(-steps, steps + 1):
        value = index * tick
        y = y_of(value)
        colour = INK if index == 0 else LINE
        elements.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" '
            f'stroke="{colour}" stroke-width="{2 if index == 0 else 1}"/>'
        )
        elements.append(_svg_text(left - 12, y + 4, f"{value:+.2f}", fill=MUTED, font_size=12,
                                  text_anchor="end", font_family=FONT))
    for position in range(1, 21):
        elements.append(_svg_text(x_of(position), bottom + 22, str(position), fill=MUTED,
                                  font_size=12, text_anchor="middle", font_family=FONT))
    elements.append(_svg_text((left + right) / 2, bottom + 48,
                              "spacer position (1 = PAM-distal 5' end, 20 = PAM-proximal)",
                              fill=INK, font_size=14, text_anchor="middle", font_family=FONT))
    elements.append(_svg_text(45, 128, "Spearman rho", fill=INK, font_size=13,
                              font_weight=700, font_family=FONT))

    for offset, (context, label, colour) in enumerate(contexts):
        entries = block["positions"][context]
        shift = (offset - 0.5) * 8
        points = []
        for entry in entries:
            x = x_of(entry["position"]) + shift
            y = y_of(entry["spearman"])
            points.append(f"{x:.1f},{y:.1f}")
            if entry["ci_low"] == entry["ci_low"]:
                elements.append(
                    f'<line x1="{x:.1f}" y1="{y_of(entry["ci_low"]):.1f}" x2="{x:.1f}" '
                    f'y2="{y_of(entry["ci_high"]):.1f}" stroke="{colour}" '
                    f'stroke-width="2" opacity="0.45"/>'
                )
        elements.append(
            f'<polyline points="{" ".join(points)}" fill="none" stroke="{colour}" '
            f'stroke-width="2.5" opacity="0.85"/>'
        )
        for entry in entries:
            x = x_of(entry["position"]) + shift
            y = y_of(entry["spearman"])
            significant = entry.get("p_adjusted", 1.0) < 0.05
            elements.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{5 if significant else 3.5}" '
                f'fill="{colour if significant else PAPER}" stroke="{colour}" stroke-width="2"/>'
            )
        legend_y = 108 + offset * 22
        elements.append(f'<rect x="320" y="{legend_y - 10}" width="14" height="14" rx="3" fill="{colour}"/>')
        elements.append(_svg_text(343, legend_y + 2, label, fill=INK, font_size=13, font_family=FONT))

    elements.append(_svg_text(
        45, height - 22,
        "Filled markers are significant after Benjamini-Hochberg correction of gene-blocked "
        "permutation p-values; bars are 95% intervals from resampling whole genes.",
        fill=MUTED, font_size=12, font_family=FONT,
    ))
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def _forest(
    path: Path,
    title: str,
    subtitle: str,
    description: str,
    rows: list[tuple[str, float, float, float]],
    zero_label: str,
) -> None:
    """Estimate-with-interval plot: one row per measure, a vertical line at zero."""
    width = 1080
    row_height = 34
    top = 175
    height = top + row_height * len(rows) + 90
    left, right = 470, 1020

    values = [v for _, estimate, low, high in rows for v in (estimate, low, high) if v == v]
    low_bound, high_bound = _nice_bounds(values)

    def x_of(value: float) -> float:
        return left + (value - low_bound) / (high_bound - low_bound) * (right - left)

    elements = _header(width, height, title, subtitle, description)

    tick = 0.05 if high_bound <= 0.2 else 0.1
    steps = int(round(high_bound / tick))
    for index in range(-steps, steps + 1):
        value = index * tick
        x = x_of(value)
        elements.append(
            f'<line x1="{x:.1f}" y1="{top - 18}" x2="{x:.1f}" y2="{top + row_height * len(rows)}" '
            f'stroke="{INK if index == 0 else LINE}" stroke-width="{2 if index == 0 else 1}"/>'
        )
        elements.append(_svg_text(x, top - 26, f"{value:+.2f}", fill=MUTED, font_size=12,
                                  text_anchor="middle", font_family=FONT))
    elements.append(_svg_text(x_of(0.0), height - 30, zero_label, fill=MUTED, font_size=12,
                              text_anchor="middle", font_family=FONT))

    for index, (label, estimate, low, high) in enumerate(rows):
        y = top + index * row_height + row_height / 2
        if index % 2 == 0:
            elements.append(
                f'<rect x="45" y="{y - row_height / 2:.1f}" width="{right - 45 + 20}" '
                f'height="{row_height}" fill="#ffffff" opacity="0.55"/>'
            )
        elements.append(_svg_text(55, y + 4, label, fill=INK, font_size=13, font_family=FONT))
        crosses_zero = not (low == low) or (low <= 0.0 <= high)
        colour = MUTED if crosses_zero else GREEN
        if low == low:
            elements.append(
                f'<line x1="{x_of(low):.1f}" y1="{y:.1f}" x2="{x_of(high):.1f}" y2="{y:.1f}" '
                f'stroke="{colour}" stroke-width="3" opacity="0.55" stroke-linecap="round"/>'
            )
        elements.append(
            f'<circle cx="{x_of(estimate):.1f}" cy="{y:.1f}" r="6" fill="{colour}"/>'
        )
        elements.append(_svg_text(right + 16, y + 4, f"{estimate:+.3f}", fill=colour,
                                  font_size=12, font_weight=700, font_family=FONT))
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_measures_figure(path: Path, results: dict) -> None:
    rows = []
    for dataset in ("doench", "crisprscan"):
        block = results["datasets"][dataset]
        tag = "human cells" if dataset == "doench" else "zebrafish"
        for entry in block["correlations"]:
            rows.append((
                f"{entry['label']}  ({tag})",
                entry["spearman"], entry["ci_low"], entry["ci_high"],
            ))
    _forest(
        path,
        "How you measure structure barely changes the answer",
        "Spearman correlation with measured editing efficiency. Bars are 95% intervals from "
        "resampling whole genes; grey means the interval includes zero.",
        "Forest plot of Spearman correlations between several accessibility measures and "
        "editing efficiency across two screens.",
        rows,
        "no correlation",
    )


def write_incremental_figure(path: Path, results: dict) -> None:
    rows = []
    for dataset in ("doench", "crisprscan"):
        block = results["datasets"][dataset]
        if "incremental" not in block:
            continue
        tag = "human cells" if dataset == "doench" else "zebrafish"
        for model in block["incremental"]["models"]:
            if model["name"] == "baseline":
                continue
            rows.append((
                f"{model['name'].replace('structure_', '+ ')}  ({tag})",
                model.get("delta", float("nan")),
                model.get("delta_ci_low", float("nan")),
                model.get("delta_ci_high", float("nan")),
            ))
    if not rows:
        return
    baselines = " · ".join(
        f"{'human' if d == 'doench' else 'zebrafish'} baseline rho = "
        f"{results['datasets'][d]['incremental']['baseline_spearman']:.3f}"
        for d in ("doench", "crisprscan")
        if "incremental" in results["datasets"][d]
    )
    _forest(
        path,
        "Does structure add anything to a sequence-only model?",
        f"Change in held-out Spearman when folding features are added to a position + "
        f"dinucleotide + G/C baseline. {baselines}.",
        "Forest plot of the change in held-out Spearman correlation when structure features "
        "are added to a sequence-only baseline.",
        rows,
        "no improvement",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Regenerate the repository's SVG figures.")
    parser.add_argument("--features", type=Path, default=Path("analysis_outputs/crispr_guide_nussinov_features.csv"))
    parser.add_argument("--study", type=Path, default=Path("analysis_outputs/study_results.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("figures"))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_workflow_figure(args.output_dir / "nussinov_workflow.svg")
    write_feature_figure(args.output_dir / "guide_feature_summary.svg", read_features(args.features))

    if args.study.exists():
        results = json.loads(args.study.read_text(encoding="utf-8"))
        write_position_figure(args.output_dir / "position_accessibility.svg", results, "doench")
        write_measures_figure(args.output_dir / "accessibility_measures.svg", results)
        write_incremental_figure(args.output_dir / "incremental_value.svg", results)
        print(f"Wrote study figures from {args.study}")
    else:
        print(f"{args.study} not found; run run_study.py for the study figures.")
    print(f"Wrote SVG figures to {args.output_dir}")


if __name__ == "__main__":
    main()
