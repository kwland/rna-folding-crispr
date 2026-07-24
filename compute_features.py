"""Fold every guide in a screen and cache the structure features to CSV.

Folding is the expensive part of the study: the partition function is O(n^3) and
the spacer-plus-scaffold molecule is 96 nt, so a full screen takes minutes rather
than seconds. Results are written once and re-read by ``run_study.py``, so the
analysis can be re-run and re-plotted without re-folding anything.

    python compute_features.py --dataset doench
    python compute_features.py --dataset crisprscan
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from datasets import Guide, load_crisprscan, load_doench_pooled
from guide_features import SGRNA_SCAFFOLD, profile_guide


LOADERS = {
    "doench": load_doench_pooled,
    "crisprscan": load_crisprscan,
}

DEFAULT_OUTPUT = {
    "doench": Path("analysis_outputs/features_doench_pooled.csv"),
    "crisprscan": Path("analysis_outputs/features_crisprscan.csv"),
}


def _worker(args: tuple[str, int, str, int | None]) -> dict[str, str]:
    """Fold one spacer. Top-level so it can be sent to a worker process."""
    spacer, seed_length, scaffold, max_span = args
    profile = profile_guide(
        spacer, seed_length=seed_length, scaffold=scaffold, max_span=max_span
    )
    return profile.as_row()


def compute(
    guides: list[Guide],
    output: Path,
    seed_length: int = 8,
    max_span: int | None = None,
    workers: int | None = None,
) -> None:
    workers = workers or max(1, (os.cpu_count() or 2) - 1)
    payload = [(g.spacer, seed_length, SGRNA_SCAFFOLD, max_span) for g in guides]

    started = time.perf_counter()
    rows: list[dict[str, str]] = []
    if workers == 1:
        for index, item in enumerate(payload):
            rows.append(_worker(item))
            _progress(index + 1, len(payload), started)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for index, row in enumerate(pool.map(_worker, payload, chunksize=8)):
                rows.append(row)
                _progress(index + 1, len(payload), started)
    sys.stderr.write("\n")

    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "guide_id", "gene", "screen", "dataset", "activity", "percentile",
        *rows[0].keys(),
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for guide, row in zip(guides, rows):
            writer.writerow({
                "guide_id": guide.guide_id,
                "gene": guide.group,
                "screen": guide.screen,
                "dataset": guide.dataset,
                "activity": f"{guide.activity:.6f}",
                "percentile": f"{guide.percentile:.6f}",
                **row,
            })
    elapsed = time.perf_counter() - started
    print(f"Wrote {output} ({len(rows)} guides, {elapsed / 60:.1f} min, {workers} workers)")


def _progress(done: int, total: int, started: float) -> None:
    if done % 50 and done != total:
        return
    elapsed = time.perf_counter() - started
    rate = done / elapsed if elapsed else 0.0
    remaining = (total - done) / rate if rate else 0.0
    sys.stderr.write(
        f"\r  folded {done}/{total} guides  "
        f"({rate:.1f}/s, ~{remaining / 60:.1f} min left)   "
    )
    sys.stderr.flush()


def main() -> None:
    parser = argparse.ArgumentParser(description="Cache folding features for a guide screen.")
    parser.add_argument("--dataset", choices=sorted(LOADERS), default="doench")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--seed-length", type=int, default=8)
    parser.add_argument(
        "--max-span", type=int, default=None,
        help="Forbid pairs spanning more than this many bases (RNAplfold-style local folding).",
    )
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--limit", type=int, default=None, help="Fold only the first N guides.")
    args = parser.parse_args()

    guides = LOADERS[args.dataset]()
    if args.limit:
        guides = guides[: args.limit]
    output = args.output or DEFAULT_OUTPUT[args.dataset]
    print(f"{args.dataset}: {len(guides)} guides, seed = last {args.seed_length} nt")
    compute(
        guides,
        output,
        seed_length=args.seed_length,
        max_span=args.max_span,
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
