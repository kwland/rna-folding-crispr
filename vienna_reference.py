"""Fold every guide with ViennaRNA, for validation against the custom model.

Why this exists
---------------
``energy_model.py`` uses the published Watson-Crick stacking values but
*simplified* wobble, loop-initiation, and multiloop terms. Those simplifications
are documented, but documenting an approximation is not the same as measuring
what it costs. This script recomputes the same quantities with ViennaRNA and the
standard Turner 2004 parameters, so every headline number in the study can be
checked against a reference implementation rather than trusted.

The intended reading of the two implementations is:

* **ViennaRNA is the physical model.** Any quantitative claim about RNA
  thermodynamics should rest on it.
* **The custom implementation is the methodological contribution.** It is
  transparent, dependency-free, verified against exhaustive enumeration, and it
  runs in a browser. Its value is that every step can be inspected, not that it
  is a better energy model.

Dependency note
---------------
ViennaRNA is an optional *validation* dependency. The rest of the project uses
only the standard library, and nothing in the analysis pipeline imports this
module: it writes a CSV that the pipeline reads. On this machine the interpreter
that has ViennaRNA installed may differ from the one that runs the project, so
run it explicitly, for example::

    python vienna_reference.py --dataset doench
    python vienna_reference.py --dataset doench --temperature 25
    python vienna_reference.py --dataset doench --max-span 40   # RNAplfold-style

Install with ``python -m pip install ViennaRNA``.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

try:
    import RNA  # type: ignore
except ImportError:  # pragma: no cover - exercised only without ViennaRNA
    RNA = None

from datasets import load_crisprscan, load_doench_pooled
from guide_features import SGRNA_SCAFFOLD


LOADERS = {"doench": load_doench_pooled, "crisprscan": load_crisprscan}
SEED_LENGTH = 8


def _require_vienna() -> None:
    if RNA is None:
        raise SystemExit(
            "ViennaRNA is not importable from this interpreter.\n"
            f"  interpreter: {sys.executable}\n"
            "  install it with: python -m pip install ViennaRNA\n"
            "This module is optional; the rest of the project does not need it."
        )


def fold_one(sequence: str, temperature: float, max_span: int | None) -> dict:
    """MFE, ensemble free energy, and per-base unpaired probability.

    The unpaired probability of a base is one minus the total probability that
    it pairs with anything, summed over the base-pair probability matrix. That
    is the same definition used in ``mccaskill.py``, so the two are directly
    comparable.
    """
    model = RNA.md()
    model.temperature = temperature
    if max_span:
        model.max_bp_span = max_span

    compound = RNA.fold_compound(sequence, model)
    structure, mfe = compound.mfe()
    # Rescaling by the MFE keeps the partition function inside double precision.
    compound.exp_params_rescale(mfe)
    _, ensemble_energy = compound.pf()

    n = len(sequence)
    unpaired = [1.0] * n
    probabilities = compound.bpp()
    for i in range(1, n + 1):
        row = probabilities[i]
        for j in range(i + 1, n + 1):
            p = row[j]
            if p:
                unpaired[i - 1] -= p
                unpaired[j - 1] -= p
    unpaired = [min(1.0, max(0.0, v)) for v in unpaired]

    return {
        "structure": structure,
        "mfe": mfe,
        "ensemble_energy": ensemble_energy,
        "unpaired": unpaired,
    }


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def profile_guide(
    spacer_dna: str,
    temperature: float,
    max_span: int | None,
    scaffold: str = SGRNA_SCAFFOLD,
) -> dict[str, str]:
    """Fold one spacer alone and inside the full sgRNA, as the study does."""
    spacer = spacer_dna.upper().replace("T", "U")
    n = len(spacer)
    seed_start = max(0, n - SEED_LENGTH)

    alone = fold_one(spacer, temperature, max_span)
    full = fold_one(spacer + scaffold, temperature, max_span)

    def mfe_flags(structure: str) -> list[float]:
        return [1.0 if c == "." else 0.0 for c in structure[:n]]

    row = {
        "vienna_mfe_energy_spacer": f"{alone['mfe']:.4f}",
        "vienna_ensemble_energy_spacer": f"{alone['ensemble_energy']:.4f}",
        "vienna_seed_mfe_spacer": f"{_mean(mfe_flags(alone['structure'])[seed_start:n]):.6f}",
        "vienna_seed_ensemble_spacer": f"{_mean(alone['unpaired'][seed_start:n]):.6f}",
        "vienna_mean_unpaired_spacer": f"{_mean(alone['unpaired'][:n]):.6f}",
        "vienna_mfe_energy_full": f"{full['mfe']:.4f}",
        "vienna_ensemble_energy_full": f"{full['ensemble_energy']:.4f}",
        "vienna_seed_mfe_full": f"{_mean(mfe_flags(full['structure'])[seed_start:n]):.6f}",
        "vienna_seed_ensemble_full": f"{_mean(full['unpaired'][seed_start:n]):.6f}",
        "vienna_mean_unpaired_full": f"{_mean(full['unpaired'][:n]):.6f}",
    }
    for index in range(n):
        row[f"vienna_unpaired_spacer_{index + 1:02d}"] = f"{alone['unpaired'][index]:.6f}"
        row[f"vienna_unpaired_full_{index + 1:02d}"] = f"{full['unpaired'][index]:.6f}"
    return row


def default_output(dataset: str, temperature: float, max_span: int | None) -> Path:
    stem = f"vienna_{dataset}"
    if abs(temperature - 37.0) > 1e-9:
        stem += f"_t{temperature:g}"
    if max_span:
        stem += f"_span{max_span}"
    return Path("analysis_outputs") / f"{stem}.csv"


def main() -> None:
    parser = argparse.ArgumentParser(description="Fold a screen with ViennaRNA.")
    parser.add_argument("--dataset", choices=sorted(LOADERS), default="doench")
    parser.add_argument("--temperature", type=float, default=37.0,
                        help="Folding temperature in Celsius (ViennaRNA default is 37).")
    parser.add_argument("--max-span", type=int, default=None,
                        help="Maximum base-pair span, the RNAplfold-style local folding limit.")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    _require_vienna()
    guides = LOADERS[args.dataset]()
    if args.limit:
        guides = guides[: args.limit]
    output = args.output or default_output(args.dataset, args.temperature, args.max_span)

    print(f"ViennaRNA {RNA.__version__} | {args.dataset}: {len(guides)} guides "
          f"| {args.temperature} C" + (f" | max span {args.max_span}" if args.max_span else ""))

    started = time.perf_counter()
    rows = []
    for index, guide in enumerate(guides, start=1):
        record = profile_guide(guide.spacer, args.temperature, args.max_span)
        rows.append({"guide_id": guide.guide_id, **record})
        if index % 250 == 0 or index == len(guides):
            rate = index / (time.perf_counter() - started)
            sys.stderr.write(f"\r  {index}/{len(guides)} ({rate:.0f}/s)   ")
            sys.stderr.flush()
    sys.stderr.write("\n")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output} ({len(rows)} guides, {time.perf_counter() - started:.1f}s)")


if __name__ == "__main__":
    main()
