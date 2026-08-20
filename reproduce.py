"""Regenerate every table, figure, and data file in the repository.

One command, ordered so each step's inputs exist before it runs:

    python reproduce.py

Steps
-----
1. verify datasets     checksums of the committed data against ``fetch_data.py``
2. tests               the full unit suite, including brute-force folding checks
3. fold both screens   cached per-guide features (the slow step, about 8 minutes)
4. run the study       the four experiments, tables, and ``study_results.json``
5. export site data    the compact file the website reads
6. train site model    ridge weights for the browser
7. figures             every SVG in ``figures/``
8. ViennaRNA reference optional, needs ViennaRNA installed
9. ViennaRNA validation optional, needs step 8

Steps 8 and 9 are skipped with a clear message if ViennaRNA is not importable,
so the pipeline still completes on a machine without it. Use ``--vienna-python``
when ViennaRNA lives in a different interpreter than the one running this script,
which is common on Windows.

    python reproduce.py --skip-folding      # reuse cached features
    python reproduce.py --quick             # small bootstrap counts, for a smoke test
    python reproduce.py --list              # show the steps without running them
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path


def vienna_available(interpreter: str) -> bool:
    result = subprocess.run(
        [interpreter, "-c", "import RNA"],
        capture_output=True, text=True,
    )
    return result.returncode == 0


def build_steps(args) -> list[tuple[str, list[str], bool]]:
    """Return (label, command, optional) for each step, in order."""
    python = sys.executable
    vienna_python = args.vienna_python or python
    boot = ["--bootstrap", "200", "--permutations", "200"] if args.quick else []

    steps: list[tuple[str, list[str], bool]] = [
        ("Verify datasets", [python, "fetch_data.py", "--verify"], False),
        ("Unit tests", [python, "-m", "unittest", "discover", "-q"], False),
    ]

    if not args.skip_folding:
        for dataset in ("doench", "crisprscan"):
            steps.append((
                f"Fold {dataset}",
                [python, "compute_features.py", "--dataset", dataset,
                 "--workers", str(args.workers)],
                False,
            ))

    steps += [
        ("Run the study", [python, "run_study.py", "--workers", str(args.workers), *boot], False),
        ("Export site data", [python, "export_site_data.py"], False),
        ("Train the site model", [python, "export_model.py"], False),
        ("Figures", [python, "make_figures.py"], False),
    ]

    if args.skip_vienna:
        return steps

    if not vienna_available(vienna_python):
        print(
            f"note: ViennaRNA is not importable from {vienna_python}, so the "
            "validation steps will be skipped.\n"
            "      install it with: python -m pip install ViennaRNA\n"
            "      or point --vienna-python at an interpreter that has it.\n"
        )
        return steps

    if not args.skip_folding:
        for dataset in ("doench", "crisprscan"):
            steps.append((f"ViennaRNA reference: {dataset}",
                          [vienna_python, "vienna_reference.py", "--dataset", dataset], True))
            for temperature in ("25", "42"):
                steps.append((f"ViennaRNA {temperature} C: {dataset}",
                              [vienna_python, "vienna_reference.py", "--dataset", dataset,
                               "--temperature", temperature], True))
            steps.append((f"ViennaRNA local folding: {dataset}",
                          [vienna_python, "vienna_reference.py", "--dataset", dataset,
                           "--max-span", "40"], True))

    steps.append(("ViennaRNA validation",
                  [python, "validate_vienna.py", *(["--bootstrap", "200"] if args.quick else [])],
                  True))
    return steps


def main() -> None:
    parser = argparse.ArgumentParser(description="Regenerate every result in the repository.")
    parser.add_argument("--skip-folding", action="store_true",
                        help="Reuse the cached per-guide features instead of re-folding.")
    parser.add_argument("--skip-vienna", action="store_true",
                        help="Skip the ViennaRNA validation even if it is installed.")
    parser.add_argument("--vienna-python", default=None,
                        help="Interpreter that has ViennaRNA, if not this one.")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--quick", action="store_true",
                        help="Small resampling counts: a smoke test, not a result.")
    parser.add_argument("--list", action="store_true", help="Print the steps and exit.")
    args = parser.parse_args()

    steps = build_steps(args)

    if args.list:
        for index, (label, command, optional) in enumerate(steps, start=1):
            tag = " (optional)" if optional else ""
            print(f"{index:2d}. {label}{tag}\n    {' '.join(command)}")
        return

    if args.quick:
        print("QUICK MODE: resampling counts are reduced. The numbers are not "
              "publication quality.\n")

    started = time.perf_counter()
    failures = []
    for index, (label, command, optional) in enumerate(steps, start=1):
        print(f"\n{'=' * 70}\n[{index}/{len(steps)}] {label}\n{'=' * 70}", flush=True)
        step_started = time.perf_counter()
        result = subprocess.run(command)
        elapsed = time.perf_counter() - step_started
        if result.returncode != 0:
            if optional:
                print(f"-- optional step failed after {elapsed:.1f}s, continuing")
                failures.append((label, True))
                continue
            raise SystemExit(f"\nStep failed: {label}\n  {' '.join(command)}")
        print(f"-- {label} finished in {elapsed:.1f}s")

    total = time.perf_counter() - started
    print(f"\n{'=' * 70}")
    print(f"Finished {len(steps)} steps in {total / 60:.1f} minutes.")
    if failures:
        print("Optional steps that did not complete:")
        for label, _ in failures:
            print(f"  - {label}")
    print("\nOutputs:")
    for path in sorted(Path("analysis_outputs").glob("*")) + sorted(Path("figures").glob("*")):
        print(f"  {path.as_posix()}")


if __name__ == "__main__":
    main()
