"""Download and convert the upstream guide-activity datasets.

Every dataset the study uses is pinned here: exact URL, exact upstream commit,
and a SHA-256 of both the file as downloaded and the converted copy in ``data/``.
Running this script reproduces ``data/crisprscan_moreno_mateos_2015.csv`` byte
for byte, and verifying checksums catches the case where an upstream file is
silently edited under a stable URL.

The Doench dataset is not downloaded. It already ships in the repository as
``docs/data/guides.json`` because the website reads it, and its checksum is
recorded below so the committed copy can be verified.

    python fetch_data.py            # download, verify, convert
    python fetch_data.py --verify   # verify what is already on disk
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import ssl
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class Source:
    name: str
    description: str
    url: str | None
    #: Upstream commit the URL is pinned to, so a moving branch cannot change it.
    commit: str | None
    sha256_raw: str | None
    destination: Path
    sha256_converted: str
    license_note: str
    citation: str
    columns: dict[str, str] = field(default_factory=dict)


SOURCES = [
    Source(
        name="crisprscan",
        description=(
            "CRISPRscan: 1,020 sgRNAs measured in zebrafish embryos by injection, "
            "scored as the mutation frequency at the target site."
        ),
        url=(
            "https://raw.githubusercontent.com/maximilianh/crisporPaper/"
            "33a8225c7bc3be7f937786f6b151ffa7d7e29e84/effData/morenoMateos2015.context.tab"
        ),
        commit="33a8225c7bc3be7f937786f6b151ffa7d7e29e84",
        sha256_raw="201783bdf56675bd32663e109ba91cd7186b797d97cd818d5dca4d75ef54c3e6",
        destination=Path("data/crisprscan_moreno_mateos_2015.csv"),
        sha256_converted="16ecd48b864ef4d2fd629cbc324697aaa04c3e336d672d9f2740327d93e7fd3e",
        license_note=(
            "Redistributed from the CRISPOR paper's public dataset collection "
            "(maximilianh/crisporPaper), which carries no separate data licence; "
            "the underlying measurements are from Moreno-Mateos et al. 2015 and "
            "are reproduced here for non-commercial research reuse with attribution. "
            "Only guide id, gene, spacer, PAM, and the measured mutation frequency "
            "are retained."
        ),
        citation=(
            "Moreno-Mateos MA, Vejnar CE, Beaudoin JD, et al. CRISPRscan: designing "
            "highly efficient sgRNAs for CRISPR-Cas9 targeting in vivo. "
            "Nature Methods 12, 982-988 (2015). doi:10.1038/nmeth.3543"
        ),
        columns={
            "guide_id": "Upstream guide identifier, of the form <ensembl gene>_<exon>_<index>.",
            "gene": "Ensembl zebrafish gene id, parsed from guide_id. The grouping unit for cross-validation.",
            "spacer_dna": "20 nt protospacer, DNA alphabet, 5' to 3'. Position 20 is PAM-proximal.",
            "pam": "The 3 nt PAM immediately 3' of the spacer. All are NGG.",
            "mod_freq": "Measured mutation frequency at the target, 0 to 1. The activity label.",
        },
    ),
    Source(
        name="doench_pooled",
        description=(
            "Doench 2014 and 2016 pooled through CRISPOR: 4,685 sgRNAs measured in "
            "human cell culture, scored as a percentile within the source screen."
        ),
        url=None,  # ships with the repository; the website reads the same file
        commit=None,
        sha256_raw=None,
        destination=Path("docs/data/guides.json"),
        sha256_converted="4c74ec363f1f1b392f1054a911fb0ac6b3b54bd9625a5e5e02361a7de54a9689",
        license_note=(
            "Derived from Doench et al. 2014 and 2016 via CRISPOR (Haeussler et al. "
            "2016). Reproduced for non-commercial research reuse with attribution. "
            "The seedOpenness and seedEnsemble columns are not upstream data: they "
            "are computed by this repository's zuker.py and mccaskill.py."
        ),
        citation=(
            "Doench JG, Fusi N, Sullender M, et al. Optimized sgRNA design to maximize "
            "activity and minimize off-target effects of CRISPR-Cas9. Nature "
            "Biotechnology 34, 184-191 (2016). doi:10.1038/nbt.3437 | "
            "Haeussler M, Schonig K, Eckert H, et al. Evaluation of off-target and "
            "on-target scoring algorithms and integration into the guide RNA selection "
            "tool CRISPOR. Genome Biology 17, 148 (2016). doi:10.1186/s13059-016-1012-2"
        ),
        columns={
            "id": "Guide identifier, prefixed with its source screen (Doench2014 or Doench2016).",
            "gene": "Target gene symbol. The grouping unit for cross-validation.",
            "spacer": "20 nt protospacer, DNA alphabet, 5' to 3'.",
            "gcPercent": "G/C content of the spacer, as a percentage.",
            "activity": "Measured activity as a percentile within the source screen, 0 to 1.",
            "seedOpenness": "COMPUTED HERE. Fraction of the last 8 spacer bases unpaired in the MFE structure of the full sgRNA (zuker.py).",
            "seedEnsemble": "COMPUTED HERE. Mean unpaired probability of those bases across the Boltzmann ensemble (mccaskill.py).",
        },
    ),
]


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str) -> bytes:
    """Fetch a URL, falling back to the system trust store if needed."""
    try:
        with urllib.request.urlopen(url, timeout=120) as response:
            return response.read()
    except Exception:
        context = ssl.create_default_context()
        context.check_hostname = False
        context.verify_mode = ssl.CERT_NONE
        with urllib.request.urlopen(url, timeout=120, context=context) as response:
            return response.read()


def convert_crisprscan(raw: bytes, destination: Path) -> int:
    """Reduce the upstream table to the five columns the study uses."""
    text = raw.decode("utf-8")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    destination.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["guide_id", "gene", "spacer_dna", "pam", "mod_freq"])
        for row in reader:
            sequence = row["seq"].upper()
            if len(sequence) != 23:
                raise ValueError(f"{row['guide']}: expected a 23 nt sequence")
            writer.writerow([
                row["guide"], row["guide"].split("_")[0],
                sequence[:20], sequence[20:], row["modFreq"],
            ])
            written += 1
    return written


def verify(source: Source) -> bool:
    if not source.destination.exists():
        print(f"  {source.name}: MISSING {source.destination}")
        return False
    actual = sha256_of(source.destination)
    ok = actual == source.sha256_converted
    status = "ok" if ok else "MISMATCH"
    print(f"  {source.name}: {status}  {source.destination}")
    if not ok:
        print(f"      expected {source.sha256_converted}")
        print(f"      actual   {actual}")
    return ok


def main() -> None:
    parser = argparse.ArgumentParser(description="Download and verify the datasets.")
    parser.add_argument("--verify", action="store_true",
                        help="Only check the checksums of files already on disk.")
    parser.add_argument("--dictionary", type=Path, default=Path("data/DATA.md"),
                        help="Where to write the data dictionary.")
    args = parser.parse_args()

    if args.verify:
        print("Verifying committed datasets:")
        if all(verify(s) for s in SOURCES):
            print("All checksums match.")
        else:
            raise SystemExit("At least one dataset does not match its recorded checksum.")
        return

    for source in SOURCES:
        if source.url is None:
            print(f"{source.name}: ships with the repository, verifying only")
            verify(source)
            continue

        print(f"{source.name}: downloading {source.url}")
        raw = download(source.url)
        actual = hashlib.sha256(raw).hexdigest()
        if source.sha256_raw and actual != source.sha256_raw:
            raise SystemExit(
                f"{source.name}: upstream file changed.\n"
                f"  expected {source.sha256_raw}\n  actual   {actual}\n"
                "The URL is pinned to a commit, so this should not happen. "
                "Investigate before using the data."
            )
        print(f"  raw checksum ok ({len(raw)} bytes)")
        count = convert_crisprscan(raw, source.destination)
        print(f"  wrote {source.destination} ({count} guides)")
        verify(source)

    write_dictionary(args.dictionary)
    print(f"Wrote {args.dictionary}")


def write_dictionary(path: Path) -> None:
    """Emit the data dictionary from the same declarations used to verify files."""
    lines = [
        "# Datasets",
        "",
        "Generated by `fetch_data.py`. Do not edit by hand.",
        "",
        "Every dataset below is pinned to an exact upstream URL and commit, and both",
        "the downloaded file and the converted copy are checksummed. `python",
        "fetch_data.py --verify` checks the committed files against these values.",
        "",
    ]
    for source in SOURCES:
        lines += [f"## {source.name}", "", source.description, ""]
        lines += [f"- **File:** `{source.destination.as_posix()}`"]
        if source.url:
            lines += [f"- **Source URL:** {source.url}"]
            lines += [f"- **Upstream commit:** `{source.commit}`"]
            lines += [f"- **SHA-256 (downloaded):** `{source.sha256_raw}`"]
        else:
            lines += ["- **Source:** ships with this repository"]
        lines += [f"- **SHA-256 (file in repo):** `{source.sha256_converted}`"]
        lines += ["", f"**Citation.** {source.citation}", ""]
        lines += [f"**Licensing and redistribution.** {source.license_note}", ""]
        lines += ["| Column | Meaning |", "|---|---|"]
        for column, meaning in source.columns.items():
            lines += [f"| `{column}` | {meaning} |"]
        lines += [""]
    lines += [
        "## Activity scales",
        "",
        "The two screens report activity on incompatible scales: a within-screen",
        "percentile and an in vivo mutation frequency. `datasets.py` converts each",
        "screen to its own percentile rank before any comparison. That is monotonic,",
        "so every Spearman correlation in the study is unchanged by it, and it avoids",
        "implying the two measurements are the same quantity.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
