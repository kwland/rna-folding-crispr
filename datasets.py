"""Loaders for the guide-activity screens used in the study.

Two screens are used, deliberately chosen to be as unlike each other as
possible while still measuring SpCas9 on-target activity:

**Doench pooled** (4,685 guides, 18 genes, human cell culture). Doench et al.
2014 and 2016, pooled through CRISPOR / Haeussler et al. 2016. Activity is a
percentile within its own source screen. Lives in ``docs/data/guides.json``
because the browser demonstration reads the same file.

**CRISPRscan** (1,020 guides, 111 genes, zebrafish embryos, in vivo). Moreno-
Mateos et al. 2015, Nature Methods 12:982. Activity is the measured mutation
frequency. Retrieved from the CRISPOR paper's public dataset collection
(``maximilianh/crisporPaper``, ``effData/morenoMateos2015.context.tab``) and
reduced to guide, gene, spacer, PAM, and activity in
``data/crisprscan_moreno_mateos_2015.csv``.

Why two screens matter: a result that holds only in the screen it was found in
is not a result. Different organism, different delivery, different readout, and
different laboratory make CRISPRscan a real out-of-distribution test rather than
another slice of the same experiment.

One caveat worth stating in any write-up: CRISPRscan guides are transcribed
in vitro from a T7/SP6 promoter, so **every** spacer begins with GG. Positions 1
and 2 therefore carry no information in that screen, and any model that leans on
them will transfer badly. That is a property of the biology of the assay, not a
data-cleaning error.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path


__all__ = [
    "Guide",
    "load_doench_pooled",
    "load_crisprscan",
    "percentile_ranks",
    "assign_percentiles",
]


@dataclass(frozen=True)
class Guide:
    """One measured guide."""

    guide_id: str
    group: str          # gene or locus; the unit that cross-validation holds out
    spacer: str         # 20 nt, DNA alphabet
    activity: float     # raw measured value, whatever the screen reported
    screen: str         # the individual screen, e.g. "Doench2016"
    dataset: str        # the collection, e.g. "doench_pooled"
    percentile: float = 0.0  # activity rank within its own screen, filled in later


def percentile_ranks(values: list[float]) -> list[float]:
    """Rank values into [0, 1], averaging ranks within ties.

    Screens report activity on incompatible scales - a survival log-fold-change
    and an in vivo mutation frequency are not the same quantity. Converting each
    screen to its own percentile puts them on one axis without pretending the
    underlying measurements are interchangeable. It is also monotonic, so every
    Spearman correlation in the study is unchanged by it.
    """
    count = len(values)
    if count == 0:
        return []
    if count == 1:
        return [0.5]
    order = sorted(range(count), key=lambda i: values[i])
    ranks = [0.0] * count
    position = 0
    while position < count:
        end = position
        while end + 1 < count and values[order[end + 1]] == values[order[position]]:
            end += 1
        shared = (position + end) / 2.0
        for index in range(position, end + 1):
            ranks[order[index]] = shared / (count - 1)
        position = end + 1
    return ranks


def assign_percentiles(guides: list[Guide]) -> list[Guide]:
    """Return the guides with ``percentile`` filled in, computed per screen."""
    by_screen: dict[str, list[int]] = {}
    for index, guide in enumerate(guides):
        by_screen.setdefault(guide.screen, []).append(index)

    output = list(guides)
    for indices in by_screen.values():
        ranks = percentile_ranks([guides[i].activity for i in indices])
        for index, rank in zip(indices, ranks):
            output[index] = Guide(**{**vars(guides[index]), "percentile": rank})
    return output


def _validate(spacer: str, guide_id: str) -> str:
    spacer = spacer.strip().upper()
    if len(spacer) != 20:
        raise ValueError(f"Guide {guide_id}: expected a 20 nt spacer, got {len(spacer)}")
    unexpected = sorted(set(spacer) - set("ACGT"))
    if unexpected:
        raise ValueError(f"Guide {guide_id}: unexpected bases {', '.join(unexpected)}")
    return spacer


def load_doench_pooled(path: Path | str = Path("docs/data/guides.json")) -> list[Guide]:
    """Load the pooled Doench 2014 + 2016 guides.

    The screen label is taken from the guide id prefix, which is what makes the
    2014-versus-2016 hold-out possible: they are separate experiments, so
    training on one and testing on the other is a real transfer test rather than
    a random split of one experiment.
    """
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    guides = []
    for record in payload["guides"]:
        guide_id = record["id"]
        screen = guide_id.split("-")[0]
        guides.append(
            Guide(
                guide_id=guide_id,
                group=record["gene"],
                spacer=_validate(record["spacer"], guide_id),
                activity=float(record["activity"]),
                screen=screen,
                dataset="doench_pooled",
            )
        )
    if not guides:
        raise ValueError(f"No guides found in {path}")
    return assign_percentiles(guides)


def load_crisprscan(
    path: Path | str = Path("data/crisprscan_moreno_mateos_2015.csv"),
) -> list[Guide]:
    """Load the CRISPRscan zebrafish guides."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. See the module docstring for its provenance."
        )
    guides = []
    with path.open(newline="", encoding="utf-8") as handle:
        for record in csv.DictReader(handle):
            guide_id = record["guide_id"]
            guides.append(
                Guide(
                    guide_id=guide_id,
                    group=record["gene"],
                    spacer=_validate(record["spacer_dna"], guide_id),
                    activity=float(record["mod_freq"]),
                    screen="MorenoMateos2015",
                    dataset="crisprscan",
                )
            )
    if not guides:
        raise ValueError(f"No guides found in {path}")
    return assign_percentiles(guides)
