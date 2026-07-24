"""Structure features for CRISPR guide spacers, computed four different ways.

The point of this module is that "seed accessibility" is not one number. It
depends on two choices that earlier passes of this project made silently:

1. **What counts as the structure.** One optimal structure (Nussinov or Zuker)
   gives every position a hard 0 or 1. The Boltzmann ensemble (McCaskill) gives
   a probability. A base that is unpaired in the MFE structure but paired in 45%
   of the ensemble is not really "accessible", and the single-structure view
   cannot say so.

2. **What molecule is folded.** A spacer does not float around on its own. In a
   real sgRNA it is followed by a 76 nt scaffold that it can pair with, so the
   spacer folded alone and the spacer folded in context are different molecules
   and can give different answers.

Every feature below is computed under all four combinations, so the study can
ask which choice - if either - changes the conclusion.

Position convention: spacer position 1 is the PAM-distal 5' end and position 20
is the PAM-proximal 3' end, the end that sits next to the scaffold. The seed is
the last ``seed_length`` positions, 13-20 by default.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from mccaskill import partition_fold
from nussinov import nussinov_fold
from zuker import zuker_fold


__all__ = ["SGRNA_SCAFFOLD", "GuideProfile", "profile_guide", "gc_fraction"]


#: Standard SpCas9 sgRNA constant region (76 nt) that follows the 20 nt spacer.
SGRNA_SCAFFOLD = (
    "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC"
)

DEFAULT_SEED_LENGTH = 8


@dataclass
class GuideProfile:
    """Structure features for one spacer under all four measurement choices."""

    spacer: str
    length: int
    gc_percent: float

    # --- accessibility of the seed, by (structure view, folded molecule)
    seed_nussinov: float          # single structure, pair-count model, spacer alone
    seed_mfe_spacer: float        # single structure, energy model, spacer alone
    seed_ensemble_spacer: float   # Boltzmann ensemble, spacer alone
    seed_mfe_full: float          # single structure, energy model, spacer + scaffold
    seed_ensemble_full: float     # Boltzmann ensemble, spacer + scaffold

    # --- whole-spacer summaries
    mean_unpaired_spacer: float
    mean_unpaired_full: float
    mfe_energy_spacer: float
    mfe_energy_full: float
    ensemble_energy_spacer: float
    ensemble_energy_full: float
    nussinov_pairs: int

    # --- per-position unpaired probability across the 20 spacer positions
    unpaired_spacer: list[float] = field(default_factory=list)
    unpaired_full: list[float] = field(default_factory=list)

    def as_row(self) -> dict[str, str]:
        """Flatten to CSV-ready strings, including the position-resolved columns."""
        row: dict[str, str] = {
            "spacer_rna": self.spacer,
            "length": str(self.length),
            "gc_percent": f"{self.gc_percent:.1f}",
            "nussinov_pairs": str(self.nussinov_pairs),
            "seed_nussinov": f"{self.seed_nussinov:.6f}",
            "seed_mfe_spacer": f"{self.seed_mfe_spacer:.6f}",
            "seed_ensemble_spacer": f"{self.seed_ensemble_spacer:.6f}",
            "seed_mfe_full": f"{self.seed_mfe_full:.6f}",
            "seed_ensemble_full": f"{self.seed_ensemble_full:.6f}",
            "mean_unpaired_spacer": f"{self.mean_unpaired_spacer:.6f}",
            "mean_unpaired_full": f"{self.mean_unpaired_full:.6f}",
            "mfe_energy_spacer": f"{self.mfe_energy_spacer:.4f}",
            "mfe_energy_full": f"{self.mfe_energy_full:.4f}",
            "ensemble_energy_spacer": f"{self.ensemble_energy_spacer:.4f}",
            "ensemble_energy_full": f"{self.ensemble_energy_full:.4f}",
        }
        for index, value in enumerate(self.unpaired_spacer, start=1):
            row[f"unpaired_spacer_{index:02d}"] = f"{value:.6f}"
        for index, value in enumerate(self.unpaired_full, start=1):
            row[f"unpaired_full_{index:02d}"] = f"{value:.6f}"
        return row


def gc_fraction(sequence: str) -> float:
    seq = sequence.upper().replace("U", "T")
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def profile_guide(
    spacer: str,
    seed_length: int = DEFAULT_SEED_LENGTH,
    scaffold: str = SGRNA_SCAFFOLD,
    max_span: int | None = None,
) -> GuideProfile:
    """Fold one spacer alone and in sgRNA context, and summarise both.

    ``max_span`` is passed through to the partition function; setting it models
    local, co-transcriptional folding in the style of RNAplfold.
    """
    if seed_length <= 0:
        raise ValueError("seed_length must be a positive integer")

    nussinov = nussinov_fold(spacer)
    rna = nussinov.sequence
    n = len(rna)
    if n == 0:
        raise ValueError("Spacer is empty")
    seed_start = max(0, n - seed_length)

    # --- spacer folded on its own
    mfe_spacer = zuker_fold(rna)
    ens_spacer = partition_fold(rna, max_span=max_span)

    nussinov_paired = {index for pair in nussinov.pairs for index in pair}
    seed_positions = range(seed_start, n)
    seed_nussinov = _mean([0.0 if p in nussinov_paired else 1.0 for p in seed_positions])

    mfe_flags = mfe_spacer.unpaired_flags()

    # --- spacer folded as part of the full sgRNA
    full = rna + scaffold
    mfe_full = zuker_fold(full)
    ens_full = partition_fold(full, max_span=max_span)
    mfe_full_flags = mfe_full.unpaired_flags()

    return GuideProfile(
        spacer=rna,
        length=n,
        gc_percent=gc_fraction(rna) * 100,
        seed_nussinov=seed_nussinov,
        seed_mfe_spacer=_mean(mfe_flags[seed_start:n]),
        seed_ensemble_spacer=ens_spacer.mean_unpaired(seed_start, n),
        seed_mfe_full=_mean(mfe_full_flags[seed_start:n]),
        seed_ensemble_full=ens_full.mean_unpaired(seed_start, n),
        mean_unpaired_spacer=ens_spacer.mean_unpaired(0, n),
        mean_unpaired_full=ens_full.mean_unpaired(0, n),
        mfe_energy_spacer=mfe_spacer.energy,
        mfe_energy_full=mfe_full.energy,
        ensemble_energy_spacer=ens_spacer.ensemble_free_energy,
        ensemble_energy_full=ens_full.ensemble_free_energy,
        nussinov_pairs=nussinov.score,
        unpaired_spacer=list(ens_spacer.unpaired[:n]),
        unpaired_full=list(ens_full.unpaired[:n]),
    )


# The list of scalar structure features offered to the predictive models lives in
# model_features.STRUCTURE_FEATURES, not here, so there is only one place for it
# to be edited and no way for two copies to drift apart.
