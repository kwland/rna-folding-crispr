# RNA folding for CRISPR guide exploration

[![tests](https://github.com/kwland/nussinov-zuker-crispr/actions/workflows/tests.yml/badge.svg)](https://github.com/kwland/nussinov-zuker-crispr/actions/workflows/tests.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-2f7d62)](https://www.python.org/)
[![MIT License](https://img.shields.io/badge/license-MIT-d6a84a)](LICENSE)

**Live interactive demo:** <https://kwland.github.io/nussinov-zuker-crispr/>

This project asks a simple biological question: if a CRISPR guide RNA folds back on itself, does that make its targeting sequence less available to bind DNA?

The short answer, tested four different ways across two independent screens, is **no** — and the interesting part is how much work it took to be confident about that.

![Four-step Nussinov workflow](figures/nussinov_workflow.svg)

## The problem with the first answer

An early version of this project folded each 20 nt spacer on its own, took the single best structure, counted how many of the last 8 bases were unpaired, and correlated that with editing efficiency. The correlation was about 0.01 — nothing.

That result had three obvious holes, and a reviewer would find all of them:

1. **One structure is not the molecule.** RNA does not sit in a single conformation. Reading accessibility off one optimal structure gives every base a hard 0 or 1 that can flip on a tenth of a kcal/mol.
2. **A bare spacer is not the molecule either.** In a real sgRNA the spacer is followed by a 76 nt scaffold it can pair with.
3. **"The last 8 bases" was an assumption**, not something the data was asked about.

And even if all three were fixed, correlating one feature against activity is the wrong question. The question a guide designer cares about is whether structure adds anything to what sequence already tells you.

This repository closes all four gaps.

## What was added

### 1. Ensemble accessibility, not single-structure accessibility

`mccaskill.py` implements the McCaskill (1990) partition function: the same O(n³) dynamic-programming family as Zuker, but with Boltzmann weights `exp(-E/RT)` instead of a minimum, plus the outside recursion that turns those into base-pair probabilities `p_ij`. Seed accessibility becomes the mean probability that each seed base is unpaired across the entire ensemble, which is what RNAplfold-style accessibility actually means.

The multiloop term in the outside pass is the awkward part — written naively it is O(n⁴). The implementation factorises the enclosing-pair sum into two accumulator tables so the whole outside pass stays O(n³).

It is verified, not just written. For sequences short enough to enumerate every pseudoknot-free structure, the partition function, every pair probability, and every unpaired probability match brute force to ~1e-15. Because the shipped parameters make multiloops carry under 0.1% of the ensemble weight, the tests also re-run the comparison with multiloops made artificially cheap, pushing them past 50% of the weight — the recursions still match exactly. Stochastic Boltzmann sampling, which uses only the inside matrices, converges to the probabilities the outside pass computes independently.

The change of instrument is real: across 4,685 guides, single-structure seed accessibility takes only **9 distinct values** (it can only be a multiple of 1/8), while the ensemble measure takes **4,668**.

![Accessibility measures versus activity](figures/accessibility_measures.svg)

### 2. Folding in real molecular context

`guide_features.py` folds every spacer twice: alone, and joined to the 76 nt sgRNA scaffold. Both are reported side by side for all 5,705 guides.

This is not a detail. Mean seed accessibility falls from 0.84 folded alone to 0.42 folded in the sgRNA, and the two views correlate at only **r = 0.17** — folding a spacer on its own tells you almost nothing about the molecule that actually exists. It also exposes a problem with the original hypothesis: a bare 20 nt spacer has an ensemble free energy of just −0.78 kcal/mol, and **22% of them are essentially unstructured**. There was barely any structure there to correlate with anything.

### 3. Position-resolved analysis

Per-base unpaired probabilities make it possible to correlate accessibility with efficiency at *every* spacer position rather than collapsing positions 13–20 into one number, with Benjamini–Hochberg correction across the 20 tests.

![Position-resolved accessibility](figures/position_accessibility.svg)

### 4. Incremental value, cross-validated and replicated

`run_study.py` builds a sequence-only baseline (position-specific bases, position-specific dinucleotides, G/C), evaluates it with **nested cross-validation that holds out whole genes**, then measures how much held-out Spearman changes when folding features are added. Confidence intervals come from resampling genes, not guides, because guides targeting the same gene are not independent.

Everything is then re-run on an independent screen: **CRISPRscan** (Moreno-Mateos et al. 2015) — 1,020 guides, 111 genes, zebrafish embryos, in vivo, a different lab and a different readout.

![Does structure add anything?](figures/incremental_value.svg)

## Findings

See [ANALYSIS.md](ANALYSIS.md) for the full account and [analysis_outputs/study_summary.md](analysis_outputs/study_summary.md) for every table.

## What is included

Folding:

- `nussinov.py` — Nussinov base-pair maximization, dependency-free.
- `energy_model.py` — nearest-neighbour free energies, plus an independent loop-by-loop structure energy evaluator and an exhaustive structure enumerator used for testing.
- `zuker.py` — minimum-free-energy folding.
- `mccaskill.py` — partition function, base-pair probabilities, Boltzmann sampling.

Analysis:

- `crispr_nussinov_analysis.py` — per-guide features under all three models (`--context sgrna` adds the scaffold).
- `guide_features.py` — spacer-alone and spacer-plus-scaffold feature extraction.
- `datasets.py` — loaders for both screens.
- `stats.py` — Spearman, ridge regression, grouped cross-validation, cluster bootstrap.
- `model_features.py` — the sequence baseline and structure feature blocks.
- `compute_features.py` — folds a whole screen in parallel and caches the result.
- `run_study.py` — the four experiments.
- `export_model.py` — fits the browser's efficiency model.

Everything above uses **only the Python standard library**.

## Installation

Python 3.10 or newer. No dependencies.

```bash
git clone https://github.com/kwland/nussinov-zuker-crispr.git
cd nussinov-zuker-crispr
python nussinov.py --sequence GGGAAACCC
```

An editable install adds the `rna-fold`, `crispr-guide-features`, `crispr-fold-screen`, and `crispr-run-study` commands:

```bash
python -m pip install -e .
```

## Reproduce the analysis

Fold a few example guides under all three models:

```bash
python crispr_nussinov_analysis.py --input examples/guides.csv --output analysis_outputs/example_guide_features.csv
```

Add the scaffold — slower, because the molecule becomes 96 nt instead of 20:

```bash
python crispr_nussinov_analysis.py --input examples/guides.csv --context sgrna --output analysis_outputs/example_sgrna_features.csv
```

Reproduce the whole study from scratch. The folding step is the expensive one (about 8 minutes for both screens on 11 cores) and is cached, so the analysis can be re-run without re-folding:

```bash
python compute_features.py --dataset doench
python compute_features.py --dataset crisprscan
python run_study.py
python make_figures.py
```

Refit the model used by the website:

```bash
python export_model.py
```

## Run the tests

```bash
python -m unittest discover -v
```

The suite covers three layers:

- **Folding correctness by brute force** — the partition function, every pair probability, and the MFE are checked against exhaustive enumeration of every pseudoknot-free structure on short sequences, including under stressed multiloop parameters. The Zuker folder reproduces the accepted yeast tRNA-Phe cloverleaf exactly.
- **Statistics** — ranking with ties, Benjamini–Hochberg monotonicity, permutation-test calibration under the null, Cholesky against known systems, and that grouped folds never split a gene.
- **The guard against fooling yourself** — nested cross-validation given pure noise must report no signal.

## How the three models differ

For an interval `(i, j)`:

- **Nussinov** stores the largest number of non-crossing pairs. Every allowed pair is worth +1. Transparent, but not physical: it over-pairs, and it cannot tell a stable helix from a weak one with the same pair count.
- **Zuker** stores the lowest free energy, using stacking, loop, and multiloop terms. One structure out.
- **McCaskill** stores the sum of `exp(-E/RT)` over all structures. Running the recursion outward as well as inward gives the probability of every pair, and hence of every base being unpaired.

The first two answer "what is the best structure?". Only the third answers "how often is this base free?", which is the question guide accessibility actually asks.

## Limitations

- **The energy parameters are partly simplified.** The ten Watson–Crick stacking values are the published Turner/Xia numbers. The G–U wobble terms, loop initiation tables, and multiloop model are documented approximations, so absolute kcal/mol will not match ViennaRNA exactly. Every approximation is flagged in `energy_model.py`.
- **Only pseudoknot-free secondary structure is modelled.** Tertiary contacts, kinetics, co-transcriptional folding, and Cas9 protein binding are all outside it. Cas9 actively unwinds and reshapes the guide, which may be why in-solution accessibility predicts so little.
- **Guide activity is not structure.** Chromatin, target context, repair pathway, expression, off-target binding, and nuclease behaviour all contribute and none are modelled.
- **The pooled Doench labels come from two screens** normalised within each source. That helps comparison but does not make the experiments identical.
- **Only 18 genes in the main screen.** Holding out whole genes is the right thing to do, but it leaves few independent units, which is why the confidence intervals are wide.
- **Every CRISPRscan spacer starts with GG** because of the T7/SP6 promoter, so positions 1–2 carry no information in that screen.
- **The 3D view is schematic**, a layout of the predicted secondary structure rather than a molecular model.
- **The committed Python example guides are synthetic** and carry no measured activity.
- **Correlation is not causation**, and a null result is evidence of absence only within the range these screens actually cover.

## Repository layout

```text
.
├── nussinov.py                  energy_model.py    zuker.py      mccaskill.py
├── guide_features.py            datasets.py        stats.py      model_features.py
├── crispr_nussinov_analysis.py  compute_features.py
├── run_study.py                 export_model.py    make_figures.py
├── test_nussinov.py             test_mccaskill.py  test_study.py
├── data/
│   ├── crispr_guide_examples.csv
│   └── crisprscan_moreno_mateos_2015.csv
├── analysis_outputs/            cached features, study_results.json, study_summary.md
├── figures/                     generated SVGs
├── docs/                        the static interactive site
└── notes/
```

## Data sources

Guide activity measurements come from published screens:

- Doench et al. (2014, 2016), pooled through CRISPOR / Haeussler et al. (2016) — `docs/data/guides.json`.
- Moreno-Mateos et al. (2015), CRISPRscan — `data/crisprscan_moreno_mateos_2015.csv`, reduced from the public CRISPOR paper dataset collection.

The ViennaRNA comparison values in the site's chapter 05 are output of the ViennaRNA package, precomputed because ViennaRNA does not run in a browser. Every folding result, accessibility value, energy, and model weight in this repository is computed by the code in it.

Key references:

- Nussinov & Jacobson (1980), *PNAS* — base-pair maximization.
- Zuker & Stiegler (1981), *Nucleic Acids Research* — minimum-free-energy folding.
- McCaskill (1990), *Biopolymers* — the partition function and base-pair probabilities.
- Ding & Lawrence (2003), *Nucleic Acids Research* — statistical sampling of RNA structures.
- Xia et al. (1998), *Biochemistry* — nearest-neighbour parameters.
- Doench et al. (2016), *Nature Biotechnology* — Rule Set 2 guide activity.
- Moreno-Mateos et al. (2015), *Nature Methods* — CRISPRscan.
- Haeussler et al. (2016), *Genome Biology* — CRISPOR.
- [ViennaRNA Package](https://www.tbi.univie.ac.at/RNA/).

## Author

Linus Tan.

## License

MIT. See [LICENSE](LICENSE).

---

Built as an explainable research prototype: small enough to inspect, verified against brute force where that is possible, and honest about a negative result. Questions, corrections, and collaboration are welcome.
