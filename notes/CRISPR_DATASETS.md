# CRISPR Guide RNA Datasets and Project Plan

## Best Dataset Targets

1. **Doench et al. 2016 / Rule Set 2**
   - Best fit for this project because it contains measured SpCas9 guide activity and is widely used for on-target guide-scoring models.
   - Use it to ask: do guides with more Nussinov-predicted spacer self-pairing, especially in the PAM-proximal seed, tend to have lower measured activity?
   - Paper: https://www.nature.com/articles/nbt.3437

2. **WU-CRISPR / Wong, Liu, and Wang 2015**
   - Good biology background source because it discusses sequence features of functional CRISPR/Cas9 guide RNAs.
   - Use it to justify features like GC content, seed sequence, and self-folding/accessibility.
   - Paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0

3. **CRISPRscan / Moreno-Mateos et al. 2015**
   - Useful comparison dataset/tool for in vivo guide efficiency, especially if your group wants a second dataset beyond human-cell screens.
   - Paper/tool: https://www.nature.com/articles/nmeth.3543

4. **ViennaRNA RNAfold**
   - Not a dataset, but the best comparison tool for your code. Nussinov counts base pairs; RNAfold uses thermodynamic free energy.
   - Use it to compare dot-bracket structures for tRNA, 5S rRNA, and guide RNAs.
   - Tool/package: https://www.tbi.univie.ac.at/RNA/

## What the Checkpoint Advice Means

The advice says not to start by downloading a huge CRISPR screen and trying to solve everything at once. Start with a small comparison set:

- **1-2 known guide RNAs**: guides from a paper, class dataset, or lab source where the activity is already known. These are the positive controls.
- **1-2 random sequences of the same length**: 20 nt sequences that look guide-sized but are not chosen because they are known to work. These are neutral controls.
- **1-2 bad or stress-test sequences**: sequences with obvious design problems, such as extreme GC content, a long T/U run, or strong predicted self-pairing. These are negative controls.

This gives the project a clean first experiment: run the exact same Nussinov code on all groups and compare `nussinov_pairs`, `self_pair_fraction`, and `seed_accessibility`. If the known guides have more accessible seed regions than the bad controls, that supports the biological idea. If not, that is still useful because it shows the limits of a simple base-pair-counting model.

The committed rows in `data/crispr_guide_examples.csv` are now fully runnable synthetic demonstrations. They intentionally avoid `known_good` claims because they do not carry measured activity labels. A future biological analysis should use exact spacer sequences and activity values from a documented public dataset rather than relabeling these examples.

## Step-by-Step Coding Plan

1. **Nussinov baseline**
   - Implement dynamic programming.
   - Print maximum base-pair count.
   - Trace back one optimal structure.
   - Represent the structure in dot-bracket notation.

2. **Validation sequences**
   - Start with short toy RNAs where the answer is easy to inspect.
   - Then test a tRNA and 5S rRNA sequence.
   - Compare with ViennaRNA RNAfold, explaining that the methods optimize different goals.

3. **CRISPR guide features**
   - For every guide spacer, convert DNA to RNA.
   - Fold the 20-nt spacer with Nussinov.
   - Record total predicted self-pairs.
   - Record whether the PAM-proximal seed region, here the last 8 nt of the spacer, is paired or unpaired.
   - Add basic GC percent.
   - Flag simple design warnings: non-20-nt spacer, low/high GC, TTTT/poly-U risk, and seed mostly paired.

4. **Dataset analysis**
   - Download a public guide-activity table.
   - Keep columns for guide sequence and measured activity.
   - Run `crispr_nussinov_analysis.py` to add folding features.
   - Plot measured activity vs. Nussinov pairs and seed accessibility.

5. **Conclusion**
   - If highly folded guides show lower activity, argue that self-structure may hide the spacer from the DNA target.
   - If the relationship is weak, explain that real guide activity also depends on chromatin, repair context, target sequence, off-targets, and Cas9/tracrRNA interactions.

## Energy Model Notes

Nussinov is deliberately simple: every allowed base pair contributes +1, and the algorithm finds the non-crossing structure with the largest number of pairs. Zuker-style folding is more physical. Instead of asking "how many pairs can I make?", it asks "which structure has the lowest free energy?"

The main terms to understand are:

- **Stacking energy**: adjacent paired bases stabilize stems. This is often the biggest stabilizing term.
- **Hairpin loop penalty**: closing a small loop costs energy because the RNA backbone has to bend.
- **Internal loop and bulge penalties**: mismatched or unpaired bases inside a stem destabilize the structure.
- **Multiloop penalty**: junctions with several branches are harder to model and use fitted parameters.
- **Initiation and terminal mismatch terms**: smaller corrections learned from experiments.

For the project, the clean story is: implement Nussinov first because it is understandable dynamic programming, then explain that ViennaRNA/RNAfold improves realism by using these experimentally fitted free-energy terms.

## Datasets Actually Used

Two screens are in the repository, chosen to be as unlike each other as possible while
still measuring SpCas9 on-target activity.

| | Doench pooled | CRISPRscan |
|---|---|---|
| Source | Doench 2014 + 2016 via CRISPOR | Moreno-Mateos et al. 2015 |
| Guides | 4,685 | 1,020 |
| Genes | 18 | 111 |
| System | human cell culture | zebrafish embryos, in vivo |
| Readout | screen percentile | mutation frequency |
| File | `docs/data/guides.json` | `data/crisprscan_moreno_mateos_2015.csv` |

The CRISPRscan table was reduced from `effData/morenoMateos2015.context.tab` in the
public CRISPOR paper dataset collection (`maximilianh/crisporPaper`) to guide id, gene,
spacer, PAM, and measured mutation frequency.

Two properties of these screens matter when reading any result:

- The two Doench screens share **no genes** (2016 covers 15, 2014 covers CD13/CD15/CD33),
  so training on one and testing on the other is a real transfer test rather than a
  reshuffle of one experiment.
- Every CRISPRscan spacer begins with **GG**, because the guides are transcribed from a
  T7/SP6 promoter. Positions 1 and 2 therefore carry no information in that screen. This
  is a property of the assay, not a data-cleaning mistake, and it is one reason a model
  trained on human-cell data can transfer poorly.

## Current Local Files

Folding:

- `nussinov.py`: dependency-free Nussinov base-pair maximization.
- `energy_model.py`: nearest-neighbour free energies, shared by every folder below.
- `zuker.py`: minimum-free-energy folding.
- `mccaskill.py`: partition function, base-pair probabilities, Boltzmann sampling.

Analysis:

- `crispr_nussinov_analysis.py`: per-guide features under all three folding models.
- `guide_features.py`: spacer-alone and spacer-plus-scaffold feature extraction.
- `datasets.py`: loaders for both screens.
- `stats.py`: Spearman, grouped cross-validation, ridge regression, cluster bootstrap.
- `model_features.py`: the sequence baseline and structure feature blocks.
- `compute_features.py`: folds a whole screen and caches the result.
- `run_study.py`: the four experiments and the results tables.

Data and outputs:

- `data/crispr_guide_examples.csv`: runnable synthetic demo input, not an experiment.
- `data/crisprscan_moreno_mateos_2015.csv`: the independent replication screen.
- `analysis_outputs/`: cached folding features, `study_results.json`, `study_summary.md`.
