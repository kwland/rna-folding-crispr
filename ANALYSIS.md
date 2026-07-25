# Does RNA secondary structure predict CRISPR guide activity?

Linus Tan · July 2026 · **exploratory research prototype**

## The question

A CRISPR guide has to expose its 20 nt spacer to find a DNA target. If the guide folds back on itself — especially across the PAM-proximal seed, the part of the spacer that matters most for recognition — that sequence is tied up in base pairs and might be less available. The hypothesis:

> Guides whose seed is more accessible should tend to have higher measured editing activity.

If true, accessibility would be a cheap, sequence-only design filter.

## The answer, and why it took four experiments to trust it

The answer is that accessibility does not predict activity, in either of two independent screens, under any of five ways of measuring it.

An earlier pass of this project reached a similar conclusion from a single number: fold each spacer alone with Nussinov, count unpaired seed bases in one optimal structure, correlate with efficiency, get ρ ≈ 0.01. That result was right but not yet believable, because three obvious objections were unanswered — and a fourth, deeper one made the whole framing wrong.

1. *You measured structure wrong.* One optimal structure is not an RNA molecule.
2. *You folded the wrong molecule.* A bare spacer is not an sgRNA.
3. *You assumed the seed.* "The last 8 bases" was an assumption, never tested.
4. *You asked the wrong question.* Whether a feature correlates with activity on its own is not the question a guide designer has. The question is whether it adds anything to what sequence already tells you.

Each experiment below closes one of these.

## Methods

### Three folding models over one energy model

| | optimizes | output per base | in this repo |
|---|---|---|---|
| **Nussinov** | most base pairs | paired / unpaired | `nussinov.py` |
| **Zuker** | lowest free energy | paired / unpaired | `zuker.py` |
| **McCaskill** | the whole Boltzmann ensemble | probability of being unpaired | `mccaskill.py` |

Zuker and McCaskill share `energy_model.py` exactly. That matters: comparing Nussinov with McCaskill would change two things at once — the objective *and* the single-structure/ensemble distinction. Comparing Zuker with McCaskill changes only the second, which is the variable being tested.

The partition function follows McCaskill (1990). The inside recursions accumulate Boltzmann weights `exp(-E/RT)` over hairpins, stacks, bulges, internal loops, and multiloops; the outside recursions turn those into the probability `p_ij` of every pair, so the probability base *i* is unpaired is `1 - Σ_j p_ij`. Seed accessibility becomes a mean probability rather than a count.

The multiloop term in the outside pass is the one piece that resists a naive implementation: written directly it sums over every enclosing pair, which is O(n⁴). Factorising that sum into two accumulator tables, updated once per finalised pair, brings the whole pass back to O(n³).

### How the folding code was checked

Correctness here is not a matter of opinion. For sequences short enough to enumerate every pseudoknot-free structure, the exact answer can be computed directly and compared:

- Partition function, every base-pair probability, and every unpaired probability match brute-force enumeration to ~1e-15 across 15 test sequences.
- The shipped parameters make multiloops carry under 0.1% of the ensemble weight, so those recursions would barely be exercised. The tests therefore re-run the same comparison with multiloops made artificially cheap and unpaired multiloop bases made costly, pushing multiloop weight above 50% and activating the `z^k` factors that are otherwise all 1. The recursions still match exactly.
- Stochastic Boltzmann sampling, which uses only the *inside* matrices, converges to the probabilities the *outside* pass computes independently: max deviation 0.010 over 4,000 samples of a 96 nt guide.
- Zuker's MFE equals the minimum over all enumerated structures, and its traceback structure re-scores to the reported energy. On yeast tRNA-Phe it recovers the accepted cloverleaf exactly (pair F1 = 1.00).
- The Python and JavaScript folders agree to the last decimal on the full 96 nt sgRNA — same MFE, same seed openness, same model prediction.

### Data

| | Doench pooled | CRISPRscan |
|---|---|---|
| source | Doench 2014 + 2016 via CRISPOR | Moreno-Mateos et al. 2015 |
| guides | 4,685 | 1,020 |
| genes | 18 | 111 |
| system | human cell culture | zebrafish embryos, in vivo |
| readout | screen percentile | mutation frequency |

The two screens differ in organism, delivery, readout, and laboratory, which is what makes the second one a real test rather than a second slice of the first. The two Doench screens themselves share **no genes**, so training on one and testing on the other is also a genuine transfer.

### How results are scored

Two choices do most of the work in keeping the numbers honest.

**Cross-validation holds out whole genes.** Guides targeting the same gene share chromatin state, expression, and local composition. A random split puts near-duplicates on both sides and inflates everything.

**Accuracy is measured within genes, not pooled across them.** This one was discovered the hard way. A model that has never seen a gene cannot know how editable that gene is — in this data the mean activity percentile of a held-out gene ranges from 0.32 to 0.77, while the model's predictions for those genes barely move (0.45 to 0.54). Pooling out-of-fold predictions across genes therefore measures mostly whether the model guessed each gene's baseline level, which it cannot. The pooled correlation collapses to about zero even when the model ranks guides *within* every gene at ρ = 0.10 to 0.30. Ranking guides within a target is also the question a designer actually asks: they have one gene and need the best guide for it. Both numbers are reported throughout so the gap is visible rather than hidden.

Confidence intervals come from resampling genes, not guides, for the same reason.

## Experiment 1 — measurement

Five definitions of seed accessibility, correlated with activity.

The ensemble measure is not a cosmetic upgrade. Across 4,685 guides, single-structure seed accessibility takes only **9 distinct values** (it can only be a multiple of 1/8); the ensemble measure takes **4,668**. The two agree in the aggregate (r = 0.87 in sgRNA context) but disagree by more than 0.25 for 6.1% of guides. So this is a genuine change of instrument.

It changes very little. Doench, with 95% intervals from resampling genes:

| measure | pooled ρ | within-gene ρ |
|---|---:|---:|
| Nussinov seed, spacer alone | +0.011 | +0.000 [−0.037, 0.036] |
| MFE seed, spacer alone | +0.021 | +0.033 [−0.025, 0.089] |
| Ensemble seed, spacer alone | −0.023 | +0.005 [−0.056, 0.059] |
| MFE seed, in sgRNA | +0.031 | +0.056 [0.008, 0.114] |
| **Ensemble seed, in sgRNA** | +0.038 | **+0.069 [0.010, 0.130]** |
| G/C percent | +0.134 | −0.016 [−0.093, 0.071] |

The sharpest measure — ensemble accessibility of the seed in the intact sgRNA — does produce the largest number, in the direction the hypothesis predicts, with an interval that excludes zero. It is worth being precise rather than triumphal about that: ρ = 0.069 means accessibility explains under half a percent of the variance in activity, and the ordering across measures (Nussinov 0.00 → MFE 0.06 → ensemble 0.07) is consistent with a better instrument recovering a real but very faint effect.

Experiment 4 settles what it is worth. Adding that exact feature to a sequence-only model changes held-out accuracy by 0.000 [−0.004, 0.004] — whatever the correlation reflects, the sequence baseline already had it.

CRISPRscan does not reproduce even this. There, seed accessibility in context is +0.024 [−0.040, 0.091], and the ensemble measure for the *bare* spacer runs the other way at −0.081 [−0.145, −0.020].

**The "you measured structure wrong" objection is answered.** The sharper instrument finds, at most, a faint association that adds nothing once sequence is accounted for — and it does not replicate.

## Experiment 2 — molecular context

Folding the spacer inside the full sgRNA changes the measurement drastically:

| | spacer alone | in the full sgRNA |
|---|---:|---:|
| mean seed accessibility | 0.84 | 0.42 |
| mean unpaired probability | 0.83 | 0.45 |
| ensemble free energy | −0.78 kcal/mol | −19.95 kcal/mol |

The correlation *between* the two views is only r = +0.17. Folding a spacer on its own tells you almost nothing about how accessible it is in the molecule that actually exists.

There is a more basic point buried in that table. A bare 20 nt spacer has an ensemble free energy of −0.78 kcal/mol; **22% of them are essentially unstructured** (mean unpaired probability above 0.99). Most spacers simply do not fold on their own — there is hardly any structure there to correlate with anything. The original hypothesis was, in effect, testing a quantity that barely varies.

In context there is plenty of structure, and it varies. It still does not predict activity.

## Experiment 3 — position

With a per-base unpaired probability, accessibility can be correlated with activity at each of the 20 spacer positions instead of averaging over an assumed window.

| | strongest position | ρ there | positions significant (BH) |
|---|---:|---:|---:|
| Doench, spacer alone | 19 | −0.054 | 3 / 20 |
| Doench, in sgRNA | 20 | +0.061 | 1 / 20 |
| CRISPRscan, spacer alone | 1 | +0.055 | 0 / 20 |
| CRISPRscan, in sgRNA | 6 | −0.121 | 12 / 20 |

The canonical seed does not survive as the right window. In the Doench data folded in context, the strongest positions are 19 and 20 — the extreme PAM-proximal end, narrower than the assumed 13–20 (mean |ρ| 0.026 inside the seed against 0.011 outside). In CRISPRscan the pattern is the other way round: the strongest positions are 3 through 13, the PAM-**distal** half, with mean |ρ| 0.091 outside the seed against 0.035 inside.

Two screens that localise an effect to opposite ends of the molecule are not describing the same mechanism. If seed occlusion were real, both should point to the same place.

*A caveat on the significance counts.* The Benjamini–Hochberg correction is applied to permutation p-values that shuffle individual guides. That shuffle breaks the gene structure, so the effective sample size is nearer the number of genes than the number of guides and those p-values are anti-conservative. The "12 of 20" for CRISPRscan should be read as "the profile is smooth and non-random-looking", not as 12 independent discoveries. The cluster-bootstrap intervals, which resample whole genes, are the trustworthy uncertainty statement, and they are what the conclusions below rest on.

## Experiment 4 — incremental value and replication

The real question: does folding add anything to a sequence-only model? Baseline is position-specific bases, position-specific dinucleotides, and G/C, evaluated by nested cross-validation holding out whole genes and scored as mean within-gene Spearman.

| | baseline ρ | + all structure features | change | 95% CI |
|---|---:|---:|---:|---|
| Doench (human cells) | 0.195 [0.112, 0.276] | 0.197 | **+0.003** | [−0.013, 0.018] |
| CRISPRscan (zebrafish) | 0.391 [0.327, 0.450] | 0.407 | **+0.016** | [0.004, 0.030] |

**In human cells, structure adds nothing.** Every structure subset — Nussinov, MFE, ensemble, spacer-only, in-context, and the full 20-position profile — moves held-out accuracy by less than 0.004, with intervals straddling zero.

**In zebrafish, structure adds a small but detectable amount**: +0.016 on a baseline of 0.391, about a 4% relative gain, with an interval that excludes zero. The ensemble-in-context features alone give +0.011 [0.003, 0.021]. This is a real result and it should not be swept into the null.

But it does not support the hypothesis, for three reasons.

1. **It points the wrong way.** In CRISPRscan the ensemble free energy of the full sgRNA correlates with activity at within-gene ρ = −0.181 [−0.250, −0.105] and mean unpaired probability at −0.134 [−0.190, −0.073]. Both say *more* structure goes with *higher* activity — the opposite of "an occluded seed blocks targeting".
2. **Its position profile contradicts the mechanism**, sitting in the PAM-distal half rather than the seed.
3. **It does not replicate.** Nothing like it appears in the human-cell data, where the same features give +0.000.

### Cross-screen transfer

| train → test | baseline ρ | + structure | change | 95% CI |
|---|---:|---:|---:|---|
| Doench2016 → Doench2014 | 0.234 | 0.284 | +0.050 | [0.014, 0.115] |
| Doench2014 → Doench2016 | 0.153 | 0.155 | +0.002 | [−0.006, 0.010] |
| Doench pooled → CRISPRscan | 0.085 | 0.103 | +0.017 | [−0.015, 0.051] |
| CRISPRscan → Doench pooled | 0.018 | 0.019 | +0.001 | [−0.006, 0.007] |

Transfer is weak in every direction, and essentially absent from zebrafish to human cells (ρ = 0.018 — a model trained on CRISPRscan knows almost nothing about Doench guides). The one apparently large structure gain, Doench2016 → Doench2014 at +0.050, rests on a test set of just **three genes**, so its interval is wide and it should not be leaned on.

### An aside that surprised me

G/C content correlates with Doench activity at pooled ρ = +0.134 — the number the earlier version of this project reported as the one feature that beat structure. Within genes it is **−0.016 [−0.093, 0.071]**: nothing at all. Its apparent usefulness there is entirely a between-gene effect, i.e. genes with higher-G/C guide sets happen to be more editable. In CRISPRscan, by contrast, within-gene G/C is a strong ρ = +0.307.

The lesson generalises past this project: a feature that looks predictive pooled across targets can carry no information for the choice a designer actually makes.

## The confound worth naming

In CRISPRscan, mean unpaired probability in sgRNA context correlates with activity at ρ = −0.13 — the largest structure-related number anywhere in this study, and pointing the *opposite* way to the hypothesis (less accessible, more active).

It is not a structure effect. Unpaired probability correlates with G/C content at ρ = −0.25, and G/C correlates with activity at ρ = +0.28 in that screen. Removing G/C by partial correlation collapses the structure signal to −0.06. Doing the same in the Doench data moves its correlation from +0.02 to +0.07 — the opposite direction. Two screens that disagree in sign after controlling for one obvious confound are describing noise, not biology.

This is exactly why the incremental test in Experiment 4 is the one that matters. A raw correlation cannot distinguish "structure predicts activity" from "structure is a proxy for G/C, and G/C predicts activity".

## Interpretation

**The original hypothesis is not supported.** Seed accessibility does not predict editing efficiency, under five ways of measuring it, in two screens, and it adds nothing to a sequence-only model in human cells. The three obvious objections to the earlier ρ ≈ 0.01 have each been answered rather than argued away: a sharper instrument finds the same nothing, the biologically correct molecule finds the same nothing, and no spacer position rescues it.

**There is a small structure-related signal in zebrafish, and honesty requires reporting it.** Folding features add +0.016 [0.004, 0.030] there. It is small, it does not replicate to human cells, its sign is backwards for the hypothesis, and its position profile is in the wrong half of the spacer.

The most likely explanation is that it is not accessibility at all. The ensemble free energy of a 96 nt sgRNA is largely a nonlinear summary of the spacer's nucleotide composition, and zebrafish activity depends on composition far more strongly than human-cell activity does (within-gene G/C ρ = 0.307 against −0.016). A feature that summarises composition in a way a linear G/C term cannot will pick up some of that. Calling it "structure" would be reading the label rather than the mechanism.

**Why a genuine null is plausible here.** A bare spacer barely folds at all — ensemble free energy −0.78 kcal/mol, with 22% essentially unstructured. In the intact sgRNA there is plenty of structure, but the sgRNA does not float free: Cas9 binds it and imposes a conformation, holding the spacer in a pre-ordered A-form helix ready for target search. An equilibrium ensemble computed in solution may simply not describe the state the guide is in when it matters. That is a mechanistic reason to expect exactly what was measured, and it is more interesting than "the model was too crude".

**What the project is now good for.** Not guide ranking. The verified partition-function implementation, the exhaustive-enumeration tests, and the grouped-CV harness are reusable for any question of the form "does this RNA feature add anything?", and the pipeline makes the next such question cheap to ask.

## Limitations

- **The energy parameters are partly simplified.** The ten Watson–Crick stacking values are the published Turner/Xia numbers; the G–U wobble terms, loop initiation tables, and multiloop model are documented approximations. Absolute kcal/mol will not match ViennaRNA. Every approximation is flagged in `energy_model.py`. The recursions are exact *given* these parameters — brute-force agreement tests the algorithm, not the physics.
- **Only pseudoknot-free secondary structure is modelled.** No pseudoknots, tertiary contacts, folding kinetics, or alternative conformations.
- **Equilibrium folding may be the wrong physical picture.** Cas9 binds the sgRNA and actively imposes a conformation; the guide is loaded into a protein, not floating free in solution. An in-solution ensemble may simply not be the relevant state, which is a plausible reason for a genuine null.
- **Guide activity is not structure.** Chromatin, target context, repair pathway, expression, off-target binding, and nuclease behaviour all contribute and none are modelled.
- **The main screen has only 18 genes.** Holding out whole genes is right but leaves few independent units, so intervals are wide.
- **Every CRISPRscan spacer begins with GG** because of the T7/SP6 promoter, so positions 1–2 carry no information there and models that lean on them transfer badly.
- **Pooled Doench labels come from two screens**, normalised within each. That aids comparison without making the experiments identical.
- **A null result is bounded by what these screens cover.** These are SpCas9, NGG PAM, 20 nt spacers, in two systems. It says nothing about other nucleases, chemically modified guides, or truncated spacers.
- **Correlation is not causation**, in either direction.

## What would change the conclusion

- A screen designed to vary predicted structure while holding G/C and position composition fixed, rather than one that happens to contain some structured guides.
- Accessibility measured on the Cas9-bound conformation rather than the free RNA.
- Direct measurement — chemical probing such as SHAPE on the actual sgRNAs — instead of prediction.

---

This is a careful negative result, not a finished predictive system. The code is written to make the assumptions inspectable and the next experiment easier to design.
