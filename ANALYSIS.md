# Predicted sgRNA accessibility adds little to sequence-based CRISPR–Cas9 activity models across two screens

Linus Tan · July 2026 · **exploratory research prototype**

## The question

A CRISPR guide has to expose its 20 nt spacer to find a DNA target. If the guide folds back on itself — especially across the PAM-proximal seed, the part of the spacer that matters most for recognition — that sequence is tied up in base pairs and might be less available. The hypothesis:

> Guides whose seed is more accessible should tend to have higher measured editing activity.

If true, accessibility would be a cheap, sequence-only design filter.

### This is a reassessment, not a first look

The idea is not new, and this work is not the first investigation of sgRNA structure. Wong, Liu and Wang (2015) built **WU-CRISPR** on exactly this premise and reported that accessibility in the PAM-proximal region, specifically around **positions 18–20**, was among their more informative features. Several later design tools carried some version of that feature forward.

What has been thin is independent, well-controlled reassessment. The specific gaps this analysis targets:

- Accessibility read off a **single predicted structure** rather than the Boltzmann ensemble.
- Folding the **bare spacer** rather than the sgRNA that exists in a cell.
- **Assuming** the predictive window instead of measuring every position.
- Reporting a **raw correlation** rather than incremental value over a sequence baseline.
- Evaluating on **randomly split** guides, which leaks between guides targeting the same gene.
- **One screen**, with no test of whether the effect transfers.

So the framing is: given a better-controlled test, how much of the WU-CRISPR accessibility claim survives? The answer below is "a small amount in one screen, and not in the direction or the position the original claim would predict".

## The answer, and why it took five experiments to trust it

Accessibility carries almost no usable signal about activity in either screen, under any of five ways of measuring it, and it adds nothing to a sequence-only model in human cells.

An earlier pass of this project reached a similar conclusion from a single number: fold each spacer alone with Nussinov, count unpaired seed bases in one optimal structure, correlate with efficiency, get ρ ≈ 0.01. That result was right but not yet believable, because four obvious objections were unanswered — and a fifth, about the model itself, would have been raised by any reviewer.

1. *You measured structure wrong.* One optimal structure is not an RNA molecule.
2. *You folded the wrong molecule.* A bare spacer is not an sgRNA.
3. *You assumed the seed.* "The last 8 bases" was an assumption, never tested.
4. *You asked the wrong question.* Whether a feature correlates with activity on its own is not the question a guide designer has. The question is whether it adds anything to what sequence already tells you.
5. *Your energy model is not the standard one.* Simplified loop terms could produce a null that a reference implementation would not.

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

### Validation against ViennaRNA

The energy parameters here use the published Watson–Crick stacking values but simplified wobble, loop-initiation, and multiloop terms. Documenting an approximation is not the same as measuring what it costs, so every headline quantity was recomputed with **ViennaRNA 2.7.2** under the standard Turner parameters (`vienna_reference.py`, `validate_vienna.py`).

The intended division of labour between the two implementations:

- **ViennaRNA is the physical model.** Any quantitative claim about RNA thermodynamics should rest on it.
- **The implementation here is the methodological contribution.** Its value is that it is transparent, dependency-free, verified against exhaustive enumeration, and small enough to run live in a browser — not that it is a better energy model.

**Agreement** (Doench, 4,685 guides, spacer folded in the full sgRNA):

| Measure | Pearson r | mean abs difference |
|---|---:|---:|
| Ensemble seed accessibility | 0.807 | 0.081 |
| Mean unpaired probability | 0.770 | 0.066 |
| MFE seed accessibility | 0.573 | 0.136 |
| Ensemble free energy | 0.833 | 4.36 kcal/mol |

Two things stand out. Free energies carry a systematic offset of about **+4.4 kcal/mol**: the simplified loop terms are not stabilising enough, exactly as the parameter documentation warns, so absolute energies from this model should never be quoted as thermodynamic values. And **ensemble measures agree substantially better than single-structure ones** (r = 0.81 against 0.57). That is expected rather than surprising: a single optimal structure can flip between two near-tied alternatives when a parameter moves by a tenth of a kcal/mol, while an average over all structures moves smoothly. It is also an argument for the ensemble measure independent of anything biological.

**Accessibility correlations barely move.** Recomputing every accessibility-to-activity correlation from ViennaRNA features shifts no value by more than 0.033, and the two implementations agree on which measures are near zero:

| Measure | custom | ViennaRNA |
|---|---:|---:|
| Ensemble seed accessibility (in sgRNA) | +0.069 | +0.075 |
| MFE seed accessibility (in sgRNA) | +0.056 | +0.066 |
| Mean unpaired probability (in sgRNA) | +0.047 | +0.066 |
| Ensemble seed accessibility (spacer alone) | +0.005 | −0.007 |

**But the free energies do not, and that changes one result.** The simplified loop terms cost real accuracy in the *energy*, and in the zebrafish screen that matters. ViennaRNA's ensemble free energy correlates with activity at within-gene ρ = −0.293, against −0.181 for the custom model, and it survives adjustment for G/C (−0.209).

The consequence shows up in the incremental test:

| Screen | custom features | ViennaRNA features |
|---|---:|---:|
| Doench (human cells) | +0.009 [−0.006, 0.025] | +0.002 [−0.004, 0.007] |
| CRISPRscan (zebrafish) | +0.009 [−0.009, 0.026] | **+0.040 [0.012, 0.072]** |

In human cells both models agree on nothing. In zebrafish, **the simplified parameters were hiding a real effect**: with proper thermodynamics, folding features add about 0.04 to a sequence-only baseline of 0.387, and the interval excludes zero. Repeating over three gene-to-fold assignments gives +0.040, +0.032, and +0.046 — a spread of 0.014, comfortably smaller than the effect, so this one is not a fold-assignment artefact of the kind that removed the earlier claim.

This is the validation earning its keep. It did not merely confirm the custom model; it found a place where the approximations changed a conclusion, and the direction of the correction is against the project's own prior result.

**Temperature and folding context** (ViennaRNA, Doench, ensemble seed accessibility):

| Setting | ρ |
|---|---:|
| 37 °C, global folding | +0.075 |
| 25 °C, global folding | +0.082 |
| 42 °C, global folding | +0.072 |
| 37 °C, local folding, max span 40 | +0.046 |

The result is stable across a 17-degree range. Restricting to local folding, which approximates a molecule folding as it is transcribed rather than equilibrating as a whole, *weakens* the association rather than strengthening it. No setting turns accessibility into a usable predictor.

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

Four choices do most of the work in keeping the numbers honest.

**Cross-validation holds out whole genes.** Guides targeting the same gene share chromatin state, expression, and local composition. A random split puts near-duplicates on both sides and inflates everything.

**Accuracy is measured within genes, not pooled across them.** This one was discovered the hard way. A model that has never seen a gene cannot know how editable that gene is — in this data the mean activity percentile of a held-out gene ranges from 0.32 to 0.77, while the model's predictions for those genes barely move (0.45 to 0.54). Pooling out-of-fold predictions across genes therefore measures mostly whether the model guessed each gene's baseline level, which it cannot. The pooled correlation collapses to about zero even when the model ranks guides *within* every gene at ρ = 0.10 to 0.30. Ranking guides within a target is also the question a designer actually asks: they have one gene and need the best guide for it. Both numbers are reported throughout so the gap is visible rather than hidden.

**Significance tests permute inside genes, not across them.** The null worth testing is "within a gene, this feature says nothing about activity". Shuffling activity freely across the whole dataset instead tests "gene identity does not matter either", a null so easily rejected that the resulting p-values are far too small. An earlier version of this analysis used the free shuffle, acknowledged in the text that the p-values were anti-conservative, and then reported significance counts from them anyway. That was not defensible, and those counts have been replaced rather than annotated. `test_study.py` measures both procedures on the same simulated data and requires the gene-blocked one to stay near its nominal error rate while the free shuffle does not.

**The ridge penalty is selected on the metric being reported.** Tuning on pooled Spearman and reporting mean within-gene Spearman would optimise a different objective than the one being judged. The human screen additionally uses **leave-one-gene-out** validation, which removes any dependence on how genes were assigned to folds. The zebrafish screen has 111 genes, where that would mean 111 nested fits per feature set, so it uses grouped five-fold repeated over three assignments and the spread across those assignments is reported alongside the estimate.

Confidence intervals come from resampling genes, not guides, for the same reason.

## Experiment 1 — measurement

Five definitions of seed accessibility, correlated with activity.

The ensemble measure is not a cosmetic upgrade. Across 4,685 guides, single-structure seed accessibility takes only **9 distinct values** (it can only be a multiple of 1/8); the ensemble measure takes **4,668**. The two agree in the aggregate (r = 0.87 in sgRNA context) but disagree by more than 0.25 for 6.1% of guides. So this is a genuine change of instrument.

It changes very little. Doench, with 95% intervals from resampling genes and p-values from shuffling activity inside genes:

| measure | pooled ρ | within-gene ρ | partial, G/C removed | gene-blocked p |
|---|---:|---:|---:|---:|
| Nussinov seed, spacer alone | +0.011 | +0.000 [−0.037, 0.036] | +0.017 | 1.00 |
| MFE seed, spacer alone | +0.021 | +0.033 [−0.025, 0.089] | +0.070 | 0.09 |
| Ensemble seed, spacer alone | −0.023 | +0.005 [−0.056, 0.059] | +0.032 | 0.82 |
| MFE seed, in sgRNA | +0.031 | +0.056 [0.008, 0.114] | +0.042 | 0.003 |
| **Ensemble seed, in sgRNA** | +0.038 | **+0.069 [0.010, 0.130]** | +0.059 | 0.001 |
| G/C percent | +0.134 | −0.016 [−0.093, 0.071] | — | 0.42 |

The sharpest measure, ensemble accessibility of the seed in the intact sgRNA, does produce the largest value, in the direction the hypothesis predicts, with an interval excluding zero and a gene-blocked p of 0.001. It survives adjustment for G/C (+0.059). Being precise rather than triumphal about it: ρ = 0.069 means accessibility accounts for under half a percent of the variance in activity, and the ordering across measures (Nussinov 0.00 → MFE 0.06 → ensemble 0.07) is what a better instrument recovering a real but very faint effect would look like.

Experiment 4 settles what it is worth. Adding that exact feature to a sequence-only model changes held-out accuracy by +0.001 [−0.003, 0.006] — whatever the correlation reflects, the sequence baseline already had it.

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

With a per-base unpaired probability, accessibility can be correlated with activity at each of the 20 spacer positions instead of averaging over an assumed window. Each correlation is the mean within gene, and each p-value comes from shuffling activity inside genes, then Benjamini–Hochberg across the 20 positions.

| | strongest position | ρ there | positions significant (BH) |
|---|---:|---:|---:|
| Doench, spacer alone | 19 | −0.042 | **0 / 20** |
| Doench, in sgRNA | 14 | +0.051 | **0 / 20** |
| CRISPRscan, spacer alone | 5 | −0.081 | **0 / 20** |
| CRISPRscan, in sgRNA | 6 | −0.131 | 8 / 20 |

**The WU-CRISPR position 18–20 claim is not reproduced.** In the human-cell screen no position survives correction in either folding context. Under the earlier, invalid free-shuffle test the same data appeared to yield three significant positions folded alone and one in context; those counts were an artefact of ignoring that guides cluster within genes.

Where positions do survive, in the zebrafish screen folded in context, they are **3 through 13**, the PAM-*distal* half. That is the opposite end of the spacer from the region the original claim identifies, and it is not what a seed-occlusion mechanism predicts.

## Experiment 4 — incremental value and replication

The question that matters: does folding add anything to a sequence-only model? Baseline is position-specific bases, position-specific dinucleotides, and G/C. The human screen uses leave-one-gene-out over its 18 genes; the zebrafish screen uses grouped five-fold repeated over three gene-to-fold assignments. In both, the ridge penalty is chosen inside each training fold on the same within-gene metric that is reported.

| | baseline ρ | + all structure features | change | 95% CI |
|---|---:|---:|---:|---|
| Doench (human cells) | 0.196 [0.112, 0.280] | 0.205 | **+0.009** | [−0.006, 0.025] |
| CRISPRscan (zebrafish) | 0.387 [0.321, 0.446] | 0.396 | **+0.009** | [−0.009, 0.026] |

**Folding adds nothing detectable in either screen.** Every structure subset — Nussinov, MFE, ensemble, spacer-only, in-context, and the full 20-position profile — moves held-out accuracy by less than 0.01, with every interval spanning zero.

This is a change from an earlier version of this analysis, and the reason is worth recording. That version reported a small but apparently real gain in zebrafish, +0.016 with an interval excluding zero. Two corrections removed it:

1. **The ridge penalty was being chosen on pooled Spearman** while the reported metric was within-gene, so the selected model was optimised for a different objective than the one being judged.
2. **A single gene-to-fold assignment was used.** Repeating over three assignments, the zebrafish baseline alone scores 0.387, 0.394, and 0.379. The spread from the split — about **0.015** — is larger than the effect being claimed. A single split can manufacture a finding of this size.

That is exactly the failure mode repeated cross-validation exists to catch, and it is why the corrected result is reported instead of the more interesting one.

### Cross-screen transfer

| train → test | baseline ρ | + structure | change | 95% CI |
|---|---:|---:|---:|---|
| Doench2016 → Doench2014 | 0.234 | 0.326 | +0.092 | [0.016, 0.199] |
| Doench2014 → Doench2016 | 0.153 | 0.155 | +0.002 | [−0.006, 0.010] |
| Doench pooled → CRISPRscan | 0.085 | 0.104 | +0.019 | [−0.010, 0.048] |
| CRISPRscan → Doench pooled | 0.018 | 0.019 | +0.001 | [−0.006, 0.007] |

Transfer is weak in every direction and essentially absent from zebrafish to human cells (ρ = 0.018). The one large apparent gain, Doench2016 → Doench2014 at +0.092, rests on a test set of **three genes**; its interval is correspondingly wide and it should not be leaned on. Reversing the same comparison gives +0.002.

### An aside that surprised me

G/C content correlates with Doench activity at pooled ρ = +0.134 — the number an earlier version of this project reported as the one feature that beat structure. Within genes it is **−0.016 [−0.093, 0.071]**: nothing at all. Its apparent usefulness there is entirely a between-gene effect, meaning genes whose guide sets happen to be higher in G/C are more editable. In CRISPRscan, by contrast, within-gene G/C is a strong ρ = +0.307.

The lesson generalises past this project: a feature that looks predictive pooled across targets can carry no information for the choice a designer actually makes.

## The confound worth naming

In CRISPRscan, mean unpaired probability in sgRNA context correlates with activity at ρ ≈ −0.13 — the largest structure-related number anywhere in this study, and pointing the *opposite* way to the hypothesis (less accessible, more active).

Much of it is composition, not structure. Unpaired probability tracks G/C content closely, and G/C is itself strongly associated with activity in that screen. `stats.partial_spearman` removes G/C by regressing both variables on its ranks and correlating the residuals; it is computed for every measure in `run_study.py` and reported in the `partial (G/C removed)` column of [analysis_outputs/study_summary.md](analysis_outputs/study_summary.md), so the values below are pipeline output rather than a side calculation.

Controlling for G/C shrinks the zebrafish association substantially, and moves the human-cell one in the *opposite* direction. Two screens that disagree in sign after adjusting for one obvious confound are describing noise, not a mechanism.

This is exactly why the incremental test in Experiment 4 is the one that matters. A raw correlation cannot distinguish "structure predicts activity" from "structure is a proxy for G/C, and G/C predicts activity". The incremental test answers that directly, because the sequence baseline already contains G/C.

## Interpretation

### What is actually supported

Two claims, and the distinction between them is the result:

> **1. Predicted seed accessibility adds no measurable value for ranking ordinary SpCas9 guides in either screen.**
>
> **2. Overall predicted folding stability adds a small amount in zebrafish (about +0.04 on a 0.387 baseline, with standard Turner parameters) and nothing in human cells. Its direction is opposite to the accessibility hypothesis: more stable folding accompanies *higher* activity.**

Every qualifier is doing work.

- **Predicted** — from a thermodynamic folding model, not measured by chemical probing.
- **Seed accessibility** specifically, as distinct from whole-molecule folding stability. The
  first carries nothing; the second carries a little, in one screen.
- **Ordinary** — standard 20 nt spacers with the standard scaffold, no unusual designs.
- **These two screens** — human cell culture and zebrafish embryos, which disagree.

The zebrafish effect deserves neither dismissal nor promotion. It is real under repeated
cross-validation and survives adjustment for G/C, but it is small, it appears only with a
reference energy model, it does not transfer to human cells, and it is carried by ensemble
free energy rather than by anything about the seed. The most defensible reading is that
*something* about overall sgRNA stability matters in that in vivo system — plausibly guide
persistence or loading rather than target interrogation — and that it is not the mechanism
the accessibility hypothesis proposes.

### What is *not* supported

This work does **not** show that RNA structure is unimportant to CRISPR, and that broader claim would contradict direct experimental evidence. Misfolding is known to matter for particular refractory guides and for engineered scaffolds: Riesenberg et al. (2022, *Nature Communications* 13:489) showed that specific sgRNAs fail through predictable misfolding and that redesigning them rescues activity. That is a real structural effect. It is compatible with everything here, because a strong effect confined to a minority of guides can be invisible in an average taken over thousands.

The right reading is about *ranking utility*, not about biology: predicted seed accessibility is not a useful general-purpose term to add to a guide-scoring model, even though structure can be decisive for individual problem guides. The zebrafish result cuts the same way from the other side — overall folding stability *does* carry a little ranking information there, so this is not a blanket claim that folding predictions are useless either.

### On the original hypothesis

Seed accessibility does not usefully predict editing efficiency here. Under five ways of measuring it, in two screens, with a reference implementation agreeing, the largest within-gene correlation found anywhere is +0.069, and adding any *accessibility* feature to a sequence-only model changes held-out accuracy by less than 0.01 with every interval spanning zero.

The one place folding does help, the zebrafish screen with ViennaRNA parameters, is driven by ensemble free energy rather than by accessibility, so it does not rescue the hypothesis.

The five objections to the earlier ρ ≈ 0.01 have each been answered rather than argued away: a sharper instrument finds the same nothing, the biologically correct molecule finds the same nothing, no spacer position survives a valid significance test, the incremental question gives nothing, and ViennaRNA under standard Turner parameters agrees throughout.

The WU-CRISPR position 18–20 claim is not reproduced as stated. No position in the human-cell screen survives gene-blocked correction in either folding context. The positions that do survive, in zebrafish, sit in the PAM-distal half, which is inconsistent with the proposed mechanism.

### What changed, and why that matters

An earlier version of this analysis reported a small but apparently real gain in the zebrafish screen: +0.016 with an interval excluding zero. It did not survive two methodological corrections — selecting the ridge penalty on the metric actually being reported, and repeating over several gene-to-fold assignments rather than trusting one. The spread across assignments (0.015) turned out to be larger than the claimed effect (0.016).

That is recorded here rather than quietly replaced, because it is the most transferable lesson in the project: with few independent groups, a single cross-validation split can manufacture an effect of exactly the size researchers find publishable.

## Limitations

- **The energy parameters are partly simplified.** The ten Watson–Crick stacking values are the published Turner/Xia numbers; the G–U wobble terms, loop initiation tables, and multiloop model are documented approximations, and the measured cost is a systematic +4.4 kcal/mol offset against ViennaRNA. Absolute kcal/mol from this model are not thermodynamic values. Every approximation is flagged in `energy_model.py`; `validate_vienna.py` quantifies the consequences and finds no conclusion depends on them. The recursions are exact *given* these parameters — brute-force agreement tests the algorithm, not the physics.
- **Only pseudoknot-free secondary structure is modelled.** No pseudoknots, tertiary contacts, folding kinetics, or alternative conformations.
- **Equilibrium folding may be the wrong physical picture.** Cas9 binds the sgRNA and actively imposes a conformation; the guide is loaded into a protein, not floating free in solution. An in-solution ensemble may simply not be the relevant state, which is a plausible reason for a genuine null.
- **Guide activity is not structure.** Chromatin, target context, repair pathway, expression, off-target binding, and nuclease behaviour all contribute and none are modelled.
- **The main screen has only 18 genes.** Holding out whole genes is right but leaves few independent units, so intervals are wide.
- **Every CRISPRscan spacer begins with GG** because of the T7/SP6 promoter, so positions 1–2 carry no information there and models that lean on them transfer badly.
- **Pooled Doench labels come from two screens**, normalised within each. That aids comparison without making the experiments identical.
- **A null result is bounded by what these screens cover.** These are SpCas9, NGG PAM, 20 nt spacers, in two systems. It says nothing about other nucleases, chemically modified guides, or truncated spacers.
- **Correlation is not causation**, in either direction.
- **The reported effects are averages over thousands of guides.** A feature can be worthless on average and still be decisive for a particular guide, which is what the misfolding literature documents. Nothing here bears on that case.
- **Significance depends on the permutation scheme.** The reported tests shuffle activity inside each gene. A free shuffle across the whole dataset produces far smaller p-values and, on this data, would report several positions as significant that the correct test does not.
- **The zebrafish screen uses grouped five-fold rather than leave-one-gene-out**, because 111 genes would mean 111 nested fits per feature set. The spread across three fold assignments is reported so the cost of that choice is visible.
- **ViennaRNA validation covers the folding model, not the biology.** It shows the conclusions do not depend on the simplified parameters. It cannot show that equilibrium folding is the right physical picture in the first place.

## What would change the conclusion

- A screen designed to vary predicted structure while holding G/C and position composition fixed, rather than one that happens to contain some structured guides.
- Accessibility measured on the Cas9-bound conformation rather than the free RNA.
- Direct measurement — chemical probing such as SHAPE on the actual sgRNAs — instead of prediction.
- A focused test on guides *predicted* to misfold badly, rather than an average over all guides. If structure matters for a minority, that minority is where to look, and it is the design that the misfolding literature supports.
- More genes. With 18 in the human screen, the gene-resampled intervals are wide enough that a real effect of the size seen here could not be distinguished from zero.

---

This is a careful negative result, not a finished predictive system. The code is written to make the assumptions inspectable and the next experiment easier to design.
