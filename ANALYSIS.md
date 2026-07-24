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

It changes nothing. In the Doench data, seed accessibility correlates with activity at ρ = +0.01 (Nussinov), +0.02 (MFE, spacer alone), −0.02 (ensemble, spacer alone), +0.03 (MFE, in sgRNA), and +0.04 (ensemble, in sgRNA). G/C content, for scale, reaches ρ = +0.13. CRISPRscan gives the same picture.

**The "you measured structure wrong" objection is answered: the sharper instrument finds the same nothing.**

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

*(filled in from the study run)*

## Experiment 4 — incremental value and replication

*(filled in from the study run)*

## The confound worth naming

In CRISPRscan, mean unpaired probability in sgRNA context correlates with activity at ρ = −0.13 — the largest structure-related number anywhere in this study, and pointing the *opposite* way to the hypothesis (less accessible, more active).

It is not a structure effect. Unpaired probability correlates with G/C content at ρ = −0.25, and G/C correlates with activity at ρ = +0.28 in that screen. Removing G/C by partial correlation collapses the structure signal to −0.06. Doing the same in the Doench data moves its correlation from +0.02 to +0.07 — the opposite direction. Two screens that disagree in sign after controlling for one obvious confound are describing noise, not biology.

This is exactly why the incremental test in Experiment 4 is the one that matters. A raw correlation cannot distinguish "structure predicts activity" from "structure is a proxy for G/C, and G/C predicts activity".

## Interpretation

*(filled in from the study run)*

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
