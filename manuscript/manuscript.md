---
title: "Folding Context Determines Measured RNA Accessibility in CRISPR Guides, and Bounds What Structure Can Add to Guide Activity Prediction"
author:
  - Linus Tan
  - Thomas Yu
date: "August 2026"
abstract: |
  Guide RNA secondary structure is frequently proposed as a determinant of CRISPR-Cas9 editing efficiency, on the reasoning that a guide sequestering its own spacer binds its DNA target less effectively. This study evaluates that proposition directly, using independent implementations of the Nussinov, Zuker, and McCaskill algorithms validated against ViennaRNA, applied to 5,705 guides with measured activity across two published screens.

  A methodological finding dominates the analysis. Accessibility computed by folding a spacer on its own bears little relation to accessibility computed by folding that same spacer inside the complete single guide RNA (Spearman correlation between the two measures, 0.17 for minimum free energy accessibility and 0.25 for ensemble accessibility). Folded alone, 22.4 percent of spacers are essentially unstructured, and mean seed accessibility falls from 0.84 to 0.42 once the 76 nt scaffold is included. Reported structural quantities therefore depend more on the folding context chosen than on the guide sequence.

  With scaffold context included, small within-gene associations between accessibility and activity are detectable in human cells (ensemble seed accessibility, within-gene Spearman 0.069, 95 percent confidence interval 0.010 to 0.130, gene-blocked permutation p = 0.001). These associations mostly do not survive as predictive value. Against a held-out sequence-only baseline evaluated by leave-one-gene-out validation, the addition of all structural features changed performance by 0.009 (95 percent confidence interval -0.006 to 0.025) in human cells, and features computed by ViennaRNA gave 0.002 (-0.004 to 0.007). One exception emerged. In zebrafish, ViennaRNA features improved held-out performance by 0.040 (0.012 to 0.072), an interval excluding zero. The gain is carried by ensemble free energy rather than by seed accessibility, and it runs opposite to the accessibility hypothesis, since more stable folding accompanies higher activity.

  Effect signs were not preserved across screens. Predicted seed accessibility does not rank ordinary guides usefully in either screen. Overall folding stability carries a small signal in one screen with reference energy parameters, which the simplified model missed.
keywords: "CRISPR-Cas9, guide RNA design, RNA secondary structure, partition function, McCaskill algorithm, cross-validation, negative result"
geometry: margin=1in
fontsize: 11pt
papersize: letter
numbersections: true
colorlinks: true
linkcolor: black
urlcolor: black
---

# Introduction

CRISPR-Cas9 locates a genomic target using a short RNA. The guide carries a spacer of approximately 20 nucleotides that base-pairs with one DNA strand, and Cas9 cleaves where that pairing occurs adjacent to a protospacer adjacent motif, which for *Streptococcus pyogenes* Cas9 is NGG [1]. Guides sharing this architecture nevertheless differ considerably in effectiveness, and two perfectly matched guides directed at the same gene may vary several-fold in editing efficiency.

Predicting that variation from sequence alone has drawn sustained attention. The most widely used models are trained on large screens in which thousands of guides were measured in parallel [2,3], and they perform well enough to serve as standard instruments in guide selection [4].

A structural proposition recurs alongside these sequence models. Because a guide RNA is single-stranded, it can fold back and pair with itself. If the spacer is locked up in its own secondary structure, the reasoning goes, less of it stays available to find and pair with DNA. The region most often identified is the seed, the PAM-proximal segment that leads target recognition and within which mismatches are least tolerated. Accessibility-based reasoning of this kind appears in guide design discussion and within feature sets constructed for scoring models [9,10].

The idea is plausible and testable, but three questions are rarely addressed together. The first concerns what should be folded. A functional guide is not a bare spacer; it carries a 76 nt scaffold that binds Cas9 and that can pair with the spacer. The second concerns how accessibility should be measured. A minimum free energy structure is a single point estimate drawn from an ensemble of conformations the molecule samples, and ensemble methods based on the partition function [11] give a different quantity. The third concerns whether any association translates into predictive value once sequence composition is accounted for, which requires held-out evaluation rather than correlation alone.

This study addresses those three questions. The folding algorithms were implemented directly instead of called from an existing package. That allowed inspection at every stage, and it allowed the same analysis to be repeated with ViennaRNA [8] as an external reference. Where a null result arises, that design distinguishes a limitation of the model from a limitation of the hypothesis.

Four contributions are reported. Independent implementations of the Nussinov, Zuker, and McCaskill algorithms are validated against ViennaRNA across complete datasets, not selected examples. Folding context is shown to dominate measured accessibility. Small within-gene associations between accessibility and activity are quantified with confidence intervals and gene-blocked permutation tests. Finally, the incremental predictive value of structural features is evaluated against a held-out sequence baseline, where seed accessibility adds nothing in either screen while overall folding stability adds a small amount in one of the two.

# Methods

## Folding algorithms

Three secondary structure algorithms were implemented from first principles, with no folding libraries.

The **Nussinov algorithm** [5] identifies the pseudoknot-free structure containing the greatest number of permitted base pairs. It populates a dynamic programming table across all subsequences, in which the entry for interval (i, j) represents the optimum among four alternatives: leaving base i unpaired, leaving base j unpaired, pairing i with j where the bases are compatible and the enclosed loop is of sufficient length, or bifurcating the interval and combining the two halves.

The **Zuker algorithm** [6] optimizes free energy rather than pair count, using three coupled matrices: V(i, j) for the energy given that i pairs with j, W(i, j) for the unconstrained energy of the interval, and WM(i, j) for intervals situated within a multiloop.

The **McCaskill algorithm** [11] computes the partition function over all pseudoknot-free structures, weighting each by its Boltzmann factor. Inside and outside recursions give base-pairing probabilities. The unpaired probability at each position is then one minus the summed pairing probability. Ensemble accessibility is defined from these probabilities. It measures something different from accessibility read off a single structure.

Watson-Crick nearest-neighbor stacking energies are the published Turner and Xia values at 37 degrees Celsius [7,12]. Loop initiation, wobble, asymmetry, and multiloop terms are simplified relative to the complete Turner parameter set, and the consequences of that simplification are quantified in Section 3.1 by direct comparison with ViennaRNA.

## Accessibility measures and folding context

The **seed** is defined as the final eight bases of the spacer. **Seed accessibility** is the fraction of those positions left unpaired, computed either from a single predicted structure or as a mean unpaired probability across the ensemble. **Mean unpaired probability** across the whole spacer was computed as an alternative that makes no assumption regarding which region matters.

Each measure was computed under two folding contexts: the spacer folded in isolation, and the spacer folded as part of the complete single guide RNA including the 76 nt scaffold. Comparing the two contexts is one of the main analyses reported here.

## Datasets

Two independent screens were analyzed. The first comprises 4,685 guides across 18 genes, pooled from two human cell screens [2,3] and obtained through CRISPOR [4], with activity expressed as a within-dataset percentile. The second comprises 1,020 guides across 111 genes from CRISPRscan [10], measured in zebrafish embryos in vivo. The two screens differ in organism, delivery, and assay, so the second is a genuine replication test and not a repetition.

## Statistical analysis

Guides aimed at the same gene share a baseline editability that has nothing to do with guide sequence. Pooled correlations across genes therefore conflate within-guide and between-gene variation. The **within-gene Spearman correlation**, computed per gene and averaged, is reported as the primary statistic throughout. Confidence intervals were obtained by bootstrapping over genes. Significance was assessed by permuting activity labels within each gene, never across genes, so that between-gene structure is preserved under the null hypothesis.

G/C content is associated with both folding and activity, so partial correlations removing G/C are reported next to the raw values.

Position-resolved analysis correlated activity against accessibility at each of the 20 spacer positions independently, with Benjamini-Hochberg correction applied across positions.

## Predictive evaluation

A sequence-only baseline was constructed from position and dinucleotide identity and G/C content, fitted by ridge regression. For the human screens, evaluation used leave-one-gene-out validation; for CRISPRscan, grouped five-fold validation repeated across three gene-to-fold assignments. The ridge penalty was chosen by within-gene Spearman correlation inside each training fold, never on the fold being predicted. Structural feature sets were then added to the baseline and the change in held-out performance recorded with bootstrap confidence intervals.

Cross-screen transfer was tested by training on one screen and applying the fitted model to another without retuning.

## External validation

The complete analysis was repeated using ViennaRNA [8] in place of the implementations presented here, across the same guides, contexts, and evaluation procedure. Agreement between the two was quantified, and the incremental value of ViennaRNA structural features was computed against the identical sequence baseline. Sensitivity to folding temperature (25, 37, and 42 degrees Celsius) and to local versus global folding was also examined.

## Availability

Code, data, and an interactive implementation are available at [repository URL]. The analysis pipeline depends only on the Python standard library and is reproducible from a single entry point.

# Results

## The implementations agree with ViennaRNA

Table 1 reports agreement between the implementations presented here and ViennaRNA across all 4,685 guides of the pooled human dataset.

: Agreement with ViennaRNA across the pooled human dataset (n = 4,685).

| Measure                                          | Pearson r | Spearman | Mean abs. diff. |
|:-------------------------------------------------|:---------:|:--------:|:---------------:|
| Ensemble seed accessibility (spacer alone)        | 0.841     | 0.892    | 0.129           |
| Ensemble seed accessibility (spacer + scaffold)   | 0.807     | 0.804    | 0.081           |
| Mean unpaired probability (spacer alone)          | 0.841     | 0.912    | 0.124           |
| Mean unpaired probability (in sgRNA)              | 0.770     | 0.753    | 0.066           |
| Ensemble free energy (spacer + scaffold)          | 0.833     | 0.819    | 4.367           |

Ensemble measures agree substantially more closely with ViennaRNA than minimum free energy measures, which reached Pearson correlations of 0.645 and 0.573 for the two folding contexts. Predicted free energies are systematically less negative, with a mean absolute difference of 4.37 kcal/mol, consistent with the simplified loop and wobble terms described in Section 2.1. The ordering of guides, which is what the subsequent analysis depends upon, is preserved.

## Folding context determines measured accessibility

What gets folded has a larger effect on reported accessibility than the guide sequence itself. Table 2 presents both contexts across the same guides.

: Effect of folding context upon accessibility measures (pooled human dataset).

| Measure                          | Spacer alone | In the sgRNA | Correlation between contexts |
|:---------------------------------|:------------:|:------------:|:----------------------------:|
| MFE seed accessibility            | 0.843        | 0.416        | 0.168                        |
| Ensemble seed accessibility       | 0.815        | 0.433        | 0.248                        |
| Mean unpaired probability         | 0.826        | 0.453        | 0.388                        |
| Ensemble free energy (kcal/mol)   | -0.780       | -19.951      | 0.311                        |

![Mean accessibility of the same 4,685 spacers, folded alone and folded inside the complete single guide RNA. The value r is the correlation between the two contexts across guides.](figures/fig1_folding_context.pdf){width=100%}

Mean seed accessibility approximately halves once the scaffold is included, and the two measures of the same underlying quantity correlate at only 0.168 for minimum free energy accessibility. A guide ranked as accessible in one context is therefore close to uninformative regarding its ranking in the other.

The isolated-spacer context is additionally degenerate. Folded alone, 1,049 of 4,685 spacers (22.4 percent) are essentially unstructured, defined as a mean unpaired probability exceeding 0.99. The equivalent figure in CRISPRscan is 166 of 1,020 (16.3 percent). A measure that gives the same maximal value to more than a fifth of the dataset cannot discriminate well.

## Associations with activity are small and context-dependent

Table 3 reports within-gene correlations with measured activity for the pooled human dataset.

: Within-gene Spearman correlations with activity, pooled human dataset (n = 4,685, 18 genes).

| Measure                                        | Within-gene rho | 95% CI            | Partial (G/C removed) | p       |
|:-----------------------------------------------|:---------------:|:-----------------:|:---------------------:|:-------:|
| Nussinov seed accessibility (spacer alone)      | 0.000           | [-0.037, 0.036]   | 0.017                 | 0.999   |
| MFE seed accessibility (spacer alone)           | 0.033           | [-0.025, 0.089]   | 0.070                 | 0.086   |
| Ensemble seed accessibility (spacer alone)      | 0.005           | [-0.056, 0.059]   | 0.032                 | 0.818   |
| MFE seed accessibility (spacer + scaffold)      | 0.056           | [0.008, 0.114]    | 0.042                 | 0.003   |
| Ensemble seed accessibility (spacer + scaffold) | 0.069           | [0.010, 0.130]    | 0.059                 | 0.001   |
| Mean unpaired probability (in sgRNA)            | 0.047           | [0.001, 0.089]    | 0.073                 | 0.014   |
| Ensemble free energy (spacer + scaffold)        | 0.055           | [-0.001, 0.113]   | 0.099                 | 0.001   |
| G/C percent                                     | -0.016          | [-0.093, 0.071]   | n/a                   | 0.416   |

![Within-gene Spearman correlation between each accessibility measure and measured activity, with 95 percent confidence intervals from resampling genes. Note that G/C content, the classic design heuristic, separates the two screens completely.](figures/fig2_correlations.pdf){width=100%}

Three observations follow. Accessibility computed on the isolated spacer shows no association with activity. The Nussinov measure returns a within-gene correlation of exactly zero. Accessibility computed with scaffold context shows a small positive association that survives permutation testing and partial correlation for G/C. Ensemble measures outperform their minimum free energy counterparts in every comparison, consistent with the agreement results in Section 3.1.

The magnitude nevertheless remains small. The strongest association, ensemble seed accessibility with scaffold context, corresponds to a within-gene correlation of 0.069.

G/C content shows the problem clearly. Its pooled correlation with activity is 0.134, but its within-gene correlation is -0.016 and not significant. The pooled figure reflects differences between genes, not a property of guides. This is why the within-gene statistic was adopted as primary.

Position-resolved analysis identified no position surviving multiple-testing correction in the human dataset under either folding context.

## Structure does not improve held-out prediction

The sequence-only baseline achieved a held-out within-gene correlation of 0.196 (95 percent confidence interval 0.112 to 0.280) under leave-one-gene-out validation. Table 4 reports the change upon adding structural features.

: Incremental value of structural features over a sequence baseline (pooled human dataset).

| Added features             | Held-out rho | Change | 95% CI of change  |
|:---------------------------|:------------:|:------:|:-----------------:|
| All structural features     | 0.205        | 0.009  | [-0.006, 0.025]   |
| Ensemble, spacer alone      | 0.202        | 0.006  | [-0.007, 0.019]   |
| MFE, spacer alone           | 0.200        | 0.005  | [-0.010, 0.018]   |
| Ensemble, full sgRNA        | 0.197        | 0.001  | [-0.003, 0.006]   |
| Nussinov                    | 0.196        | 0.000  | [-0.001, 0.001]   |
| Positional profile, spacer  | 0.196        | -0.000 | [-0.008, 0.008]   |

Every confidence interval includes zero. The largest point estimate, obtained by adding all structural features, is 0.009.

![Change in held-out Spearman correlation when structural features are added to the sequence-only baseline. Intervals crossing zero are shown in grey. Only the zebrafish result computed with reference Turner parameters excludes zero.](figures/fig3_incremental.pdf){width=100%}

The same evaluation performed with ViennaRNA features produced a change of 0.002 (95 percent confidence interval -0.004 to 0.007). In human cells the absence of incremental value is therefore not attributable to the simplified energy model presented here, since the reference implementation behaves equivalently.

Zebrafish behaved differently, and the difference is instructive. Against the CRISPRscan baseline of 0.387, the implementations presented here changed held-out performance by 0.009 (95 percent confidence interval -0.009 to 0.026), consistent with the human result. ViennaRNA features on the same guides changed it by 0.040 (0.012 to 0.072), an interval excluding zero, with ensemble features alone giving 0.042 (0.015 to 0.069). The gain is carried by ensemble free energy rather than by seed accessibility, and its direction is opposite to the accessibility hypothesis: more stable predicted folding accompanies higher measured activity. The simplified loop terms used here were sufficient to obscure that effect, which is the practical argument for validating a from-scratch model against a reference implementation rather than trusting agreement on selected examples.

Sensitivity analysis found the result stable across folding temperature. Within-gene correlation for ViennaRNA ensemble seed accessibility was 0.075 at 37 degrees, 0.082 at 25 degrees, and 0.072 at 42 degrees. Restricting folding to a local span of 40 nucleotides reduced it to 0.046.

## Effects do not replicate across screens

Table 5 compares the two screens directly.

: Within-gene correlations with activity across two independent screens.

| Measure                                     | Human cells (n = 4,685) | Zebrafish in vivo (n = 1,020) |
|:--------------------------------------------|:-----------------------:|:-----------------------------:|
| Ensemble seed accessibility (spacer alone)   | 0.005                   | -0.081                        |
| Ensemble seed accessibility (with scaffold)  | 0.069                   | 0.024                         |
| Mean unpaired probability (in sgRNA)         | 0.047                   | -0.134                        |
| Ensemble free energy (with scaffold)         | 0.055                   | -0.181                        |
| G/C percent                                  | -0.016                  | 0.307                         |

Several measures reverse sign between screens. Mean unpaired probability within the sgRNA is positively associated with activity in human cells and negatively associated in zebrafish, with both intervals excluding zero. G/C content is the clearest case: negligible within genes in human cells, and strongly positive in zebrafish (0.307, 95 percent confidence interval 0.244 to 0.371).

In CRISPRscan the position-resolved analysis did identify signal, with 8 of 20 positions significant after correction under scaffold context, the strongest at position 6 (rho = -0.131). No such structure appeared in the human data.

![Correlation with activity at each of the 20 spacer positions, spacer folded inside the sgRNA. Filled bars survive Benjamini-Hochberg correction. The shaded band marks positions 18 to 20, where WU-CRISPR reported accessibility to be predictive.](figures/fig4_positions.pdf){width=100%}

The incremental evaluation nevertheless returned the same answer in both screens. Against a CRISPRscan baseline of 0.387, adding all structural features changed held-out performance by 0.009 (95 percent confidence interval -0.009 to 0.026).

Cross-screen transfer was inconsistent. Training on Doench 2016 and testing on Doench 2014 improved performance by 0.092 (95 percent confidence interval 0.016 to 0.199) when structure was added, the single case in which structural features helped by a margin excluding zero. The reverse direction gave 0.002, and transfer between the human and zebrafish screens gave 0.019 and 0.001 respectively.

# Discussion

The central finding splits in two. Predicted seed accessibility does not improve prediction of CRISPR guide activity beyond what guide sequence already supplies, in either screen, whether computed by the implementations presented here or by ViennaRNA. Overall predicted folding stability is different. It adds a small but measurable amount in zebrafish when reference energy parameters are used, in the direction opposite to the accessibility hypothesis, and nothing in human cells.

The validation against ViennaRNA does most of the work here. A null result from a simplified energy model invites one obvious objection: that the model was not good enough. Repeating the evaluation with the reference implementation, which agrees with the present implementation at Spearman 0.80 to 0.91 for ensemble measures, answers it in both directions. In human cells both implementations return intervals including zero, so the null is not an artifact of simplified parameters. In zebrafish only the reference implementation finds a signal, so in that screen the simplified parameters were indeed hiding one. A study that had implemented folding and stopped there would have reported a clean null and been wrong about one of its two screens.

The finding about folding context may matter more than the negative result. Accessibility measured on an isolated spacer correlates with accessibility measured on the same spacer inside the complete guide at 0.168, and more than a fifth of spacers are unstructured when folded alone. Analyses that fold the spacer alone are therefore measuring something substantially different from analyses that include the scaffold, and the two should not be compared directly. The associations with activity appear only under scaffold context and disappear without it, so the choice changes the result.

The small associations that do appear need careful description. Under scaffold context in human cells, ensemble seed accessibility correlates with activity at 0.069 within genes, surviving permutation testing and adjustment for G/C. This is a real association by conventional standards. It is also far too small to filter guides with, and it does not translate into held-out predictive value. The sequence baseline appears to capture the relevant information already.

The failure to replicate across screens deserves emphasis. Signs reverse for several measures, and the two screens disagree sharply regarding G/C content. Differences in organism, delivery, expression, and assay all plausibly contribute. Whatever the cause, a design rule taken from one screen should not be assumed to transfer. These data give a concrete case where it does not.

The single case in which structure helped, transfer from Doench 2016 to Doench 2014, is reported for completeness. It is one of four transfer comparisons, and the other three returned intervals including zero, so it is better read as variability than as evidence.

Two possibilities remain open. Thermodynamic folding predictions describe equilibrium conformations, whereas a guide RNA is loaded into Cas9 co-transcriptionally, so kinetic accessibility during loading may differ from equilibrium accessibility. The models employed here also exclude protein interaction, and the scaffold in a functional complex is bound by Cas9 rather than free in solution.

# Limitations

Loop initiation, wobble, asymmetry, and multiloop terms are simplified relative to the complete Turner parameter set, and absolute free energies are correspondingly less negative than ViennaRNA values by approximately 4.4 kcal/mol. Agreement analysis in Section 3.1 and replication of the principal evaluation using ViennaRNA both address this limitation directly.

The models predict pseudoknot-free secondary structure only. Pseudoknots, tertiary contacts, folding kinetics, alternative conformations, and interactions with Cas9 protein lie outside their scope.

Activity labels are pooled across screens conducted under differing conditions and normalized within source datasets. Within-dataset percentiles render values comparable without rendering the experiments identical.

The human dataset contains 18 genes. Leave-one-gene-out validation over 18 groups gives wide confidence intervals, and the baseline interval of 0.112 to 0.280 shows it.

Analysis is correlational throughout and cannot establish that structure causes a change in editing activity.

Finally, absence of incremental value against this particular baseline does not exclude the possibility that structure contributes within a different model class or feature representation.

# Conclusion

Independent implementations of the Nussinov, Zuker, and McCaskill algorithms were validated against ViennaRNA and applied to 5,705 guides across two screens. Measured accessibility depends substantially on whether the spacer is folded alone or inside the complete guide, to the extent that the two contexts yield weakly correlated quantities and the isolated context is degenerate for approximately a fifth of guides. With scaffold context, small within-gene associations between accessibility and activity are detectable in human cells. Those associations do not improve on a held-out sequence baseline. Features computed by ViennaRNA do not improve on it either.

Predicted seed accessibility, the quantity the original hypothesis names, does not rank ordinary guides usefully in either screen. Overall folding stability carries a small signal in one screen of the two, only with reference energy parameters, and pointing the opposite way. Neither result says RNA structure is unimportant for guide function, since an effect confined to a minority of guides would be invisible in an average over thousands. Reports of structural effects should state the folding context and the energy model they were computed under, because this study shows both change the answer.

# Acknowledgments

The authors thank [teacher / mentor names] for guidance, and [institution] for supporting this work. The authors further thank Henry Zhang and Brady Lambert for early discussions concerning comparison methodology and the biological motivation for the structural hypothesis.

# Author contributions

L.T. implemented the Nussinov, Zuker, and McCaskill algorithms, constructed the feature extraction and statistical pipeline, performed the ViennaRNA validation, developed the interactive implementation, and drafted the manuscript. T.Y. assembled the pooled guide dataset, generated ViennaRNA reference values, and contributed to the predictive modeling. Both authors analyzed the results and approved the final manuscript.

# Data availability

Code, data, and an interactive implementation are available at [repository URL]. Activity measurements originate from Doench et al. 2014 and 2016, obtained through CRISPOR, and from Moreno-Mateos et al. 2015.

# Figure legends

**Figure 1.** Schematic representation of the three folding algorithms, presenting the Nussinov recursion cases, the Zuker decomposition, and the McCaskill inside and outside recursions from which base-pairing probabilities are obtained.

**Figure 2.** Agreement with ViennaRNA. Scatter plots of ensemble seed accessibility and ensemble free energy against ViennaRNA values across the pooled human dataset, with identity lines indicated.

**Figure 3.** Effect of folding context. Distribution of seed accessibility for spacers folded alone and within the complete single guide RNA, with the correlation between contexts indicated.

**Figure 4.** Within-gene correlations with activity for each accessibility measure and folding context, with bootstrap confidence intervals, presented for both screens.

**Figure 5.** Incremental value of structural features over the sequence baseline, presenting the change in held-out within-gene correlation with confidence intervals, for both the present implementation and ViennaRNA.

# References

1. Jinek M, Chylinski K, Fonfara I, Hauer M, Doudna JA, Charpentier E. A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity. *Science*. 2012;337(6096):816-821.

2. Doench JG, Fusi N, Sullender M, et al. Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nature Biotechnology*. 2016;34(2):184-191.

3. Doench JG, Hartenian E, Graham DB, et al. Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation. *Nature Biotechnology*. 2014;32(12):1262-1267.

4. Haeussler M, Schonig K, Eckert H, et al. Evaluation of off-target and on-target scoring algorithms and integration into the guide RNA selection tool CRISPOR. *Genome Biology*. 2016;17(1):148.

5. Nussinov R, Jacobson AB. Fast algorithm for predicting the secondary structure of single-stranded RNA. *Proceedings of the National Academy of Sciences*. 1980;77(11):6309-6313.

6. Zuker M, Stiegler P. Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information. *Nucleic Acids Research*. 1981;9(1):133-148.

7. Xia T, SantaLucia J, Burkard ME, et al. Thermodynamic parameters for an expanded nearest-neighbor model for formation of RNA duplexes with Watson-Crick base pairs. *Biochemistry*. 1998;37(42):14719-14735.

8. Lorenz R, Bernhart SH, Honer zu Siederdissen C, et al. ViennaRNA Package 2.0. *Algorithms for Molecular Biology*. 2011;6:26.

9. Wong N, Liu W, Wang X. WU-CRISPR: characteristics of functional guide RNAs for the CRISPR/Cas9 system. *Genome Biology*. 2015;16:218.

10. Moreno-Mateos MA, Vejnar CE, Beaudoin JD, et al. CRISPRscan: designing highly efficient sgRNAs for CRISPR-Cas9 targeting in vivo. *Nature Methods*. 2015;12(10):982-988.

11. McCaskill JS. The equilibrium partition function and base pair binding probabilities for RNA secondary structure. *Biopolymers*. 1990;29(6-7):1105-1119.

12. Turner DH, Mathews DH. NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure. *Nucleic Acids Research*. 2010;38(Database issue):D280-D282.

13. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57(1):289-300.
