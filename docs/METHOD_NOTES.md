# Method notes

- PPI network: STRING v12 export, confidence >= 0.4 combined score (see manuscript Methods for exact
  thresholds); hub consensus: cytoHubba (MCC, MNC, Degree, Closeness, Radiality, EPC), top-10 per
  algorithm; cross-validated against a custom Python implementation (records of the comparison
  are available from the authors).
- Machine learning: three selector families (LASSO logistic regression, random forest, SVM-RFE),
  intersection-based signature; external testing in GSE48452 (and per-cohort splits as described).
- Docking: Vina, exhaustiveness 32 (see script); receptor preparation as documented in `code/07_docking`.
- MD: OpenMM 8.6.1; Amber14SB + OBC2; OpenFF 2.2.0 ("Sage") for PFOA (MMFF94 charges, net charge -1);
  hydrogen-mass repartitioning (3 amu); 4-fs timestep; hydrogen bonds constrained; 300 K Langevin
  (1 ps^-1); 100 ns production per complex; frames every 50 ps.
- MM-GBSA: single-trajectory approach; 25 snapshots per complex (50-100 ns); user is referred to the
  manuscript for the interpretation caveats (anionic ligand; relative ranking).

## Paths and environment

Scripts retain the original project-root variables (e.g., `BASE`, `base`, `RUNS`); update them to your local paths before rerunning.

## NHANES anchor (revision round 2026-09-20)

Survey-weighted linear models (R `survey`; two-cycle PFAS subsample weights, strata/PSU design) of log-transformed ALT/AST/GGT and FLI on log2 serum PFOA (sum of linear and branched isomers), adjusted for age, sex, race/ethnicity, BMI, cotinine and cycle (NHANES 2015-2018; n = 3,212 adults). Outputs in `data/nhanes/`.

## cytoHubba EPC note (2026-09-20)

A naive fixed-threshold percolation variant of EPC degenerates on this dense network (mass ties) and is excluded; the final analysis uses the official cytoHubba plugin EPC, reproduced with an adaptive-threshold implementation (`run_09`; Spearman rho = 0.9965; top-10 overlap 8/10 within Monte-Carlo noise). Full score exports are provided in `data/`.

## Round-2 revision (2026-09-20, v1.5)
- NHANES: FLI reported on the point scale (-0.46 points per doubling); per-model n in nhanes_results.csv (scale/unit/est columns); sensitivity analysis added (nhanes_sensitivity.R; heavy alcohol use, diabetes, lipid-lowering medication); ALT quartile analysis exported.
- Specificity: run_10_specificity_controls.py added (code/02_target_merge); inputs in data/specificity/; disease_in_universe field added to specificity_controls.json.
- Table S6 repacked from the corrected hub-score file (official EPC export; superseded fixed-threshold variant kept as a labelled column).
- GSE126848 panels relabelled (NW/NAFL/NASH; formula titles removed); hub-panel labels decluttered; signature genes highlighted in the feature-selection figure.
- requirements.txt pinned to the exact tested versions; LICENSE-DATA (CC BY 4.0) added; per-sample signature scores deposited under data/cohort/.


## Round-3 revision (2026-09-20, v1.6)

- Authors: the Author Contributions statement, CITATION.cff and LICENSE were updated to match the final author list.
- Section 2.1: TargetNet/PharmMapper counts now distinguish raw tool outputs (623/297) from
  post-harmonisation retention in the compound-side pool (495/284); verified bijectively against
  TargetNet_targets_symbols.tsv / PharmMapper_targets_symbols.tsv / PFOA_target_pool_v1.tsv.
- Section 2.13: one sentence added noting the analytic hypergeometric expectation (full disease
  pool) is a conservative upper bound relative to the array-universe-conditioned permutation null.
- Table S15 (new): expression and provenance of the 15 consensus hub genes and 5 signature genes
  (GSE89632 limma; GSE126848 Welch t + BH). Added to Supplementary docx, datasets zip (CSV) and
  back matter; hub/signature DE evidence paragraph added at the end of Section 3.3.
- Discussion: 'inflammatory backbone' reframed as 'inflammatory module' with explicit note that
  IL6/IL1B were down-regulated in NASH liver and not significant in the independent cohort;
  framed as compound-side-annotated network convergence rather than direct induction evidence.
- Limitations: hub prioritisation is topological, not a set of disease-response genes (Table S15).
- Abstract: 'inflammatory arm' -> 'inflammatory module' (word count unchanged, 198).
- Figure 4 regenerated with collision-free greedy label placement (all six panels verified
  overlap-free by rendered-bbox geometry check).
- Table S14 sensitivity analysis re-verified end-to-end from the six source XPT files
  (ALQ/DIQ/RXQ_RX _I/_J); output byte-identical to TableS14_nhanes_sensitivity.csv.

## Round-4 revision (2026-09-21, v1.8)

- Final author list and corresponding-author contact applied across CITATION.cff, LICENSE and README.
- Repository documentation consolidated; source comments and console output standardized to English.
