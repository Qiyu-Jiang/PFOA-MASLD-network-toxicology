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
## NHANES anchor

Survey-weighted linear models (R `survey`; two-cycle PFAS subsample weights, strata/PSU design) of log-transformed ALT/AST/GGT and FLI on log2 serum PFOA (sum of linear and branched isomers), adjusted for age, sex, race/ethnicity, BMI, cotinine and cycle (NHANES 2015-2018; n = 3,212 adults). Outputs in `data/nhanes/`.

## cytoHubba EPC note

A naive fixed-threshold percolation variant of EPC degenerates on this dense network (mass ties) and is excluded; the final analysis uses the official cytoHubba plugin EPC, reproduced with an adaptive-threshold implementation (`run_09`; Spearman rho = 0.9965; top-10 overlap 8/10 within Monte-Carlo noise). Full score exports are provided in `data/`.

## Changelog

- NHANES anchor reworked to survey-weighted models with FLI on the point scale, per-model n, ALT-quartile and covariate sensitivity analyses (`data/nhanes/`).
- Specificity controls script (`run_10_specificity_controls.py`) and inputs (`data/specificity/`) added; hypergeometric expectation documented as a conservative upper bound.
- Hub-score export repacked from the corrected official-EPC file (superseded fixed-threshold variant kept as a labelled column).
- TargetNet/PharmMapper counts documented to distinguish raw tool outputs from post-harmonisation retention.
- Table S15 (hub/signature expression and provenance) added to the supplementary package.
- Author list, CITATION.cff, LICENSE and README aligned to the final version; requirements.txt pinned to tested versions; LICENSE-DATA (CC BY 4.0) added.
- Repository documentation consolidated; source comments and console output standardized to English.
- Docking note clarified for MMP9 (1GKC): the catalytic zinc is removed during receptor preparation (metal-site limitation); docking summary wording matches the manuscript.
- Figures regenerated: collision-free Figure 4 label placement; Figure 13 with explicit A/B/C panels and corrected axis symbols; Figure 6 colour scale switched to -log10(p.adjust); Figure 2 caption wording aligned.
