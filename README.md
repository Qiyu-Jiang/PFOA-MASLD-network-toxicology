# PFOA-MASLD network toxicology: integrated multi-omics, machine learning, docking and MD analyses

Code and processed data accompanying the manuscript:

> *Integrated Network Toxicology, Machine Learning and Molecular Docking Prioritize Candidate Targets Linking
> PFOA Exposure to Non-Alcoholic Fatty Liver Disease: A Computational Study* (submitted to **Toxics**).

The study integrates (i) network toxicology (PFOA targets from CTD, SwissTargetPrediction, TargetNet, PharmMapper and SEA;
disease genes from the GSE89632 discovery cohort and the DisGeNET, GeneCards and CTD databases), (ii) a cross-validated
six-algorithm cytoHubba consensus on a STRING PPI network, (iii) functional enrichment, (iv) a three-selector
machine-learning signature with external testing and nested cross-validation, (v) molecular docking of PFOA against
eight selected network targets, (vi) 100-ns all-atom molecular dynamics (MD, implicit solvent) with MM-GBSA estimates for
the three strongest complexes, (vii) independent public-data validation in a liver-biopsy cohort (GSE126848) and
single-cell data (GSE344087; 22,674 cells), and (viii) a human exposure anchor in NHANES 2015-2018 (serum PFOA vs
liver enzymes; see `code/08_public_validation/nhanes_analysis.R`).

## Repository layout

```
code/        analysis scripts in execution order (R + Python)
data/        processed result tables (see below)
md_inputs/   MD system inputs (OpenMM system.xml, minimized states, start structures)
docs/        method notes and repository documentation
```

## Key processed data

| File | Content |
|---|---|
| `data/intersection_genes_479.txt` | 479 shared PFOA-NAFLD genes |
| `data/hub_scores_all_nodes.csv` | per-node scores for the six algorithms (official EPC export included; the superseded fixed-threshold EPC variant is kept as a clearly labelled column) |
| `data/hub_consensus_top10_six_algorithms.csv` | top-10 per algorithm (official plugin + custom implementation check) |
| `data/cytohubba_official_scores.csv` | full official Cytoscape cytoHubba export (457 nodes x 12 methods) |
| `data/cytohubba_epc_replication.csv`, `data/cytohubba_epc_comparison_meta.json` | EPC adaptive-threshold replication vs official export (rho = 0.9965; top-10 overlap 8/10) |
| `data/five_algorithm_core_targets.csv` | five-algorithm robust core (strict/extended) |
| `data/ml/` | feature matrices, five-gene signature, selector lists, nested-CV and ROC data |
| `data/cohort/gse126848_score_per_sample.csv` | per-sample five-gene signature scores (GSE126848) |
| `data/specificity_controls.json` | hypergeometric / permutation / stratified specificity analyses (script: `code/02_target_merge/run_10_specificity_controls.py`; inputs in `data/specificity/`) |
| `data/nhanes/` | NHANES 2015-2018 analysis dataset, results, ALT-quartile and covariate-sensitivity analyses |
| `data/docking/redocking_validation.csv` | native-ligand redocking validation (pose-recovery heavy-atom RMSD for 8 co-crystallized complexes; Table S17) |
| `data/deg/`, `data/md/`, `data/mmgbsa/`, `data/docking/`, `data/enrichment/` | differential expression, MD time series, MM-GBSA, docking and enrichment outputs |

## Reproduction

1. **Environment.** Python >= 3.10 with the packages in `requirements.txt`; R >= 4.3 with `limma`, `clusterProfiler`,
   `org.Hs.eg.db`, `ggplot2`, `pROC`, `pheatmap`, `survey`. MD additionally requires `openmm`, `openmmforcefields`,
   `openff-toolkit`, `pdbfixer`, `mdtraj`.
2. **Pipeline order.** `code/01` -> `code/02` -> `code/03` -> `code/04` -> `code/05` -> `code/06` -> `code/07`;
   public-data validation and the NHANES anchor in `code/08`; MD in `code/09` (systems built with `md_build2.py`;
   inputs in `md_inputs/`; analyses with `md_analyze_final.py` and `mmgbsa.py`).
3. **Inputs.** Raw datasets are public: GEO accessions GSE89632, GSE48452, GSE126848, GSE344087; CDC NHANES 2015-2018
   (https://www.cdc.gov/nchs/nhanes/); STRING v12; PDB entries 3O96, 1K7L, 2PRG and the other docked receptors.

## Software versions

See `docs/METHOD_NOTES.md` for software/pipeline parameters. Scripts retain the original project-root variables
(e.g. `BASE`, `base`, `RUNS`); update them to your local paths before rerunning.

## License

Code: MIT (see `LICENSE`). Processed data: CC BY 4.0 (see `LICENSE-DATA`).

## Contact

Qiyu Jiang (corresponding author), email: jiangqiyu@hbmu.edu.cn
