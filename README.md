# GSE188872 source-paper-matched placental analysis

This repository contains a focused reanalysis of first-trimester human fetal placenta bulk RNA-seq from [GEO GSE188872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188872). It retains the complete input/QC record and one inferential workflow: the analysis matched as closely as possible to Figure 10 of the source paper.

## Retained analysis

The paper-matched workflow uses:

- the deposited raw integer-count columns;
- the authors' reported sex labels: 18 female and 16 male placentas;
- an unadjusted DESeq2 model with male minus female as the effect;
- nominal gene p < 0.05 only;
- the paper-style absolute log2 fold-change threshold of 0.6 where explicitly stated.

It recovers 41 of the 50 genes displayed in Figure 10 at nominal p < 0.05 and all five genes named in the paper text. The general GO reconstruction is retained with its full numerical table. Inflammasome GO testing is retained in both strict paper-rule and nominal-only forms; no inflammasome term reaches nominal one-sided Fisher p < 0.05. The circular figures therefore show nominally significant genes within selected terms, not significant term-level enrichment.

An inflammasome-only preranked GSEA uses the complete DESeq2 Wald-statistic ranking. Of 41 predefined inflammasome GO terms, 11 contain at least five ranked genes and are tested. Three reach nominal GSEA p < 0.05, all toward the female-higher end of the ranking: NLRP3 inflammasome complex, positive regulation of NLRP3 inflammasome complex assembly, and canonical inflammasome complex. These related GO terms overlap substantially and are not independent discoveries. No BH threshold is used.

## QC retained

The QC record includes count-column validation, metadata/ID mapping, reported-versus-expression sex concordance, library size and normalization summaries, mean-variance diagnostics, PCA, and the H-19760 sensitivity audit. Four samples reported as male have female-like expression, but reported sex is deliberately used in the final analysis to match the paper.

## Provenance and raw-data policy

Raw GEO payloads are fetched into `data/raw/`, which is gitignored because GEO is authoritative and the files are re-fetchable. [`R/00_fetch_data.R`](R/00_fetch_data.R) downloads missing files and fails on any mismatch against [`data/checksums.md5`](data/checksums.md5). The supplied TPM archive is retained only as an input/provenance check; the final analysis uses raw counts.

The author metadata come from [`shiny_app/human_seq_metadata.csv`](https://github.com/bendevlin18/human-fetal-RNASeq/blob/main/shiny_app/human_seq_metadata.csv). The fixed inflammasome GO term and membership snapshots used by the paper-matched analysis are stored under `data/metadata/`.

## Reproduction and R isolation

All writable R state is confined to this repository. [`scripts/run_r_local.ps1`](scripts/run_r_local.ps1) refuses to run outside the fixed project path and directs `R_USER`, the package library, cache, temporary files, and `renv` state into ignored folders under the project root.

Run the workflow in order from PowerShell:

```powershell
.\scripts\run_r_local.ps1 scripts\bootstrap_renv.R
.\scripts\run_r_local.ps1 R\00_fetch_data.R
.\scripts\run_r_local.ps1 R\01_qc_sex_assignment.R
.\scripts\run_r_local.ps1 R\02_filter_normalize.R
.\scripts\run_r_local.ps1 R\03_paper_matched_reanalysis.R
.\scripts\run_r_local.ps1 R\04_inflammasome_gsea.R
```

## Repository map

- `notebook/` — QC records and the paper-matched analysis record
- `R/` — data acquisition, QC, and the single retained analysis script
- `data/metadata/` — committed metadata and fixed GO membership inputs
- `results/` — QC tables and full-precision paper-matched numerical results
- `figures/` — QC figures and paper-matched figures
- `logs/` — stage-specific execution and R session information
