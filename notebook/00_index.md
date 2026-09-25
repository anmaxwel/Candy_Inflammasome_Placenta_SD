# Electronic lab notebook index

## Entries

| Date | Stage | Entry | Status |
|---|---:|---|---|
| 2026-09-24 | 01 | [ID mapping and transcriptomic sex QC](2026-09-24_01_qc.md) | complete |
| 2026-09-24 | 02 | [Filtering, normalization, and global QC](2026-09-24_02_normalization.md) | complete |
| 2026-09-24 | 03 | [Source-paper-matched reanalysis](2026-09-24_03_paper_matched_reanalysis.md) | complete; awaiting review |
| 2026-09-24 | 04 | [Inflammasome-only preranked GSEA](2026-09-24_04_inflammasome_gsea.md) | complete |

## Current decisions

| Date | Stage | Decision | Status |
|---|---:|---|---|
| 2026-09-24 | 01 | D01 — Use the integer `.x` columns in the GEO raw-count CSV; reject the paired CPM-like `.y` columns. | standing |
| 2026-09-24 | 01 | D02 — Retain expression-based sex assignment as QC and report the four annotation/expression discrepancies. | standing |
| 2026-09-24 | 01 | D03 — Use the authors' reported sex labels in the final analysis because the objective is source-paper reproduction. | standing |
| 2026-09-24 | 01 | D04 — Treat `trig` as already log10-scaled and do not transform it again. | standing |
| 2026-09-24 | 02 | D05 — Retain all 16,399 deposited genes because the expression filter removed none. | standing |
| 2026-09-24 | 02 | D06 — Retain normalization and PCA outputs as QC; they do not define the final DESeq2 model. | standing |
| 2026-09-24 | 03 | D07 — Use an unadjusted reported-sex DESeq2 model as the single retained inferential analysis. | standing |
| 2026-09-24 | 03 | D08 — Use nominal p only; the strict paper-style DEG rule is p < 0.05 and absolute log2 fold change > 0.6. | standing |
| 2026-09-24 | 03 | D09 — Treat recovery of 41/50 Figure 10 genes and all five text-named genes as strong but incomplete reproduction. | standing |
| 2026-09-24 | 03 | D10 — Report nominally significant inflammasome genes, but do not claim term-level enrichment because no inflammasome GO term has Fisher p < 0.05. | standing |
| 2026-09-24 | 04 | D11 — Restrict GSEA testing and figures to the predefined inflammasome GO collection. | standing |
| 2026-09-24 | 04 | D12 — Use the complete paper-matched DESeq2 Wald-statistic ranking and nominal GSEA p < 0.05, without a BH threshold. | standing |
| 2026-09-24 | 04 | D13 — Interpret the three female-higher nominal terms as an overlapping NLRP3/canonical-inflammasome signal rather than three independent discoveries. | standing |
