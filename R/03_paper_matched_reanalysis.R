options(stringsAsFactors = FALSE)
set.seed(42)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "utils", "theme_mobile.R"))
source(file.path(root, "R", "utils", "log_session.R"))
assert_project_local_r(root)

required_packages <- c("AnnotationDbi", "DESeq2", "GO.db", "ggplot2", "org.Hs.eg.db")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing project-local packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

paper50 <- c(
  "KRT5", "GABRP", "MUC16", "GRID1", "CACNB4", "ENOX1", "TMC3", "CPZ", "FCN1", "LYZ",
  "RN7SL4P", "FTH1", "FTH1P8", "HBZ", "SLC4A1", "EPB42", "NFE2", "MYL4", "KLF1", "HBM",
  "AHSP", "C17orf99", "SPTA1", "ANK1", "GATA1", "TFR2", "MT1G", "HMGCS2", "APOM", "LIPC",
  "HAL", "F2", "PLG", "AKR1C2", "CPS1", "PIEZO2", "RHCE", "SLC38A5", "DDX3P1", "PCDH11Y",
  "NLGN4Y", "DDX3Y", "PRKY", "IGFBP1", "RBP1", "CFH", "APOC1", "HAMP", "MT2A", "TXN"
)
paper_named <- c("GATA1", "APOM", "CPS1", "CFH", "KRT5")

input_paths <- c(
  counts = file.path(root, "data", "raw", "GSE188872_placenta_raw_counts.csv.gz"),
  metadata = file.path(root, "results", "sex_concordance.tsv"),
  inflammasome_terms = file.path(root, "data", "metadata", "inflammasome_go_terms.tsv"),
  inflammasome_membership = file.path(
    root,
    "data",
    "metadata",
    "inflammasome_go_term_gene_membership.tsv"
  )
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing required input(s): ", paste(basename(missing_inputs), collapse = ", "), call. = FALSE)
}

raw <- read.csv(input_paths[["counts"]], check.names = FALSE)
x_cols <- grep("\\.x$", names(raw), value = TRUE)
if (!length(x_cols)) stop("No integer .x count columns were found.", call. = FALSE)
counts <- as.matrix(raw[, x_cols, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- as.character(raw$gene)
if (anyDuplicated(rownames(counts))) stop("Count-matrix gene symbols are not unique.", call. = FALSE)
if (any(counts < 0) || any(abs(counts - round(counts)) > 1e-8)) {
  stop("Paper reproduction requires the non-negative integer .x columns.", call. = FALSE)
}

count_columns <- sub("\\.x$", "", x_cols)
metadata <- read.delim(input_paths[["metadata"]], check.names = FALSE)
metadata <- metadata[match(count_columns, metadata$count_column), , drop = FALSE]
if (anyNA(metadata$sample_id)) stop("Count columns did not map to Stage 01 metadata.", call. = FALSE)
if (!identical(metadata$count_column, count_columns)) stop("Sample order changed during mapping.", call. = FALSE)
metadata$reported_sex <- factor(metadata$annotated_sex, levels = c("Female", "Male"))
if (!identical(as.integer(table(metadata$reported_sex)), c(18L, 16L))) {
  stop("The paper-reproduction cohort is not 18 reported females and 16 reported males.", call. = FALSE)
}
colnames(counts) <- metadata$sample_id

sample_data <- data.frame(row.names = metadata$sample_id, sex = metadata$reported_sex)
dds <- DESeq2::DESeqDataSetFromMatrix(counts, sample_data, ~ sex)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
deseq <- as.data.frame(DESeq2::results(dds, contrast = c("sex", "Male", "Female")))
deseq$gene <- rownames(deseq)

chromosome <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = deseq$gene,
  column = "CHR",
  keytype = "SYMBOL",
  multiVals = "first"
)
chromosome <- unname(chromosome[deseq$gene])

deg <- data.frame(
  gene = deseq$gene,
  chromosome = chromosome,
  sex_chromosome = chromosome %in% c("X", "Y"),
  base_mean = deseq$baseMean,
  log2_fold_change_male_minus_female = deseq$log2FoldChange,
  standard_error = deseq$lfcSE,
  wald_statistic = deseq$stat,
  nominal_p = deseq$pvalue,
  stringsAsFactors = FALSE
)
deg$nominal_p_lt_0_05 <- !is.na(deg$nominal_p) & deg$nominal_p < 0.05
deg$paper_deg_rule <- deg$nominal_p_lt_0_05 &
  !is.na(deg$log2_fold_change_male_minus_female) &
  abs(deg$log2_fold_change_male_minus_female) > 0.6
deg$direction <- ifelse(
  deg$log2_fold_change_male_minus_female > 0,
  "Male higher",
  "Female higher"
)
deg <- deg[order(deg$nominal_p, deg$gene, na.last = TRUE), , drop = FALSE]
write_tsv(deg, file.path(root, "results", "paper_reproduction_deg.tsv"))
write_tsv(
  deg[deg$paper_deg_rule, , drop = FALSE],
  file.path(root, "results", "paper_reproduction_deg_nominal_p_lfc.tsv")
)

## QC sensitivity audit for the source-paper genes named in the article text.
without_h19760 <- metadata$sample_id != "H-19760"
sample_data_no_h19760 <- data.frame(
  row.names = metadata$sample_id[without_h19760],
  sex = droplevels(metadata$reported_sex[without_h19760])
)
dds_no_h19760 <- DESeq2::DESeqDataSetFromMatrix(
  counts[, without_h19760, drop = FALSE],
  sample_data_no_h19760,
  ~ sex
)
dds_no_h19760 <- DESeq2::DESeq(dds_no_h19760, quiet = TRUE)
deseq_no_h19760 <- as.data.frame(DESeq2::results(
  dds_no_h19760,
  contrast = c("sex", "Male", "Female")
))
deseq_no_h19760$gene <- rownames(deseq_no_h19760)
named_audit <- deg[match(paper_named, deg$gene), c(
  "gene", "log2_fold_change_male_minus_female", "nominal_p"
), drop = FALSE]
named_without_h19760 <- deseq_no_h19760[match(paper_named, deseq_no_h19760$gene), c(
  "gene", "log2FoldChange", "pvalue"
), drop = FALSE]
names(named_without_h19760) <- c(
  "gene",
  "log2_fold_change_male_minus_female_excluding_H19760",
  "nominal_p_excluding_H19760"
)
named_audit <- merge(named_audit, named_without_h19760, by = "gene", sort = FALSE)
named_audit <- named_audit[match(paper_named, named_audit$gene), , drop = FALSE]
write_tsv(
  named_audit,
  file.path(root, "results", "paper_reproduction_named_genes_H19760_audit.tsv")
)

paper50_stats <- deg[match(paper50, deg$gene), , drop = FALSE]
paper50_stats$paper_heatmap_order <- seq_along(paper50)
paper50_stats <- paper50_stats[, c(
  "paper_heatmap_order", setdiff(names(paper50_stats), "paper_heatmap_order")
)]
write_tsv(paper50_stats, file.path(root, "results", "paper_reproduction_figure10_gene_stats.tsv"))

## Figure 10A-style heatmap using the paper's displayed 50 genes.
vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
vst_matrix <- SummarizedExperiment::assay(vsd)
heat <- vst_matrix[paper50, , drop = FALSE]
heat <- t(scale(t(heat)))
heat[!is.finite(heat)] <- 0
heat <- pmax(pmin(heat, 2.5), -2.5)

order_within_sex <- function(indices) {
  if (length(indices) < 2L) return(indices)
  indices[stats::hclust(stats::dist(t(heat[, indices, drop = FALSE])))$order]
}
female_idx <- which(metadata$reported_sex == "Female")
male_idx <- which(metadata$reported_sex == "Male")
sample_order <- c(order_within_sex(female_idx), order_within_sex(male_idx))

row_order <- stats::hclust(stats::dist(heat))$order
gene_order <- rownames(heat)[row_order]
heat_long <- data.frame(
  gene = rep(rownames(heat), times = ncol(heat)),
  sample_id = rep(colnames(heat), each = nrow(heat)),
  z_score = as.vector(heat),
  stringsAsFactors = FALSE
)
heat_long$x <- match(heat_long$sample_id, metadata$sample_id[sample_order])
heat_long$y <- match(heat_long$gene, rev(gene_order))
sex_bar <- data.frame(
  x = seq_along(sample_order),
  reported_sex = metadata$reported_sex[sample_order],
  stringsAsFactors = FALSE
)

p_heat <- ggplot2::ggplot(heat_long, ggplot2::aes(x, y, fill = z_score)) +
  ggplot2::geom_tile() +
  ggplot2::geom_tile(
    data = sex_bar[sex_bar$reported_sex == "Female", , drop = FALSE],
    ggplot2::aes(x = x, y = 51.4),
    inherit.aes = FALSE,
    width = 1,
    height = 0.8,
    fill = okabe_ito[["sky_blue"]]
  ) +
  ggplot2::geom_tile(
    data = sex_bar[sex_bar$reported_sex == "Male", , drop = FALSE],
    ggplot2::aes(x = x, y = 51.4),
    inherit.aes = FALSE,
    width = 1,
    height = 0.8,
    fill = okabe_ito[["vermillion"]]
  ) +
  ggplot2::annotate(
    "text",
    x = c(mean(seq_along(female_idx)), length(female_idx) + mean(seq_along(male_idx))),
    y = 52.2,
    label = c("Reported female", "Reported male"),
    fontface = "bold",
    size = 3.8
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    name = "Row z-score"
  ) +
  ggplot2::scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  ggplot2::scale_y_continuous(
    breaks = seq_along(gene_order),
    labels = rev(gene_order),
    limits = c(0.5, 52.7),
    expand = c(0, 0),
    position = "right"
  ) +
  ggplot2::labs(
    title = "Source-paper reproduction: 50 reported sex-associated genes",
    subtitle = "18 reported female and 16 reported male placentas; VST expression scaled within gene",
    x = NULL,
    y = NULL,
    caption = "Gene set is transcribed from Figure 10A of the source paper; samples are clustered within reported sex."
  ) +
  theme_mobile(11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 7.8),
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.margin = ggplot2::margin(12, 20, 12, 12)
  )
ggplot2::ggsave(
  file.path(root, "figures", "paper_01_figure10_heatmap.png"),
  p_heat,
  width = 8,
  height = 10,
  dpi = 220,
  units = "in",
  bg = "white"
)

tested <- deg[!is.na(deg$nominal_p), , drop = FALSE]
universe <- tested$gene
selected_paper_rule <- tested$gene[tested$paper_deg_rule]
selected_nominal_only <- tested$gene[tested$nominal_p_lt_0_05]

## Genome-wide biological-process ORA using propagated GO annotations.
go_columns <- AnnotationDbi::columns(org.Hs.eg.db::org.Hs.eg.db)
if (!all(c("GOALL", "ONTOLOGYALL") %in% go_columns)) {
  stop("org.Hs.eg.db does not expose propagated GOALL annotations.", call. = FALSE)
}
go_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = universe,
  keytype = "SYMBOL",
  columns = c("GOALL", "ONTOLOGYALL")
))
go_map <- unique(go_map[
  !is.na(go_map$GOALL) & go_map$ONTOLOGYALL == "BP",
  c("SYMBOL", "GOALL"),
  drop = FALSE
])
names(go_map) <- c("gene", "go_id")
go_map <- go_map[go_map$gene %in% universe, , drop = FALSE]
term_names <- AnnotationDbi::Term(GO.db::GOTERM)

ora_from_map <- function(
  map,
  term_ids,
  scope_label,
  selected_genes,
  minimum_term_size = 2L
) {
  rows <- lapply(term_ids, function(go_id) {
    members <- unique(map$gene[map$go_id == go_id])
    members <- intersect(members, universe)
    hits <- intersect(members, selected_genes)
    a <- length(hits)
    n <- length(members)
    k <- length(selected_genes)
    total <- length(universe)
    expected <- n * k / total
    fisher <- if (n >= minimum_term_size && a > 0 && k > 0) {
      stats::fisher.test(
        matrix(c(a, n - a, k - a, total - n - k + a), nrow = 2, byrow = TRUE),
        alternative = "greater"
      )
    } else {
      NULL
    }
    hit_rows <- deg[match(hits, deg$gene), , drop = FALSE]
    data.frame(
      analysis_scope = scope_label,
      go_id = go_id,
      ontology = "BP",
      term_name = unname(term_names[[go_id]]),
      universe_genes = total,
      selected_deg_genes = k,
      measured_term_genes = n,
      selected_deg_in_term = a,
      expected_selected_in_term = expected,
      fold_enrichment = ifelse(expected > 0, a / expected, NA_real_),
      fisher_nominal_p = if (is.null(fisher)) NA_real_ else fisher$p.value,
      male_higher_genes = sum(hit_rows$log2_fold_change_male_minus_female > 0, na.rm = TRUE),
      female_higher_genes = sum(hit_rows$log2_fold_change_male_minus_female < 0, na.rm = TRUE),
      selected_gene_symbols = paste(sort(hits), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$nominally_enriched_p_lt_0_05 <- !is.na(out$fisher_nominal_p) & out$fisher_nominal_p < 0.05
  out[order(out$fisher_nominal_p, -out$selected_deg_in_term, out$term_name), , drop = FALSE]
}

general_go <- ora_from_map(
  go_map,
  sort(unique(go_map$go_id)),
  "PAPER_REPRODUCTION_GENOMEWIDE",
  selected_genes = selected_paper_rule,
  minimum_term_size = 5L
)
write_tsv(general_go, file.path(root, "results", "paper_reproduction_general_go_nominal_ora.tsv"))

source_terms <- c("immune system process", "humoral immune response", "immune response")
source_term_audit <- general_go[tolower(general_go$term_name) %in% source_terms, , drop = FALSE]
source_term_audit$rank_by_nominal_p <- match(
  source_term_audit$go_id,
  general_go$go_id[!is.na(general_go$fisher_nominal_p)]
)
write_tsv(source_term_audit, file.path(root, "results", "paper_reproduction_source_go_term_audit.tsv"))

general_plot <- general_go[
  general_go$measured_term_genes >= 5 & general_go$selected_deg_in_term > 0 &
    general_go$nominally_enriched_p_lt_0_05,
  ,
  drop = FALSE
]
general_plot <- general_plot[seq_len(min(15L, nrow(general_plot))), , drop = FALSE]
if (!nrow(general_plot)) stop("No nominally enriched general BP terms were found.", call. = FALSE)
general_plot$dominant_direction <- ifelse(
  general_plot$male_higher_genes > general_plot$female_higher_genes,
  "More male-higher genes",
  ifelse(
    general_plot$female_higher_genes > general_plot$male_higher_genes,
    "More female-higher genes",
    "Direction tied"
  )
)
general_plot$term_label <- vapply(
  general_plot$term_name,
  function(x) paste(strwrap(x, width = 38), collapse = "\n"),
  character(1)
)
general_plot$term_label <- factor(general_plot$term_label, levels = rev(general_plot$term_label))
general_plot$p_label <- paste0(
  general_plot$selected_deg_in_term,
  " genes; P=",
  formatC(general_plot$fisher_nominal_p, digits = 2, format = "g")
)

direction_colours <- c(
  "More male-higher genes" = okabe_ito[["vermillion"]],
  "More female-higher genes" = okabe_ito[["blue"]],
  "Direction tied" = okabe_ito[["bluish_green"]]
)
p_general <- ggplot2::ggplot(
  general_plot,
  ggplot2::aes(x = -log10(fisher_nominal_p), y = term_label)
) +
  ggplot2::geom_col(ggplot2::aes(fill = dominant_direction), width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = p_label),
    hjust = -0.08,
    size = 3.6
  ) +
  ggplot2::scale_fill_manual(values = direction_colours, name = NULL) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.35))) +
  ggplot2::labs(
    title = "Paper-reproduction biological-process enrichment",
    subtitle = paste0(
      length(selected_paper_rule),
      " DEGs: nominal gene P < 0.05 and |log2FC| > 0.6; top 15 nominal GO terms"
    ),
    x = expression(-log[10]("one-sided Fisher P")),
    y = NULL,
    caption = "No BH threshold is used. Full nominal results are retained in the accompanying table."
  ) +
  theme_mobile(12) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 10, lineheight = 0.9),
    legend.position = "bottom",
    panel.grid.major.y = ggplot2::element_blank(),
    plot.title.position = "plot",
    plot.caption.position = "plot"
  )
ggplot2::ggsave(
  file.path(root, "figures", "paper_02_general_go.png"),
  p_general,
  width = 8.2,
  height = 8.2,
  dpi = 200,
  units = "in",
  bg = "white"
)

## Inflammasome-focused ORA, reusing the audited Stage 04 GO branch.
inflammasome_terms <- read.delim(input_paths[["inflammasome_terms"]], check.names = FALSE)
inflammasome_membership <- read.delim(
  input_paths[["inflammasome_membership"]],
  check.names = FALSE
)
inflammasome_membership <- unique(inflammasome_membership[
  inflammasome_membership$retained_filterByExpr &
    inflammasome_membership$gene %in% universe,
  c("go_id", "gene"),
  drop = FALSE
])
inflammasome_go <- ora_from_map(
  inflammasome_membership,
  inflammasome_terms$go_id,
  "PAPER_REPRODUCTION_INFLAMMASOME_P_LFC",
  selected_genes = selected_paper_rule,
  minimum_term_size = 2L
)
inflammasome_go$branch_scope <- ifelse(
  inflammasome_terms$in_core_branch[match(inflammasome_go$go_id, inflammasome_terms$go_id)] &
    inflammasome_terms$in_downstream_branch[match(inflammasome_go$go_id, inflammasome_terms$go_id)],
  "Reached from both search tiers",
  ifelse(
    inflammasome_terms$in_core_branch[match(inflammasome_go$go_id, inflammasome_terms$go_id)],
    "Assembly / activation core",
    "Downstream cytokine tier"
  )
)
write_tsv(
  inflammasome_go,
  file.path(root, "results", "paper_reproduction_inflammasome_go_p_lfc_ora.tsv")
)

inflammasome_nominal_only <- ora_from_map(
  inflammasome_membership,
  inflammasome_terms$go_id,
  "PAPER_REPRODUCTION_INFLAMMASOME_NOMINAL_ONLY",
  selected_genes = selected_nominal_only,
  minimum_term_size = 2L
)
inflammasome_nominal_only$branch_scope <- inflammasome_go$branch_scope[
  match(inflammasome_nominal_only$go_id, inflammasome_go$go_id)
]
write_tsv(
  inflammasome_nominal_only,
  file.path(root, "results", "paper_reproduction_inflammasome_go_nominal_only_ora.tsv")
)

inflammasome_plot <- inflammasome_nominal_only[
  inflammasome_nominal_only$selected_deg_in_term > 0,
  ,
  drop = FALSE
]
inflammasome_enriched <- inflammasome_plot[
  inflammasome_plot$nominally_enriched_p_lt_0_05,
  ,
  drop = FALSE
]
if (nrow(inflammasome_enriched)) {
  inflammasome_plot <- inflammasome_enriched[seq_len(min(15L, nrow(inflammasome_enriched))), , drop = FALSE]
  inflammasome_title <- "Paper-reproduction inflammasome GO enrichment"
  inflammasome_subtitle <- paste0(
    nrow(inflammasome_enriched),
    " term(s) with nominal one-sided Fisher P < 0.05"
  )
} else {
  inflammasome_plot <- inflammasome_plot[seq_len(min(15L, nrow(inflammasome_plot))), , drop = FALSE]
  inflammasome_title <- "No nominal inflammasome GO enrichment"
  inflammasome_subtitle <- "Top overlaps shown; all one-sided Fisher P values are at least 0.05"
}
inflammasome_plot$dominant_direction <- ifelse(
  inflammasome_plot$male_higher_genes > inflammasome_plot$female_higher_genes,
  "More male-higher genes",
  ifelse(
    inflammasome_plot$female_higher_genes > inflammasome_plot$male_higher_genes,
    "More female-higher genes",
    "Direction tied"
  )
)
inflammasome_plot$term_label <- vapply(
  paste0(inflammasome_plot$term_name, " (", inflammasome_plot$go_id, ")"),
  function(x) paste(strwrap(x, width = 42), collapse = "\n"),
  character(1)
)
inflammasome_plot$term_label <- factor(
  inflammasome_plot$term_label,
  levels = rev(inflammasome_plot$term_label)
)
inflammasome_plot$point_label <- paste0(
  inflammasome_plot$selected_deg_in_term,
  "/",
  inflammasome_plot$measured_term_genes,
  "; P=",
  formatC(inflammasome_plot$fisher_nominal_p, digits = 2, format = "g")
)

p_inflammasome <- ggplot2::ggplot(
  inflammasome_plot,
  ggplot2::aes(x = fold_enrichment, y = term_label)
) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_segment(
    ggplot2::aes(x = 1, xend = fold_enrichment, yend = term_label),
    colour = "grey70",
    linewidth = 1
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = dominant_direction, size = selected_deg_in_term),
    shape = 21,
    colour = "black",
    stroke = 0.8
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = point_label),
    nudge_x = 0.12,
    hjust = 0,
    size = 3.7,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(values = direction_colours, name = NULL) +
  ggplot2::scale_size_continuous(range = c(4, 9), guide = "none") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.3))) +
  ggplot2::labs(
    title = paste0(inflammasome_title, ": nominal-P-only genes"),
    subtitle = inflammasome_subtitle,
    x = "Fold enrichment of nominal-P DEGs",
    y = NULL,
    caption = paste(
      strwrap(
        "Gene selection is nominal DESeq2 P < 0.05 with no fold-change or BH threshold. Labels show selected DEGs / measured term genes and nominal one-sided Fisher P.",
        width = 72
      ),
      collapse = "\n"
    )
  ) +
  theme_mobile(12) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 10.5, lineheight = 0.9),
    legend.position = "bottom",
    panel.grid.major.y = ggplot2::element_blank(),
    plot.title.position = "plot",
    plot.caption.position = "plot"
  )
ggplot2::ggsave(
  file.path(root, "figures", "paper_03_inflammasome_go.png"),
  p_inflammasome,
  width = 8,
  height = 7.8,
  dpi = 200,
  units = "in",
  bg = "white"
)

circle_pool <- if (nrow(inflammasome_enriched)) {
  inflammasome_enriched
} else {
  inflammasome_nominal_only
}
circle_pool <- circle_pool[circle_pool$selected_deg_in_term >= 2, , drop = FALSE]
circle_pool <- circle_pool[order(
  -circle_pool$selected_deg_in_term,
  circle_pool$fisher_nominal_p,
  circle_pool$term_name
), , drop = FALSE]
circle_terms <- circle_pool[seq_len(min(3L, nrow(circle_pool))), , drop = FALSE]
circle_rows <- list()
circle_prefixes <- c("paper_04", "paper_05", "paper_06")[seq_len(nrow(circle_terms))]
if (nrow(circle_terms)) {
  circle_gene_union <- unique(unlist(strsplit(
    circle_terms$selected_gene_symbols,
    "; ",
    fixed = TRUE
  )))
  circle_gene_union <- circle_gene_union[nzchar(circle_gene_union)]
  shared_limit <- max(abs(deg$log2_fold_change_male_minus_female[
    deg$gene %in% circle_gene_union
  ]), na.rm = TRUE)
  if (!is.finite(shared_limit) || shared_limit == 0) shared_limit <- 1

  for (i in seq_len(nrow(circle_terms))) {
    term <- circle_terms[i, , drop = FALSE]
    genes <- strsplit(term$selected_gene_symbols, "; ", fixed = TRUE)[[1]]
    genes <- genes[nzchar(genes)]
    plot_data <- deg[match(genes, deg$gene), c(
      "gene", "chromosome", "log2_fold_change_male_minus_female", "nominal_p"
    ), drop = FALSE]
    plot_data <- plot_data[order(
      plot_data$log2_fold_change_male_minus_female,
      plot_data$gene
    ), , drop = FALSE]
    plot_data$plot_id <- seq_len(nrow(plot_data))
    plot_data$label <- paste0(
      plot_data$gene,
      "\nP=",
      formatC(plot_data$nominal_p, digits = 2, format = "g")
    )
    plot_data$label_angle <- 90 - 360 * (plot_data$plot_id - 0.5) / nrow(plot_data)
    plot_data$label_hjust <- ifelse(plot_data$label_angle < -90, 1, 0)
    plot_data$label_angle <- ifelse(
      plot_data$label_angle < -90,
      plot_data$label_angle + 180,
      plot_data$label_angle
    )

    p_circle <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = plot_id, y = 1, fill = log2_fold_change_male_minus_female)
    ) +
      ggplot2::geom_col(width = 0.88, colour = "white", linewidth = 1) +
      ggplot2::geom_text(
        ggplot2::aes(
          y = 1.18,
          label = label,
          angle = label_angle,
          hjust = label_hjust
        ),
        size = 4.2,
        lineheight = 0.9,
        colour = "black",
        show.legend = FALSE
      ) +
      ggplot2::coord_polar(theta = "x", start = 0, clip = "off") +
      ggplot2::scale_fill_gradient2(
        low = okabe_ito[["blue"]],
        mid = "white",
        high = okabe_ito[["vermillion"]],
        midpoint = 0,
        limits = c(-shared_limit, shared_limit),
        name = "Male - female\nlog2 fold change"
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0.5, nrow(plot_data) + 0.5),
        breaks = NULL
      ) +
      ggplot2::scale_y_continuous(
        limits = c(-0.85, 1.85),
        breaks = NULL,
        expand = c(0, 0)
      ) +
      ggplot2::labs(
        title = paste(strwrap(term$term_name, width = 44), collapse = "\n"),
        subtitle = paste0(
          term$go_id,
          " | ", term$selected_deg_in_term,
          " nominal-P DEGs among ", term$measured_term_genes,
          " measured members\nNominal one-sided Fisher P = ",
          formatC(term$fisher_nominal_p, digits = 3, format = "g")
        ),
        x = NULL,
        y = NULL,
        caption = "Only nominal DESeq2 P < 0.05 genes are displayed; no fold-change or BH threshold is used."
      ) +
      theme_mobile(12) +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, lineheight = 0.95),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5,
          size = 11,
          lineheight = 1.05,
          margin = ggplot2::margin(b = 56)
        ),
        plot.caption = ggplot2::element_text(size = 10.5),
        plot.margin = ggplot2::margin(24, 60, 18, 60)
      )

    file_name <- paste0(
      circle_prefixes[[i]],
      "_circle_",
      gsub(":", "_", term$go_id),
      ".png"
    )
    ggplot2::ggsave(
      file.path(root, "figures", file_name),
      p_circle,
      width = 8,
      height = 8.5,
      dpi = 200,
      units = "in",
      bg = "white"
    )
    circle_rows[[i]] <- data.frame(
      rank_by_selected_gene_count = i,
      go_id = term$go_id,
      term_name = term$term_name,
      selected_deg_in_term = term$selected_deg_in_term,
      measured_term_genes = term$measured_term_genes,
      fisher_nominal_p = term$fisher_nominal_p,
      nominally_enriched_p_lt_0_05 = term$nominally_enriched_p_lt_0_05,
      selected_gene_symbols = term$selected_gene_symbols,
      figure = file.path("figures", file_name),
      stringsAsFactors = FALSE
    )
  }
}
circle_summary <- if (length(circle_rows)) do.call(rbind, circle_rows) else data.frame()
write_tsv(
  circle_summary,
  file.path(root, "results", "paper_reproduction_inflammasome_circle_terms.tsv")
)

summary <- data.frame(
  metric = c(
    "reported_female_samples", "reported_male_samples", "tested_genes",
    "nominal_p_lt_0_05_genes", "paper_rule_degs", "paper_rule_male_higher",
    "paper_rule_female_higher", "paper_figure10_genes_nominal_p_lt_0_05",
    "general_bp_terms_nominal_p_lt_0_05",
    "inflammasome_terms_p_lfc_nominal_fisher_lt_0_05",
    "inflammasome_terms_nominal_only_fisher_lt_0_05"
  ),
  value = c(
    sum(metadata$reported_sex == "Female"),
    sum(metadata$reported_sex == "Male"),
    nrow(tested),
    sum(tested$nominal_p_lt_0_05),
    length(selected_paper_rule),
    sum(tested$paper_deg_rule & tested$log2_fold_change_male_minus_female > 0),
    sum(tested$paper_deg_rule & tested$log2_fold_change_male_minus_female < 0),
    sum(paper50_stats$nominal_p_lt_0_05, na.rm = TRUE),
    sum(general_go$nominally_enriched_p_lt_0_05, na.rm = TRUE),
    sum(inflammasome_go$nominally_enriched_p_lt_0_05, na.rm = TRUE),
    sum(inflammasome_nominal_only$nominally_enriched_p_lt_0_05, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, file.path(root, "results", "paper_reproduction_summary.tsv"))

cat("\nPAPER-REPRODUCTION SUMMARY\n")
print(summary, row.names = FALSE)
cat("\nSOURCE-PAPER GO TERM AUDIT\n")
print(source_term_audit, row.names = FALSE)
cat("\nTOP GENERAL BP TERMS\n")
print(general_plot[, c(
  "go_id", "term_name", "measured_term_genes", "selected_deg_in_term",
  "fold_enrichment", "fisher_nominal_p", "selected_gene_symbols"
)], row.names = FALSE)
cat("\nINFLAMMASOME GO TERMS WITH NOMINAL FISHER P < 0.05\n")
print(inflammasome_enriched, row.names = FALSE)
cat("\nINFLAMMASOME CIRCLE TERMS\n")
print(circle_summary, row.names = FALSE)

log_session("03", root)
