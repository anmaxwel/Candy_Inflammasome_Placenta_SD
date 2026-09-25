options(stringsAsFactors = FALSE)
set.seed(42)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "utils", "theme_mobile.R"))
source(file.path(root, "R", "utils", "log_session.R"))
assert_project_local_r(root)

required_packages <- c("fgsea", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing project-local packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

format_p <- function(x) {
  ifelse(x < 0.001, formatC(x, digits = 2, format = "e"), formatC(x, digits = 3, format = "f"))
}

input_paths <- c(
  deg = file.path(root, "results", "paper_reproduction_deg.tsv"),
  terms = file.path(root, "data", "metadata", "inflammasome_go_terms.tsv"),
  membership = file.path(root, "data", "metadata", "inflammasome_go_term_gene_membership.tsv")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing required input(s): ", paste(basename(missing_inputs), collapse = ", "), call. = FALSE)
}

deg <- read.delim(input_paths[["deg"]], check.names = FALSE)
required_deg <- c("gene", "wald_statistic", "log2_fold_change_male_minus_female", "nominal_p")
if (!all(required_deg %in% names(deg))) {
  stop("Paper-reproduction DEG table has an unexpected schema.", call. = FALSE)
}
rank_data <- deg[
  !is.na(deg$gene) & nzchar(deg$gene) & is.finite(deg$wald_statistic),
  required_deg,
  drop = FALSE
]
rank_data <- rank_data[!duplicated(rank_data$gene), , drop = FALSE]
rank_data <- rank_data[order(rank_data$wald_statistic, decreasing = TRUE), , drop = FALSE]
ranks <- rank_data$wald_statistic
names(ranks) <- rank_data$gene
if (!all(diff(ranks) <= 0)) stop("GSEA ranking is not decreasing.", call. = FALSE)

terms <- read.delim(input_paths[["terms"]], check.names = FALSE)
membership <- read.delim(input_paths[["membership"]], check.names = FALSE)
required_terms <- c(
  "go_id", "ontology", "term_name", "in_core_branch", "in_downstream_branch", "is_anchor"
)
if (!all(required_terms %in% names(terms))) {
  stop("Inflammasome term table has an unexpected schema.", call. = FALSE)
}
required_membership <- c("go_id", "gene", "retained_filterByExpr")
if (!all(required_membership %in% names(membership))) {
  stop("Inflammasome membership table has an unexpected schema.", call. = FALSE)
}

membership <- unique(membership[
  membership$retained_filterByExpr & membership$gene %in% names(ranks),
  c("go_id", "gene"),
  drop = FALSE
])
pathways_all <- split(membership$gene, membership$go_id)
pathways_all <- lapply(pathways_all, unique)
measured_size <- lengths(pathways_all)

min_size <- 5L
max_size <- 500L
term_audit <- terms[, required_terms, drop = FALSE]
term_audit$measured_ranked_genes <- unname(measured_size[term_audit$go_id])
term_audit$measured_ranked_genes[is.na(term_audit$measured_ranked_genes)] <- 0L
term_audit$gsea_eligible <- term_audit$measured_ranked_genes >= min_size &
  term_audit$measured_ranked_genes <= max_size
term_audit$exclusion_reason <- ifelse(
  term_audit$gsea_eligible,
  "",
  ifelse(term_audit$measured_ranked_genes < min_size, "fewer_than_5_ranked_genes", "more_than_500_ranked_genes")
)
term_audit <- term_audit[order(!term_audit$gsea_eligible, -term_audit$measured_ranked_genes, term_audit$go_id), ]
write_tsv(
  term_audit,
  file.path(root, "results", "paper_reproduction_inflammasome_gsea_term_audit.tsv")
)

eligible_ids <- term_audit$go_id[term_audit$gsea_eligible]
pathways <- pathways_all[intersect(eligible_ids, names(pathways_all))]
if (!length(pathways)) stop("No inflammasome GO terms met the GSEA size boundary.", call. = FALSE)

gsea_raw <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = ranks,
  minSize = min_size,
  maxSize = max_size,
  eps = 0,
  scoreType = "std",
  nproc = 1
)
gsea <- as.data.frame(gsea_raw)
names(gsea)[names(gsea) == "pathway"] <- "go_id"
gsea$leading_edge_genes <- vapply(gsea$leadingEdge, paste, collapse = "; ", character(1))
gsea$leading_edge_count <- lengths(gsea$leadingEdge)
gsea$leadingEdge <- NULL
gsea$term_name <- terms$term_name[match(gsea$go_id, terms$go_id)]
gsea$ontology <- terms$ontology[match(gsea$go_id, terms$go_id)]
gsea$in_core_branch <- terms$in_core_branch[match(gsea$go_id, terms$go_id)]
gsea$in_downstream_branch <- terms$in_downstream_branch[match(gsea$go_id, terms$go_id)]
gsea$is_anchor <- terms$is_anchor[match(gsea$go_id, terms$go_id)]
gsea$direction <- ifelse(gsea$NES > 0, "Male higher", "Female higher")
gsea$nominal_p_lt_0_05 <- gsea$pval < 0.05
names(gsea)[names(gsea) == "pval"] <- "nominal_p"
names(gsea)[names(gsea) == "size"] <- "pathway_size"
gsea$padj <- NULL
gsea <- gsea[order(gsea$nominal_p, -abs(gsea$NES), gsea$go_id), c(
  "go_id", "ontology", "term_name", "ES", "NES", "nominal_p", "log2err", "pathway_size",
  "direction", "nominal_p_lt_0_05", "leading_edge_count", "leading_edge_genes",
  "in_core_branch", "in_downstream_branch", "is_anchor"
)]
write_tsv(
  gsea,
  file.path(root, "results", "paper_reproduction_inflammasome_gsea.tsv")
)

summary <- data.frame(
  metric = c(
    "ranked_genes", "inflammasome_terms_total", "inflammasome_terms_gsea_eligible",
    "inflammasome_terms_nominal_p_lt_0_05", "male_higher_nominal_terms",
    "female_higher_nominal_terms"
  ),
  value = c(
    length(ranks), nrow(term_audit), sum(term_audit$gsea_eligible),
    sum(gsea$nominal_p_lt_0_05),
    sum(gsea$nominal_p_lt_0_05 & gsea$NES > 0),
    sum(gsea$nominal_p_lt_0_05 & gsea$NES < 0)
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  summary,
  file.path(root, "results", "paper_reproduction_inflammasome_gsea_summary.tsv")
)

okabe_ito <- c(
  "Male higher" = "#0072B2",
  "Female higher" = "#D55E00"
)

## Figure 1: inflammasome-only GSEA overview.
display_n <- min(15L, nrow(gsea))
overview <- gsea[seq_len(display_n), , drop = FALSE]
overview$term_label <- paste0(
  vapply(overview$term_name, function(x) paste(strwrap(x, 43), collapse = "\n"), character(1)),
  " (", overview$go_id, ")"
)
overview$term_label <- factor(overview$term_label, levels = rev(overview$term_label))
overview$point_size <- -log10(pmax(overview$nominal_p, .Machine$double.xmin))
overview$p_label <- paste0("P=", format_p(overview$nominal_p))

p_overview <- ggplot2::ggplot(
  overview,
  ggplot2::aes(x = NES, y = term_label, colour = direction)
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey55") +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = NES, yend = term_label),
    linewidth = 0.8,
    colour = "grey72"
  ) +
  ggplot2::geom_point(ggplot2::aes(size = point_size), alpha = 0.95) +
  ggplot2::geom_text(
    ggplot2::aes(label = p_label),
    colour = "black",
    size = 3.6,
    hjust = ifelse(overview$NES >= 0, -0.12, 1.12),
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = okabe_ito) +
  ggplot2::scale_size_continuous(name = expression(-log[10](P)), range = c(3, 8)) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.23, 0.23))) +
  ggplot2::labs(
    title = "Inflammasome-only preranked GSEA",
    subtitle = paste0(
      "Full DESeq2 Wald ranking\n", sum(gsea$nominal_p_lt_0_05),
      " of ", nrow(gsea), " eligible GO terms have nominal P < 0.05"
    ),
    x = "Normalized enrichment score (male minus female ranking)",
    y = NULL,
    colour = NULL,
    caption = "Only predefined inflammasome GO terms are tested and displayed.\nNo BH threshold is used."
  ) +
  theme_mobile(12) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(18, 70, 18, 18)
  )
ggplot2::ggsave(
  file.path(root, "figures", "paper_gsea_01_inflammasome_overview.png"),
  p_overview,
  width = 10.5,
  height = 9.5,
  dpi = 220,
  units = "in",
  bg = "white"
)

## Figure 2: running enrichment curves for the three lowest nominal-P terms.
curve_terms <- gsea$go_id[seq_len(min(3L, nrow(gsea)))]
curve_rows <- list()
hit_rows <- list()
for (i in seq_along(curve_terms)) {
  go_id <- curve_terms[[i]]
  result_row <- gsea[gsea$go_id == go_id, , drop = FALSE]
  genes <- pathways[[go_id]]
  hits <- names(ranks) %in% genes
  hit_weights <- abs(ranks) * hits
  hit_total <- sum(hit_weights)
  miss_total <- sum(!hits)
  running <- cumsum(ifelse(hits, hit_weights / hit_total, -1 / miss_total))
  facet_label <- paste0(
    result_row$term_name, " (", go_id, ")\nNES=",
    formatC(result_row$NES, digits = 3, format = "f"),
    "; nominal P=", format_p(result_row$nominal_p)
  )
  curve_rows[[i]] <- data.frame(
    rank_position = seq_along(ranks),
    running_score = running,
    direction = result_row$direction,
    facet_label = facet_label,
    stringsAsFactors = FALSE
  )
  hit_rows[[i]] <- data.frame(
    rank_position = which(hits),
    facet_label = facet_label,
    stringsAsFactors = FALSE
  )
}
curve_data <- do.call(rbind, curve_rows)
hit_data <- do.call(rbind, hit_rows)

p_curves <- ggplot2::ggplot(
  curve_data,
  ggplot2::aes(x = rank_position, y = running_score, colour = direction)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.45) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_rug(
    data = hit_data,
    ggplot2::aes(x = rank_position),
    inherit.aes = FALSE,
    sides = "b",
    colour = "black",
    alpha = 0.38,
    linewidth = 0.25
  ) +
  ggplot2::facet_wrap(~ facet_label, ncol = 1, scales = "free_y") +
  ggplot2::scale_colour_manual(values = okabe_ito) +
  ggplot2::labs(
    title = "Top inflammasome GSEA running-enrichment curves",
    subtitle = "Terms selected by lowest nominal GSEA P value",
    x = "Position in DESeq2 Wald ranking (male higher to female higher)",
    y = "Running enrichment score",
    colour = NULL,
    caption = "Vertical ticks mark term genes in the full ranked list. Only inflammasome GO terms are shown."
  ) +
  theme_mobile(11.5) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 10.5),
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(16, 18, 16, 18)
  )
ggplot2::ggsave(
  file.path(root, "figures", "paper_gsea_02_inflammasome_curves.png"),
  p_curves,
  width = 9,
  height = 10,
  dpi = 220,
  units = "in",
  bg = "white"
)

## Figure 3: leading-edge genes for the lowest nominal-P inflammasome term.
top_result <- gsea[1, , drop = FALSE]
top_genes <- strsplit(top_result$leading_edge_genes, "; ", fixed = TRUE)[[1]]
leading <- rank_data[match(top_genes, rank_data$gene), , drop = FALSE]
leading <- leading[order(abs(leading$wald_statistic), decreasing = TRUE), , drop = FALSE]
leading <- leading[seq_len(min(20L, nrow(leading))), , drop = FALSE]
leading$gene <- factor(leading$gene, levels = rev(leading$gene))
leading$direction <- ifelse(leading$wald_statistic > 0, "Male higher", "Female higher")
leading$p_label <- paste0("P=", format_p(leading$nominal_p))

p_leading <- ggplot2::ggplot(
  leading,
  ggplot2::aes(x = wald_statistic, y = gene, colour = direction)
) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.5) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = wald_statistic, yend = gene),
    colour = "grey72",
    linewidth = 0.8
  ) +
  ggplot2::geom_point(size = 4) +
  ggplot2::geom_text(
    ggplot2::aes(label = p_label),
    colour = "black",
    size = 3.5,
    hjust = ifelse(leading$wald_statistic >= 0, -0.15, 1.15),
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = okabe_ito) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.25, 0.25))) +
  ggplot2::labs(
    title = paste(strwrap(top_result$term_name, width = 52), collapse = "\n"),
    subtitle = paste0(
      top_result$go_id, "; leading-edge genes; NES=",
      formatC(top_result$NES, digits = 3, format = "f"),
      "; nominal GSEA P=", format_p(top_result$nominal_p)
    ),
    x = "DESeq2 Wald statistic (male minus female)",
    y = NULL,
    colour = NULL,
    caption = "Up to 20 leading-edge genes are displayed; the complete leading edge is retained in the result table."
  ) +
  theme_mobile(12) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.major.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(18, 45, 18, 18)
  )
ggplot2::ggsave(
  file.path(root, "figures", "paper_gsea_03_inflammasome_leading_edge.png"),
  p_leading,
  width = 8,
  height = 8,
  dpi = 220,
  units = "in",
  bg = "white"
)

cat("\nINFLAMMASOME GSEA SUMMARY\n")
print(summary, row.names = FALSE)
cat("\nINFLAMMASOME GSEA RESULTS\n")
print(gsea, row.names = FALSE)

log_session("04", root)
