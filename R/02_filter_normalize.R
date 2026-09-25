options(stringsAsFactors = FALSE)
set.seed(42)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "utils", "theme_mobile.R"))
source(file.path(root, "R", "utils", "log_session.R"))
assert_project_local_r(root)

for (dir in c("results", "figures", "logs")) {
  dir.create(file.path(root, dir), recursive = TRUE, showWarnings = FALSE)
}

source(file.path(root, "R", "00_fetch_data.R"))

if (!requireNamespace("edgeR", quietly = TRUE) || !requireNamespace("limma", quietly = TRUE)) {
  stop("Stage 02 requires the project-local edgeR and limma packages.", call. = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

save_plot_grid <- function(plots, path, nrow, ncol, width_px = 1100, height_px = 1000, dpi = 150) {
  grDevices::png(path, width = width_px, height = height_px, res = dpi, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  layout <- grid::grid.layout(nrow, ncol)
  grid::pushViewport(grid::viewport(layout = layout))
  for (i in seq_along(plots)) {
    row <- ((i - 1L) %/% ncol) + 1L
    col <- ((i - 1L) %% ncol) + 1L
    print(plots[[i]], vp = grid::viewport(layout.pos.row = row, layout.pos.col = col))
  }
  invisible(path)
}

robust_z <- function(x) {
  scale_value <- stats::mad(x, center = stats::median(x), constant = 1.4826)
  if (!is.finite(scale_value) || scale_value == 0) return(rep(0, length(x)))
  (x - stats::median(x)) / scale_value
}

counts_path <- file.path(root, "data", "raw", "GSE188872_placenta_raw_counts.csv.gz")
sex_path <- file.path(root, "results", "sex_concordance.tsv")
if (!file.exists(sex_path)) stop("Stage 01 output is missing: results/sex_concordance.tsv", call. = FALSE)

raw <- read.csv(counts_path, check.names = FALSE)
x_cols <- grep("\\.x$", names(raw), value = TRUE)
count_columns <- sub("\\.x$", "", x_cols)
counts <- as.matrix(raw[, x_cols, drop = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- raw$gene

stage1 <- read.delim(sex_path, check.names = FALSE)
match_index <- match(count_columns, stage1$count_column)
if (anyNA(match_index)) stop("Stage 02 count columns did not all match Stage 01 frozen metadata.", call. = FALSE)
metadata <- stage1[match_index, , drop = FALSE]
if (!identical(metadata$count_column, count_columns)) stop("Sample ordering mismatch.", call. = FALSE)
if (!all(metadata$primary_include)) stop("An ambiguous Stage 01 sample entered PRIMARY.", call. = FALSE)

metadata$sex <- factor(metadata$inferred_sex, levels = c("Female", "Male"))
metadata$trig_z <- as.numeric(scale(metadata$trig_as_reported))
metadata$age_z <- as.numeric(scale(metadata$age_dpc))
design <- stats::model.matrix(~ sex + trig_z + age_z, data = metadata)
rownames(design) <- metadata$sample_id
if (qr(design)$rank != ncol(design)) stop("Stage 03 design is rank deficient.", call. = FALSE)

dge_unfiltered <- edgeR::DGEList(counts = counts, genes = data.frame(gene = rownames(counts)))
keep <- edgeR::filterByExpr(dge_unfiltered, group = metadata$sex)
dge <- dge_unfiltered[keep, , keep.lib.sizes = FALSE]
dge <- edgeR::normLibSizes(dge, method = "TMM")

filter_status <- data.frame(
  gene = rownames(counts),
  present_in_deposited_matrix = TRUE,
  retained_filterByExpr = as.logical(keep),
  removal_reason = ifelse(keep, "retained", "low_expression_filterByExpr")
)
write_tsv(filter_status, file.path(root, "results", "gene_filter_status.tsv"))

filter_summary <- data.frame(
  input_genes = nrow(counts),
  retained_genes = sum(keep),
  removed_genes = sum(!keep),
  retained_percent = 100 * mean(keep),
  filter = "edgeR::filterByExpr",
  group = "Stage 01 inferred sex (Female/Male)"
)
write_tsv(filter_summary, file.path(root, "results", "filtering_summary.tsv"))

raw_library_size <- colSums(counts)
filtered_library_size <- dge$samples$lib.size
norm_factor <- dge$samples$norm.factors
effective_library_size <- filtered_library_size * norm_factor
library_z <- robust_z(log2(raw_library_size))
norm_z <- robust_z(log2(norm_factor))
library_outlier <- abs(library_z) > 3
norm_factor_outlier <- abs(norm_z) > 3

normalization <- data.frame(
  sample_id = metadata$sample_id,
  inferred_sex = metadata$inferred_sex,
  raw_library_size = as.numeric(raw_library_size),
  filtered_library_size = as.numeric(filtered_library_size),
  tmm_norm_factor = as.numeric(norm_factor),
  effective_library_size = as.numeric(effective_library_size),
  robust_z_log2_raw_library = as.numeric(library_z),
  robust_z_log2_norm_factor = as.numeric(norm_z),
  library_size_outlier_3mad = library_outlier,
  norm_factor_outlier_3mad = norm_factor_outlier,
  author_note = metadata$author_notes
)
write_tsv(normalization, file.path(root, "results", "library_normalization.tsv"))

voom_object <- limma::voom(dge, design = design, plot = FALSE, save.plot = TRUE)

log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
pca <- stats::prcomp(t(log_cpm), center = TRUE, scale. = FALSE)
variance_percent <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
pca_scores <- data.frame(
  sample_id = metadata$sample_id,
  inferred_sex = metadata$inferred_sex,
  age_dpc = metadata$age_dpc,
  trig_as_reported = metadata$trig_as_reported,
  author_exclude = metadata$author_exclude,
  author_notes = metadata$author_notes,
  pca$x[, seq_len(min(4L, ncol(pca$x))), drop = FALSE],
  check.names = FALSE
)
names(pca_scores)[grepl("^PC", names(pca_scores))] <- paste0("PC", seq_len(sum(grepl("^PC", names(pca_scores)))))
write_tsv(pca_scores, file.path(root, "results", "pca_scores.tsv"))

pca_variance <- data.frame(
  component = paste0("PC", seq_len(min(4L, length(variance_percent)))),
  variance_percent = variance_percent[seq_len(min(4L, length(variance_percent)))]
)
write_tsv(pca_variance, file.path(root, "results", "pca_variance.tsv"))

association_rows <- list()
for (pc in paste0("PC", 1:4)) {
  score <- pca_scores[[pc]]
  sex_model <- stats::lm(score ~ metadata$sex)
  sex_anova <- stats::anova(sex_model)
  association_rows[[length(association_rows) + 1L]] <- data.frame(
    component = pc,
    covariate = "inferred_sex",
    association_metric = "male_minus_female_mean",
    estimate = mean(score[metadata$sex == "Male"]) - mean(score[metadata$sex == "Female"]),
    r_squared = summary(sex_model)$r.squared,
    p_value = sex_anova$`Pr(>F)`[[1]],
    status = "tested"
  )
  numeric_covariates <- list(age_dpc = metadata$age_dpc, trig_as_reported = metadata$trig_as_reported)
  for (covariate_name in names(numeric_covariates)) {
    values <- numeric_covariates[[covariate_name]]
    test <- stats::cor.test(score, values, method = "pearson")
    association_rows[[length(association_rows) + 1L]] <- data.frame(
      component = pc,
      covariate = covariate_name,
      association_metric = "pearson_r",
      estimate = unname(test$estimate),
      r_squared = unname(test$estimate)^2,
      p_value = test$p.value,
      status = "tested"
    )
  }
  association_rows[[length(association_rows) + 1L]] <- data.frame(
    component = pc,
    covariate = "author_exclude",
    association_metric = "not_testable",
    estimate = NA_real_,
    r_squared = NA_real_,
    p_value = NA_real_,
    status = "single level: all included samples have exclude=n"
  )
}
pca_associations <- do.call(rbind, association_rows)
pca_associations$p_adjust_bh <- NA_real_
tested_rows <- pca_associations$status == "tested"
pca_associations$p_adjust_bh[tested_rows] <- stats::p.adjust(
  pca_associations$p_value[tested_rows],
  method = "BH"
)
write_tsv(pca_associations, file.path(root, "results", "pca_associations.tsv"))

fit <- limma::lmFit(voom_object, design)
fit <- limma::eBayes(fit, robust = TRUE)
sex_coefficient <- which(colnames(design) == "sexMale")
if (length(sex_coefficient) != 1L) stop("Could not identify the male-vs-female design coefficient.", call. = FALSE)
sex_t <- fit$t[, sex_coefficient]
mean_sex_t <- mean(sex_t, na.rm = TRUE)
sanity <- data.frame(
  coefficient = "sexMale (male vs female)",
  n_filtered_genes = length(sex_t),
  mean_moderated_t = mean_sex_t,
  median_moderated_t = stats::median(sex_t, na.rm = TRUE),
  sd_moderated_t = stats::sd(sex_t, na.rm = TRUE),
  fraction_positive = mean(sex_t > 0, na.rm = TRUE),
  threshold_abs_mean_t = 0.1,
  pass = abs(mean_sex_t) <= 0.1,
  normalization = "TMM + voom",
  design = "~ inferred_sex + scale(trig_as_reported_log10) + scale(age_dpc)"
)
write_tsv(sanity, file.path(root, "results", "normalization_sanity.tsv"))

analysis_parameters <- data.frame(
  parameter = c(
    "random_seed", "primary_samples", "female_samples", "male_samples",
    "trig_scale", "filter", "normalization", "voom_design", "outlier_rule"
  ),
  value = c(
    "42", nrow(metadata), sum(metadata$sex == "Female"), sum(metadata$sex == "Male"),
    "already log10 as deposited; no second log applied", "edgeR::filterByExpr(group=inferred sex)",
    "edgeR TMM", "~ inferred_sex + scale(trig_as_reported) + scale(age_dpc)",
    "absolute robust z > 3 using median/MAD on log2 metric"
  )
)
write_tsv(analysis_parameters, file.path(root, "results", "stage02_parameters.tsv"))

voom_points <- data.frame(
  mean_log_count = voom_object$voom.xy$x,
  sqrt_residual_sd = voom_object$voom.xy$y
)
voom_line <- data.frame(
  mean_log_count = voom_object$voom.line$x,
  sqrt_residual_sd = voom_object$voom.line$y
)
p02 <- ggplot2::ggplot(voom_points, ggplot2::aes(mean_log_count, sqrt_residual_sd)) +
  ggplot2::geom_point(colour = okabe_ito[["sky_blue"]], alpha = 0.25, size = 0.8) +
  ggplot2::geom_line(
    data = voom_line,
    ggplot2::aes(mean_log_count, sqrt_residual_sd),
    colour = okabe_ito[["vermillion"]],
    linewidth = 1.2
  ) +
  ggplot2::labs(
    title = "Voom mean–variance trend after TMM normalization",
    subtitle = paste0(sum(keep), " expressed genes across ", nrow(metadata), " placentas"),
    x = "Mean log2 count",
    y = "Square-root residual standard deviation",
    caption = "Red line: fitted mean–variance trend used to derive precision weights."
  ) +
  theme_mobile(12)
ggplot2::ggsave(
  file.path(root, "figures", "02_voom_mean_variance.png"),
  p02,
  width = 7.33,
  height = 5.33,
  dpi = 150,
  units = "in",
  bg = "white"
)

pca_plot_data <- pca_scores
pca_plot_data$flag_label <- ifelse(pca_plot_data$sample_id == "H-19760", "H-19760", "")
base_pca <- function(title, subtitle = NULL) {
  ggplot2::ggplot(pca_plot_data, ggplot2::aes(PC1, PC2)) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = sprintf("PC1 (%.1f%%)", variance_percent[[1]]),
      y = sprintf("PC2 (%.1f%%)", variance_percent[[2]])
    ) +
    theme_mobile(12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12),
      plot.subtitle = ggplot2::element_text(size = 10),
      legend.position = "bottom"
    )
}
pca_sex <- base_pca("Inferred sex") +
  ggplot2::geom_point(ggplot2::aes(colour = inferred_sex, shape = inferred_sex), size = 3) +
  ggplot2::geom_text(ggplot2::aes(label = flag_label), nudge_y = 0.8, colour = "black", size = 3.4) +
  ggplot2::scale_colour_manual(
    values = c(Female = okabe_ito[["reddish_purple"]], Male = okabe_ito[["blue"]]),
    name = "Inferred sex"
  ) +
  ggplot2::scale_shape_manual(values = c(Female = 16, Male = 17), name = "Inferred sex")
pca_age <- base_pca("Gestational age", "days post conception") +
  ggplot2::geom_point(ggplot2::aes(colour = age_dpc), size = 3) +
  ggplot2::geom_text(ggplot2::aes(label = flag_label), nudge_y = 0.8, colour = "black", size = 3.4) +
  ggplot2::scale_colour_gradientn(
    colours = c(okabe_ito[["sky_blue"]], okabe_ito[["yellow"]], okabe_ito[["vermillion"]]),
    name = "Age (d.p.c.)"
  )
pca_trig <- base_pca("Decidual triglyceride proxy", "reported log10 values") +
  ggplot2::geom_point(ggplot2::aes(colour = trig_as_reported), size = 3) +
  ggplot2::geom_text(ggplot2::aes(label = flag_label), nudge_y = 0.8, colour = "black", size = 3.4) +
  ggplot2::scale_colour_gradientn(
    colours = c(okabe_ito[["sky_blue"]], okabe_ito[["yellow"]], okabe_ito[["vermillion"]]),
    name = "Reported log10 trig",
    breaks = c(0.5, 1.0, 1.5),
    labels = c("0.5", "1.0", "1.5")
  )
pca_exclude <- base_pca("Author exclusion flag", "all matrix samples: exclude = n") +
  ggplot2::geom_point(colour = okabe_ito[["bluish_green"]], shape = 16, size = 3) +
  ggplot2::geom_text(ggplot2::aes(label = flag_label), nudge_y = 0.8, colour = "black", size = 3.4)
save_plot_grid(
  list(pca_sex, pca_age, pca_trig, pca_exclude),
  file.path(root, "figures", "03_pca_panels.png"),
  nrow = 2,
  ncol = 2,
  width_px = 1100,
  height_px = 1050
)

normalization$sample_factor <- factor(normalization$sample_id, levels = rev(metadata$sample_id))
normalization$any_outlier <- normalization$library_size_outlier_3mad | normalization$norm_factor_outlier_3mad
point_colours <- c(`FALSE` = okabe_ito[["blue"]], `TRUE` = okabe_ito[["vermillion"]])
lib_limits <- stats::median(log2(raw_library_size)) + c(-3, 3) * stats::mad(log2(raw_library_size))
norm_limits <- stats::median(log2(norm_factor)) + c(-3, 3) * stats::mad(log2(norm_factor))

p_lib <- ggplot2::ggplot(normalization, ggplot2::aes(log2(raw_library_size), sample_factor)) +
  ggplot2::geom_vline(xintercept = lib_limits, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(ggplot2::aes(colour = library_size_outlier_3mad, shape = inferred_sex), size = 3) +
  ggplot2::scale_colour_manual(values = point_colours, name = ">3 MAD") +
  ggplot2::scale_shape_manual(values = c(Female = 16, Male = 17), name = "Inferred sex") +
  ggplot2::labs(
    title = "Raw library size",
    subtitle = "orange: >3 MAD | circle: F | triangle: M",
    x = "log2 total counts",
    y = "Sample"
  ) +
  theme_mobile(12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 12),
    plot.subtitle = ggplot2::element_text(size = 9),
    legend.position = "none"
  )
p_norm <- ggplot2::ggplot(normalization, ggplot2::aes(log2(tmm_norm_factor), sample_factor)) +
  ggplot2::geom_vline(xintercept = norm_limits, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_vline(xintercept = 0, colour = "grey75") +
  ggplot2::geom_point(ggplot2::aes(colour = norm_factor_outlier_3mad, shape = inferred_sex), size = 3) +
  ggplot2::scale_colour_manual(values = point_colours, name = ">3 MAD") +
  ggplot2::scale_shape_manual(values = c(Female = 16, Male = 17), name = "Inferred sex") +
  ggplot2::labs(
    title = "TMM normalization factor",
    subtitle = "orange: >3 MAD | circle: F | triangle: M",
    x = "log2 TMM factor",
    y = NULL
  ) +
  theme_mobile(12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 12),
    plot.subtitle = ggplot2::element_text(size = 9),
    legend.position = "none"
  )
save_plot_grid(
  list(p_lib, p_norm),
  file.path(root, "figures", "04_library_size_norm_factors.png"),
  nrow = 1,
  ncol = 2,
  width_px = 1100,
  height_px = 1100
)

cat("\nFILTERING SUMMARY\n")
print(filter_summary, row.names = FALSE)
cat("\nLIBRARY / NORMALIZATION OUTLIERS (>3 robust MAD)\n")
print(
  normalization[
    normalization$library_size_outlier_3mad | normalization$norm_factor_outlier_3mad,
    c("sample_id", "inferred_sex", "raw_library_size", "tmm_norm_factor", "library_size_outlier_3mad", "norm_factor_outlier_3mad", "author_note")
  ],
  row.names = FALSE
)
cat("\nPCA VARIANCE\n")
print(pca_variance, row.names = FALSE)
cat("\nPCA ASSOCIATIONS (PC1-PC4)\n")
print(pca_associations, row.names = FALSE)
cat("\nMANDATORY GLOBAL-SHIFT SANITY CHECK\n")
print(sanity, row.names = FALSE)

log_session("02", root)

if (!sanity$pass) {
  stop(
    sprintf(
      "Stage 02 global-shift check failed: abs(mean moderated t) = %.4f > 0.1. Diagnose before Stage 03.",
      abs(mean_sex_t)
    ),
    call. = FALSE
  )
}
