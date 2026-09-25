options(stringsAsFactors = FALSE)
set.seed(42)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "utils", "theme_mobile.R"))
source(file.path(root, "R", "utils", "log_session.R"))
assert_project_local_r(root)

for (dir in c("results", "figures", "logs", file.path("data", "metadata"))) {
  dir.create(file.path(root, dir), recursive = TRUE, showWarnings = FALSE)
}

source(file.path(root, "R", "00_fetch_data.R"))

strip_quotes <- function(x) sub('^"(.*)"$', '\\1', x)
parse_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  get_field <- function(prefix) {
    line <- lines[startsWith(lines, prefix)]
    if (length(line) != 1L) stop("Expected one ", prefix, " line.", call. = FALSE)
    strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
  }
  titles <- get_field("!Sample_title")
  geo <- get_field("!Sample_geo_accession")
  source_name <- get_field("!Sample_source_name_ch1")
  characteristic_lines <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
  characteristic_values <- lapply(characteristic_lines, function(line) {
    strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
  })
  out <- data.frame(
    geo_accession = geo,
    title = titles,
    source_name = source_name,
    stringsAsFactors = FALSE
  )
  for (values in characteristic_values) {
    key <- sub(":.*$", "", values[[1]])
    key <- gsub("[^A-Za-z0-9]+", "_", tolower(key))
    out[[key]] <- sub("^[^:]+:\\s*", "", values)
  }
  out$sample_id <- sub("\\s+[BP]$", "", out$title)
  out$tissue_code <- sub("^.*\\s+([BP])$", "\\1", out$title)
  out
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

counts_path <- file.path(root, "data", "raw", "GSE188872_placenta_raw_counts.csv.gz")
tpm_path <- file.path(root, "data", "raw", "GSE188872_placenta_tpm.csv.gz")
metadata_path <- file.path(root, "data", "metadata", "human_seq_metadata.csv")
series_path <- file.path(root, "data", "raw", "GSE188872_series_matrix.txt.gz")

raw <- read.csv(counts_path, check.names = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE)
geo <- parse_series_matrix(series_path)
geo_placenta <- geo[geo$tissue_code == "P", , drop = FALSE]
write_tsv(geo_placenta, file.path(root, "data", "metadata", "geo_pdata.tsv"))

x_cols <- grep("\\.x$", names(raw), value = TRUE)
y_cols <- grep("\\.y$", names(raw), value = TRUE)
if (!length(x_cols) || length(x_cols) != length(y_cols)) {
  stop("Expected paired .x/.y sample columns in GEO raw-count file.", call. = FALSE)
}
x_base <- sub("\\.x$", "", x_cols)
y_base <- sub("\\.y$", "", y_cols)
if (!identical(x_base, y_base)) stop(".x/.y sample column order differs.", call. = FALSE)

counts <- as.matrix(raw[, x_cols, drop = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- raw$gene
scaled <- as.matrix(raw[, y_cols, drop = FALSE])
storage.mode(scaled) <- "numeric"
if (any(counts < 0) || any(abs(counts - round(counts)) > 1e-8)) {
  stop("The selected .x columns are not non-negative integer counts.", call. = FALSE)
}

column_audit <- do.call(rbind, lapply(seq_along(x_cols), function(i) {
  nonzero <- counts[, i] > 0
  ratios <- scaled[nonzero, i] / counts[nonzero, i]
  data.frame(
    sample_column = x_base[[i]],
    x_role = "integer_raw_counts_selected",
    y_role = "per_sample_scaled_values_rejected",
    x_library_sum_retained_genes = sum(counts[, i]),
    y_sum_retained_genes = sum(scaled[, i]),
    y_over_x_scale = stats::median(ratios),
    y_over_x_scale_sd = stats::sd(ratios),
    stringsAsFactors = FALSE
  )
}))
write_tsv(column_audit, file.path(root, "results", "input_column_audit.tsv"))

sample_id <- sub("^X([0-9]+)_.*$", "H-\\1", x_base)
sample_map <- data.frame(sample_id = sample_id, count_column = x_base)
sample_map <- merge(sample_map, metadata, by.x = "sample_id", by.y = "sampleID", all.x = TRUE)
if (anyNA(sample_map$sex)) stop("At least one count column did not map to metadata.", call. = FALSE)
sample_map <- sample_map[match(sample_id, sample_map$sample_id), , drop = FALSE]

metadata_unmatched <- setdiff(metadata$sampleID, sample_id)
counts_unmatched <- setdiff(sample_id, metadata$sampleID)
geo_unmatched_from_counts <- setdiff(geo_placenta$sample_id, sample_id)
counts_unmatched_from_geo <- setdiff(sample_id, geo_placenta$sample_id)
mapping_summary <- data.frame(
  comparison = c(
    "count_columns_to_GitHub_metadata",
    "GitHub_metadata_not_in_counts",
    "count_samples_not_in_GitHub_metadata",
    "GEO_placenta_samples_not_in_counts",
    "count_samples_not_in_GEO_placenta"
  ),
  n = c(
    sum(sample_id %in% metadata$sampleID),
    length(metadata_unmatched),
    length(counts_unmatched),
    length(geo_unmatched_from_counts),
    length(counts_unmatched_from_geo)
  ),
  sample_ids = c(
    paste(sample_id[sample_id %in% metadata$sampleID], collapse = ";"),
    paste(metadata_unmatched, collapse = ";"),
    paste(counts_unmatched, collapse = ";"),
    paste(geo_unmatched_from_counts, collapse = ";"),
    paste(counts_unmatched_from_geo, collapse = ";")
  )
)
write_tsv(mapping_summary, file.path(root, "results", "id_mapping_summary.tsv"))

geo_overlap <- merge(
  metadata,
  geo_placenta,
  by.x = "sampleID",
  by.y = "sample_id",
  suffixes = c("_github", "_geo")
)
geo_age <- as.numeric(sub("\\s*d\\.p\\.c\\.$", "", geo_overlap$age_geo))
geo_trig <- as.numeric(geo_overlap$log10_maternal_decidual_triglycerides)
disagreement <- data.frame(
  sample_id = geo_overlap$sampleID,
  sex_disagrees = geo_overlap$sex != geo_overlap$fetal_sex,
  age_disagrees = as.numeric(geo_overlap$age_github) != geo_age,
  trig_disagrees = abs(as.numeric(geo_overlap$trig) - geo_trig) > 1e-6,
  github_sex = geo_overlap$sex,
  geo_sex = geo_overlap$fetal_sex,
  github_age = as.numeric(geo_overlap$age_github),
  geo_age = geo_age,
  github_trig = as.numeric(geo_overlap$trig),
  geo_log10_trig = geo_trig
)
write_tsv(
  disagreement[disagreement$sex_disagrees | disagreement$age_disagrees | disagreement$trig_disagrees, ],
  file.path(root, "results", "metadata_disagreements.tsv")
)

markers <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY", "PCDH11Y")
tpm <- read.csv(tpm_path, check.names = FALSE)
tpm_genes <- as.character(tpm[[1]])
library_sizes <- colSums(counts)
cpm <- sweep(counts, 2, library_sizes / 1e6, "/")

marker_audit <- do.call(rbind, lapply(markers, function(gene) {
  present <- gene %in% rownames(counts)
  data.frame(
    gene = gene,
    present_raw_counts = present,
    present_tpm = gene %in% tpm_genes,
    raw_count_min = if (present) min(counts[gene, ]) else NA_real_,
    raw_count_max = if (present) max(counts[gene, ]) else NA_real_,
    cpm_min = if (present) min(cpm[gene, ]) else NA_real_,
    cpm_max = if (present) max(cpm[gene, ]) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write_tsv(marker_audit, file.path(root, "results", "sex_marker_presence.tsv"))

y_genes <- intersect(setdiff(markers, "XIST"), rownames(counts))
if (!"XIST" %in% rownames(counts) || length(y_genes) < 1L) {
  stop("Insufficient X/Y markers for transcriptomic sex assignment.", call. = FALSE)
}

xist_cpm <- cpm["XIST", ]
y_mean_cpm <- colMeans(cpm[y_genes, , drop = FALSE])
log2_xist <- log2(xist_cpm + 0.5)
log2_y <- log2(y_mean_cpm + 0.5)
largest_gap_midpoint <- function(x) {
  ordered <- sort(x)
  gaps <- diff(ordered)
  i <- which.max(gaps)
  c(cutoff = mean(ordered[c(i, i + 1L)]), gap_low = ordered[[i]], gap_high = ordered[[i + 1L]])
}
sex_score <- log2_y - log2_xist
score_cut <- largest_gap_midpoint(sex_score)
in_empirical_gap <- sex_score > score_cut[["gap_low"]] & sex_score < score_cut[["gap_high"]]
inferred <- ifelse(sex_score > score_cut[["cutoff"]], "Male", "Female")

sex_concordance <- data.frame(
  sample_id = sample_id,
  count_column = x_base,
  annotated_sex = sample_map$sex,
  inferred_sex = inferred,
  discordant = sample_map$sex != inferred,
  xist_cpm = as.numeric(xist_cpm),
  y_mean_cpm = as.numeric(y_mean_cpm),
  log2_xist_cpm_plus_0_5 = as.numeric(log2_xist),
  log2_y_mean_cpm_plus_0_5 = as.numeric(log2_y),
  sex_score_log2_y_minus_log2_xist = as.numeric(sex_score),
  in_empirical_ambiguous_gap = in_empirical_gap,
  age_dpc = as.numeric(sample_map$age),
  trig_as_reported = as.numeric(sample_map$trig),
  author_exclude = sample_map$exclude,
  author_notes = sample_map$Notes,
  primary_include = inferred != "Ambiguous",
  sensitivity_include = inferred != "Ambiguous" & sample_map$sex == inferred,
  stringsAsFactors = FALSE
)
for (gene in y_genes) sex_concordance[[paste0(tolower(gene), "_cpm")]] <- as.numeric(cpm[gene, ])
write_tsv(sex_concordance, file.path(root, "results", "sex_concordance.tsv"))

cutoffs <- data.frame(
  metric = paste0(
    "log2(mean(", paste(y_genes, collapse = ", "),
    ") CPM + 0.5) - log2(XIST CPM + 0.5)"
  ),
  cutoff = score_cut[["cutoff"]],
  gap_low = score_cut[["gap_low"]],
  gap_high = score_cut[["gap_high"]],
  n_inside_open_gap = sum(in_empirical_gap),
  rule_for_male = "score > cutoff"
)
write_tsv(cutoffs, file.path(root, "results", "sex_call_cutoffs.tsv"))

primary <- sex_concordance[sex_concordance$primary_include, , drop = FALSE]
sensitivity <- sex_concordance[sex_concordance$sensitivity_include, , drop = FALSE]
sample_sets <- rbind(
  data.frame(set = "PRIMARY", sample_id = primary$sample_id, inferred_sex = primary$inferred_sex),
  data.frame(set = "SENSITIVITY", sample_id = sensitivity$sample_id, inferred_sex = sensitivity$inferred_sex)
)
write_tsv(sample_sets, file.path(root, "results", "frozen_sample_sets.tsv"))

balance_one <- function(variable, label) {
  male <- primary[primary$inferred_sex == "Male", variable]
  female <- primary[primary$inferred_sex == "Female", variable]
  test <- wilcox.test(male, female, exact = FALSE)
  data.frame(
    variable = label,
    male_n = length(male),
    male_median = median(male),
    male_iqr = IQR(male),
    female_n = length(female),
    female_median = median(female),
    female_iqr = IQR(female),
    test = "Wilcoxon rank-sum",
    p_value = unname(test$p.value)
  )
}
covariate_balance <- rbind(
  balance_one("age_dpc", "gestational_age_dpc"),
  balance_one("trig_as_reported", "trig_as_reported_GEO_labels_log10")
)
write_tsv(covariate_balance, file.path(root, "results", "covariate_balance.tsv"))

plot_data <- sex_concordance
plot_data$label <- ifelse(plot_data$discordant, plot_data$sample_id, "")
label_data <- plot_data[plot_data$discordant, , drop = FALSE]
label_data <- label_data[order(label_data$log2_xist_cpm_plus_0_5, decreasing = TRUE), , drop = FALSE]
label_data$label_x <- 4.75
label_data$label_y <- seq(11.25, 10.05, length.out = nrow(label_data))
shape_values <- c(Female = 16, Male = 17, Ambiguous = 4)
colour_values <- c(Female = okabe_ito[["reddish_purple"]], Male = okabe_ito[["blue"]])

p <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = log2_y_mean_cpm_plus_0_5,
    y = log2_xist_cpm_plus_0_5,
    colour = annotated_sex,
    shape = inferred_sex
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = -score_cut[["cutoff"]],
    linetype = "dashed",
    colour = "grey40"
  ) +
  ggplot2::geom_point(size = 3.5, stroke = 1) +
  ggplot2::geom_segment(
    data = label_data,
    ggplot2::aes(
      x = log2_y_mean_cpm_plus_0_5,
      y = log2_xist_cpm_plus_0_5,
      xend = label_x - 0.08,
      yend = label_y
    ),
    inherit.aes = FALSE,
    colour = "grey45",
    linewidth = 0.4,
    show.legend = FALSE
  ) +
  ggplot2::geom_label(
    data = label_data,
    ggplot2::aes(x = label_x, y = label_y, label = sample_id),
    inherit.aes = FALSE,
    colour = "black",
    fill = "white",
    linewidth = 0.2,
    size = 3.4,
    hjust = 0,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = colour_values, name = "Annotated sex") +
  ggplot2::scale_shape_manual(values = shape_values, name = "Inferred sex") +
  ggplot2::labs(
    title = "Expression-based sex calls: 12 male, 22 female",
    subtitle = "Four reported SRY-male samples cluster with expression-inferred females",
    x = paste0("log2(mean ", paste(y_genes, collapse = "/"), " CPM + 0.5)"),
    y = "log2(XIST CPM + 0.5)",
    caption = paste(
      strwrap(
        paste0(
          "n = ", nrow(plot_data),
          "; colour shows reported SRY-PCR sex, shape shows expression-inferred sex. ",
          "The dashed diagonal is the empirical sex-score cutoff; labels mark discordant samples."
        ),
        width = 82
      ),
      collapse = "\n"
    )
  ) +
  theme_mobile(12)

ggplot2::ggsave(
  file.path(root, "figures", "01_xist_vs_y.png"),
  p,
  width = 7.33,
  height = 5.8,
  dpi = 150,
  units = "in",
  bg = "white"
)

cat("\nID MAPPING SUMMARY\n")
print(mapping_summary, row.names = FALSE)
cat("\nSEX MARKER PRESENCE\n")
print(marker_audit, row.names = FALSE)
cat("\nANNOTATED x INFERRED SEX\n")
print(table(Annotated = sex_concordance$annotated_sex, Inferred = sex_concordance$inferred_sex))
cat("\nDISCORDANT SAMPLES\n")
print(sex_concordance[sex_concordance$discordant, c("sample_id", "annotated_sex", "inferred_sex")], row.names = FALSE)
cat("\nFROZEN SAMPLE SETS\n")
print(table(Set = sample_sets$set, Inferred = sample_sets$inferred_sex))
cat("\nCOVARIATE BALANCE (PRIMARY)\n")
print(covariate_balance, row.names = FALSE)
cat("\nAUTHOR EXCLUDE FLAGS IN COUNT MATRIX\n")
print(table(sex_concordance$author_exclude, useNA = "ifany"))
cat("\nSEX CALL CUTOFFS\n")
print(cutoffs, row.names = FALSE)

log_session("01", root)
