root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "utils", "log_session.R"))
assert_project_local_r(root)

dir.create(file.path(root, "data", "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "data", "metadata"), recursive = TRUE, showWarnings = FALSE)

sources <- c(
  GSE188872_placenta_raw_counts.csv.gz = paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188872/suppl/",
    "GSE188872_placenta_raw_counts.csv.gz"
  ),
  GSE188872_placenta_tpm.csv.gz = paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188872/suppl/",
    "GSE188872_placenta_tpm.csv.gz"
  ),
  GSE188872_series_matrix.txt.gz = paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188872/matrix/",
    "GSE188872_series_matrix.txt.gz"
  )
)

for (name in names(sources)) {
  destination <- file.path(root, "data", "raw", name)
  if (!file.exists(destination)) {
    download.file(sources[[name]], destination, mode = "wb", quiet = FALSE)
  }
}

metadata_path <- file.path(root, "data", "metadata", "human_seq_metadata.csv")
if (!file.exists(metadata_path)) {
  download.file(
    paste0(
      "https://raw.githubusercontent.com/bendevlin18/human-fetal-RNASeq/",
      "main/shiny_app/human_seq_metadata.csv"
    ),
    metadata_path,
    mode = "wb",
    quiet = FALSE
  )
}

checksum_file <- file.path(root, "data", "checksums.md5")
if (!file.exists(checksum_file)) {
  stop("Missing committed checksum manifest: data/checksums.md5", call. = FALSE)
}
manifest <- read.table(
  checksum_file,
  col.names = c("md5", "path"),
  stringsAsFactors = FALSE,
  comment.char = "#"
)
for (i in seq_len(nrow(manifest))) {
  path <- file.path(root, manifest$path[[i]])
  if (!file.exists(path)) stop("Missing input: ", manifest$path[[i]], call. = FALSE)
  observed <- unname(tools::md5sum(path))
  if (!identical(tolower(observed), tolower(manifest$md5[[i]]))) {
    stop(
      "Checksum mismatch for ", manifest$path[[i]],
      "; expected ", manifest$md5[[i]], " but observed ", observed,
      call. = FALSE
    )
  }
}

cat("All input checksums verified.\n")

