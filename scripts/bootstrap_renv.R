root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
expected <- normalizePath(
  "C:/Users/antho/Desktop/git_repos/Candy_Inflammasome_Placenta_SD",
  winslash = "/",
  mustWork = TRUE
)
stopifnot(identical(tolower(root), tolower(expected)))

local_lib <- file.path(root, ".r-local", "library")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .Library.site, .Library))

if (!requireNamespace("renv", quietly = TRUE, lib.loc = local_lib)) {
  install.packages("renv", lib = local_lib, repos = "https://cloud.r-project.org")
}

renv::init(bare = TRUE, restart = FALSE, force = TRUE)
renv::install(c("ggplot2", "BiocManager"), prompt = FALSE)
BiocManager::install(version = "3.23", ask = FALSE, update = FALSE)
original_r_profile_user <- Sys.getenv("R_PROFILE_USER", unset = "")
Sys.setenv(R_PROFILE_USER = "")
BiocManager::install(
  c("edgeR", "limma", "DESeq2", "fgsea", "KEGGREST", "AnnotationDbi", "GO.db", "org.Hs.eg.db"),
  ask = FALSE,
  update = FALSE
)
Sys.setenv(R_PROFILE_USER = original_r_profile_user)
renv::snapshot(type = "explicit", prompt = FALSE)

cat("Project-local R bootstrap complete.\n")
cat("R:", R.version.string, "\n")
cat("Bioconductor:", as.character(BiocManager::version()), "\n")
cat("Libraries:\n")
cat(paste0("  ", .libPaths(), collapse = "\n"), "\n")
