source("renv/activate.R")
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
local_library <- file.path(project_root, ".r-local", "library")
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_library, .Library.site, .Library))

if (file.exists(file.path(project_root, "renv", "activate.R"))) {
  source(file.path(project_root, "renv", "activate.R"))
}

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  stringsAsFactors = FALSE
)

