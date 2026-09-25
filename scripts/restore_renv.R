root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
expected <- normalizePath(
  "C:/Users/antho/Desktop/git_repos/Candy_Inflammasome_Placenta_SD",
  winslash = "/",
  mustWork = TRUE
)
stopifnot(identical(tolower(root), tolower(expected)))
if (!requireNamespace("renv", quietly = TRUE)) {
  stop("renv is not available in the project-local library; run scripts/bootstrap_renv.R first.")
}
renv::restore(prompt = FALSE)

