status <- renv::status()
if (!isTRUE(status$synchronized)) {
  stop("renv lockfile and project library are not synchronized.", call. = FALSE)
}
cat("renv_synchronized\n")

