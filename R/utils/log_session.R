assert_project_local_r <- function(root = getwd()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  writable_paths <- c(
    Sys.getenv("R_USER"),
    Sys.getenv("R_LIBS_USER"),
    Sys.getenv("R_CACHE_DIR"),
    Sys.getenv("XDG_CACHE_HOME"),
    Sys.getenv("TEMP"),
    Sys.getenv("TMP"),
    Sys.getenv("TMPDIR"),
    Sys.getenv("RENV_PATHS_CACHE"),
    Sys.getenv("RENV_PATHS_LIBRARY")
  )
  writable_paths <- writable_paths[nzchar(writable_paths)]
  resolved <- normalizePath(writable_paths, winslash = "/", mustWork = FALSE)
  outside <- resolved[!startsWith(tolower(resolved), paste0(tolower(root), "/"))]
  if (length(outside)) {
    stop(
      "R writable path(s) escaped the project root: ",
      paste(outside, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

git_sha <- function(root = getwd()) {
  sha <- tryCatch(
    system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (!length(sha)) "UNBORN" else sha[[1]]
}

log_session <- function(stage, root = getwd()) {
  assert_project_local_r(root)
  log_dir <- file.path(root, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(log_dir, sprintf("stage_%s_sessionInfo.txt", stage))
  con <- file(out, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    sprintf("timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("git_sha_at_execution: %s", git_sha(root)),
    sprintf("project_root: %s", normalizePath(root, winslash = "/")),
    sprintf("R_USER: %s", Sys.getenv("R_USER")),
    sprintf("R_LIBS_USER: %s", Sys.getenv("R_LIBS_USER")),
    sprintf("TEMP: %s", Sys.getenv("TEMP")),
    "",
    capture.output(sessionInfo())
  ), con)
  invisible(out)
}

