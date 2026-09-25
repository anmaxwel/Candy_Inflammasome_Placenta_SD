root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
r_files <- list.files(
  file.path(root, "R"),
  pattern = "\\.[Rr]$",
  recursive = TRUE,
  full.names = TRUE
)
script_files <- list.files(
  file.path(root, "scripts"),
  pattern = "\\.[Rr]$",
  recursive = FALSE,
  full.names = TRUE
)
files <- c(r_files, script_files)
for (file in files) {
  parse(file = file)
  cat("parse_ok", sub(paste0("^", root, "/?"), "", normalizePath(file, winslash = "/")), "\n")
}

