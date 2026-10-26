# Install the bundled package into a private study library, leaving other libraries intact.
here <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    p <- get0("ofile", envir = sys.frame(i), inherits = FALSE)
    if (is.character(p) && length(p) == 1L && basename(p) == "install_arext12_package.R")
      return(dirname(normalizePath(p, winslash = "/")))
  }
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) return(dirname(normalizePath(sub("^--file=", "", arg[1]), winslash = "/")))
  getwd()
})
lib <- Sys.getenv("AREXT12_R_LIBRARY", file.path(Sys.getenv("USERPROFILE", path.expand("~")), "arext12_library"))
bundle <- file.path(here, "package", "bbq_0.1.2.tar.gz")
stopifnot(file.exists(bundle))
if ("bbq" %in% loadedNamespaces()) stop("Restart R before installing the bundled bbq package.")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
install.packages(bundle, repos = NULL, type = "source", lib = lib)
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
stopifnot(packageVersion("bbq", lib.loc = lib) == "0.1.2")
cat("Installed bbq 0.1.2 in", normalizePath(lib, winslash = "/"), "\n")
cat("Next: source run_arext12.R from this study folder.\n")
