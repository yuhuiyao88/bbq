# Use the private study package when it has been installed.
.study_library <- Sys.getenv("AREXT12_R_LIBRARY", file.path(Sys.getenv("USERPROFILE", path.expand("~")), "arext12_library"))
if (dir.exists(.study_library)) .libPaths(c(.study_library, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

# Package and study integration checks. Model fits belong to smoke_driver_arext12.R.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]), winslash = "/")) else getwd()
source(file.path(here, "verify_config.R"))
tmp <- tempfile("arext12_verify_"); dir.create(tmp)
Sys.setenv(TV11_SCRIPT = file.path(here, "simulation_arext12.R"),
           TV11_WORK_DIR = file.path(tmp, "work"), TV11_BOX_DIR = file.path(tmp, "box"),
           AREXT12_DEFINITIONS_ONLY = "1")
source(file.path(here, "simulation_arext12.R"))
stopifnot(identical(TAUS, c(.25, .5, .75)), packageVersion("bbq") >= "0.1.2")
ex <- parse(file.path(here, "simulation_arext12.R"), keep.source = FALSE)
obj <- vapply(Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]), ex),
              function(x) as.character(x[[2]]), "")
src <- deparse(body(run_batch_simulation), width.cutoff = 500L)
i <- grep("export_names <- c\\(", src); j <- i + which(grepl("\\)\\s*$", src[i:(i+60)]))[1] - 1
exported <- eval(parse(text = paste(src[i:j], collapse = "\n"))[[1]][[3]])
seen <- character(); needed <- character(); queue <- "run_single_simulation"
while (length(queue)) {
  nm <- queue[1]; queue <- queue[-1]
  if (nm %in% seen) next
  seen <- c(seen, nm)
  g <- codetools::findGlobals(get(nm), merge = FALSE)
  loc <- intersect(unique(c(g$functions, g$variables)), obj)
  needed <- union(needed, loc)
  queue <- c(queue, loc[vapply(loc, function(n) is.function(get(n)), TRUE)])
}
stopifnot(length(setdiff(needed, exported)) == 0L)
cat("PASS: all worker-side study globals exported\n")
sup <- new.env(parent = globalenv())
source(file.path(here, "run_arext12_supervised.R"), local = sup)
stopifnot(identical(sup$.tv11_basis_sig(), BASIS_SIG))
cat("PASS: driver/supervisor fingerprints agree\n")
changed <- file.path(tmp, "changed.R")
txt <- readLines(file.path(here, "simulation_arext12.R"))
txt <- sub("c(0.25, 0.5, 0.75)", "c(0.2, 0.5, 0.8)", txt, fixed = TRUE)
writeLines(txt, changed)
stopifnot(!identical(.tv11_basis_sig_from(changed, TAUS), BASIS_SIG))
cat("PASS: changed quantile grid invalidates provenance\n")
stopifnot(sup$arext12_preflight(cores = 1L, target = 100L, sweep_stray = FALSE))
cat("AR_EXT12 INTEGRATION CHECKS PASSED\n")
