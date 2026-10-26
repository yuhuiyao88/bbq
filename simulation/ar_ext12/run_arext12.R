# ar_ext12: full ar_ext7 grid, three quantiles and location/scale inference.
# See README.md. Default workstation size 80; AREXT12_CORES overrides it.
## ---- knobs ------------------------------------------------------------------
CORES        <- as.integer(Sys.getenv("AREXT12_CORES", "80"))    # workers
TARGET       <- 100L   # replications per config (the ladder fills evenly toward this)
SPEED_FACTOR <- 1      # how many times SLOWER this box is than the reference desktop.
                       # Scales the per-model TIMEOUT only -- it cannot change a fit
                       # that succeeds, it decides how long to wait before giving up.
                       # Set too LOW on a slow machine and healthy fits get killed and
                       # retried 6x each, which is what stranded cells in the old study.

## ---- paths ------------------------------------------------------------------
.here <- local({
  # When this file is sourced, --file= can name an outer launcher instead.
  for (i in rev(seq_len(sys.nframe()))) {
    of <- get0("ofile", envir = sys.frame(i), inherits = FALSE)
    if (is.character(of) && length(of) == 1L &&
        basename(of) == "run_arext12.R" && file.exists(of))
      return(dirname(normalizePath(of, winslash = "/")))
  }
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1L) {
    script <- sub("^--file=", "", a)
    if (basename(script) == "run_arext12.R" && file.exists(script))
      return(dirname(normalizePath(script, winslash = "/")))
  }
  if (file.exists(file.path(getwd(), "run_arext12.R")))
    return(normalizePath(getwd(), winslash = "/"))
  stop("Cannot locate run_arext12.R. Source the complete script by its full path; do not run selected lines.")
})

USER_HOME <- Sys.getenv("USERPROFILE", unset = path.expand("~"))
# Prefer the study's private installation and propagate its search path to workers.
study_library <- Sys.getenv("AREXT12_R_LIBRARY", file.path(USER_HOME, "arext12_library"))
if (dir.exists(study_library)) .libPaths(c(study_library, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
if ("bbq" %in% loadedNamespaces() && getNamespaceVersion("bbq") < "0.1.2")
  stop("An older bbq is already loaded. Restart R, then source run_arext12.R.")
# Use a study-specific override so an earlier ar_ext8 session cannot reuse its work dir.
WORK <- Sys.getenv("AREXT12_WORK_DIR", unset = file.path(USER_HOME, "arext12_run"))

# A smoke script may have been sourced earlier in this R session. Do not carry
# its single-config restriction or verbose test settings into the full study.
Sys.unsetenv(c("TV11_ONLY_CONFIG", "TV11_FORCE_CLUSTER", "TV11_VERBOSE", "AREXT12_DEFINITIONS_ONLY"))

Sys.setenv(TV11_SCRIPT         = file.path(.here, "simulation_arext12.R"),
           TV11_BOX_DIR        = .here,
           TV11_WORK_DIR       = WORK,
           TV11_CORES          = as.character(CORES),
           TV11_NSIM           = "120",          # POOL of sim_ids; must EXCEED TARGET
           TV11_INIT_BATCH     = "16",
           TV11_INIT_DELAY     = "180",
           TV11_TIMEOUT_SCALE  = as.character(SPEED_FACTOR))

cat("ar_ext12\n  script :", Sys.getenv("TV11_SCRIPT"),
    "\n  work   :", WORK, "\n  box    :", .here,
    "\n  cores  :", CORES, "  target:", TARGET, "  timeout scale:", SPEED_FACTOR, "\n\n")

source(file.path(.here, "run_arext12_supervised.R"))

## ---- 1. preflight -----------------------------------------------------------
# arext12_preflight() returns a bare logical via invisible(), not a list.
pf_ok <- isTRUE(arext12_preflight(target = TARGET, cores = CORES, sweep_stray = FALSE))
if (!pf_ok) {
  cat("\nPREFLIGHT FAILED -- not starting. Fix the items above.\n")
} else {

  ## ---- 2. where we stand ----------------------------------------------------
  cat("\n--- deficit before launch ---\n")
  arext12_status(target = TARGET)

  ## ---- 3. run ---------------------------------------------------------------
  # recycle_hours 24: the child is relaunched on that cadence so a slow leak or a
  # wedged worker cannot hold the run forever. In-flight tasks are lost at each
  # recycle, so retain the established interval until full-grid timing is measured.
  cat("\n--- launching the supervised ladder ---\n")
  run_arext12_supervised(target = TARGET, cores = CORES, recycle_hours = 24, sweep_stray = FALSE)
}

# -----------------------------------------------------------------------------
# To stop:  press Esc / Ctrl-C, then  run_arext12_kill()
# To check: arext12_status(target = 100)
# To aggregate when done (or any time -- it reads whatever exists):
#   Rscript build_app_csv_from_tasks.R
# -----------------------------------------------------------------------------
