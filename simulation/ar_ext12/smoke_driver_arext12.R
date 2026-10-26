# Use the private study package when it has been installed.
.study_library <- Sys.getenv("AREXT12_R_LIBRARY", file.path(Sys.getenv("USERPROFILE", path.expand("~")), "arext12_library"))
if (dir.exists(.study_library)) .libPaths(c(.study_library, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

# smoke_driver_arext12.R -- run the REAL study driver end to end: 1 replicate, 1 config, 1 worker.
#
# Author's rule (2026-09-05): before the simulation code is handed over, the actual driver
# (simulation_arext12.R, not a definitions-only extract) must run with 1 replicate on the
# quickest configuration with 1 thread. TV11_FORCE_CLUSTER makes the driver create a
# 1-worker PSOCK cluster, so clusterExport() and the worker-side code path are exercised;
# a plain 1-core run takes the serial path and would miss a missing export, which is what
# killed every task of the 2026-09-04 run ("object 'CV_IQ_EM_MAX_ITER' not found").
# ar_ext12 (2026-09-25): NO EM -- lambda_iq is CV-tuned on IQ_LAMBDA_GRID; the checks below reflect that.
#
# Usage:  Rscript "<folder>/smoke_driver_arext12.R"      or      source("<folder>/smoke_driver_arext12.R")
#         SMOKE_CONFIG=<id> selects one of the 33 manuscript cells (default 19 = T09 null).
#         SMOKE_OUT=<dir>   output root (default <TEMP>/arext12_driver_smoke); kept after the run.
# Writes only under that root (a work dir and a stand-in "box" dir); nothing touches the study folder.

# ---- locate the study folder (same logic as smoke_arext12.R) ----
.here <- local({
  d <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    of <- tryCatch(get("ofile", envir = sys.frame(i), inherits = FALSE), error = function(e) NULL)
    if (is.character(of) && length(of) == 1L && nzchar(of)) { d <- dirname(normalizePath(of, winslash = "/", mustWork = FALSE)); break }
  }
  if (is.null(d)) {
    a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(a) == 1L) d <- dirname(normalizePath(sub("^--file=", "", a), winslash = "/", mustWork = FALSE))
  }
  if (is.null(d) && nzchar(Sys.getenv("TV11_BOX_DIR"))) d <- Sys.getenv("TV11_BOX_DIR")
  if (is.null(d)) d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  d
})
driver <- file.path(.here, "simulation_arext12.R")
if (!file.exists(driver)) stop("simulation_arext12.R not found next to this script: ", .here)

# ---- persistent output root (tempdir() would vanish with the session), one task, one worker ----
root <- Sys.getenv("SMOKE_OUT", file.path(Sys.getenv("TMPDIR", unset = Sys.getenv("TEMP", unset = tempdir())), "arext12_driver_smoke"))
tmp  <- file.path(root, format(Sys.time(), "%Y%m%d_%H%M%S"))
work <- file.path(tmp, "work"); box <- file.path(tmp, "box")
dir.create(work, recursive = TRUE, showWarnings = FALSE); dir.create(box, recursive = TRUE, showWarnings = FALSE)
cfg_id <- Sys.getenv("SMOKE_CONFIG", "19")
Sys.unsetenv("AREXT12_DEFINITIONS_ONLY")
Sys.setenv(TV11_SCRIPT = driver, TV11_BOX_DIR = box, TV11_WORK_DIR = work,
           TV11_NSIM = "1", TV11_TARGET = "1", TV11_CORES = "1", TV11_FORCE_CLUSTER = "1",
           TV11_ONLY_CONFIG = cfg_id, TV11_INIT_BATCH = "1", TV11_INIT_DELAY = "0",
           TV11_VERBOSE = "1")
cat("=== driver smoke: 1 replicate | config", cfg_id, "| 1 worker (forced PSOCK cluster) ===\n")
cat("  driver :", driver, "\n  work   :", work, "\n  box    :", box, "\n")
cat("  bbq    :", as.character(utils::packageVersion("bbq")), "\n\n")
t0 <- Sys.time()
source(driver)     # into the global env: run_batch_simulation() exports from .GlobalEnv; the bottom block runs the study
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# ---- verdict ----
fail <- 0L
ok <- function(cond, msg, extra = "") {
  cat(sprintf("  [%s] %-60s %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg, extra))
  if (!isTRUE(cond)) fail <<- fail + 1L
}
cat(sprintf("\n=== verdict (%.1f min including Stan compilation) ===\n", mins))
task <- file.path(work, sprintf("arext12_task_config_%s_sim_1.rds", cfg_id))
ok(file.exists(task), "task checkpoint written by the worker", basename(task))
if (file.exists(task)) {
  o <- readRDS(task); r <- o$rows
  ok(is.data.frame(r), "checkpoint carries a rows data.frame", sprintf("%d rows x %d cols", nrow(r), ncol(r)))
  ok(isTRUE(nrow(r) == 376L), "376 rows = 8 cpm + 8 arima_cpm + 2 x (6 priors x 30 methods)", nrow(r))
  expected <- make_configs()[[as.integer(cfg_id)]]
  expected_blocks <- nrow(block_bounds(expected$n_obs, expected$warmup_period, expected$block_length))
  ok(all(r$n_blocks == expected_blocks & r$block_length == expected$block_length &
         r$shift_length == expected$shift_length & r$setting == expected$setting),
     "driver checkpoint records the selected manuscript configuration")
  ok(identical(as.character(o$basis_sig), BASIS_SIG),
     "checkpoint fingerprint matches the loaded three-quantile study", as.character(o$basis_sig))
  ok(all(c("n_retries", "fit_seed") %in% names(r)), "rows carry n_retries and fit_seed")
  bq <- r[r$method %in% c("BQQ", "arima_BQQ") | grepl("_w$", as.character(r$adjustment)), ]
  ok(nrow(bq) == 360L, "360 BQQ rows (2 streams x 6 priors x 30 methods)", nrow(bq))
  ok(!any(bq$n_errors == 1L, na.rm = TRUE), "no failed BQQ model (n_errors == 1 marks a timeout or error)",
     sprintf("%d of %d rows failed", sum(bq$n_errors == 1L, na.rm = TRUE), nrow(bq)))
  ok(!any(grepl("maxent", as.character(r$adjustment), fixed = TRUE)), "no maxent rows (basis retired)")
  ok(all(is.finite(bq$val_score)) && all(is.finite(bq$val_pinball)), "finite val_score and val_pinball on every BQQ row",
     sprintf("%d NA", sum(!is.finite(bq$val_score))))
  ok(all(is.finite(bq$best_lambda_iq2)) && all(sqrt(bq$best_lambda_iq2) %in% c(0, 1, 10, 100, 1000)), "lambda_iq2 on every BQQ row is a CV winner from IQ_LAMBDA_GRID^2 (ar_ext12: no EM)")
  ok(all(is.na(bq$iq_em_iters)) && all(is.na(bq$iq_em_converged)), "iq_em_iters / iq_em_converged are NA on every BQQ row (no EM ran)")
  ok(all(bq$iq_mode == "CV"), "iq_mode label is 'CV'")
  ok(length(unique(bq$prior_gamma)) == 6L, "all 6 priors present", paste(sort(unique(as.character(bq$prior_gamma))), collapse = " "))
  ok(length(unique(bq$adjustment)) == 30L, "30 distinct BQQ methods", length(unique(bq$adjustment)))
  ok(all(!is.na(bq$fit_seed)), "fit_seed recorded on every BQQ row",
     sprintf("retries per model: %s", paste(unique(bq$n_retries), collapse = ",")))
  cat("\n  per prior (stream BQQ): lambda_iq (CV winner), val_score, best spike_sd\n")
  one <- bq[bq$method == bq$method[1] & bq$adjustment == bq$adjustment[1], ]
  if (nrow(one)) print(data.frame(prior = one$prior_gamma, lambda_iq = sqrt(one$best_lambda_iq2),
                                  val_score = round(one$val_score, 4), spike_sd = one$best_spike_sd,
                                  lambda_nc = one$best_lambda_nc), row.names = FALSE)
  cat("\n  flags by method family (rows flagging any block):\n")
  fam <- sub("_(ui|t2)_.*$", "", as.character(bq$adjustment))
  print(table(family = fam, any_signal = bq$any_signal))
}
ok(file.exists(file.path(box, basename(task))), "checkpoint synced to the box dir")
wl <- Sys.glob(file.path(work, "arext12_progress.log.worker_*.log"))
ok(length(wl) >= 1L, "worker log written", if (length(wl)) basename(wl[1]) else "")
if (length(wl)) {
  lines <- readLines(wl[1], warn = FALSE)
  ok(!any(grepl("error", lines, ignore.case = TRUE) & !grepl("n_errors|0 error", lines)), "no error lines in the worker log")
  cat("\n  worker log, model timing lines:\n")
  ok(any(grepl("CV start|CV done|final fit", lines)), "verbose progress lines captured in the worker log")
  tl <- grep("CV start|CV done|final fit", lines, value = TRUE)
  cat(paste0("    ", if (length(tl) > 40) c(head(tl, 20), "    ...", tail(tl, 20)) else tl), sep = "\n")
}
cat(sprintf("\n%s  (%d failure%s)  outputs kept under %s\n", if (fail == 0L) "DRIVER SMOKE PASSED" else "DRIVER SMOKE FAILED",
            fail, if (fail == 1L) "" else "s", tmp))
if (fail > 0L && !interactive()) quit(status = 1L)
