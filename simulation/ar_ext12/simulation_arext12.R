# ar_ext12 -- Section 4 grid of the BBQ manuscript on the bbq package.
# taus = (0.25, 0.5, 0.75). QSS and L-moment inference use location/scale only.
# Requires bbq >= 0.1.2; checkpoints from other studies must never be pooled here.
# 33 configurations = 11 settings (the ten of ar_ext10 plus HK of ar_ext11) x
# {null, mean, variance}. Six priors, score CV, raw/prewhitened streams and
# benchmarks are inherited from ar_ext10. No EM. See README.md.
# Use run_arext12.R for production or smoke_driver_arext12.R for the real smoke.

# ---------------------------------------------------------------------------
# SECTION 0: console capture -- FIRST, before anything that can fail
# ---------------------------------------------------------------------------
# Child-side console capture. system(cmd, wait = FALSE) on Windows bypasses the shell,
# so the supervisor's '> tv11_console.log 2>&1' redirect silently did nothing: the log
# was never created and every crash/stall message from this process was lost. The child
# now captures its OWN stdout+stderr whenever the supervisor provides TV11_CONSOLE.
#
# This MUST come before library(bbq) below. It used to sit after it, which meant a
# package that failed to load -- a locked/updated rstan DLL, a half-written library
# tree -- killed the child with the console log not yet open, leaving a zero-byte file
# and no error text anywhere. That is indistinguishable, from the supervisor's side,
# from "Rscript never started", and it is the most likely reading of the 2026-08-04
# outage: 18 launches, no memory ever allocated, no manifest, nothing logged.
if (nzchar(Sys.getenv("TV11_CONSOLE", ""))) {
  .tv11_console_con <- file(Sys.getenv("TV11_CONSOLE"), open = "wt")
  sink(.tv11_console_con)
  sink(.tv11_console_con, type = "message")
}
cat("[", format(Sys.time(), "%H:%M:%S"), "] child started; R", R.version.string, "\n")
cat("[", format(Sys.time(), "%H:%M:%S"), "] loading packages (parallel, cpm, bbq, PearsonDS, forecast)...\n")
flush.console()

# ---------------------------------------------------------------------------
# SECTION 1: Packages + paths
# ---------------------------------------------------------------------------
suppressMessages({
  library(parallel)
  library(cpm)
  library(bbq)
  library(PearsonDS)
  library(forecast)          # [ar_ext E5] reference R2 = cpm on auto.arima residuals
})
cat("[", format(Sys.time(), "%H:%M:%S"), "] packages loaded; bbq",
    as.character(utils::packageVersion("bbq")), "\n")
flush.console()

# Requires bqq >= 0.4.5: unit-information intercept prior as the beta0_scale
# default (Kass & Wasserman 1995), prior-center MAP initialization, and tight
# optimizer tolerances. Refuse to run against an older package, where the
# defaults are a flat intercept prior and a random-init/loose-tolerance MAP.
stopifnot(utils::packageVersion("bbq") >= "0.1.2")   # one CV function, complete-data stop only, maxent retired
# 0.6.7 REQUIRED (2026-09-05, author rulings). Estimation is one chain per fit: the pilot
# initializes the iterations once, each EM step maximizes Eq. (17) from the previous
# solution, the update is the fixed point (C.9), and the iterations stop when the
# relative gain of the complete-data log posterior (Eq. 17 with r(m-1) log lambda_iq)
# is below iq_em_lp_tol = 1%. That is the only tolerance; iq_em_tol / iq_em_mc_tol /
# iq_em_switch_tol and the "hybrid" step are gone. The E-step is the closed-form
# folded-normal mean under the Laplace approximation (iq_em_estep = "closed").
# Cross-validation has ONE function, cv_copss_grid(y, taus, H, X, w, grid, base_args,
# loss, seed, verbose): every model argument goes through base_args to getModel(), so
# the tuning fits run the same chain as the final fit by construction. The maxent
# basis is retired; detectChangepoints_gamma(basis = "all") means quantile + qss + lmom.
# Priors aligned with Appendix B: spike-and-slab LASSO in Rockova-George rate form,
# adaptive LASSO local rates lambda_{q,j}^2 ~ Gamma(a_lambda, b_lambda) (Leng et al. 2014).
# 0.6.0 REQUIRED, not merely preferred: cv_copss_grid() gained loss = "score"
# (the exact unconstrained score-likelihood criterion) and now returns BOTH
# val_score and val_pinball. Under 0.5.x the loss argument does not exist --
# do.call would error -- and the both-loss columns this study records would be
# absent.
# 0.6.4 REQUIRED (2026-09-04): EM tolerances 0.1 (tuning fits) / 0.01 (final fit),
# lambda_iq2 = NULL pilot-scale start, map_init_values, cv_winner_init(), map_init_used.
# 0.6.3 (2026-09-03): cv_copss_grid() no longer forces adaptive_iq = FALSE,
# so the EM for lambda_iq2 runs INSIDE every CV fit as well as in the final fit
# (the manuscript's single estimation procedure), and getModel() gained the hybrid
# update iq_em_step = "hybrid" with iq_em_switch_tol (the former Appendix C Eq. (C.10), dropped in 0.6.6; (C.10) now labels the initial value).


# Durable (Box) destination vs. LOCAL work dir for hot I/O.
# ROOT CAUSE of the "stuck" runs: every per-task .rds + the shared progress log were
# written to the Box-synced folder. Box's sync client locks/scans files as 80 workers
# write them, so throughput collapses to ~single-thread AND, when the mount briefly goes
# unresponsive, every worker blocks inside a native file syscall that setTimeLimit CANNOT
# interrupt -> total freeze with no error. Fix: keep all hot per-task I/O on a LOCAL disk
# and copy only the aggregated outputs back to Box at the end. Both overridable by env var.
#
# PORTABILITY: no machine-specific paths. The durable (Box) destination defaults to
# THIS SCRIPT'S OWN FOLDER, so the study runs unchanged on any workstation as long as
# the folder is synced; the local work dir defaults to <home>/tv11_run. Override
# either with TV11_BOX_DIR / TV11_WORK_DIR.

# Treat an env var that is set-but-empty as unset.
.tv11_env <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

# Locate this script's own directory. source() is checked FIRST (its `ofile` is the
# more specific context): if this file is sourced from a script that Rscript is
# running, --file= names the OUTER script, not this one. Rscript is the fallback,
# then the working directory.
.tv11_script_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    of <- tryCatch(get("ofile", envir = sys.frame(i), inherits = FALSE),
                   error = function(e) NULL)
    if (is.character(of) && length(of) == 1L && nzchar(of))
      return(dirname(normalizePath(of, winslash = "/", mustWork = FALSE)))
  }
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1L)
    return(dirname(normalizePath(sub("^--file=", "", a), winslash = "/", mustWork = FALSE)))
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

.tv11_home <- function() Sys.getenv("USERPROFILE", unset = path.expand("~"))

box_dir  <- .tv11_env("TV11_BOX_DIR",  .tv11_script_dir())
work_dir <- .tv11_env("TV11_WORK_DIR", file.path(.tv11_home(), "arext12_run"))
for (d in c(box_dir, work_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
setwd(work_dir)          # keep all scratch (incl. any relative writes) on local disk
output_dir <- work_dir   # back-compat: downstream file.path(output_dir, "arext") -> local
cat("arext paths:\n  box_dir  (durable) :", box_dir,
    "\n  work_dir (hot I/O) :", work_dir, "\n")

# (console capture moved to SECTION 0, above the package block)

# EARLY liveness breadcrumb -- write the master PID before ANY slow work.
# The full manifest (master + every worker PID) is written after makeCluster, but on 80
# workers that is 10-30 min into startup, and the supervisor's "no PID manifest ->
# failed to start" rule fires at 15 min. On 2026-08-04 that killed a child which had
# not yet finished starting, 18 times in a row, and the study sat idle for 4.5 h with
# the watchdog log reporting nothing but "starting up". Writing the PID here makes
# "this process exists" observable from the first second, so a slow start can never
# again be mistaken for a failed one. The post-cluster write below overwrites this
# file with the full list; the master stays at index 1 either way, so kill semantics
# (taskkill /F /T on the master first) are unchanged.
try(writeLines(as.character(Sys.getpid()), file.path(work_dir, "arext12_run.pids")),
    silent = TRUE)
cat("startup: master pid", Sys.getpid(), "at", format(Sys.time()), "\n")

# Watchdog heartbeat: touch a per-worker file at every R-level safe point (task start +
# each model attempt). Checkpoint .rds files are the WRONG liveness signal for a
# watchdog -- one replication takes 60-140 min, so "no new checkpoint for 20 min" fires
# on a perfectly healthy run. Beats are minutes-granular, so fleet-wide beat silence
# really does mean frozen. The supervisor watches tv11_hb_* mtimes.
HB_DIR <- work_dir
tv11_beat <- function() {
  try(cat(format(Sys.time(), "%H:%M:%S"),
          file = file.path(HB_DIR, paste0("tv11_hb_", Sys.getpid()))), silent = TRUE)
  invisible()
}

# Atomic RDS write: write to a temp file then rename it into place, so every
# per-task / per-config checkpoint on disk is complete and the run is resumable.
save_rds_atomic <- function(obj, path) {
  tmp <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

# ---------------------------------------------------------------------------
# SECTION 2: shift-arm constants + DGP building blocks
# ---------------------------------------------------------------------------
# DESIGN POINT SET 2026-08-13 (user): every arm sits at TV = 0.383, which is
# Cohen's d = 1.0002 on the normal mean arm (TV = 2*Phi(d/2) - 1).  The previous
# design used TV = 0.38 for mean/variance/skewness and a CAPPED 0.11 / 0.162 for
# kurtosis, so those three arms barely move; the REAL change is that the KURTOSIS
# arm now reaches the same TV as everything else.
#
# d = 1.0 was chosen over d = 1.5 deliberately.  At d = 1.5 the mean arm carries a
# standardized signal of 6.56 and the scale arm 12.16 -- power ~1.00 for every
# method, so the panel stops discriminating.  At d = 1.0 the mean arm scores 4.37,
# high power but not saturated, and the comparison stays informative.
#
# *** THIS INVALIDATES EVERY CHECKPOINT WRITTEN BEFORE 2026-08-13. ***
# All four arms changed, so old and new replications describe different
# alternatives and CANNOT be pooled.  BASIS_SIG does not fingerprint the DGP any
# more than it fingerprints CV grids, so nothing in the code will catch a mixed
# pool -- the pre-2026-08-13 checkpoints must be moved out of this folder.
#
# WHY THE KURTOSIS ARM CHANGED FAMILY.  The old kurtosis arm was Pearson, and a
# Pearson kurtosis shift saturates near TV 0.16 however extreme the kurtosis:
# kurtosis is a TAIL property while TV is a MASS property, and Pearson type VII
# puts its kurtosis where there is no mass to move.  Measured consequence at
# block length 30: the old arm displaced the five fitted quantiles by at most
# 0.141 and carried a standardized QSS signal of 0.98 -- against 4.37 for the
# mean arm.  That, not the monitoring basis or the fusion penalty, is why
# kurtosis detection sat at the nominal line.
#
# The replacement is a SCALE MIXTURE OF NORMALS (contaminated normal).  It raises
# kurtosis by restructuring the BULK, so TV moves with it, and it holds
# mean = 0, variance = 1, skewness = 0 EXACTLY by construction rather than by
# numerical moment-matching:
#   X = [(1-e) N(0,1) + e N(0,c^2)] / sqrt((1-e) + e c^2)
# Both components are centred, so all odd moments vanish; the divisor pins the
# variance.  Only the 4th moment moves.  Standardized signal 2.66 vs 0.98.
#
# CAVEAT to state in the write-up: at FIXED variance, raising kurtosis must shrink
# the IQR (mass leaves the shoulders for the core and the tails).  This arm's IQR
# falls 73%, so its detectability arrives mainly through the QUARTILE contrast,
# not the 4th-order one.  It is a pure fourth-MOMENT change, not a pure
# fourth-CONTRAST change.  The old Pearson arm had the same property, weakly.
#
# S9 (skewed baseline) kurtosis arm: a scale mixture AMPLIFIES skewness by
# E[S^3]/E[S^2]^1.5 >= 1 (Lyapunov), so the mixed Pearson base uses skew 0.2085
# to land the mixture on 0.6.  Verified on 6e5 draws: mean 0.000, var 1.000,
# kurt 20.9, TV 0.383.
TV_ALL_ARMS <- 0.383        # = 2*Phi(d/2) - 1 with d = 1.0002
TV_BY_ARM <- c(null = 0, mean = TV_ALL_ARMS, variance = TV_ALL_ARMS,
               skewness = TV_ALL_ARMS, kurtosis = TV_ALL_ARMS)
TV_DEVIATED      <- TV_ALL_ARMS   # deviated mean, variance, skewness
TV_DEVIATED_KURT <- TV_ALL_ARMS   # no longer capped: the mixture reaches the target
POST_MEAN <- 1.0002    # normal mean arm     (TV = 0.383; Cohen's d = 1.0)
POST_SD   <- 2.3094    # normal variance arm (sd; variance = 5.3334; TV = 0.383)
HK_POST_MEAN <- 0.923235749420706  # HK location arm, from ar_ext11 (innovation TV = 0.383)
HK_POST_SD   <- 2.408834595863566  # HK scale multiplier, from ar_ext11 (innovation TV = 0.383)
TAUS      <- c(0.25, 0.5, 0.75)   # ar_ext12: location/scale grid

# ---------------------------------------------------------------------------
# BASIS FINGERPRINT -- guards against pooling results computed under different
# contrast weights.
# ---------------------------------------------------------------------------
# bqq 0.5.1 redefined the "lmom" basis: it integrates lambda over the WHOLE unit
# interval instead of the truncated grid, so the weights changed (the first
# L-skewness weight moved 0.095133 -> 0.079479) and detections changed with them.
# A 0.5.0 checkpoint and a 0.5.1 checkpoint have the SAME 248 rows and the SAME
# method names, so nothing in the row schema can tell them apart. Every
# checkpoint therefore carries a fingerprint of the actual weight matrices, and
# the supervisor refuses to import a file whose fingerprint differs from the
# running one. This is generic: any future change to any basis invalidates old
# checkpoints automatically, without anyone remembering to bump a constant.
# Study-schema hash (2026-09-05, Codex full-test defect 3): the fingerprint also covers the
# study DESIGN, so a change to the DGP, the config grid, the CV grid, the priors, the EM
# controls or the detection code invalidates old checkpoints without anyone bumping a
# constant. It is computed from the PARSED source of this file (comments and formatting do
# not count) for the top-level objects named below. The supervisor evaluates these same
# two definitions from the driver text, so there is exactly one formula.
SCHEMA_HASH_NAMES <- c(
  "make_configs", "create_config", "get_default_config", "build_cv_grid", "DEFAULT_PRIORS",
  "simulate_data", "draw_baseline", "draw_post", "BASE_PARAMS", "POST_PARAMS", "POST_MEAN",
  "POST_SD", "HK_POST_MEAN", "HK_POST_SD", "PEARSON_SKEW", "KURT_MIX_NORMAL", "KURT_MIX_SKEWED", "TAUS",
  "TV_ALL_ARMS", "TV_BY_ARM", "TV_DEVIATED", "TV_DEVIATED_KURT",
  "IQ_LAMBDA_GRID",
  "CV_LAPLACE_DRAWS", "MAP_OPTS", "BQQ_BASE_FIELDS", "LMOM_METHODS", "CPM_TYPES", "ARIMA_ARGS",
  "run_bqq_detection", "run_single_simulation", "shape_flags", "block_confusion",
  "compute_signal_blocks", "CONFIG_MAP", ".assert_config_map")
.tv11_schema_hash <- function(script, names = SCHEMA_HASH_NAMES) {
  ex <- parse(script, keep.source = FALSE)
  got <- list()
  for (e in ex)
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% names)
      got[[as.character(e[[2]])]] <- paste(deparse(e[[3]], width.cutoff = 500L), collapse = "\n")
  miss <- setdiff(names, names(got))
  if (length(miss)) stop("schema hash: not defined at top level of ", script, ": ", paste(miss, collapse = ", "))
  txt <- paste(names, vapply(names, function(n) got[[n]], ""), sep = " <- ", collapse = "\n\n")
  tf <- tempfile(fileext = ".R"); on.exit(unlink(tf))
  con <- file(tf, open = "wb"); writeLines(txt, con, sep = "\n", useBytes = TRUE); close(con)
  substr(unname(tools::md5sum(tf)), 1L, 8L)
}
# bqq version + L-moment weights (0.6.7: maxent retired) + study-schema hash.
.tv11_basis_sig_from <- function(script, taus) {
  w <- tryCatch(as.numeric(bbq:::.bbq_lmom_weights(taus)), error = function(e) NA_real_)
  if (anyNA(w)) stop("cannot compute the basis fingerprint -- is bbq installed?")
  chk <- sum(round(w, 9) * seq_along(w))          # position-weighted, so a swap shows
  funs <- c("detectChangepoints_gamma", ".bbq_lmom_weights", ".bbq_block_tests",
            ".bbq_adj_family", ".bbq_cell_layer")
  code <- unlist(lapply(funs, function(n) c(n, deparse(getFromNamespace(n, "bbq")))))
  tf <- tempfile(); on.exit(unlink(tf), add = TRUE)
  writeLines(code, tf, useBytes = TRUE)
  engine <- substr(unname(tools::md5sum(tf)), 1L, 8L)
  sprintf("bbq%s-%s-%s-%s", utils::packageVersion("bbq"),
          formatC(chk, format = "f", digits = 6), .tv11_schema_hash(script), engine)
}
# This file's own path: TV11_SCRIPT if set (the supervisor sets it for the child, the smoke
# harness for a definitions-only load), else the script directory.
.tv11_this_script <- local({
  p <- .tv11_env("TV11_SCRIPT", file.path(.tv11_script_dir(), "simulation_arext12.R"))
  if (!file.exists(p) && nzchar(Sys.getenv("TV11_BOX_DIR", "")))
    p <- file.path(Sys.getenv("TV11_BOX_DIR"), "simulation_arext12.R")
  if (!file.exists(p)) stop("cannot locate simulation_arext12.R for the schema hash (set TV11_SCRIPT): ", p)
  normalizePath(p, winslash = "/", mustWork = TRUE)
})
BASIS_SIG <- .tv11_basis_sig_from(.tv11_this_script, TAUS)
cat("  basis fingerprint:", BASIS_SIG, "\n")
writeLines(capture.output(sessionInfo()), file.path(work_dir, "arext12_sessionInfo.txt"))

# Pearson moment -> parameter helper (moments use RAW kurtosis = excess + 3).
.pf <- function(m, v, s, k)
  PearsonDS::pearsonFitM(moments = c(mean = m, variance = v, skewness = s, kurtosis = k))

# Normal-baseline shape arms (named constants; reused by the verification block).
PEARSON_SKEW <- .pf(0, 1, 1.0485, 3)       # skewness arm  (TV = 0.383, kurtosis held at 3)
# PEARSON_KURT retired 2026-08-13: a Pearson kurtosis shift cannot reach TV 0.383
# (it saturates near 0.16). The kurtosis arm is now a scale mixture -- see
# KURT_MIX_* below and the design-point note above.

# Scale-mixture kurtosis arms.  X = S*W / sqrt(E[S^2]) with S in {1, c} w.p. (1-e, e).
#   normal baseline : W ~ N(0,1)                 -> mean 0, var 1, skew 0 EXACTLY
#   skewed baseline : W ~ Pearson(0,1,skew_b,3)  -> skew_b chosen so the MIXTURE
#                     has skew 0.6, since a scale mixture multiplies skewness by
#                     E[S^3]/E[S^2]^1.5 >= 1 (Lyapunov).
KURT_MIX_NORMAL <- list(e = 0.10, c = 7.5165)                     # kurtosis 22.38
KURT_MIX_SKEWED <- list(e = 0.10, c = 6.6696, skew_base = 0.2428)  # kurtosis 20.85

# Pre-change (baseline) distribution per underlying-distribution factor level:
#   normal       : N(0,1)                              -> drawn via rnorm
#   right_skewed : Pearson(0, 1, 0.6, 3.0)  (excess kurt 0;   right-skewed)
# The `platykurtic` baseline was REMOVED 2026-08-13 (user). It was never referenced
# by make_configs() -- the settings are normal x {indep, AR, T, l, isolated} plus
# right_skewed --
# and its arms were never recalibrated to the new design point, so keeping it would
# have left a stale, silently-wrong scenario in the file.
BASE_PARAMS <- list(
  normal       = NULL,
  right_skewed = .pf(0, 1, 0.6, 3.0),
  high_kurtosis = .pf(0, 1, 0, 4)            # HK (ar_ext11)
)

# Post-change (shifted) arm per (baseline, moment). kind = "norm" (mean/sd) or
# "pearson" (Pearson params). Deviated arms calibrated to TV = 0.20 (right-skew
# kurtosis capped at excess 29.99, TV = 0.162).
POST_PARAMS <- list(
  normal = list(
    mean     = list(kind = "norm",    mean = POST_MEAN, sd = 1),
    variance = list(kind = "norm",    mean = 0,         sd = POST_SD),
    skewness = list(kind = "pearson", p = PEARSON_SKEW),
    kurtosis = list(kind = "scalemix", mix = KURT_MIX_NORMAL, base = NULL)
  ),
  right_skewed = list(                            # ALL arms at TV = 0.383
    mean     = list(kind = "pearson", p = .pf(1.002300, 1.000000, 0.600000, 3.000000)),
    variance = list(kind = "pearson", p = .pf(0.000000, 4.420900, 0.600000, 3.000000)),
    skewness = list(kind = "pearson", p = .pf(0.000000, 1.000000, 1.093400, 3.000000)),
    kurtosis = list(kind = "scalemix", mix = KURT_MIX_SKEWED,
                    base = .pf(0, 1, 0.2428, 3))
  ),
  high_kurtosis = list(                           # HK (ar_ext11), TV = 0.383
    mean     = list(kind = "pearson", p = .pf(HK_POST_MEAN, 1, 0, 4)),
    variance = list(kind = "pearson", p = .pf(0, HK_POST_SD^2, 0, 4))
  )
)

# Per-(baseline, shift) TV for the shift_size_tv output column. AR cells carry the
# normal (innovation-level) TV; the marginal TV is diluted by the AR(1) filter.
scenario_tv <- function(distribution, shift_type) {
  if (shift_type == "null") return(0)
  switch(distribution,
    normal       = as.numeric(TV_BY_ARM[shift_type]),
    right_skewed = if (shift_type == "kurtosis") TV_DEVIATED_KURT else TV_DEVIATED,
    high_kurtosis = if (shift_type %in% c("mean", "variance")) TV_ALL_ARMS else NA_real_,
    NA_real_)
}

# ---------------------------------------------------------------------------
# [ar_ext E9] MARGINAL total-variation distance, per (baseline, arm, phi).
# ---------------------------------------------------------------------------
# scenario_tv() above returns the INNOVATION-level TV -- 0.383 in every AR cell,
# because the arms are calibrated on the innovation stream and held fixed across phi
# (GRID.md section 4, closed 2026-08-18). That is the design target. It is NOT what
# the chart sees: the AR(1) filter amplifies a mean shift, leaves the normal-baseline
# variance arm invariant, and DILUTES both shape arms toward Gaussian. At phi = 0.8
# the normal skewness arm is observed at TV 0.072, not 0.383.
#
# Without this column any downstream table reading `shift_size_tv` labels the P08
# skewness arm "TV = 0.383" when the observed distance is a fifth of that. Both
# columns are written: `shift_size_tv` (design) and `shift_size_tv_marginal`
# (observed).
#
# Measured by ar_ext/calibrate_marginal_tv.R -- histogram estimate over 3e6 draws per
# cell, using this file's own draw_baseline()/draw_post() and the same AR filter and
# standardisation as simulate_data(). The normal rows reproduce GRID.md section 4.
#
# SCOPE: this is the STATIONARY post-change marginal for a SUSTAINED shift. For the
# isolated pattern (setting ISO) the process does not fully reach that law inside the
# shifted window, so the value is an UPPER BOUND there. phi = 0 rows equal the
# innovation TV, as they must.
TV_MARGINAL <- data.frame(
  distribution = c(rep("normal", 16), rep("right_skewed", 16)),
  shift_type   = rep(rep(c("mean", "variance", "skewness", "kurtosis"), each = 4), 2),
  ar_phi       = rep(c(0.0, 0.2, 0.5, 0.8), 8),
  tv_marginal  = c(0.3828, 0.4601, 0.6131, 0.8671,    # normal / mean
                   0.3827, 0.3833, 0.3833, 0.3830,    # normal / variance
                   0.3824, 0.3153, 0.1771, 0.0719,    # normal / skewness
                   0.3828, 0.3567, 0.3003, 0.1895,    # normal / kurtosis
                   0.3835, 0.4618, 0.6147, 0.8686,    # right_skewed / mean
                   0.3835, 0.3775, 0.3606, 0.3493,    # right_skewed / variance
                   0.3820, 0.2874, 0.1341, 0.0419,    # right_skewed / skewness
                   0.3826, 0.3554, 0.2898, 0.1766),   # right_skewed / kurtosis
  stringsAsFactors = FALSE)

# phi = 0 is encoded as dependence == "independent" (see the encoding trap in
# GRID.md section 0), so read the level from `dependence` first and only then ar_phi.
scenario_tv_marginal <- function(distribution, shift_type, dependence, ar_phi) {
  if (identical(shift_type, "null")) return(0)
  phi <- if (identical(dependence, "independent")) 0 else as.numeric(ar_phi)
  i <- which(TV_MARGINAL$distribution == distribution &
             TV_MARGINAL$shift_type   == shift_type &
             abs(TV_MARGINAL$ar_phi - phi) < 1e-9)
  if (length(i) != 1L) return(NA_real_)
  TV_MARGINAL$tv_marginal[i]
}

# Draw n observations from the pre (baseline) distribution.
draw_baseline <- function(distribution, n) {
  if (n <= 0) return(numeric(0))
  if (distribution == "normal") return(rnorm(n, 0, 1))
  PearsonDS::rpearson(n, BASE_PARAMS[[distribution]])
}

# Draw n observations from the post (shifted) arm; null == baseline.
draw_post <- function(distribution, shift_type, n) {
  if (n <= 0) return(numeric(0))
  if (shift_type == "null") return(draw_baseline(distribution, n))
  spec <- POST_PARAMS[[distribution]][[shift_type]]
  if (is.null(spec)) stop("No post arm for ", distribution, "/", shift_type)
  if (spec$kind == "norm") rnorm(n, spec$mean, spec$sd)
  else if (spec$kind == "scalemix") .draw_scalemix(spec, n)
  else PearsonDS::rpearson(n, spec$p)
}

# Scale mixture of a centred base:  X = S*W / sqrt(E[S^2]),  S in {1, c} wp (1-e, e).
# Both components are centred, so the mean is 0 and (for a symmetric base) every odd
# moment vanishes; the divisor pins the variance at 1. Only the shape moves.
.scalemix_s2 <- function(mix) sqrt((1 - mix$e) + mix$e * mix$c^2)
.draw_scalemix <- function(spec, n) {
  mix <- spec$mix
  w <- if (is.null(spec$base)) rnorm(n) else PearsonDS::rpearson(n, spec$base)
  S <- ifelse(stats::runif(n) < mix$e, mix$c, 1)
  S * w / .scalemix_s2(mix)
}
.scalemix_density <- function(spec) {
  mix <- spec$mix; s2 <- .scalemix_s2(mix)
  fb <- if (is.null(spec$base)) function(x) stats::dnorm(x) else
        function(x) PearsonDS::dpearson(x, params = spec$base)
  function(x) ((1 - mix$e) * fb(x * s2) + mix$e * fb(x * s2 / mix$c) / mix$c) * s2
}

# Analytic pre / post densities (for the overlap-MC verification block).
pre_density <- function(distribution) {
  if (distribution == "normal") return(function(x) dnorm(x, 0, 1))
  p <- BASE_PARAMS[[distribution]]; function(x) PearsonDS::dpearson(x, params = p)
}
post_density <- function(distribution, shift_type) {
  if (shift_type == "null") return(pre_density(distribution))
  spec <- POST_PARAMS[[distribution]][[shift_type]]
  if (spec$kind == "norm") function(x) dnorm(x, spec$mean, spec$sd)
  else if (spec$kind == "scalemix") .scalemix_density(spec)
  else function(x) PearsonDS::dpearson(x, params = spec$p)
}

# ---------------------------------------------------------------------------
# SECTION 2A: startup DGP verification block
# ---------------------------------------------------------------------------
# Draws a large sample from each arm (empirical moments) and computes the TV vs
# N(0,1) via the overlap estimator TV = 1 - E_{X~post}[min(1, f_pre(X)/f_post(X))],
# which is robust to the Pearson density's boundary singularity (the old fine-grid
# integrator returned Inf for the higher-skew arm). Returns (invisibly) the table.
verify_dgp <- function(N = 3e5, seed = 1) {
  set.seed(seed)
  emp_moments <- function(z) {
    m <- mean(z); v <- mean((z - m)^2)
    c(mean = m, variance = v,
      skewness = mean((z - m)^3) / v^1.5,
      exkurt   = mean((z - m)^4) / v^2 - 3)
  }
  # overlap-MC TV between an arm's post and its baseline pre (X ~ pre, so smooth).
  tv_arm <- function(distribution, shift_type) {
    fpre <- pre_density(distribution); fpost <- post_density(distribution, shift_type)
    X  <- draw_baseline(distribution, N)
    fp <- fpre(X); ok <- is.finite(fp) & fp > 0
    1 - mean(pmin(1, fpost(X[ok]) / fp[ok]))
  }
  cat("================================================================\n")
  cat("DGP VERIFICATION  (per-arm target TV vs overlap-MC estimate)\n")
  cat("================================================================\n")
  cat(sprintf("%-13s %-9s %8s %8s %8s %9s %8s %8s\n",
              "baseline", "shift", "mean", "var", "skew", "exkurt", "TV_tgt", "TV_est"))
  shifts <- c("null", "mean", "variance", "skewness", "kurtosis")
  out <- list()
  for (dist in c("normal", "right_skewed")) {
    for (st in shifts) {
      z  <- draw_post(dist, st, N)
      mo <- emp_moments(z)
      tt <- scenario_tv(dist, st)
      te <- if (st == "null") 0 else tv_arm(dist, st)
      cat(sprintf("%-13s %-9s %8.4f %8.4f %8.4f %9.3f %8.3f %8.4f\n",
                  dist, st, mo["mean"], mo["variance"], mo["skewness"], mo["exkurt"], tt, te))
      out[[paste(dist, st, sep = "/")]] <- c(mo, TV_target = tt, TV = te)
    }
  }
  cat("----------------------------------------------------------------\n")
  cat("all baselines: mean/var/skew TV=0.38 | kurtosis capped: normal 0.11, deviated 0.162\n")
  cat("        (a pure kurtosis shift saturates and cannot reach TV=0.38).\n")
  cat("NOTE: sample excess kurtosis is severely downward-biased for very heavy tails;\n")
  cat("      the overlap-MC TV is the authoritative calibration check.\n")
  cat("================================================================\n\n")
  invisible(out)
}

# ---------------------------------------------------------------------------
# SECTION 3: Data generation -- distribution x dependence x shift
# ---------------------------------------------------------------------------
# distribution in {normal, right_skewed} sets the pre/post shape;
# dependence in {independent, ar} sets the temporal structure. Under AR the (pre,
# post) draws are the AR(1) innovations, so the shift enters the innovation stream.
#' @return list(y, shift_point, shift_end, params). shift_type == "null" has no shift
#'   region (shift_end = NA) and the whole series is the baseline distribution.
simulate_data <- function(n, distribution = "normal", dependence = "independent",
                          shift_type, shift_start = 0.8, shift_pattern = "sustained",
                          shift_length = 30, seed = NULL, ar_phi = 0.5) {
  if (!is.null(seed)) set.seed(seed)
  shift_point <- round(shift_start * n)

  if (shift_type == "null") {
    idx_post  <- integer(0)
    shift_end <- NA_integer_
  } else if (shift_pattern == "isolated") {
    shift_end <- min(shift_point + shift_length - 1L, n)
    idx_post  <- shift_point:shift_end
  } else {                                  # sustained
    shift_end <- n
    idx_post  <- shift_point:n
  }
  idx_pre <- setdiff(seq_len(n), idx_post)

  # Driving series s: pre = baseline, post = shifted arm (both per `distribution`).
  s <- numeric(n)
  s[idx_pre]  <- draw_baseline(distribution, length(idx_pre))
  s[idx_post] <- draw_post(distribution, shift_type, length(idx_post))

  if (dependence == "independent") {
    y <- s
  } else if (dependence == "ar") {
    # AR(1): y_t = phi*y_{t-1} + s_t; stationary init; standardized to ~unit pre variance.
    ar_sd <- 1 / sqrt(1 - ar_phi^2)
    y <- numeric(n)
    y[1] <- rnorm(1, 0, ar_sd)
    for (t in 2:n) y[t] <- ar_phi * y[t - 1] + s[t]
    y <- y / ar_sd
  } else {
    stop("Unknown dependence: ", dependence, ". Use 'independent' or 'ar'.")
  }

  list(y = y, shift_point = shift_point, shift_end = shift_end,
       params = list(n = n, distribution = distribution, dependence = dependence,
                     shift_type = shift_type, shift_start = shift_start,
                     shift_pattern = shift_pattern, shift_point = shift_point,
                     shift_end = shift_end))
}

# ---------------------------------------------------------------------------
# SECTION 4: Block geometry + signal / confusion helpers (CHANGE 4)
# ---------------------------------------------------------------------------
# Block b covers observations [w + (b-1)*l + 1, min(w + b*l, n)], with the last
# block extending to n. r = ceil((n - w)/l). Matches getSustainedShift/Isolated
# and the package's h_to_obs() mapping in detectChangepoints_gamma.
block_bounds <- function(n, w, l) {
  r <- as.integer((n - w) / l) + as.integer(((n - w) %% l) > 0L)
  b <- seq_len(r)
  data.frame(block = b,
             obs_start = w + (b - 1L) * l + 1L,
             obs_end   = pmin(w + b * l, n))
}

# A block is a SIGNAL block if any of its observations lies within `tolerance`
# of the ONSET, i.e. inside the FIXED window [shift_point - tol, shift_point + tol];
# none for the null arm.
#
# CHANGED 2026-08-04.  This previously used `hi <- shift_end + tolerance`, which
# made the truth window depend on the block length:
#   * sustained -> shift_end = n, so every block from the onset to the end of the
#     series was truth-positive.  At tau = 292, n = 365, w = 30 that is 4 signal
#     blocks at l = 30, 6 at l = 14, 11 at l = 7.
#   * isolated  -> shift_end = shift_point + shift_length - 1, and shift_length is
#     set to block_length, so that window scaled with l as well.
# Because a cpm run contributes a single change point, cpm_block_flags() can flag
# at most ONE block, so cpm sensitivity was capped at 1/#signal blocks =
# 0.250 / 0.167 / 0.091 for l = 30 / 14 / 7 -- an arithmetic penalty that grew as l
# shrank and had nothing to do with detection ability.  The window is now anchored
# to the onset and is 2*tolerance + 1 = 15 observations wide for every setting.
compute_signal_blocks <- function(bounds, shift_point, shift_end, shift_type, tolerance) {
  r <- nrow(bounds)
  if (shift_type == "null") return(rep(FALSE, r))
  lo <- shift_point - tolerance
  hi <- shift_point + tolerance      # was shift_end + tolerance (l-dependent)
  (bounds$obs_end >= lo) & (bounds$obs_start <= hi)
}

# Block-level confusion + onset localization for a per-block flag vector.
#   is_signal  : truth-positive blocks (shift region +/- tolerance; region-based).
#   tolerance  : half-width of the ONSET window [shift_point +/- tolerance] used by
#                hit_within_tol and (indirectly) detection rate.
# closest_signal_time = obs_start of the flagged block nearest the onset tau (the change
#   estimate); detection_delay = closest_signal_time - tau (signed); hit_within_tol = 1
#   if any flagged block overlaps the onset window. For the null arm there is no onset,
#   so closest_signal_time / detection_delay / hit_within_tol are NA.
block_confusion <- function(flags, is_signal, bounds, shift_point, shift_type, tolerance) {
  flags   <- as.logical(flags)
  is_null <- shift_type == "null"
  flagged <- which(flags)
  if (length(flagged) > 0 && !is_null) {
    # distance of each block to tau (0 if tau lies within the block)
    d_to_tau <- pmax(0, bounds$obs_start - shift_point, shift_point - bounds$obs_end)
    nearest  <- flagged[which.min(d_to_tau[flagged])]
    closest_signal_time <- bounds$obs_start[nearest]
    detection_delay     <- closest_signal_time - shift_point
    overlaps_win <- (bounds$obs_end   >= shift_point - tolerance) &
                    (bounds$obs_start <= shift_point + tolerance)
    hit_within_tol <- as.integer(any(flags & overlaps_win))
  } else {
    closest_signal_time <- NA_real_
    detection_delay     <- NA_real_
    hit_within_tol      <- if (is_null) NA_integer_ else 0L
  }
  list(
    TP = sum(flags & is_signal),
    FP = sum(flags & !is_signal),
    FN = sum(!flags & is_signal),
    TN = sum(!flags & !is_signal),
    any_signal          = as.integer(any(flags)),
    n_flagged           = sum(flags),
    closest_signal_time = closest_signal_time,
    detection_delay     = detection_delay,
    hit_within_tol      = hit_within_tol,
    # ADDED 2026-08-04.  The WHICH, not just the how-many: "|"-separated indices of
    # every flagged block ("" when nothing flagged).  Until now only counts were
    # persisted, so any change to the truth window (compute_signal_blocks) forced a
    # full re-run of the fits.  With the flag set stored, TP/FP/TN/FN and every
    # metric derived from them can be RE-SCORED offline from the task files:
    #   flags <- seq_len(n_blocks) %in% as.integer(strsplit(flagged_blocks,"|",fixed=TRUE)[[1]])
    # and block geometry is reproducible from (n_obs, warmup, block_length).
    flagged_blocks      = if (length(flagged)) paste(flagged, collapse = "|") else ""
  )
}

# ---------------------------------------------------------------------------
# SECTION 5: cpm baselines (REUSED verbatim from v3)
# ---------------------------------------------------------------------------
CPM_TYPES <- c("Student", "Bartlett", "GLR",
               "Mann-Whitney", "Mood", "Lepage",
               "Kolmogorov-Smirnov", "Cramer-von-Mises")

run_cpm_batch <- function(y, shift_point, tolerance, cpmType,
                          alpha = 0.05, lambda = NA) {
  variant_label <- paste0("cpm-", cpmType)
  tryCatch({
    res <- cpm::detectChangePointBatch(x = y, cpmType = cpmType,
                                       alpha = alpha, lambda = lambda)
    cp_locations <- if (isTRUE(res$changeDetected) &&
                        !is.null(res$changePoint) &&
                        is.finite(res$changePoint)) {
      as.integer(res$changePoint)
    } else integer(0)
    list(method = variant_label, detected_points = cp_locations,
         n_detections = length(cp_locations), error = NULL)
  }, error = function(e) {
    list(method = variant_label, detected_points = integer(0),
         n_detections = 0, error = conditionMessage(e))
  })
}

run_cpm_variants <- function(y, shift_point, tolerance) {
  results <- list()
  for (t in CPM_TYPES) {
    results[[t]] <- run_cpm_batch(y, shift_point, tolerance,
                                  cpmType = t, alpha = 0.05, lambda = NA)
  }
  results
}

# Map a cpm changepoint to the block containing it (comparable block scheme).
# A cp inside the warm-up (<= w, no block) flags nothing.
cpm_block_flags <- function(cp, bounds) {
  flags <- rep(FALSE, nrow(bounds))
  if (length(cp) == 1L && is.finite(cp)) {
    b <- which(bounds$obs_start <= cp & cp <= bounds$obs_end)
    if (length(b) == 1L) flags[b] <- TRUE
  }
  flags
}

# ---------------------------------------------------------------------------
# SECTION 5B: reference R2 -- cpm on auto.arima residuals      [ar_ext E5, NEW]
# ---------------------------------------------------------------------------
# The standard "prewhiten, then chart" remedy for autocorrelation, and the steel-man
# for this study's headline claim: under AR(1) the raw-series charts blow their
# false-alarm rate, and R2 asks whether prewhitening repairs it.
#
# d = 0 IS A DIRECTIVE (user, 2026-08-18), not a default. The DGP is a stationary
# AR(1) with no unit root, so any d > 0 is the sustained level shift tripping the
# stationarity test rather than a property of the process. Left free, auto.arima
# selects (1,1,1) on the mean arm and differences the shift away. Measured over 300
# replications at phi = 0.8, standardised post-change residual mean against an
# innovation truth of 1.000:  free 0.224  ->  d = 0  0.341. The shape arms are
# unaffected -- they already selected d = 0. See ar_ext/r2_d0.R.
#
# The model is fitted to the WHOLE series, contamination included, because that is
# what a practitioner does (GRID.md OPEN-A, closed 2026-08-18). Under the null the
# whole-series and Phase-I fits coincide, so the SIZE curve -- the result this study
# exists to produce -- is unaffected by that choice.
#
# CAVEAT to carry into the write-up, not to repair: because the shape arms dilute in
# the MARGINAL but not in the innovations, prewhitening approximately UNDOES the
# dilution (post-change skewness 0.352 in y vs 0.999 in these residuals at phi = 0.8).
# On the shape arms at high phi, R2 and BQQ therefore read different signal strengths
# off the same replication; that is a comparison of streams, not of methods.
ARIMA_ARGS <- list(d = 0, seasonal = FALSE, stepwise = TRUE, approximation = FALSE)

# [ar_ext E10, 2026-08-19] The auto.arima fit is now factored out so BOTH prewhitened
# references -- cpm AND BQQ -- read the SAME residual series from ONE fit. Fitting twice
# would be wasteful and, worse, the two streams could diverge if auto.arima's stepwise
# search broke a tie differently, so the four families would no longer be paired.
fit_arima_once <- function(y) {
  fit <- tryCatch(do.call(forecast::auto.arima, c(list(y = y), ARIMA_ARGS)),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(e = NULL, err = "auto.arima failed",
                                p = NA_integer_, q = NA_integer_,
                                ar1 = NA_real_, has_mean = NA_integer_))
  e <- tryCatch(as.numeric(stats::residuals(fit)), error = function(e) NULL)
  if (is.null(e) || !all(is.finite(e)))
    return(list(e = NULL, err = "non-finite arima residuals",
                p = NA_integer_, q = NA_integer_, ar1 = NA_real_, has_mean = NA_integer_))
  ord <- tryCatch(forecast::arimaorder(fit), error = function(e) rep(NA_integer_, 3))
  cf  <- tryCatch(stats::coef(fit), error = function(e) numeric(0))
  list(e = e, err = NULL,
       p = as.integer(ord[1]), q = as.integer(ord[3]),
       ar1 = if ("ar1" %in% names(cf)) as.numeric(cf[["ar1"]]) else NA_real_,
       has_mean = as.integer(any(c("mean", "intercept", "drift") %in% names(cf))))
}

run_arima_cpm_variants <- function(y, shift_point, tolerance, af = NULL) {
  na_out <- function(msg) {
    res <- stats::setNames(lapply(CPM_TYPES, function(tt)
      list(method = paste0("arima_cpm-", tt), detected_points = integer(0),
           n_detections = 0, error = msg)), CPM_TYPES)
    list(results = res, p = NA_integer_, q = NA_integer_,
         ar1 = NA_real_, has_mean = NA_integer_)
  }
  if (is.null(af)) af <- fit_arima_once(y)
  if (!is.null(af$err)) return(na_out(af$err))
  res <- list()
  for (tt in CPM_TYPES)
    res[[tt]] <- run_cpm_batch(af$e, shift_point, tolerance,
                               cpmType = tt, alpha = 0.05, lambda = NA)
  list(results = res, p = af$p, q = af$q, ar1 = af$ar1, has_mean = af$has_mean)
}

# ---------------------------------------------------------------------------
# SECTION 6: CV grid builder -- ar_ext12: lambda_iq IS ON THE GRID, NO EM
# ---------------------------------------------------------------------------
# Prior shrinkage: spike_slab / spike_slab_lasso search spike_sd; the LASSO-type priors
# search lambda_lasso2_b (grids unchanged from ar_ext5). lambda_nc is a single value.
# lambda_iq is a FIXED penalty weight per grid row, tuned by the same score CV:
IQ_LAMBDA_GRID <- c(0, 1, 10, 100, 1000)   # lambda_iq values; 0 = no interquantile fusion.
                                           # getModel() takes lambda_iq2 = IQ_LAMBDA_GRID^2.
                                           # Part of SCHEMA_HASH_NAMES: change it and the
                                           # fingerprint changes, so old checkpoints cannot pool.
# TV11_VERBOSE=1 (set by smoke_arext12.R) prints one timestamped line per model attempt,
# per CV grid row, and at the start/end of each final fit. Off in the study.
TV11_VERBOSE <- nzchar(Sys.getenv("TV11_VERBOSE", ""))
vsay <- function(fmt, ...) {
  if (!TV11_VERBOSE) return(invisible(NULL))
  msg <- paste0(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(fmt, ...), "\n")
  cat(msg, sep = "")
  # on a PSOCK worker stdout is not captured; also append to that worker's log file
  wl <- get0(".tv11_wlog", envir = .GlobalEnv, inherits = FALSE)
  if (is.character(wl) && length(wl) == 1L && nzchar(wl)) try(cat(msg, file = wl, append = TRUE, sep = ""), silent = TRUE)
  invisible(NULL)
}
CV_LAPLACE_DRAWS <- 200L       # Laplace draws inside CV fits: nothing reads them (no EM,
                               # no detection on CV fits), so the count only bounds a
                               # throw-away cost. The final fit uses 20000.

# ---------------------------------------------------------------------------
# MAP optimizer controls -- rstan::optimizing's OWN defaults (user directive,
# 2026-08-02), replacing the deliberately tight bqq >= 0.4.5 defaults, to cut
# fit time. Applied to BOTH the CV fits and the final fit so they cannot drift.
#
#   argument           bqq default   -> set here (rstan default)
#   map_tol_rel_grad   1e2              1e7
#   map_tol_rel_obj    1e2              1e4
#   map_history_size   25               5
#   map_iter           10000            2000
#   map_tol_obj/grad/param  1e-12/1e-8/1e-8  unchanged (already rstan's values)
#
# NOTE (bqq ?getModel): the tight defaults exist because "looser relative
# tolerances (e.g., rstan's default tol_rel_obj = 1e4) can terminate ... while
# reporting convergence, leaving outer-quantile coefficients stranded far from
# the optimum", and a larger map_history_size "typically reduc[es] the iteration
# count without loosening any stopping tolerance". Replications produced under
# these settings are therefore NOT numerically identical to earlier ones.
# CALIBRATED 2026-08-02 by direct measurement (n=365 variance arm, bqq 0.4.10).
#
# (a) The FULL rstan default bundle (rel_grad 1e7, rel_obj 1e4, history 5,
#     iter 2000) is 3.6-9.0x faster but CHANGES THE RESULTS -- it stops on a
#     plateau with coefficients stranded, while reporting return_code 0:
#       adaptive_lasso    logp 357.48 ->  290.99 :   0 flags ->  50   CHANGED
#       spike_slab_lasso  logp  14.25 ->   11.59 :  43 flags ->  22   CHANGED
#     Rejected. The destructive parts are map_history_size and map_iter.
#
# (b) Loosening ONLY the relative tolerances, keeping history 25 / iter 10000,
#     is FREE: 15/15 comparisons across 3 cases and 2 independent draws came
#     back bit-identical (logp delta +0.00e+00, same max|gamma|, same flags):
#       adaptive_lasso  seed  7: 289.1s -> 109.5s (2.64x)  identical
#       adaptive_lasso  seed 11: 106.5s -> 105.1s (1.01x)  identical
#       spike_slab_lasso seed 7:  24.3s ->  24.8s (0.98x)  identical
#     The gain is case-dependent (only fits grinding against the gradient test
#     benefit), so expect a modest average speedup, not 2.6x.
#
# map_tol_rel_grad is the ONLY value changed from the bqq default.
# ---------------------------------------------------------------------------
# TIMEOUT_SCALE -- how many times SLOWER this machine is than the reference box
# the per-model timeout anchor (1900 s) was calibrated on. Purely a give-up
# deadline: it cannot alter any fit that succeeds, it only decides how long to
# wait before declaring failure. Set too low on a slow machine, healthy fits are
# killed and retried 6x each, which is what stranded the S2 adaptive_lasso cells.
#   reference desktop -> 1        4x-slower workstation -> 4
TIMEOUT_SCALE <- {
  v <- suppressWarnings(as.numeric(Sys.getenv("TV11_TIMEOUT_SCALE", "1")))
  if (is.na(v) || v <= 0) 1 else v
}
cat("  TIMEOUT_SCALE      :", TIMEOUT_SCALE,
    "(per-model timeout multiplier; set TV11_TIMEOUT_SCALE)\n")

MAP_OPTS <- list(
  map_tol_obj      = 1e-12,   # bqq + rstan default
  map_tol_grad     = 1e-8,    # bqq + rstan default
  map_tol_param    = 1e-8,    # bqq + rstan default
  map_tol_rel_grad = 1e4,     # CHANGED from 1e2 (rstan's is 1e7) -- results identical
  map_tol_rel_obj  = 1e2,     # bqq default kept; rstan's 1e4 is what strands fits
  map_history_size = 25L,     # bqq default kept; cutting to 5 was destructive
  map_iter         = 10000L   # bqq default kept; cutting to 2000 was destructive
)
build_cv_grid <- function(prior, config) {
  shrink <- if (prior %in% c("spike_slab", "spike_slab_lasso")) list(spike_sd = config$spike_sd_grid)
            else                                                list(lambda_lasso2_b = config$lambda_lasso2_b_grid)
  # ar_ext12: lambda_iq2 IS in the grid (IQ_LAMBDA_GRID^2), crossed with the shrinkage grid.
  do.call(expand.grid,
          c(list(lambda_nc = config$lambda_nc_grid), shrink, list(lambda_iq2 = IQ_LAMBDA_GRID^2)))
}

# ---------------------------------------------------------------------------
# SECTION 7: (removed) in-script four-contrast QSS calibration.
# ---------------------------------------------------------------------------
# The QSS shape test is now OWNED BY THE PACKAGE. detectChangepoints_gamma() with
# With three quantiles, basis = "qss" builds two contrasts (median, IQR), studentizes each cell,
# whitens across cells (z-tilde = Sigma^{-1/2} gammabar), and returns the
# union-intersection (UI) and Hotelling-T^2 tests with the full {raw, Bonferroni,
# Holm, BH, calibrated} adjustment family. See run_bqq_detection() below.

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# SECTION 7B: L-moment basis -- OWNED BY THE PACKAGE
# ---------------------------------------------------------------------------
# The L-moment family used to be reimplemented here. As of bqq 0.5.0 it is a package
# basis: detectChangepoints_gamma(basis = "all") builds the quantile, QSS and L-moment
# contrasts and scores all three with the same whitening (z~ = Sigma^{-1/2} cbar), the
# same UI / Hotelling-T^2 statistics and the same {raw, Bonferroni, Holm, BH, calibrated}
# family. The maximum-entropy basis was retired with bqq 0.6.7 (2026-09-05). The in-script
# replicas have been DELETED: nothing here recomputes a contrast, a weight matrix or a test.
#
# The package derivation (shifted-Legendre quadrature with the truncation projection)
# is documented in MONITORING_BASES_math.md.
#
# All that remains is a reshape: the package returns flagged block INDICES at
# det$tests[[fam]][[stat]][[adjust]], and the results schema wants a logical vector
# per method under the existing "<fam>_<stat>_<adjust>_w" naming.
shape_flags <- function(det, fam, r) {
  out <- list()
  for (st in c("ui", "t2")) {
    node <- det$tests[[fam]][[if (st == "ui") "ui" else "hotelling_t2"]]
    if (is.null(node))
      stop("detectChangepoints_gamma() returned no '", fam, "' / '", st,
           "' member; was it called with basis = \"all\" and both statistics?")
    for (ad in c("raw", "bonf", "holm", "bh", "calib"))
      out[[paste0(fam, "_", st, "_", ad, "_w")]] <- seq_len(r) %in% node[[ad]]
  }
  out
}

# ---------------------------------------------------------------------------
# SECTION 8: BQQ detection -- ONE MAP+Laplace fit; both bases x both statistics
# ---------------------------------------------------------------------------
# One BQQ model per prior. Both shrinkage and lambda_iq2 are score-CV tuned.
# Every tuning and final fit uses adaptive_iq = FALSE.
run_bqq_detection <- function(y, config, prior, sim_seed) {
  n     <- length(y)
  l     <- config$block_length
  w     <- config$warmup_period

  H <- if (config$shift_pattern == "isolated")
    getIsolatedShift(n = n, l = l, w = w)
  else
    getSustainedShift(n = n, l = l, w = w)

  # ---- cross-validation over (lambda_nc, prior shrinkage) ----
  cv_grid <- build_cv_grid(prior, config)
  # ar_ext12: lambda_iq2 IS in cv_grid (IQ_LAMBDA_GRID^2 x the shrinkage grid). NO EM:
  # every tuning fit is one MAP fit at its row's (lambda_nc, shrinkage, lambda_iq2).
  # beta0_loc / beta0_scale are left at their bqq >= 0.4.5 defaults: warm-up
  # empirical-quantile centering with unit-information scale (computed per
  # replication from each simulated series' own warm-up window).
  # MAP optimizer controls: MAP_OPTS (rstan defaults) -- see the constant.
  cv_base_args <- c(list(
    adaptive_gamma    = TRUE,
    prior_gamma       = prior,
    adaptive_iq       = FALSE,             # ar_ext12: no EM in the CV fits
    laplace_n_samples = CV_LAPLACE_DRAWS
  ), MAP_OPTS)
  # [ar_ext5 CHANGE 6] loss = "score" is passed EXPLICITLY, not left to the bqq
  # 0.6.0 default (also "score"), so the selection criterion is visible here and
  # survives any future change to the package default. One sweep computes both
  # criteria; val_score does the ranking, val_pinball rides along for reporting.
  vsay("  %s | n = %d, r = %d | CV start: %d grid points x 2 folds, no EM (lambda_iq on the grid)",
       prior, n, ncol(H), nrow(cv_grid))
  t_cv <- Sys.time()
  cv <- cv_copss_grid(y = y, taus = TAUS, H = H, X = NULL, w = w,
                      grid = cv_grid, base_args = cv_base_args,
                      loss = "score",
                      seed = sim_seed, verbose = TV11_VERBOSE)
  vsay("  %s | CV done in %.1f min", prior, as.numeric(difftime(Sys.time(), t_cv, units = "mins")))
  best <- cv[1, , drop = FALSE]            # min val_loss (returned sorted ascending)
  best_lambda_nc <- best$lambda_nc[1]
  val_loss       <- best$val_loss[1]       # == val_score: the selection criterion
  val_score      <- best$val_score[1]
  val_pinball    <- best$val_pinball[1]    # pinball AT THE SCORE-SELECTED config
  is_spike       <- prior %in% c("spike_slab", "spike_slab_lasso")
  best_spike_sd  <- if (is_spike) best$spike_sd[1] else NA_real_
  best_lasso_b   <- if (is_spike) NA_real_ else best$lambda_lasso2_b[1]
  best_lambda_iq2 <- best$lambda_iq2[1]   # ar_ext12: the CV winner (IQ_LAMBDA_GRID^2), not an EM value
  # ---- final MAP fit + Laplace posterior draws (single fit) ----
  final_args <- c(list(
    y = y, taus = TAUS, H = H, X = NULL, w = w,
    lambda_nc         = best_lambda_nc,
    adaptive_gamma    = TRUE,
    # ar_ext12: lambda_iq2 is the CV winner, held FIXED; no EM in the final fit either.
    # Coefficients start from the LASSO pilot (Sec. 2.4), as in every CV fit.
    lambda_iq2        = best_lambda_iq2,
    adaptive_iq       = FALSE,
    prior_gamma       = prior,
    fit_method        = "map",
    map_hessian       = TRUE,
    # RAISED 1000 -> 10000 (2026-08-09).  .bqq_block_tests() whitens with a Sigma
    # ESTIMATED from these draws, and Sigma^{-1/2} amplifies the noisiest SMALL
    # eigenvalues, so achieved FAP exceeds nominal even when the posterior is
    # perfectly calibrated.  Measured (ncell=4, alpha=.05, Bonferroni, nominal .05):
    #   r=12 d=48  n_iter= 1000  d/n=.048 -> FAP .068   <- what the 08-06 run used
    #   r=12 d=48  n_iter= 4000  d/n=.012 -> FAP .060
    #   r=12 d=48  n_iter=16000  d/n=.003 -> FAP .046   <- nominal
    # At 20000, d/n = .0024 -- past the 16000 row above, i.e. at nominal.
    laplace_n_samples = 20000,
    seed              = sim_seed,
    verbose           = TV11_VERBOSE
  ), MAP_OPTS)
  if (is_spike) final_args$spike_sd <- best_spike_sd
  else          final_args$lambda_lasso2_b <- best_lasso_b
  vsay("  %s | final fit start: winner %s, lambda_iq = %g (CV), 20000 draws", prior,
       if (is_spike) sprintf("spike_sd = %s", best_spike_sd) else sprintf("lambda_lasso2_b = %s", best_lasso_b), sqrt(best_lambda_iq2))
  t_fit <- Sys.time()
  fit <- do.call(getModel, final_args)
  vsay("  %s | final fit done in %.1f min: lambda_iq = %g (CV winner), init = %s", prior,
       as.numeric(difftime(Sys.time(), t_fit, units = "mins")), sqrt(best_lambda_iq2),
       if (is.null(fit$map_init_used)) "?" else fit$map_init_used)

  # ---- ar_ext12: no EM. best_lambda_iq2 is the CV winner (set above); the two EM
  #      diagnostics are NA and kept only so the row schema matches ar_ext5's. ----
  iq_em_converged  <- NA
  iq_em_iters      <- NA_integer_

  # ---- detection from the ONE shared fit: both bases (quantile / qss) x both
  #      statistics (ui / hotelling_t2). ONE analytic detection call; posterior
  #      whitening (z-tilde = Sigma^{-1/2} gammabar) is ALWAYS applied -> 20 methods.
  #      The calibrated member is the analytic single-step max. ----
  common <- list(fit, taus = TAUS, l = l, w = w, signal_position = "first", y = y,
                 alpha = 0.05, basis = "all",                # quantile + qss + lmom (0.6.7)
                 statistic = c("ui", "hotelling_t2"), seed = sim_seed)
  det <- do.call(detectChangepoints_gamma, common)
  stopifnot(identical(as.integer(dim(det$tests$quantile$z)), c(3L, ncol(H))),
            identical(as.integer(dim(det$tests$qss$z)), c(2L, ncol(H))),
            identical(as.integer(dim(det$tests$lmom$z)), c(2L, ncol(H))))
  db  <- det$detected_blocks

  # ---- L-moment family, read off the SAME detection call. No Laplace draws are
  #      touched here: the package built and scored these. ----
  lm_flags <- shape_flags(det, "lmom", r = nrow(db))

  read_flags <- function(db, suffix)
    stats::setNames(lapply(BQQ_BASE_FIELDS, function(col) as.logical(db[[col]])),
                    paste0(names(BQQ_BASE_FIELDS), suffix))

  list(
    r = nrow(db),
    flags = c(read_flags(db, "_w"), lm_flags),
    best_lambda_nc  = best_lambda_nc,
    best_spike_sd   = best_spike_sd,
    best_lambda_iq2 = best_lambda_iq2,
    iq_em_converged = iq_em_converged,
    iq_em_iters     = iq_em_iters,
    val_loss        = val_loss,
    val_score       = val_score,
    val_pinball     = val_pinball,
    error = NULL
  )
}

# Base BQQ charts: 4 charts (quantile / qss x UI / Hotelling-T^2) crossed with ALL 5
# adjustments (raw / Bonferroni / Holm / BH / calibrated single-step max) -> 20 methods,
# each a detected_blocks column. Posterior whitening (z-tilde = Sigma^{-1/2} gammabar)
# is always applied; the methods retain the "_w" suffix.
BQQ_BASE_FIELDS <- c(
  q_ui_raw     = "significant_raw",          q_ui_bonf    = "significant_bonf",
  q_ui_holm    = "significant_holm",         q_ui_bh      = "significant_bh",
  q_ui_calib   = "significant_calib",
  q_t2_raw     = "significant_wald_raw",      q_t2_bonf    = "significant_wald_bonf",
  q_t2_holm    = "significant_wald_holm",     q_t2_bh      = "significant_wald_bh",
  q_t2_calib   = "significant_wald_calib",
  qss_ui_raw   = "significant_qss_raw",       qss_ui_bonf  = "significant_qss_bonf",
  qss_ui_holm  = "significant_qss_holm",      qss_ui_bh    = "significant_qss_bh",
  qss_ui_calib = "significant_qss_calib",
  qss_t2_raw   = "significant_qss_t2_raw",    qss_t2_bonf  = "significant_qss_t2_bonf",
  qss_t2_holm  = "significant_qss_t2_holm",   qss_t2_bh    = "significant_qss_t2_bh",
  qss_t2_calib = "significant_qss_t2_calib"
)
# The 10 ADDED L-moment methods (SECTION 7B). Same naming convention; they are not
# detected_blocks columns, they come from shape_flags() reading det$tests.
LMOM_METHODS   <- as.vector(outer(c("ui", "t2"), c("raw", "bonf", "holm", "bh", "calib"),
                                  function(a, b) paste0("lmom_", a, "_", b, "_w")))
# 30 methods scored per replication: the ORIGINAL 20 plus 10 L-moment. The maxent
# family was retired with bqq 0.6.7 (author, 2026-09-05). All three families come from
# ONE detectChangepoints_gamma(basis = "all") call on the SAME single MAP+Laplace fit.
BQQ_METHODS <- c(paste0(names(BQQ_BASE_FIELDS), "_w"), LMOM_METHODS)

# ---------------------------------------------------------------------------
# SECTION 9: single replication -> per-method DETAIL rows (CHANGE 4)
# ---------------------------------------------------------------------------
# One CELL replication = one shared data draw scored by cpm (once) AND all 6 BQQ
# models (6 priors). cpm and every BQQ model see the SAME y and are
# emitted with the IDENTICAL row schema, so the methods are directly comparable.
run_single_simulation <- function(sim_id, config, ci = NA_integer_,
                                  max_retries = 3, log_file = NULL, seed = NULL) {
  log_msg <- function(msg) if (!is.null(log_file))
    cat(sprintf("[%s] cfg %s sim %d: %s\n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"), as.character(ci), sim_id, msg),
        file = log_file, append = TRUE)

  mr <- if (!is.null(config$max_retries)) config$max_retries else max_retries
  tv11_beat()   # task-start beat for the external watchdog

  # Shared task seed: supplied by the caller (which advances it across restarts so a
  # hard-hung draw is never re-hit) or the deterministic default. Recorded in meta$seed.
  task_seed <- if (!is.null(seed)) as.integer(seed) else as.integer(config$seed + sim_id)

  # ---- ONE shared data draw (retry only if the draw itself errors) ----
  attempt <- 0L; sim <- NULL
  while (is.null(sim) && attempt <= mr) {
    current_seed <- if (attempt == 0L) task_seed else
      as.integer((as.numeric(Sys.time()) * 1000 + Sys.getpid() * 131 + attempt) %% 2147483647) + 1L
    sim <- tryCatch(
      simulate_data(n = config$n_obs, distribution = config$distribution,
                    dependence = config$dependence, shift_type = config$shift_type,
                    shift_start = config$shift_start, shift_pattern = config$shift_pattern,
                    shift_length = config$shift_length, seed = current_seed,
                    ar_phi = config$ar_phi),
      error = function(e) { log_msg(paste("draw error:", conditionMessage(e))); NULL })
    if (is.null(sim)) attempt <- attempt + 1L
  }
  y <- sim$y; shift_point <- sim$shift_point; shift_end <- sim$shift_end

  # ---- cpm benchmarks: ONCE per draw (prior/IQ-independent) ----
  cpm_results <- run_cpm_variants(y, shift_point, config$tolerance)
  # ---- [ar_ext E5] reference R2: the same 8 baselines on auto.arima residuals.
  #      Also once per draw, and on the SAME y, so R1/R2/BQQ stay paired. ----
  # [ar_ext E10] ONE auto.arima fit, shared by the prewhitened cpm and the
  # prewhitened BQQ stream, so all four families are paired on the same draw.
  arima_fit <- fit_arima_once(y)
  arima_ref <- run_arima_cpm_variants(y, shift_point, config$tolerance, af = arima_fit)

  # ---- block geometry + truth (shared by cpm and every BQQ model) ----
  bounds    <- block_bounds(config$n_obs, config$warmup_period, config$block_length)
  is_signal <- compute_signal_blocks(bounds, shift_point, shift_end,
                                      config$shift_type, config$tolerance)
  r_blocks  <- nrow(bounds)
  n_sig     <- sum(is_signal)

  # ---- cell-level meta (identical on every row: cpm + all 6 BQQ models) ----
  meta <- data.frame(
    config_id = ci, sim_id = sim_id, seed = task_seed, setting = config$setting,
    distribution = config$distribution, dependence = config$dependence,
    # [ar_ext E4] with four phi levels the detail rows cannot otherwise tell P02 from
    # P08 except through config_id. Imported lmom_3cfg rows have no such column and the
    # importer adds it (0 for ids 1-5, 0.5 for ids 6-10) -- see ar_ext/REUSE.md.
    ar_phi = config$ar_phi,
    shift_type = config$shift_type, shift_pattern = config$shift_pattern,
    shift_length = config$shift_length, n_obs = config$n_obs,
    block_length = config$block_length, warmup = config$warmup_period,
    shift_start = config$shift_start, tolerance = config$tolerance,
    shift_size_tv = scenario_tv(config$distribution, config$shift_type),
    # [ar_ext E9] design TV (above) vs OBSERVED TV of the filtered marginal (below).
    shift_size_tv_marginal = scenario_tv_marginal(config$distribution, config$shift_type,
                                                  config$dependence, config$ar_phi),
    true_shift_point = shift_point, shift_end = shift_end,
    is_null = as.integer(config$shift_type == "null"),
    n_blocks = r_blocks, n_signal_blocks = n_sig,
    # ADDED 2026-08-04.  Stamps WHICH truth definition scored this row, so old and
    # new task files can never be pooled silently.  "onset_window" = signal blocks
    # are those overlapping [tau - tol, tau + tol]; rows written before this date
    # carry no such column and used the l-dependent [tau - tol, shift_end + tol].
    truth_rule = "onset_window", stringsAsFactors = FALSE
  )

  mk_row <- function(method, adjustment, prior, iq_mode, conf,
                     best_lambda_nc = NA_real_, best_spike_sd = NA_real_,
                     best_lambda_iq2 = NA_real_, iq_em_converged = NA,
                     iq_em_iters = NA_integer_,
                     # [ar_ext5 CHANGE 6] val_loss == val_score (the selection
                     # criterion); val_pinball is the pinball loss at that same
                     # score-selected configuration. NA on cpm/arima_cpm rows.
                     val_loss = NA_real_,
                     val_score = NA_real_, val_pinball = NA_real_,
                     n_errors = 0L,
                     # 2026-09-05: failed attempts before the recorded fit, and its seed
                     n_retries = 0L, fit_seed = NA_integer_,
                     cp_obs = NA_real_,
                     # [ar_ext E5] arima_cpm diagnostics; NA on cpm and BQQ rows.
                     arima_p = NA_integer_, arima_q = NA_integer_,
                     arima_ar1 = NA_real_, arima_mean = NA_integer_) {
    cbind(meta, data.frame(
      method = method, adjustment = adjustment,
      prior_gamma = prior, iq_mode = iq_mode,
      TP = conf$TP, FP = conf$FP, TN = conf$TN, FN = conf$FN,
      any_signal = conf$any_signal, n_flagged = conf$n_flagged,
      closest_signal_time = conf$closest_signal_time,
      detection_delay = conf$detection_delay, hit_within_tol = conf$hit_within_tol,
      # the flagged block set, and (cpm only) the raw change-point observation --
      # together these make every block-level metric re-scorable without refitting
      flagged_blocks = conf$flagged_blocks, cp_obs = cp_obs,
      best_lambda_nc = best_lambda_nc, best_spike_sd = best_spike_sd,
      best_lambda_iq2 = best_lambda_iq2,
      iq_em_converged = iq_em_converged, iq_em_iters = iq_em_iters,
      val_loss = val_loss, val_score = val_score, val_pinball = val_pinball,
      n_errors = n_errors, n_retries = n_retries, fit_seed = fit_seed,
      arima_p = arima_p, arima_q = arima_q,
      arima_ar1 = arima_ar1, arima_mean = arima_mean,
      stringsAsFactors = FALSE
    ))
  }

  rows <- list()

  # ---- cpm rows (mapped to the block scheme; same confusion + schema as BQQ) ----
  for (t in CPM_TYPES) {
    cr <- cpm_results[[t]]
    flags <- cpm_block_flags(cr$detected_points, bounds)
    conf  <- block_confusion(flags, is_signal, bounds, shift_point,
                             config$shift_type, config$tolerance)
    # finer cpm localization: use the exact changepoint when it lands within tolerance
    if (!is.na(conf$hit_within_tol) && conf$hit_within_tol == 1L &&
        length(cr$detected_points) == 1L) {
      conf$closest_signal_time <- cr$detected_points[1]
      conf$detection_delay     <- cr$detected_points[1] - shift_point
    }
    rows[[length(rows) + 1]] <- mk_row(
      "cpm", t, NA_character_, NA_character_, conf,
      n_errors = as.integer(!is.null(cr$error)),
      # raw change-point location: finer than the block index, so cpm can be
      # re-scored against ANY block length, not just the one it was run under
      cp_obs = if (length(cr$detected_points) == 1L &&
                   is.finite(cr$detected_points[1])) cr$detected_points[1]
               else NA_real_)
  }

  # ---- [ar_ext E5] reference R2 rows: identical schema and identical scoring to the
  #      R1 rows above, so R1 and R2 differ ONLY in the series the charts were run on.
  for (t in CPM_TYPES) {
    cr    <- arima_ref$results[[t]]
    flags <- cpm_block_flags(cr$detected_points, bounds)
    conf  <- block_confusion(flags, is_signal, bounds, shift_point,
                             config$shift_type, config$tolerance)
    if (!is.na(conf$hit_within_tol) && conf$hit_within_tol == 1L &&
        length(cr$detected_points) == 1L) {
      conf$closest_signal_time <- cr$detected_points[1]
      conf$detection_delay     <- cr$detected_points[1] - shift_point
    }
    rows[[length(rows) + 1]] <- mk_row(
      "arima_cpm", t, NA_character_, NA_character_, conf,
      n_errors = as.integer(!is.null(cr$error)),
      cp_obs = if (length(cr$detected_points) == 1L &&
                   is.finite(cr$detected_points[1])) cr$detected_points[1]
               else NA_real_,
      arima_p = arima_ref$p, arima_q = arima_ref$q,
      arima_ar1 = arima_ref$ar1, arima_mean = arima_ref$has_mean)
  }

  # ---- BQQ: 6 models = 6 priors. ar_ext12: lambda_iq2 is CV-tuned on IQ_LAMBDA_GRID^2,
  #      so the iq_mode label is "CV" (constant; kept for schema stability). ----
  iq_grids <- list(CV = NA)
  # [ar_ext E10, 2026-08-19] FOUR families now, not three: BQQ and cpm each run on the
  # raw series AND on the auto.arima residuals. The prewhitened BQQ stream is the
  # symmetric counterpart of reference R2 -- previously only the reference was
  # prewhitened, which a referee would fairly call an unequal comparison.
  #
  # COST: this DOUBLES the BQQ work, which is essentially all of the run's cost.
  #
  # CAVEAT inherited from Section 5B and now applying to BQQ too: because the shape
  # arms dilute in the MARGINAL but not in the innovations, prewhitening approximately
  # UNDOES that dilution. On the shape arms at high phi the raw and prewhitened streams
  # therefore read DIFFERENT signal strengths off the same replication -- that is a
  # comparison of streams, not of methods. It is a fair comparison on the null arm
  # (where there is no shift to dilute) and on the mean/variance arms.
  bqq_streams <- list(list(tag = "BQQ",       y = y),
                      list(tag = "arima_BQQ", y = arima_fit$e))
  for (.st in bqq_streams) {
   stream_tag <- .st$tag
   y_stream   <- .st$y
   stream_dead <- is.null(y_stream) || !all(is.finite(y_stream))
   if (stream_dead)
     log_msg(sprintf("stream %s unavailable (%s) -- emitting failed rows",
                     stream_tag, if (!is.null(arima_fit$err)) arima_fit$err else "NULL"))
  for (prior in config$priors) {
    for (iqm in names(iq_grids)) {
      grid_pts <- nrow(build_cv_grid(prior, config))
      # TIMEOUT_SCALE: the 1900 s anchor was calibrated on a REFERENCE machine, so
      # on a slower box every model legitimately takes proportionally longer while
      # the budget would not move -- the timeout then fires on healthy fits, each
      # costing 6 retries (~7 h on one prior) and stranding the whole replication.
      # Set TV11_TIMEOUT_SCALE to how many times slower this machine is.
      # Budget in FITS, not grid points: the CV is 2 folds per grid point and the final
      # fit is ONE fit (ar_ext12: no EM). The 1900 s anchor was calibrated at 45 grid
      # points + 1 final fit = 91 fits.
      n_fits   <- 2L * grid_pts + 1L
      est_sec  <- TIMEOUT_SCALE * 1900 * (n_fits / 91) *
                  ((config$n_obs * r_blocks) / (365 * 11))
      bqq_floor <- if (is.null(config$bqq_timeout_floor)) 900   else config$bqq_timeout_floor
      bqq_cap   <- if (is.null(config$bqq_timeout_cap))   21600 else config$bqq_timeout_cap
      bqq_cap   <- bqq_cap * TIMEOUT_SCALE   # else the cap re-binds on a slow box
      bqq_timeout <- min(bqq_cap, max(bqq_floor, ceiling(3 * est_sec)))

      # per-model timeout + RERUN: on any failure or TIMEOUT, retry up to max_retries.
      # 2026-09-05 (author ruling after the Codex full test): every attempt fits the SAME
      # series y_stream, so the model stays paired with cpm and the other 11 models of the
      # replication; only the model seed changes between attempts (a time-based seed, so a
      # relaunch never repeats a stuck initialization). Before this date a retry re-drew the
      # series, which silently broke the within-replication pairing. The retry count and
      # the seed of the successful attempt are recorded on every row (n_retries, fit_seed).
      # After max_retries -> mark failed (all-negative, n_errors=1).
      bqq <- NULL; m_attempt <- 0L
      repeat {
        if (stream_dead) break
        tv11_beat()   # one beat per model attempt: minutes-granular liveness signal
        vsay("stream %s | prior %s | attempt %d | per-model timeout %d s", stream_tag, prior, m_attempt + 1L, bqq_timeout)
        bqq_err <- NULL
        y_fit <- y_stream                                   # always the shared series
        model_seed <- if (m_attempt == 0L) current_seed else
          as.integer((as.numeric(Sys.time()) * 1000 + Sys.getpid() * 131 +
                      m_attempt * 17) %% 2147483647) + 1L
        if (is.null(y_fit) || !all(is.finite(y_fit))) { bqq <- NULL; bqq_err <- "no stream"
          m_attempt <- m_attempt + 1L; if (m_attempt > mr) break else next }
        bqq <- tryCatch({
          setTimeLimit(cpu = Inf, elapsed = bqq_timeout, transient = TRUE)
          res <- run_bqq_detection(y_fit, config, prior, model_seed)
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
          res
        }, error = function(e) {
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
          bqq_err <<- conditionMessage(e); NULL
        })
        if (!is.null(bqq) && is.null(bqq$error)) break
        m_attempt <- m_attempt + 1L
        log_msg(sprintf("BQQ prior=%s iq=%s attempt %d FAILED: %s", prior, iqm, m_attempt,
                        if (!is.null(bqq_err)) bqq_err else "NULL"))
        vsay("  FAILED attempt %d (%s | %s): %s", m_attempt, stream_tag, prior,
             if (!is.null(bqq_err)) bqq_err else if (!is.null(bqq$error)) bqq$error else "NULL")
        if (m_attempt > mr) break
      }
      failed <- stream_dead || is.null(bqq) || !is.null(bqq$error)
      if (failed) log_msg(sprintf("BQQ GAVE UP prior=%s iq=%s after %d attempts", prior, iqm, m_attempt))

      for (m in BQQ_METHODS) {
        flags <- if (failed) rep(FALSE, r_blocks) else bqq$flags[[m]]
        conf  <- block_confusion(flags, is_signal, bounds, shift_point,
                                 config$shift_type, config$tolerance)
        rows[[length(rows) + 1]] <- mk_row(
          stream_tag, m, prior, iqm, conf,
          arima_p    = if (stream_tag == "arima_BQQ") arima_fit$p    else NA_integer_,
          arima_q    = if (stream_tag == "arima_BQQ") arima_fit$q    else NA_integer_,
          arima_ar1  = if (stream_tag == "arima_BQQ") arima_fit$ar1  else NA_real_,
          arima_mean = if (stream_tag == "arima_BQQ") arima_fit$has_mean else NA_integer_,
          best_lambda_nc = if (failed) NA_real_ else bqq$best_lambda_nc,
          best_spike_sd  = if (failed) NA_real_ else bqq$best_spike_sd,
          best_lambda_iq2 = if (failed) NA_real_ else bqq$best_lambda_iq2,
          iq_em_converged = if (failed) NA      else bqq$iq_em_converged,
          iq_em_iters     = if (failed) NA_integer_ else bqq$iq_em_iters,
          val_loss       = if (failed) NA_real_ else bqq$val_loss,
          val_score      = if (failed) NA_real_ else bqq$val_score,
          val_pinball    = if (failed) NA_real_ else bqq$val_pinball,
          n_errors       = as.integer(failed),
          n_retries      = m_attempt,
          fit_seed       = if (failed) NA_integer_ else as.integer(model_seed))
      }
    }
  }
  }   # end stream loop [ar_ext E10]

  list(sim_id = sim_id, ci = ci, rows = do.call(rbind, rows))
}

# ---------------------------------------------------------------------------
# SECTION 10: config helpers + grid (CHANGE grid)
# ---------------------------------------------------------------------------
DEFAULT_PRIORS <- c("group_lasso", "lasso", "spike_slab",
                    "het_group_lasso", "adaptive_lasso", "spike_slab_lasso")

get_default_config <- function() {
  list(
    n_sim = 100,
    # [ar_ext E1] the DEFAULT is now D0 = normal + AR(1) phi = 0.5, not the
    # independent cell. Every other setting moves exactly one factor from this.
    setting = "D0", n_obs = 365, block_length = 30, warmup_period = 30,
    distribution = "normal", dependence = "ar",
    shift_type = "mean", shift_pattern = "sustained", shift_length = 30,
    shift_start = 0.8,
    priors = DEFAULT_PRIORS,          # 6 priors fit in every cell
    # CV grids (per BQQ model). lambda_nc (non-crossing) is de-emphasized to one value;
    # the lambda_iq grid (IQ_ON_GRID, which includes 0 = no fusion) is supplied by
    # run_single_simulation. Shrinkage: spike_sd (spike priors) / lambda_lasso2_b.
    lambda_nc_grid       = 50,
    lambda_lasso2_b_grid = c(0.01, 0.05, 0.1, 0.5, 1),
    # spike_sd grid NARROWED 2026-08-11.  The 2026-08-09 widening to
    # c(.01,.025,.05,.1,.15,.25,.4,.6) let CV reach values that BREAK SIZE, and the
    # 40-replication run showed the damage is entirely attributable to the small end:
    #
    #   spike_sd <= 0.05   71/155 null replications flagged   (0.458)
    #   spike_sd >= 0.10    1/85  null replications flagged   (0.012)
    #   Fisher exact p = 6.5e-16, odds ratio 70 (95% CI 12-2833)
    #
    # In S1 ALL FIVE of spike_slab's false alarms came from the 21 replications where
    # CV picked spike_sd <= 0.05; the other 19 gave zero.  The same mechanism explains
    # spike_slab_lasso's catastrophic 0.526 size -- it is not a different failure, it
    # just lands on a tight spike far more often (553/598 fits vs 190/598).
    # Mechanism: a tight spike makes the null component near a point mass, so any
    # noise-driven gamma is pushed into the slab and flagged.  A wider spike absorbs it.
    #
    # 0.6 is dropped too: selected 255 times for spike_slab and yields skewness power
    # 0.000, against 0.500 at spike_sd = 0.1.  That upper cut rests on only 4
    # replications, so it is the weaker half of this change.
    #
    # NOTE the objectives conflict: CV minimizes predictive loss, which a tight spike
    # improves, while size control wants the opposite.  This grid works by denying CV
    # the trade.  That is a deliberate restriction and belongs in the write-up.
    #
    # NOT a superset of the 2026-08-09 grid -- results are NOT poolable with the
    # 61-100 run, which is archived under archive_widegrid_2026-08-11/.  BASIS_SIG
    # does NOT fingerprint CV grids, so nothing in the code will catch a re-import;
    # keeping those checkpoints out of this folder is what enforces the separation.
    spike_sd_grid        = c(0.1, 0.15, 0.25, 0.4),
    tolerance = 7,
    # bqq_timeout_floor RAISED 900 -> 3600 (2026-08-02). This is a give-up deadline,
    # NOT a design parameter: it cannot change any fit that succeeds, it only decides
    # how long to wait before declaring failure. At 900 s the floor was the BINDING
    # budget for n = 180 (est_sec puts 3*est at only ~14 min there), and adaptive_lasso
    # -- measured at 262 s idle, ~9 min on the run workstation -- blew it routinely:
    # 87.5% of config 8 (S2 variance) and 47.9% of config 7 (S2 mean) adaptive_lasso
    # fits were recorded as failures, i.e. all-negative rows that bias size/power DOWN.
    # 3600 s gives ~6.7x headroom there. n = 365 cells are unaffected (3*est > 3600
    # already), so only S2/S3 budgets change.
    # bqq_timeout_cap RAISED 21600 -> 25000 (2026-08-14) for the S5 (l = 7) setting.
    # The cap is multiplied by TIMEOUT_SCALE, so at SCALE = 4 the old value gave a
    # 86,400 s ceiling while S5's 3 x est_sec wants 97,418 s (r = 47 blocks drives
    # est_sec up 4.3x over S1). The CAP would then bind instead of the 3x rule,
    # cutting S5's margin to 2.66x and re-creating the failure mode that fabricated
    # 87.5% of config 8's adaptive_lasso failures in the old study, where the FLOOR
    # bound at n = 180. 25000 * 4 = 100,000 s > 97,418, so the 3x rule governs every
    # setting again. This is a give-up deadline; it cannot change a fit that succeeds.
    bqq_timeout_floor = 3600, bqq_timeout_cap = 25000, max_retries = 5,
    n_cores = NULL, ar_phi = 0.5, seed = 2026
  )
}

create_config <- function(...) {
  config <- get_default_config()
  args <- list(...)
  for (nm in names(args)) {
    if (nm %in% names(config)) config[[nm]] <- args[[nm]]
    else warning("Unknown configuration parameter: ", nm)
  }
  config
}

# Config grid = 40 cells (EXPANDED from 15, 2026-08-14). A config is a CELL =
# (setting x shift type); every cell is scored by cpm (once) AND all 6 BQQ models
# (6 priors) on one draw. EIGHT settings, each changing at most ONE factor from the
# default (warmup w = 30 fixed) -- the OFAT grid of update_tv_shift_11 minus the
# retired platykurtic setting:
#   S1 default : T365 l30 sustained independent normal
#   S7 AR(0.5)              (dependence)
#   S9 right_skewed         (underlying distribution)
#   S2 T180  | S3 T90       (time-window length)
#   S4 l14   | S5 l7        (block length)
#   S6 isolated             (shift pattern)
# Every setting runs all 5 shift arms {null, mean, variance, skewness, kurtosis},
# all at TV = 0.383 under the 2026-08-13 design point.  -> 8 x 5 = 40 configs.
#
# *** ORDER IS APPEND-ONLY AND MUST NOT BE RESORTED. ***
# S1/S7/S9 keep config ids 1-15 exactly as the 597 checkpoints written between
# 2026-08-13 and 2026-08-14 assume. Re-sorting into OFAT order (S1,S2,...,S9)
# would silently re-point config 6 from S7-null to S2-null, and NOTHING would
# catch it: BASIS_SIG fingerprints the lmom weights, the bqq version and (since 2026-09-05)
# a hash of the study design (SCHEMA_HASH_NAMES),
# not the config grid -- the same blind spot it has for the CV grid and the DGP.
# New settings are therefore APPENDED at ids 16-40. CONFIG_MAP below is the guard.
make_configs <- function(n_sim = 100, priors = DEFAULT_PRIORS) {
  # ar_ext12: 11 settings x 3 arms = 33 configs (HK appended last; skewness and
  # kurtosis arms dropped). [ar_ext E2] OFAT: each setting moves at most
  # ONE factor from the default D0. P00 and D0 deliberately occupy ids 1-10, the same
  # ids they hold in lmom_3cfg as S1 and S7, so those 1000 banked replications import
  # with no re-keying -- see ar_ext/REUSE.md.
  D <- list(setting = "D0", n_obs = 365, block_length = 30,
            shift_pattern = "sustained", distribution = "normal",
            dependence = "ar", ar_phi = 0.5)
  settings <- list(
    # phi = 0 MUST be "independent", never ar with ar_phi = 0: the ar branch draws a
    # stationary y[1] and consumes one extra RNG draw, which desynchronises the whole
    # stream and breaks bit-exact reproduction of the banked cells (GRID.md section 0).
    modifyList(D, list(setting = "P00", dependence = "independent", ar_phi = 0)),  #  1-5  <- lmom_3cfg S1
    D,                                                                             #  6-10 <- lmom_3cfg S7  DEFAULT
    modifyList(D, list(setting = "P02", ar_phi = 0.2)),                            # 11-15
    modifyList(D, list(setting = "P08", ar_phi = 0.8)),                            # 16-20
    modifyList(D, list(setting = "K",   distribution = "right_skewed")),           # 21-25
    modifyList(D, list(setting = "T18", n_obs = 180)),                             # 26-30
    modifyList(D, list(setting = "T09", n_obs = 90)),                              # 31-35
    modifyList(D, list(setting = "L14", block_length = 14)),                       # 36-40
    modifyList(D, list(setting = "L07", block_length = 7)),                        # 41-45
    modifyList(D, list(setting = "ISO", shift_pattern = "isolated")),              # 28-30
    modifyList(D, list(setting = "HK",  distribution = "high_kurtosis"))           # 31-33 (ar_ext11)
  )
  # ar_ext12: the manuscript arms only (Section 4.1): no shift, location, scale.
  shifts <- c("null", "mean", "variance")
  configs <- list()
  for (s in settings) {
    for (st in shifts) {
      configs[[length(configs) + 1]] <- create_config(
        n_sim = n_sim, setting = s$setting, n_obs = s$n_obs,
        block_length = s$block_length, warmup_period = 30,
        distribution = s$distribution, dependence = s$dependence,
        shift_type = st, shift_pattern = s$shift_pattern,
        ar_phi = s$ar_phi,           # [ar_ext E2] WITHOUT THIS every cell silently
                                     # inherits the default 0.5 and the phi axis is a
                                     # no-op. Config 16 must come back with 0.8.
        shift_length = s$block_length, priors = priors)
    }
  }
  .assert_config_map(configs)
  configs
}

# Guard against the re-sorting hazard described above. Any change to the ORDER of
# `settings` (as opposed to appending to it) trips this immediately, instead of
# silently pooling two different settings under one config id months later.
# [ar_ext E3] append-only, never re-sort: these ids key the checkpoints, and ids 1-10
# additionally key the imported lmom_3cfg files.
CONFIG_MAP <- c(P00 = 1L, D0 = 4L, P02 = 7L, P08 = 10L, K = 13L,
                T18 = 16L, T09 = 19L, L14 = 22L, L07 = 25L, ISO = 28L, HK = 31L)
.assert_config_map <- function(configs) {
  if (length(configs) != 33L)
    stop("config grid is ", length(configs), " cells, expected 33", call. = FALSE)
  for (nm in names(CONFIG_MAP)) {
    i <- CONFIG_MAP[[nm]]
    if (!identical(configs[[i]]$setting, nm))
      stop("CONFIG ORDER DRIFT: id ", i, " is setting '", configs[[i]]$setting,
           "', expected '", nm, "'. Existing checkpoints are keyed on config id; ",
           "re-sorting invalidates them silently. Append new settings, never reorder.",
           call. = FALSE)
    if (!identical(configs[[i]]$shift_type, "null"))
      stop("CONFIG ORDER DRIFT: id ", i, " should be the 'null' arm of ", nm,
           ", got '", configs[[i]]$shift_type, "'", call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# SECTION 11: recovery helper (REUSED verbatim from v3)
# ---------------------------------------------------------------------------
recover_from_task_files <- function(filename, configs_to_run, configs) {
  all_results <- list(); missing_tasks <- list()
  for (ci in configs_to_run) {
    config <- configs[[ci]]; ci_results <- list()
    for (sim_id in 1:config$n_sim) {
      task_file <- sprintf("%s_task_config_%d_sim_%d.rds", filename, ci, sim_id)
      if (file.exists(task_file)) {
        ci_results[[length(ci_results) + 1]] <- readRDS(task_file)
      } else {
        missing_tasks[[length(missing_tasks) + 1]] <- list(ci = ci, sim_id = sim_id, config = config)
      }
    }
    all_results[[as.character(ci)]] <- ci_results
  }
  list(results = all_results, missing = missing_tasks)
}

# ---------------------------------------------------------------------------
# SECTION 12: parallel batch driver (structure REUSED from v3; row output)
# ---------------------------------------------------------------------------
run_batch_simulation <- function(configs, filename = "tv11",
                                  save_individual = TRUE, resume = TRUE,
                                  max_cores = 30, init_batch_size = 10,
                                  init_delay_sec = 60, stagger_init = TRUE,
                                  run_verify = TRUE, box_dir = NULL) {
  n_configs <- length(configs)
  cat("================================================================\n")
  cat("TV11 BATCH SIMULATION STUDY  (", n_configs, "configs )\n")
  cat("================================================================\n")
  if (run_verify) verify_dgp()

  csv_file <- paste0(filename, "_combined.csv")
  log_file <- paste0(filename, "_progress.log")

  # --- BALANCED SCHEDULING (env TV11_TARGET; SCHEDULING ONLY) ---------------
  # Levels every config up to TV11_TARGET completed replications in round-robin
  # waves, instead of shuffling all remaining tasks uniformly. NOTHING about the
  # study design changes: the config grid, the CV grid, the DGP and every
  # hyperparameter are untouched, and task (ci, sim_id) keeps exactly the seed it
  # would have had (config$seed + 200000 + sim_id). Only the ORDER in which the
  # already-defined tasks are dispatched changes, so results are bit-identical to
  # what the uniform schedule would eventually have produced for the same tasks.
  # Unset (default) -> the original uniform-shuffle behavior, unchanged.
  tv11_target <- suppressWarnings(as.integer(Sys.getenv("TV11_TARGET", "")))
  balanced    <- !is.na(tv11_target) && tv11_target > 0L

  if (!resume) {
    cat("", file = log_file, append = FALSE)
    if (file.exists(csv_file)) file.remove(csv_file)
    for (g in c(Sys.glob(paste0(filename, "_config_*.rds")),
                Sys.glob(paste0(filename, "_task_config_*.rds")),
                Sys.glob(paste0(filename, "_task_config_*.try")),
                Sys.glob(paste0(filename, "_progress.log.worker_*.log")),
                file.path(dirname(filename), "arext12_run.pids"),
                paste0(filename, ".rds"))) if (file.exists(g)) file.remove(g)
  }

  # --- prime the LOCAL work dir from the durable (Box) dir ONCE (resume only) ---
  # Import any checkpoints already produced on Box so a first local run continues the
  # work done there instead of starting over. Cheap after the first pass (skips existing).
  if (resume && !is.null(box_dir) &&
      normalizePath(box_dir, mustWork = FALSE) !=
      normalizePath(dirname(filename), mustWork = FALSE)) {
    # Timestamped so the console log shows how long the Box scan took -- this step
    # reads ~1500 names off the synced mount and is the one place a stalled Box Drive
    # can hang the child before the cluster exists (cf. the 2026-07-13 mount freeze).
    cat("[", format(Sys.time(), "%H:%M:%S"), "] scanning box_dir for importable checkpoints...\n")
    base <- basename(filename); wd <- dirname(filename)
    # Also scan box_dir/tv11_run: stale local work dirs get parked there by hand
    # (e.g. the 2026-07-13 run), and those checkpoints must survive a fresh machine.
    src_dirs <- c(box_dir, file.path(box_dir, "tv11_run"))
    imp <- unlist(lapply(src_dirs, function(d)
      c(Sys.glob(file.path(d, paste0(base, "_task_config_*.rds"))),
        Sys.glob(file.path(d, paste0(base, "_config_*.rds"))))))
    # Import is fingerprint-gated. Without this, deleting stale checkpoints from
    # work_dir does nothing: the next launch copies them straight back from box_dir,
    # and the box sync copies anything in work_dir back to box_dir, so the two
    # directories keep re-seeding each other and old results appear to respawn.
    # Skipping them here makes a delete in work_dir actually stick.
    n_imp <- 0L; n_skip <- 0L
    for (f in imp) { dst <- file.path(wd, basename(f))
      if (file.exists(dst)) next
      if (grepl("_task_config_", basename(f), fixed = TRUE)) {
        o <- tryCatch(readRDS(f), error = function(e) NULL)
        if (is.null(o) || is.null(o$basis_sig) ||
            !identical(as.character(o$basis_sig), as.character(BASIS_SIG))) {
          n_skip <- n_skip + 1L; next
        }
      }
      if (isTRUE(file.copy(f, dst, overwrite = FALSE))) n_imp <- n_imp + 1L }
    cat("[", format(Sys.time(), "%H:%M:%S"), "] box scan done:", length(imp), "candidate(s),",
        n_imp, "newly imported ->", wd, "\n")
    if (n_skip > 0L)
      cat("           ", n_skip, "skipped: different contrast basis (current stamp ",
          BASIS_SIG, ")\n", sep = "")
  }

  # rows from a per-config rds
  rows_from_rds <- function(sim_result) sim_result$rows

  # --- which configs still need running (resume) ---
  if (resume) {
    completed <- integer(0)
    for (ci in seq_len(n_configs))
      if (file.exists(paste0(filename, "_config_", ci, ".rds"))) completed <- c(completed, ci)
    configs_to_run <- setdiff(seq_len(n_configs), completed)
    cat(length(completed), "of", n_configs, "configs already complete.\n")
  } else configs_to_run <- seq_len(n_configs)
  # TV11_ONLY_CONFIG = "31" (or "31,36"): run only these config ids. Used by
  # smoke_driver_arext12.R to run the real driver on one config; ids stay the positions
  # in the full grid, so checkpoint names and seeds are the study's own.
  .only <- Sys.getenv("TV11_ONLY_CONFIG", "")
  if (nzchar(.only)) {
    configs_to_run <- intersect(configs_to_run, as.integer(strsplit(.only, ",")[[1]]))
    cat("TV11_ONLY_CONFIG: restricting to config id(s)", paste(configs_to_run, collapse = ", "), "\n")
  }

  # --- count remaining tasks BEFORE the cluster ---
  n_actual_tasks <- 0L
  for (ci in configs_to_run) for (i in 1:configs[[ci]]$n_sim)
    if (!file.exists(sprintf("%s_task_config_%d_sim_%d.rds", filename, ci, i)))
      n_actual_tasks <- n_actual_tasks + 1L
  n_cores <- min(max_cores, max(1L, n_actual_tasks))
  cat("Using n_cores =", n_cores, "(", n_actual_tasks, "remaining tasks )\n")

  batch_cl <- NULL; worker_pids <- integer(0)
  # TV11_FORCE_CLUSTER=1 makes a 1-worker PSOCK cluster instead of the serial path, so a
  # 1-thread smoke goes through clusterExport() and catches a missing export.
  # ar_ext12: PSOCK workers on every OS (ar_ext10 used them on Windows only, so macOS
  # and Linux fell back to the sequential loop below, one core for every task).
  if (n_cores > 1 || nzchar(Sys.getenv("TV11_FORCE_CLUSTER", ""))) {
    cat("[", format(Sys.time(), "%H:%M:%S"), "] creating PSOCK cluster with", n_cores, "workers...\n")
    flush.console()
    batch_cl <- parallel::makeCluster(n_cores, timeout = 86400)
    cat("[", format(Sys.time(), "%H:%M:%S"), "] cluster up.\n")
    # Capture worker PIDs and write the manifest IMMEDIATELY -- BEFORE the slow staggered
    # init -- so the external watchdog can see the run is alive throughout the ~15-min
    # initialization and does not mistake a still-initializing run for a crashed one.
    worker_pids <- unlist(parallel::clusterCall(batch_cl, Sys.getpid))
    try(writeLines(as.character(c(Sys.getpid(), worker_pids)),
                   file.path(dirname(filename), "arext12_run.pids")), silent = TRUE)
    # Workers beat from the moment they are alive, not from their first task. tv11_beat()
    # is not exported yet at this point (that happens after init), so the beat is written
    # inline against the exported HB_DIR. Without this the whole ~15-25 min ramp is
    # heartbeat-silent and indistinguishable, from outside, from a wedged startup.
    init_fun <- function(i) {
      .beat <- function() try(cat(format(Sys.time(), "%H:%M:%S"),
                                  file = file.path(.hb, paste0("tv11_hb_", Sys.getpid()))),
                              silent = TRUE)
      .beat()
      batch_id <- ceiling(i / .bs)
      if (batch_id > 1) Sys.sleep((batch_id - 1) * .ds)
      .beat()
      require(cpm); require(bbq); require(PearsonDS); require(forecast)
      .beat()
      options(mc.cores = 1)
      tryCatch(bbq::getModel(y = rnorm(20), taus = c(0.25, 0.5, 0.75),
                             H = matrix(0, 20, 1), w = 10, fit_method = "map",
                             map_hessian = FALSE, map_iter = 1, verbose = FALSE),
               error = function(e) NULL)
      gc(); list(worker = i, ready = TRUE)
    }
    if (stagger_init) {
      .bs <- init_batch_size; .ds <- init_delay_sec; .hb <- HB_DIR
      clusterExport(batch_cl, c(".bs", ".ds", ".hb"), envir = environment())
      ir <- parallel::parLapplyLB(batch_cl, seq_len(n_cores), init_fun)
      cat("Cluster ready:", sum(sapply(ir, function(r) isTRUE(r$ready))), "workers.\n")
    } else {
      parallel::clusterEvalQ(batch_cl, { require(cpm); require(bbq); require(PearsonDS); require(forecast); options(mc.cores = 1); TRUE })
    }
  }
  # PID manifest path (written above, right after the cluster came up); removed on exit.
  pidfile <- file.path(dirname(filename), "arext12_run.pids")
  on.exit({
    if (!is.null(batch_cl)) tryCatch(parallel::stopCluster(batch_cl), error = function(e) NULL)
    for (pid in worker_pids) tryCatch(tools::pskill(pid), error = function(e) NULL)
    try(unlink(pidfile), silent = TRUE)
  }, add = TRUE)

  # export functions + constants to workers (everything run_single_simulation needs)
  export_names <- c(
    "simulate_data", "draw_baseline", "draw_post",
    "BASE_PARAMS", "POST_PARAMS", "POST_MEAN", "POST_SD", "HK_POST_MEAN", "HK_POST_SD", "PEARSON_SKEW",
    # PEARSON_KURT was REMOVED 2026-08-13 (the kurtosis arm is a scale mixture now).
    # Leaving it here would make clusterExport() fail on every worker.
    "KURT_MIX_NORMAL", "KURT_MIX_SKEWED",
    ".scalemix_s2", ".draw_scalemix", ".scalemix_density",
    "TAUS", "TV_ALL_ARMS", "TV_BY_ARM", "TV_DEVIATED", "TV_DEVIATED_KURT", "scenario_tv",
    "run_cpm_batch", "run_cpm_variants", "CPM_TYPES", "cpm_block_flags",
    # [ar_ext E5 / E9]
    "run_arima_cpm_variants", "fit_arima_once", "ARIMA_ARGS", "TV_MARGINAL", "scenario_tv_marginal",
    "block_bounds", "compute_signal_blocks", "block_confusion",
    "build_cv_grid", "IQ_LAMBDA_GRID",
    # 2026-09-05: every constant the worker-side functions read must be listed here. The
    # 2026-09-04 study died on all 5000 tasks with "object 'CV_IQ_EM_MAX_ITER' not found"
    # because constants were added to the driver but not to this list; the old smoke ran
    # in one process and could not see the gap. smoke_arext12.R now has a static export
    # check (W2) and smoke_driver_arext12.R runs the driver itself through a 1-worker cluster.
    "CV_LAPLACE_DRAWS",
    "TV11_VERBOSE", "vsay",
    "MAP_OPTS", "TIMEOUT_SCALE",
    "run_bqq_detection", "BQQ_BASE_FIELDS", "BQQ_METHODS",
    "shape_flags", "LMOM_METHODS",
    "run_single_simulation", "save_rds_atomic", "tv11_beat", "HB_DIR", "BASIS_SIG")

  if (!is.null(batch_cl)) {
    clusterExport(batch_cl, export_names, envir = .GlobalEnv)

    .lf <- log_file; .fn <- filename
    run_sim_worker <- function(task) {
      # Per-worker log (one file per PID) instead of every worker appending to ONE shared
      # log on the mount -- concurrent appends to a single file were a serialization /
      # lock point. The .rds checkpoints (not the log) are the source of truth for resume.
      wlog <- sprintf("%s.worker_%d.log", .lf, Sys.getpid())
      assign(".tv11_wlog", wlog, envir = .GlobalEnv)   # vsay() progress lines land here too
      # Restart-safe seed. A hard hang (native/rstan wedge that setTimeLimit cannot
      # interrupt) is killed by the watchdog and the task is redone; without this it would
      # re-hit the SAME deterministic seed and kill-loop forever. Persist a per-task attempt
      # counter BEFORE the (possibly hanging) fit, so each relaunch advances to a NEW seed;
      # the seed used is recorded in the output for reproducibility.
      try_file <- sprintf("%s_task_config_%d_sim_%d.try", .fn, task$ci, task$sim_id)
      k <- if (file.exists(try_file)) suppressWarnings(as.integer(readLines(try_file, n = 1L)[1])) else 0L
      if (length(k) == 0L || is.na(k)) k <- 0L
      writeLines(as.character(k + 1L), try_file)                       # advance-and-persist first
      task_seed <- if (k == 0L) as.integer(task$config$seed + task$sim_id) else
        as.integer(((task$config$seed + task$sim_id) * 1000003 + k * 97 + 1) %% 2147483647) + 1L
      result <- run_single_simulation(task$sim_id, task$config, ci = task$ci,
                                      log_file = wlog, seed = task_seed)
      save_rds_atomic(list(ci = task$ci, sim_id = task$sim_id, seed = task_seed,
                           basis_sig = BASIS_SIG, rows = result$rows),
                      sprintf("%s_task_config_%d_sim_%d.rds", .fn, task$ci, task$sim_id))
      if (file.exists(try_file)) try(file.remove(try_file), silent = TRUE)  # completed -> drop marker
      cat(sprintf("[%s] config %d sim %d/%d done (seed %d)\n", format(Sys.time(), "%H:%M:%S"),
                  task$ci, task$sim_id, task$config$n_sim, task_seed), file = wlog, append = TRUE)
      gc()   # rstan workers grow all day without this (observed 53 -> 11 tasks/hr decay)
      list(ci = task$ci, sim_id = task$sim_id)
    }

    if (length(configs_to_run) > 0) {
      tasks <- list()
      for (ci in configs_to_run) {
        config <- configs[[ci]]
        for (i in 1:config$n_sim) {
          if (file.exists(sprintf("%s_task_config_%d_sim_%d.rds", filename, ci, i))) next
          rc <- config; if (resume) rc$seed <- config$seed + 200000L
          tasks[[length(tasks) + 1]] <- list(ci = ci, sim_id = i, config = rc)
        }
      }
      if (balanced) {
        # Keep only the tasks still needed to bring each config up to the target,
        # taking the lowest missing sim_ids (deterministic, seed-preserving), then
        # interleave configs round-robin so wave k advances EVERY deficit config by
        # one replication. Stop the run at any point and the panel stays balanced.
        ci_of  <- vapply(tasks, function(t) t$ci, integer(1))
        # Count USABLE replications: one where every BQQ fit converged. A replication
        # holding a timed-out fit (n_errors == 1) carries all-FALSE detection flags
        # that look like genuine non-detections downstream, so counting it would
        # leave the cell quietly under-filled. Topping up is purely ADDITIVE -- the
        # existing file stays exactly where it is; a NEW sim_id is scheduled instead.
        # Unusable ONLY if a failed fit is positively visible; a checkpoint format
        # without n_errors counts as usable (else it would reschedule forever).
        .clean <- function(f) {
          o <- tryCatch(readRDS(f), error = function(e) NULL)
          if (is.null(o) || is.null(o$rows)) return(FALSE)
          # A checkpoint computed under a DIFFERENT contrast basis must NOT count.
          # Pre-0.5.1 files carry no stamp, have the same 248 rows and the same
          # method names, and are therefore indistinguishable by schema alone --
          # counting them is what made a fresh study report "target met" and
          # schedule 0 tasks.
          if (is.null(o$basis_sig) ||
              !identical(as.character(o$basis_sig), as.character(BASIS_SIG))) return(FALSE)
          if (!"n_errors" %in% names(o$rows)) return(TRUE)
          !any(o$rows$n_errors == 1L, na.rm = TRUE)
        }
        done_n <- vapply(seq_len(n_configs), function(ci) {
          fs <- sprintf("%s_task_config_%d_sim_%d.rds",
                        filename, ci, seq_len(configs[[ci]]$n_sim))
          fs <- fs[file.exists(fs)]
          if (!length(fs)) return(0L)
          sum(vapply(fs, .clean, logical(1)))
        }, integer(1))
        # TV11_PICK: which end of each config's missing sim_ids to take.
        # "low" (default) takes the lowest ids; "high" takes the highest. Two
        # machines sharing one checkpoint pool set opposite ends, so they work
        # toward each other and only overlap once the deficit is nearly closed --
        # otherwise both would pick the SAME ids and duplicate every replication.
        pick_high <- identical(tolower(Sys.getenv("TV11_PICK", "low")), "high")
        keep <- logical(length(tasks))
        for (ci in unique(ci_of)) {
          idx  <- which(ci_of == ci)                  # already ascending in sim_id
          need <- max(0L, tv11_target - done_n[ci])
          if (need > 0L)
            keep[if (pick_high) utils::tail(idx, need) else utils::head(idx, need)] <- TRUE
        }
        cat(sprintf("  task pick end: %s (TV11_PICK)\n", if (pick_high) "HIGH sim_ids" else "low sim_ids"))
        tasks <- tasks[keep]; ci_of <- ci_of[keep]
        if (length(tasks) > 0L) {
          wave  <- stats::ave(seq_along(ci_of), ci_of, FUN = seq_along)
          tasks <- tasks[order(wave, ci_of)]
        }
        {
          all_f <- Sys.glob(paste0(filename, "_task_config_*_sim_*.rds"))
          if (length(all_f)) {
            stale <- sum(!vapply(all_f, .clean, logical(1)))
            if (stale > 0L)
              cat(sprintf("  %d of %d checkpoint(s) on disk do NOT count (failed fit, or a\n",
                          stale, length(all_f)),
                  sprintf("    different contrast basis -- current stamp %s).\n", BASIS_SIG), sep = "")
          }
        }
        cat(sprintf("BALANCED MODE: target %d/config -> %d tasks over %d deficit configs.\n",
                    tv11_target, length(tasks), length(unique(ci_of))))
        cat("  per-config deficit:",
            paste(sprintf("c%d:%d", sort(unique(ci_of)),
                          pmax(0L, tv11_target - done_n[sort(unique(ci_of))])),
                  collapse = " "), "\n")
      } else {
        set.seed(2026); tasks <- tasks[sample(length(tasks))]
      }
      cat("Distributing", length(tasks), "tasks across", n_cores, "workers...\n")
      # chunk.size = 1 -> one task per dispatch. The default (NULL) pre-splits the
      # list into 2*nnodes CONTIGUOUS chunks, which both defeats load balancing and
      # scrambles the round-robin wave order above.
      if (length(tasks) > 0L)
        tryCatch(parallel::parLapplyLB(batch_cl, tasks, run_sim_worker, chunk.size = 1),
                 error = function(e) cat("parLapplyLB error:", conditionMessage(e), "\n"))
      tryCatch(parallel::stopCluster(batch_cl), error = function(e) NULL)
      batch_cl <- NULL
      for (pid in worker_pids) tryCatch(tools::pskill(pid), error = function(e) NULL)
      gc()
    }
  } else {
    # sequential fallback
    for (ci in configs_to_run) {
      config <- configs[[ci]]
      for (i in 1:config$n_sim) {
        if (file.exists(sprintf("%s_task_config_%d_sim_%d.rds", filename, ci, i))) next
        try_file <- sprintf("%s_task_config_%d_sim_%d.try", filename, ci, i)
        k <- if (file.exists(try_file)) suppressWarnings(as.integer(readLines(try_file, n = 1L)[1])) else 0L
        if (length(k) == 0L || is.na(k)) k <- 0L
        writeLines(as.character(k + 1L), try_file)
        task_seed <- if (k == 0L) as.integer(config$seed + i) else
          as.integer(((config$seed + i) * 1000003 + k * 97 + 1) %% 2147483647) + 1L
        res <- run_single_simulation(i, config, ci = ci, log_file = log_file, seed = task_seed)
        save_rds_atomic(list(ci = ci, sim_id = i, seed = task_seed,
                             basis_sig = BASIS_SIG, rows = res$rows),
                        sprintf("%s_task_config_%d_sim_%d.rds", filename, ci, i))
        if (file.exists(try_file)) try(file.remove(try_file), silent = TRUE)
      }
    }
  }

  # --- aggregate per-task rows -> per-config rds ---
  # DESTRUCTIVE: writes tv11_config_<ci>.rds and then DELETES that config's task
  # files. Under balanced scheduling a config is only PARTIALLY filled, so doing
  # this would (a) throw away the resumable per-task checkpoints and (b) make the
  # next resume treat the partial config as COMPLETE and never top it up again.
  # In balanced mode we therefore leave the task files alone and let
  # build_app_csv_from_tasks.R do the aggregation (it reads the task files directly).
  if (balanced) {
    cat("BALANCED MODE: skipping per-config aggregation + task-file cleanup.\n")
    cat("  Aggregate with build_app_csv_from_tasks.R instead.\n")
    .clean2 <- function(f) {
      o <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(o) || is.null(o$rows)) return(FALSE)
      if (is.null(o$basis_sig) ||
          !identical(as.character(o$basis_sig), as.character(BASIS_SIG))) return(FALSE)
      if (!"n_errors" %in% names(o$rows)) return(TRUE)
      !any(o$rows$n_errors == 1L, na.rm = TRUE)
    }
    have <- vapply(seq_len(n_configs), function(ci) {
      fs <- sprintf("%s_task_config_%d_sim_%d.rds", filename, ci, seq_len(configs[[ci]]$n_sim))
      fs <- fs[file.exists(fs)]
      if (!length(fs)) return(0L)
      sum(vapply(fs, .clean2, logical(1)))
    }, integer(1))
    short <- which(have < tv11_target)
    cat(sprintf("  STATUS: %d/%d configs at >= %d reps (min %d, median %g, max %d).\n",
                sum(have >= tv11_target), n_configs, tv11_target,
                min(have), stats::median(have), max(have)))
    if (length(short) == 0L) {
      cat("  *** TARGET MET for every config -- stop the supervisor now. ***\n")
    } else {
      cat("  still short:",
          paste(sprintf("c%d(%d)", short, tv11_target - have[short]), collapse = " "), "\n")
    }
    if (!is.null(box_dir) &&
        normalizePath(box_dir, mustWork = FALSE) !=
        normalizePath(dirname(filename), mustWork = FALSE)) {
      n_ok <- 0L
      for (f in Sys.glob(paste0(filename, "_task_config_*.rds")))
        if (isTRUE(try(file.copy(f, file.path(box_dir, basename(f)), overwrite = TRUE),
                       silent = TRUE))) n_ok <- n_ok + 1L
      cat("  Synced", n_ok, "task checkpoints to Box.\n")
    }
    return(invisible(list(configs = configs, combined = NULL, csv = NA_character_)))
  }
  final <- recover_from_task_files(filename, configs_to_run, configs)
  for (ci in configs_to_run) {
    ci_data <- final$results[[as.character(ci)]]
    if (length(ci_data) == 0) next
    ci_data <- ci_data[order(sapply(ci_data, function(r) r$sim_id))]
    rows <- do.call(rbind, lapply(ci_data, function(r) r$rows))
    if (save_individual)
      save_rds_atomic(list(config = configs[[ci]], rows = rows),
                      paste0(filename, "_config_", ci, ".rds"))
    for (sim_id in 1:configs[[ci]]$n_sim) {
      tf <- sprintf("%s_task_config_%d_sim_%d.rds", filename, ci, sim_id)
      if (file.exists(tf)) file.remove(tf)
      trf <- sprintf("%s_task_config_%d_sim_%d.try", filename, ci, sim_id)
      if (file.exists(trf)) file.remove(trf)
    }
  }

  # --- rebuild combined CSV from all per-config rds ---
  all_rows <- list()
  for (ci in seq_len(n_configs)) {
    f <- paste0(filename, "_config_", ci, ".rds")
    if (file.exists(f)) all_rows[[length(all_rows) + 1]] <- rows_from_rds(readRDS(f))
  }
  combined <- if (length(all_rows) > 0) do.call(rbind, all_rows) else NULL
  if (!is.null(combined)) {
    write.csv(combined, csv_file, row.names = FALSE)
    cat("Combined CSV written:", csv_file, "(", nrow(combined), "rows )\n")
  }
  save_rds_atomic(list(configs = configs, combined = combined), paste0(filename, ".rds"))

  # --- sync durable outputs back to Box (hot per-task I/O stayed local) ---
  if (!is.null(box_dir) &&
      normalizePath(box_dir, mustWork = FALSE) !=
      normalizePath(dirname(filename), mustWork = FALSE)) {
    outs <- c(Sys.glob(paste0(filename, "_config_*.rds")), csv_file, paste0(filename, ".rds"))
    for (f in outs) if (file.exists(f))
      try(file.copy(f, file.path(box_dir, basename(f)), overwrite = TRUE), silent = TRUE)
    cat("Synced aggregated results to Box:", box_dir, "\n")
  }
  invisible(list(configs = configs, combined = combined, csv = csv_file))
}

# ---------------------------------------------------------------------------
# SECTION 13: quick run (MANDATORY end-to-end check) + gated driver
# ---------------------------------------------------------------------------
run_quick <- function() {
  quick_configs <- make_configs(n_sim = 2, priors = c("spike_slab"))
  res <- run_batch_simulation(quick_configs,
                              filename = file.path(output_dir, "arext12"),
                              max_cores = 3, init_batch_size = 3,
                              init_delay_sec = 5, resume = FALSE, run_verify = TRUE)
  if (!is.null(res$combined)) {
    cat("\n---- HEAD of arext12_combined.csv (confusion + localization columns) ----\n")
    show <- res$combined[, c("config_id", "sim_id", "setting", "shift_type", "method",
                             "adjustment", "prior_gamma", "iq_mode", "TP", "FP", "TN", "FN",
                             "any_signal", "hit_within_tol", "detection_delay", "n_errors")]
    print(utils::head(show, 24))
  }
  invisible(res)
}

#mode <- Sys.getenv("TV11_MODE", "")
#if (identical(mode, "quick")) {
#  run_quick()
#} else if (identical(mode, "full")) {
  # TV11_NSIM RAISED 100 -> 120 (2026-08-14). This is the POOL of available sim_ids
  # per config, NOT the amount of work: the ladder bounds work through TV11_TARGET
  # (balanced scheduling, see run_batch_simulation), so a larger pool costs nothing.
  # It must EXCEED the ladder target. `have` counts only USABLE replications while
  # sim_ids run 1:n_sim, so at n_sim = 100 with a target of 100 there is ZERO
  # headroom: one timed-out fit anywhere leaves have = 99, need = 1, and no free
  # sim_id to fill it -- the target becomes permanently unreachable and the ladder
  # stops on a no-op child. The risk is real for the new S4/S5 settings, which carry
  # 24 and 47 blocks per replication; the old study lost 87.5% of one config's
  # adaptive_lasso fits to exactly this kind of timeout pressure.
if (!identical(Sys.getenv("AREXT12_DEFINITIONS_ONLY"), "1")) {
  configs <- make_configs(n_sim = as.integer(Sys.getenv("TV11_NSIM", "120")))
  cat(sprintf("Generated %d configs for the ar_ext12 manuscript grid.\n", length(configs)))
  # Workstation-safe staggered init: only TV11_INIT_BATCH workers load rstan per wave,
  # TV11_INIT_DELAY seconds apart (defaults 4 / 120s). resume=TRUE -> re-run the SAME
  # command any time to resume from the per-task checkpoints.
  # "arext12", NOT "arext": this prefix names every checkpoint
  # (arext12_task_config_*.rds). The supervisor's status/sync and the aggregator all
  # glob arext12_*; "arext" here would write files they can never see AND collide
  # with the 3-family ar_ext study's checkpoint names (bug found 2026-08-19: the
  # rename sweep missed it because the prefix is built by sprintf, so the literal
  # string "arext_task_config_" appears nowhere in this file).
  run_batch_simulation(configs, filename = file.path(output_dir, "arext12"),
                       max_cores       = as.integer(Sys.getenv("TV11_CORES", "80")),
                       init_batch_size = as.integer(Sys.getenv("TV11_INIT_BATCH", "16")),
                       init_delay_sec  = as.integer(Sys.getenv("TV11_INIT_DELAY", "180")),
                       resume = TRUE, box_dir = box_dir)
#} else {
#  cat("tv11 simulation script loaded.\n")
#  cat("  Sys.setenv(TV11_MODE='quick'); source('simulation_tv11.R')  # quick check\n")
#  cat("  Sys.setenv(TV11_MODE='full');  source('simulation_tv11.R')  # full grid\n")
#}

}
