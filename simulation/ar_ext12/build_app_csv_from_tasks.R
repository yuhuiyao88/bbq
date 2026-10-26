# Use the private study package when it has been installed.
.study_library <- Sys.getenv("AREXT12_R_LIBRARY", file.path(Sys.getenv("USERPROFILE", path.expand("~")), "arext12_library"))
if (dir.exists(.study_library)) .libPaths(c(.study_library, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

# =============================================================================
# ar_ext/build_app_csv_from_tasks.R — aggregate the ar_ext per-replication task
# checkpoints into one row per (config x method x prior x IQ mode).
#
# Adapted from the LOCKED lmom_3cfg/build_app_csv_from_tasks.R. EVERY metric is
# computed by the same expression as there, so ar_ext numbers and lmom_3cfg numbers
# are directly comparable. Four things differ, all forced by the ar_ext design:
#
#   1. ar_phi is READ FROM THE ROW, not hardcoded. The lmom_3cfg aggregator has
#        ar_phi = if (f$dependence == "ar") 0.5 else NA_real_
#      which is correct there (one AR level) and WRONG here (four). ar_ext detail
#      rows carry ar_phi (change E4) and it is used.
#   2. `distribution` is no longer overwritten with "ar". lmom_3cfg wrote
#        distribution = if (f$dependence == "ar") "ar" else f$distribution
#      which folds the dependence axis into the distribution column. With four phi
#      levels that collapses P00/D0/P02/P08 into one label and loses normal vs
#      right_skewed. Here `distribution` stays faithful and the dependence axis is
#      carried by `dependence` + `ar_phi` + `setting`. SEE REUSE.md — this is the
#      one schema change a plot script must know about.
#   3. arima_cpm (reference R2) is handled everywhere cpm is: dropped when it
#      errored, and replicated across (prior, IQ) so app filters keep it visible.
#   4. shift_size_marginal (change E9) is carried alongside shift_size, and
#      pct_signal is written explicitly.
#
# Output columns are a SUPERSET of the lmom_3cfg ones; nothing was renamed or
# removed, so a script written against lmom3_app_combined.csv still runs.
#
# Run:  Rscript build_app_csv_from_tasks.R
# =============================================================================

suppressMessages(library(parallel))

arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
study_dir <- if (length(arg) == 1) dirname(normalizePath(sub("--file=", "", arg))) else normalizePath(".")
src_dirs <- c(study_dir, file.path(study_dir, "arext12_run"), file.path(study_dir, "tmp"))
local_work <- Sys.getenv("AREXT12_WORK_DIR", file.path(Sys.getenv("USERPROFILE", unset = path.expand("~")), "arext12_run"))
if (dir.exists(local_work) &&
    !(normalizePath(local_work) %in% normalizePath(src_dirs[dir.exists(src_dirs)])))
  src_dirs <- c(src_dirs, local_work)

# ---- 1. discover + dedup task files ----------------------------------------
files <- unlist(lapply(src_dirs, function(d)
  Sys.glob(file.path(d, "arext12_task_config_*_sim_*.rds"))))
m  <- regmatches(basename(files),
                 regexec("^arext12_task_config_(\\d+)_sim_(\\d+)\\.rds$", basename(files)))
ok <- lengths(m) == 3
files <- files[ok]
ci  <- as.integer(vapply(m[ok], `[`, "", 2))
sim <- as.integer(vapply(m[ok], `[`, "", 3))
key <- paste(ci, sim)
mt  <- file.info(files)$mtime
ord <- order(key, mt, decreasing = c(FALSE, TRUE), method = "radix")
files_u <- files[ord][!duplicated(key[ord])]
cat(sprintf("Found %d task files (%d unique (config, sim) pairs)\n",
            length(files), length(files_u)))
if (!length(files_u)) stop("no arext task files found in: ",
                           paste(src_dirs, collapse = ", "), call. = FALSE)

# ---- 2. read all rows --------------------------------------------------------
# Reject incompatible checkpoints before collecting their rows.
sup <- new.env(parent = globalenv())
source(file.path(study_dir, "run_arext12_supervised.R"), local = sup)
sup$TV11_SCRIPT_DEFAULT <- file.path(study_dir, "simulation_arext12.R")
read_one <- function(f) {
  if (!isTRUE(sup$.tv11_schema_ok(f)))
    stop("Incompatible ar_ext12 checkpoint: ", f, call. = FALSE)
  readRDS(f)$rows
}
.n_cores <- if (.Platform$OS.type == "windows") 1L else 8L
chunks <- mclapply(files_u, read_one, mc.cores = .n_cores)
bad <- sum(vapply(chunks, is.null, TRUE))
if (bad > 0) cat("WARNING:", bad, "task files failed to read (skipped)\n")
chunks <- chunks[!vapply(chunks, is.null, TRUE)]

all_cols <- unique(unlist(lapply(chunks, names)))
chunks <- lapply(chunks, function(x) {
  for (mm in setdiff(all_cols, names(x))) x[[mm]] <- NA
  x[all_cols]
})
dat <- do.call(rbind, chunks)
cat(sprintf("Read %d detail rows\n", nrow(dat)))

# ---- 2a. provenance guards ---------------------------------------------------
if ("truth_rule" %in% names(dat)) {
  tr <- ifelse(is.na(dat$truth_rule), "pre-2026-08-04 (shift_end window)", dat$truth_rule)
  tb <- table(tr)
  cat("truth-window generations:",
      paste(sprintf("%s = %d rows", names(tb), as.integer(tb)), collapse = " | "), "\n")
  if (length(tb) > 1L)
    cat("*** WARNING: MIXED truth definitions. TP/FP/TN/FN and everything built on\n",
        "***   them are NOT comparable across the two. power/size (any_signal) ARE.\n", sep = "")
}
# ar_phi must be present on every row: without it the phi axis is unrecoverable and
# P02/P08 are indistinguishable except through config_id. Imported lmom_3cfg rows get
# it from import_banked.R; new rows get it from the engine (E4).
if (!"ar_phi" %in% names(dat) || any(is.na(dat$ar_phi))) {
  n_na <- if ("ar_phi" %in% names(dat)) sum(is.na(dat$ar_phi)) else nrow(dat)
  stop("ar_phi missing on ", n_na, " rows. Imported lmom_3cfg checkpoints must be run ",
       "through import_banked.R, which stamps it (0 for ids 1-5, 0.5 for ids 6-10). ",
       "See REUSE.md.", call. = FALSE)
}
cat("phi levels present:",
    paste(sprintf("%.1f (%d rows)", as.numeric(names(table(dat$ar_phi))),
                  as.integer(table(dat$ar_phi))), collapse = " | "), "\n")

# ---- 2b. drop failed fits ----------------------------------------------------
# n_errors == 1 marks a fit that errored or timed out; its detections are recorded as
# all-negative, which drags size and power DOWN. lmom_3cfg dropped these for BQQ; here
# arima_cpm is included, because an auto.arima failure has exactly the same effect on
# the R2 rows. n_sim is recomputed per group below, so the reported count self-corrects.
err <- dat$method %in% c("BQQ", "arima_cpm") & !is.na(dat$n_errors) & dat$n_errors == 1
if (any(err)) {
  by_m <- table(dat$method[err])
  cat(sprintf("Dropping %d errored rows (%.2f%% of all rows): %s\n",
              sum(err), 100 * sum(err) / nrow(dat),
              paste(sprintf("%s=%d", names(by_m), as.integer(by_m)), collapse = ", ")))
  dat <- dat[!err, ]
}

per_cfg <- tapply(dat$sim_id, dat$config_id, function(s) length(unique(s)))
cat(sprintf("Sims per config: min %d / median %d / max %d over %d configs\n",
            min(per_cfg), round(median(per_cfg)), max(per_cfg), length(per_cfg)))

# ---- 3. per-replication derived metrics -------------------------------------
# IDENTICAL to lmom_3cfg. Do not "improve" these: it would break comparability.
neg <- dat$FP + dat$TN
pos <- dat$TP + dat$FN
dat$rep_fpr  <- ifelse(neg > 0, dat$FP / neg, NA_real_)
dat$rep_spec <- ifelse(neg > 0, dat$TN / neg, NA_real_)
dat$rep_sens <- ifelse(pos > 0, dat$TP / pos, NA_real_)
dat$rep_f1   <- ifelse(pos > 0, 2 * dat$TP / (2 * dat$TP + dat$FP + dat$FN), NA_real_)

# ---- 4. aggregate over sims --------------------------------------------------
pg  <- ifelse(is.na(dat$prior_gamma), "-", dat$prior_gamma)
iqm <- ifelse(is.na(dat$iq_mode),     "-", dat$iq_mode)
grp <- interaction(dat$config_id, dat$method, dat$adjustment, pg, iqm,
                   drop = TRUE, sep = "\r")
idx <- split(seq_len(nrow(dat)), grp)

mean_ <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
sd_   <- function(x) { x <- x[is.finite(x)]; if (length(x) > 1) sd(x) else NA_real_ }

agg <- do.call(rbind, lapply(idx, function(ii) {
  d <- dat[ii, ]; f <- d[1, ]
  data.frame(
    config_id = f$config_id, setting = f$setting,
    distribution = f$distribution,          # faithful; NOT overwritten with "ar"
    dependence = f$dependence, shift_type = f$shift_type,
    shift_pattern = f$shift_pattern, design_type = f$shift_pattern,
    shift_length = f$shift_length, n_obs = f$n_obs,
    block_length = f$block_length, warmup = f$warmup,
    shift_start = f$shift_start, tolerance = f$tolerance,
    shift_size          = f$shift_size_tv,            # design (innovation level)
    shift_size_marginal = f$shift_size_tv_marginal,   # observed (filtered marginal)
    gk_g = NA_real_, gk_k = NA_real_,
    ar_phi = f$ar_phi,                      # per-row, not hardcoded
    n_sim = length(unique(d$sim_id)),
    method = f$method, adjustment = f$adjustment,
    prior_gamma = if (is.na(f$prior_gamma)) NA_character_ else f$prior_gamma,
    adaptive_iq = if (is.na(f$iq_mode)) NA_integer_ else as.integer(f$iq_mode == "IQ"),
    # --- detection metrics, expression-for-expression identical to lmom_3cfg ---
    power_mean = mean_(d$any_signal),          power_sd = sd_(d$any_signal),
    detection_rate_mean = mean_(d$hit_within_tol),
    detection_rate_sd   = sd_(d$hit_within_tol),
    detection_delay_mean = mean_(d$detection_delay),
    detection_delay_sd   = sd_(d$detection_delay),
    fpr_mean = mean_(d$rep_fpr),               fpr_sd = sd_(d$rep_fpr),
    sensitivity_mean = mean_(d$rep_sens),      sensitivity_sd = sd_(d$rep_sens),
    specificity_mean = mean_(d$rep_spec),      specificity_sd = sd_(d$rep_spec),
    f1_mean = mean_(d$rep_f1),                 f1_sd = sd_(d$rep_f1),
    # --- fitting / validation measures ---
    # [ar_ext5 CHANGE 6] val_loss == val_score (the CV selection criterion, the
    # exact unconstrained score likelihood); val_pinball is the average pinball
    # loss AT the score-selected configuration. Both are reported per the
    # 2026-08-20 decision. val_loss_mean/sd kept so lmom_3cfg-era readers run.
    val_loss_mean = mean_(d$val_loss),         val_loss_sd = sd_(d$val_loss),
    val_score_mean = mean_(d$val_score),       val_score_sd = sd_(d$val_score),
    val_pinball_mean = mean_(d$val_pinball),   val_pinball_sd = sd_(d$val_pinball),
    error_rate = mean_(d$n_errors),
    # --- arima_cpm model selection (NA on cpm and BQQ rows) ---
    arima_ar1_mean = mean_(d$arima_ar1),
    arima_p_mode = if (all(is.na(d$arima_p))) NA_integer_ else
      as.integer(names(sort(table(d$arima_p), decreasing = TRUE))[1]),
    arima_q_mode = if (all(is.na(d$arima_q))) NA_integer_ else
      as.integer(names(sort(table(d$arima_q), decreasing = TRUE))[1]),
    stringsAsFactors = FALSE
  )
}))
rownames(agg) <- NULL

# LOCKED REPORTING RULE: every stats/detection table carries an explicit
# "% of signal" column -- the significance rate, i.e. size on a null arm and power
# on a shift arm. power_mean is that quantity; pct_signal makes it unambiguous.
agg$pct_signal <- 100 * agg$power_mean

# ---- 5. replicate reference rows across (prior, IQ) so app filters keep them ---
is_ref <- agg$method %in% c("cpm", "arima_cpm")
combos <- unique(agg[!is_ref, c("prior_gamma", "adaptive_iq")])
ref_rep <- do.call(rbind, lapply(seq_len(nrow(combos)), function(k) {
  b <- agg[is_ref, ]
  b$prior_gamma <- combos$prior_gamma[k]
  b$adaptive_iq <- combos$adaptive_iq[k]
  b
}))
out <- rbind(agg[!is_ref, ], ref_rep)
out <- out[order(out$config_id, out$method, out$prior_gamma,
                 out$adaptive_iq, out$adjustment), ]

# ---- 6. write ----------------------------------------------------------------
dst <- file.path(study_dir, "arext12_app_combined.csv")
write.csv(out, dst, row.names = FALSE)
cat(sprintf("Wrote %d aggregated rows -> %s\n", nrow(out), dst))
cat("NOTE: not copied into shiny_app/ -- that folder currently serves the LOCKED\n",
    "  lmom_3cfg study. Point the app at this file deliberately, see REUSE.md.\n", sep = "")

# ---- 7. coverage summary -----------------------------------------------------
cat("\nCoverage by setting:\n")
cv <- unique(out[, c("config_id", "setting", "shift_type", "ar_phi", "n_sim")])
cv <- cv[order(cv$config_id), ]
for (s in unique(cv$setting)) {
  d <- cv[cv$setting == s, ]
  cat(sprintf("  %-4s phi %.1f  ids %2d-%2d  n_sim min %3d / max %3d\n",
              s, d$ar_phi[1], min(d$config_id), max(d$config_id),
              min(d$n_sim), max(d$n_sim)))
}
