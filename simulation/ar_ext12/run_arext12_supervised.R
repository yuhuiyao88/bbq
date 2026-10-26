# =============================================================================
# run_arext12_supervised.R  -- crash/stall watchdog for the ar_ext study. (v2)
#
# WHY THIS EXISTS
#   The per-model setTimeLimit() inside the study can only interrupt R at an
#   interpreter safe-point. It CANNOT interrupt native code -- a Stan/BLAS matrix
#   op, or a file syscall that blocks because a mount froze. When that happens the
#   run sits "alive but frozen" with no error, and the only reliable cure is to
#   kill the process at the OS level. This watchdog runs in ITS OWN interactive R
#   session, launches the study as a separate Rscript process, and restarts it
#   whenever it stops showing signs of life. Per-task .rds checkpoints make every
#   restart cheap: only in-flight tasks are redone.
#
# POST-MORTEM OF THE 2026-07-13 OUTAGE (froze at 22:02, dead for 11+ hours)
#   v1 used "a new checkpoint .rds within stall_min = 20 minutes" as liveness.
#   That signal was wrong twice over:
#    (1) One replication legitimately takes 60-140 min, and throughput decayed all
#        day (53 -> 11 checkpoints/hr as the 80 rstan workers grew in memory), so a
#        natural >20-min checkpoint gap eventually HAD to occur -> false "stall".
#    (2) After the kill+relaunch, the fleet needs ~14 min of staggered init plus a
#        full task duration before the FIRST new checkpoint -- far beyond the
#        30-min grace + 20-min stall budget -> every relaunch was itself killed
#        before it could prove progress. Permanent kill loop, zero progress.
#   Two defects compounded it:
#    (3) The "failed to start" branch killed only PIDs in the manifest, so a child
#        killed BEFORE writing its manifest was never killed at all -- relaunches
#        could stack multiple 80-worker clusters onto the machine.
#    (4) system(cmd, wait = FALSE) on Windows bypasses the shell, so the
#        '> tv11_console.log 2>&1' redirect never did anything: no console log was
#        ever written and the child's own error messages were lost.
#
# WHAT v2 DOES DIFFERENTLY
#   - The study now writes per-worker HEARTBEAT files (tv11_hb_<pid>) at every
#     R-level safe point (task start + each model attempt), so liveness is
#     minutes-granular and independent of how long a full replication takes.
#   - Tiered stall rule: if >= busy_workers workers were beating when the silence
#     began, silence > stall_busy_min means frozen (a busy fleet is never all-silent
#     for 45 min). In the end-of-run tail (few workers on long fits) the limit is
#     stall_tail_min, so one legitimately long fit is never kill-looped.
#   - Relaunches terminate the tracked study process tree. The global Rscript.exe
#     sweep is disabled by default and by run_arext12.R, so other R jobs are spared.
#     An explicit sweep_stray = TRUE is a legacy dedicated-workstation option
#     and requires an interactive supervisor.
#   - Planned recycle every recycle_hours to shed worker memory growth (the cause
#     of the 53 -> 11/hr throughput decay). Checkpoint-resume makes this cheap.
#   - The child captures its own console via TV11_CONSOLE (see simulation_arext12.R),
#     one file per launch; the watchdog appends every decision + free RAM to
#     tv11_watchdog.log in the work dir. Next post-mortem will have data.
#
# HOW TO RUN (interactive R or Rscript with the default sweep_stray = FALSE)
#   source("<path>/ar_ext12/run_arext12.R")                  # or this file directly
#   arext12_preflight()                                    # validate this machine
#   arext12_status(target = 20)                            # what is missing
#   run_arext12_supervised(target = 20, recycle_hours = 24)
#
#   To stop: press Esc / Ctrl-C in the console, then run_arext12_kill() to stop the
#   child study process too.
#
# PORTABLE: every path is derived from THIS file's own location, so the folder can
# be copied/synced to any workstation and run unchanged. Overrides: TV11_SCRIPT,
# TV11_BOX_DIR, TV11_WORK_DIR, or the function arguments.
#
# BALANCED MODE (`target = k`): levels EVERY config up to k replications in
# round-robin waves instead of shuffling the whole remaining grid uniformly. It
# schedules ONLY the missing (config, sim) tasks, keeps each task's original seed,
# and never deletes or rewrites an existing checkpoint -- so stopping at any moment
# leaves a balanced panel. Omit `target` for the original full-grid behavior.
# In balanced mode the supervisor stops when all 50 configs reach the target.
# =============================================================================

# ---- portable paths --------------------------------------------------------
# NO machine-specific paths: everything is derived from THIS file's own location,
# so the study runs unchanged on any workstation. Override with TV11_SCRIPT /
# TV11_BOX_DIR / TV11_WORK_DIR, or by passing the arguments explicitly.
.tv11_env <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

# source() is checked FIRST: if this file is sourced from a script that Rscript is
# running, --file= names the OUTER script, not this one.
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

TV11_DIR            <- .tv11_script_dir()
TV11_BOX_DEFAULT    <- .tv11_env("TV11_BOX_DIR", TV11_DIR)
TV11_SCRIPT_DEFAULT <- .tv11_env("TV11_SCRIPT", file.path(TV11_DIR, "simulation_arext12.R"))
TV11_WORK_DEFAULT   <- .tv11_env("TV11_WORK_DIR", file.path(.tv11_home(), "arext12_run"))

# WORKSTATION REQUIREMENT (standing, not negotiable): 80 workers. Never reduced.
# A smaller machine must pass `cores =` explicitly -- see run_tv11_local.R, which
# sets 80% of this computer's logical cores.
TV11_WORKSTATION_CORES <- 80L
tv11_default_cores <- function() TV11_WORKSTATION_CORES

# 80% of the local machine's logical cores, for a box that is not the workstation.
tv11_cores_80pct <- function() {
  n <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)))
  if (is.na(n) || n < 1L) n <- 4L
  max(1L, as.integer(floor(0.8 * n)))
}

# Stan is loaded in exactly FOUR initialization waves, whatever the worker count.
tv11_init_batch <- function(cores) max(1L, as.integer(ceiling(cores / 4)))

# ---- helpers ---------------------------------------------------------------
.tv11_is_windows <- function() .Platform$OS.type == "windows"

.tv11_pid_alive <- function(pid) {
  if (is.na(pid)) return(FALSE)
  if (.tv11_is_windows()) {
    out <- suppressWarnings(try(system(sprintf('tasklist /FI "PID eq %d" /NH', pid),
                                       intern = TRUE, ignore.stderr = TRUE), silent = TRUE))
    if (inherits(out, "try-error")) return(FALSE)
    any(grepl(paste0("\\b", pid, "\\b"), out))
  } else !inherits(try(tools::pskill(pid, signal = 0L), silent = TRUE), "try-error")
}

.tv11_kill_pid_tree <- function(pid) {
  if (is.na(pid)) return(invisible())
  if (.tv11_is_windows())
    suppressWarnings(try(system(sprintf("taskkill /F /T /PID %d", pid),
                                intern = TRUE, ignore.stderr = TRUE), silent = TRUE))
  else try(tools::pskill(pid), silent = TRUE)
  invisible()
}

# Kill the whole study process tree using the PID manifest the study writes.
.tv11_kill_run <- function(pidfile) {
  if (file.exists(pidfile)) {
    pids <- suppressWarnings(as.integer(readLines(pidfile, warn = FALSE)))
    pids <- pids[!is.na(pids)]
    for (p in pids) .tv11_kill_pid_tree(p)   # master first (index 1) cascades via /T
    try(unlink(pidfile), silent = TRUE)
  }
  invisible()
}

# Sweep ANY Rscript.exe. Needed on every relaunch, not just at startup: a child
# killed before it wrote its PID manifest is invisible to .tv11_kill_run, and
# leaving it alive stacks 80-worker clusters until the machine dies (defect 3).
.tv11_kill_stray_workers <- function() {
  if (.tv11_is_windows())
    suppressWarnings(try(system("taskkill /F /IM Rscript.exe",
                                intern = TRUE, ignore.stderr = TRUE), silent = TRUE))
  invisible()
}

# module-level handle to the current child so run_arext12_kill() can find it
.tv11_state <- new.env(parent = emptyenv())

run_arext12_kill <- function(work_dir = TV11_WORK_DEFAULT, sweep_stray = FALSE) {
  .tv11_kill_run(file.path(work_dir, "arext12_run.pids"))
  if (isTRUE(sweep_stray)) {
    .tv11_kill_stray_workers()
    message("Killed the tv11 process tree (manifest PIDs + Rscript.exe sweep).")
  } else {
    message("Killed the tv11 process tree (manifest PIDs only; unrelated Rscript jobs left alone).")
  }
}

# ---- copy checkpoints from the LOCAL work dir to the durable (Box) folder ----
# Hot per-task I/O stays on a local disk (writing it to the synced folder was the
# root cause of the 2026-07-13 freeze), so results must be pushed to box_dir
# separately. The supervisor calls this every sync_min minutes and once more on
# completion; call it BY HAND after stopping a run early, or before moving the
# study to another workstation, to guarantee nothing is left behind.
#   since = NULL  -> copy everything (a full, belt-and-braces sync)
#   since = <time>-> copy only what changed after that time (the periodic mode)
# Copies work_dir -> box_dir only. It never deletes and never writes to work_dir.
tv11_sync <- function(work_dir = TV11_WORK_DEFAULT, box_dir = TV11_BOX_DEFAULT,
                      since = NULL, quiet = FALSE) {
  if (is.null(box_dir) || !nzchar(box_dir)) return(invisible(0L))
  bx <- tryCatch(normalizePath(box_dir,  mustWork = FALSE), error = function(e) box_dir)
  wk <- tryCatch(normalizePath(work_dir, mustWork = FALSE), error = function(e) work_dir)
  if (identical(bx, wk)) {
    if (!quiet) message("work_dir and box_dir are the same folder; nothing to sync.")
    return(invisible(0L))
  }
  if (!dir.exists(work_dir)) {
    if (!quiet) message("work_dir does not exist: ", work_dir)
    return(invisible(0L))
  }
  if (!dir.exists(box_dir)) try(dir.create(box_dir, recursive = TRUE), silent = TRUE)
  patt  <- c("arext12_task_config_*.rds", "arext12_config_*.rds",
             "arext12_combined.csv", "arext12_progress.log", "arext12.rds")
  files <- unlist(lapply(patt, function(p) Sys.glob(file.path(work_dir, p))))
  if (!length(files)) return(invisible(0L))
  if (!is.null(since)) {
    fi    <- file.info(files)
    files <- files[!is.na(fi$mtime) & fi$mtime > since]
  }
  if (!length(files)) return(invisible(0L))
  n_ok <- sum(vapply(files, function(f)
    isTRUE(tryCatch(file.copy(f, file.path(box_dir, basename(f)), overwrite = TRUE),
                    error = function(e) FALSE)), logical(1)))
  if (!quiet)
    message(sprintf("box-sync: %d/%d file(s) copied -> %s", n_ok, length(files), box_dir))
  invisible(n_ok)
}

# ---- import completed replications from OTHER tv11 study folders -------------
# The study can be run from more than one folder (e.g. a workstation folder and a
# second machine's folder). Each machine writes its own; this pulls the others'
# completed replications in so no machine ever recomputes work another has done.
# COPY ONLY: overwrite = FALSE, no delete, no move -- the source folders are never
# modified. Idempotent, so it is safe to call on every launch.

# Schema fingerprint of the CURRENT design. Older studies (e.g. tv_shift_11) use
# the SAME arext12_task_config_<ci>_sim_<s>.rds file names but a different method set
# and row schema, so importing one blindly would silently mix incompatible results
# into the pool. Every source folder is gated on this before anything is copied.
#   188 rows = 8 cpm + 6 priors x 30 methods (20 pkg + 10 lmom)   -- pre-maxent
#   248 rows = 8 cpm + 6 priors x 40 methods (+ 10 maxent)         -- lmom_3cfg era
#   256 rows = 8 cpm + 8 arima_cpm + 6 priors x 40 methods         -- ar_ext, 3 families
#   496 rows = 8 cpm + 8 arima_cpm + 2 x (6 priors x 40)           -- ar_ext5 under bqq 0.6.6
#   376 rows = 8 cpm + 8 arima_cpm + 2 x (6 priors x 30)           -- THIS STUDY (bqq 0.6.7,
#                                                                     maxent retired), the ONLY
#                                                                     count accepted here.
#   ar_ext (3-family, 256 rows) is a DIFFERENT experiment: it has no arima_BQQ rows at
#   all, so pooling would leave that family permanently short. The checkpoint prefix is
#   also distinct (arext12_), so the two studies cannot collide even in one folder.
#   278 rows = 8 cpm + 6 priors x 45 methods (+  5 PIP)            -- PIP era, RETIRED.
#
#   *** 256 WAS MISSING FROM THIS GUARD (found 2026-08-19). Every ar_ext checkpoint
#   written before the four-family change has 256 rows and was therefore REJECTED by
#   .tv11_schema_ok(), so tv11_seed_from() could not import any of them across
#   machines. The local run was unaffected because it globs the folder directly.
#   Add the row count HERE whenever the method registry changes. ***
#     PIP has been removed from the study; 278-row files are historical and are no
#     longer produced. They are still accepted so an existing pool stays readable,
#     but their 5 pip_* rows have no counterpart in new replications.
#   Since bqq 0.6.7 (2026-09-05) only 376-row files are produced and accepted; the
#   maxent family is retired, so the 496-row 0.6.6 files never enter the pool.
#   The tv11 study writes 128-row files under DIFFERENT names, but this guard is the
#   real protection: a 128-row file can never enter this pool even if renamed.
.tv11_schema_ok <- function(f) {
  o <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(o) || is.null(o$rows)) return(FALSE)
  r <- o$rows
  # A checkpoint from bqq 0.5.0 has the SAME 248 rows and the SAME method names
  # as one from 0.5.1, but a different L-moment basis underneath, so the row
  # schema alone cannot separate them. Require the fingerprint written by
  # simulation_arext12.R and require it to MATCH the one this session computes;
  # unstamped files are pre-0.5.1 and are rejected.
  sig_ok <- !is.null(o$basis_sig) &&
            identical(as.character(o$basis_sig), .tv11_basis_sig())
  sig_ok &&
    all(c("method","adjustment","prior_gamma","n_errors") %in% names(r)) &&
    # [ar_ext5 CHANGE 6] val_score/val_pinball are REQUIRED: an ar_ext4 checkpoint
    # has the same 496 rows and method names but pinball-selected hyperparameters
    # and no such columns. The basis fingerprint (bqq0.6.0- vs bqq0.5.2-) already
    # separates them today; this column check keeps doing so even if ar_ext4 were
    # ever re-run under bqq 0.6.0.
    all(c("val_score","val_pinball") %in% names(r)) &&
    nrow(r) == 376L &&
    any(r$adjustment == "q_ui_calib_w",    na.rm = TRUE) &&  # the 20 originals
    any(r$adjustment == "lmom_ui_calib_w", na.rm = TRUE)     # the 10 added
}

tv11_seed_from <- function(src_dirs, dst = TV11_BOX_DEFAULT, quiet = FALSE) {
  # also look inside each source's parked local work dirs, as build_app_csv does
  dirs <- unlist(lapply(src_dirs, function(d)
    c(d, file.path(d, "arext12_run"), file.path(d, "tmp"))))
  dirs <- unique(dirs[dir.exists(dirs)])
  dirs <- dirs[normalizePath(dirs, mustWork = FALSE) !=
                 normalizePath(dst, mustWork = FALSE)]
  if (!length(dirs)) {
    if (!quiet) message("no other study folder found to import from.")
    return(invisible(0L))
  }
  fs <- character(0)
  for (d in dirs) {
    cand <- Sys.glob(file.path(d, "arext12_task_config_*_sim_*.rds"))
    if (!length(cand)) next
    # gate the whole folder on a sampled file's schema
    if (!.tv11_schema_ok(cand[1])) {
      if (!quiet)
        message(sprintf("SKIPPED %s: %d file(s) do NOT match the current row schema ",
                        d, length(cand)),
                "(different study version -- importing them would mix incompatible results)")
      next
    }
    fs <- c(fs, cand)
  }
  if (!length(fs)) return(invisible(0L))
  bn  <- basename(fs)                       # same (config, sim) twice -> keep newest
  ord <- order(bn, file.info(fs)$mtime, decreasing = c(FALSE, TRUE), method = "radix")
  fs  <- fs[ord][!duplicated(bn[ord])]
  n <- 0L
  for (f in fs) {
    d <- file.path(dst, basename(f))
    if (!file.exists(d) && isTRUE(file.copy(f, d, overwrite = FALSE))) n <- n + 1L
  }
  if (!quiet)
    message(sprintf("imported %d new replication(s) (%d already here)",
                    n, length(fs) - n))
  invisible(n)
}

# Every tv11 study folder that is a sibling of this one, so the two machines
# automatically pick up each other's results.
tv11_sibling_dirs <- function(this_dir = TV11_BOX_DEFAULT) {
  p <- dirname(this_dir)
  d <- Sys.glob(file.path(p, "*"))
  d <- d[dir.exists(d)]
  d <- d[vapply(d, function(x) length(Sys.glob(file.path(x, "arext12_task_config_*.rds"))) > 0 ||
                  file.exists(file.path(x, "simulation_arext12.R")), logical(1))]
  d[normalizePath(d, mustWork = FALSE) != normalizePath(this_dir, mustWork = FALSE)]
}

# ---- status: how many USABLE replications does each config have? ------------
# A replication counts only if EVERY BQQ fit in it succeeded (no n_errors == 1).
# A replication holding a timed-out fit carries all-FALSE detection flags, which
# are indistinguishable from genuine non-detections downstream, so counting it
# would silently under-fill the cell. Targeting usable replications tops those
# cells up ADDITIVELY -- no existing file is ever moved, deleted or rewritten.
# Unions the LOCAL work dir and the durable box dir by (config, sim), exactly as
# build_app_csv_from_tasks.R does. Read-only -- safe to call at any time.
#
# Reading every checkpoint each poll would be wasteful, so verdicts are cached on
# (path, mtime); after the first scan only newly written files are read.
.tv11_clean_cache <- new.env(parent = emptyenv())
# The basis fingerprint, computed here so the supervisor can use it BEFORE the study
# script loads (preflight and the census both run first, when BASIS_SIG does not yet
# exist). The formula must stay identical to the one in simulation_arext12.R; the
# preflight below asserts the two agree.
.tv11_sig_cache <- new.env(parent = emptyenv())
# 2026-09-05: ONE formula. The driver defines SCHEMA_HASH_NAMES, .tv11_schema_hash() and
# .tv11_basis_sig_from(); the supervisor evaluates those three definitions from the driver's
# parsed text (plus TAUS) and calls the same function on the same file, so the two stamps
# cannot drift. The stamp is bqq<version>-<lmom weights>-<study-schema hash>.
.tv11_driver_defs <- function(script, names) {
  ex <- parse(script, keep.source = FALSE); env <- new.env(parent = globalenv())
  for (e in ex)
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% names) eval(e, env)
  miss <- setdiff(names, ls(env, all.names = TRUE))
  if (length(miss)) stop("driver does not define: ", paste(miss, collapse = ", "))
  env
}
.tv11_basis_sig <- function(script = TV11_SCRIPT_DEFAULT) {
  if (!is.null(.tv11_sig_cache$sig)) return(.tv11_sig_cache$sig)
  sig <- tryCatch({
    env <- .tv11_driver_defs(script, c("TAUS", "SCHEMA_HASH_NAMES", ".tv11_schema_hash", ".tv11_basis_sig_from"))
    env$.tv11_basis_sig_from(script, env$TAUS)
  }, error = function(e) NA_character_)
  .tv11_sig_cache$sig <- sig
  sig
}

# mtimes are passed in already batched: one vectorized file.info() over the whole
# set instead of one stat call per file (that alone was ~8 s per poll on a synced
# mount). Only files whose mtime changed -- i.e. the ones just written -- are read.
.tv11_rep_clean <- function(f, mt) {
  hit <- .tv11_clean_cache[[f]]
  if (!is.null(hit) && identical(hit$mtime, mt)) return(hit$clean)
  o  <- tryCatch(readRDS(f), error = function(e) NULL)
  # Two ways a replication fails to count.
  #  (a) a fit failed  -- n_errors == 1 somewhere.
  #  (b) it was computed under a DIFFERENT contrast basis. A pre-0.5.1 checkpoint
  #      has the same 248 rows and the same method names as a current one, so
  #      nothing in the row schema distinguishes them; only the stamp does. An
  #      UNSTAMPED file is pre-0.5.1 by definition. Counting those as usable is
  #      what made a fresh study report "600 usable, 0 needed" and do nothing.
  sig <- .tv11_basis_sig()
  cl <- if (is.null(o) || is.null(o$rows)) FALSE
        else if (is.null(o$basis_sig) || !identical(as.character(o$basis_sig), sig)) FALSE
        else if (!"n_errors" %in% names(o$rows)) TRUE
        else !any(o$rows$n_errors == 1L, na.rm = TRUE)
  assign(f, list(mtime = mt, clean = cl), envir = .tv11_clean_cache)
  cl
}

# ---- THE config count. MUST equal length(make_configs()) in simulation_arext12.R ----
# Eleven manuscript settings x three shift arms (null, mean, variance). No imported checkpoints.
#
# This number is NOT cosmetic. arext12_status() sizes its tally with
# tabulate(nbins = n_configs) and reports need over seq_len(n_configs), so a value
# SMALLER than the real grid makes every config above it invisible: they are never
# counted, never scheduled, and the ladder reports "already met" and skips. That is
# exactly what happened on the 2026-08-14 09:57 run, which expanded make_configs()
# to 40 but left this at 15 -- it saw "configs at >= target: 15/15, TASKS NEEDED: 0"
# and jumped straight to target 60 on the old three settings.
# A value LARGER than the real grid is safe but reports phantom deficits forever.
TV11_N_CONFIGS <- 33L

arext12_status <- function(target = 20L, work_dir = TV11_WORK_DEFAULT,
                        box_dir = TV11_BOX_DEFAULT, n_configs = TV11_N_CONFIGS, n_sim = 100L,
                        clean_only = TRUE, quiet = FALSE) {
  dirs <- unique(c(work_dir, box_dir)); dirs <- dirs[dir.exists(dirs)]
  fs <- unlist(lapply(dirs, function(d) Sys.glob(file.path(d, "arext12_task_config_*_sim_*.rds"))))
  bn <- basename(fs)
  m  <- regmatches(bn, regexec("^arext12_task_config_(\\d+)_sim_(\\d+)\\.rds$", bn))
  ok <- lengths(m) == 3
  fs <- fs[ok]
  ci <- as.integer(vapply(m[ok], `[`, "", 2))
  si <- as.integer(vapply(m[ok], `[`, "", 3))
  keep <- !duplicated(paste(ci, si))          # same (config, sim) in both dirs
  fs <- fs[keep]; ci <- ci[keep]
  reps <- tabulate(ci, nbins = n_configs)
  have <- if (clean_only && length(fs)) {
    mts <- file.info(fs)$mtime                       # ONE batched stat over all files
    tabulate(ci[mapply(.tv11_rep_clean, fs, mts, USE.NAMES = FALSE)], nbins = n_configs)
  } else reps
  need <- pmax(0L, as.integer(target) - have)
  if (!quiet) {
    cat(sprintf("tv11 status  (target %d USABLE reps/config, n_sim cap %d)\n", target, n_sim))
    cat(sprintf("  work_dir: %s\n  box_dir : %s\n", work_dir, box_dir))
    cat(sprintf("  %d replications on disk, %d usable (all fits converged)\n",
                sum(reps), sum(have)))
    # Say WHY a file on disk does not count. Without this the census silently
    # reports 0 usable next to 600 on disk and looks like data loss.
    if (sum(reps) > sum(have)) {
      stale <- 0L
      if (length(fs)) {
        sig <- .tv11_basis_sig()
        stale <- sum(vapply(fs, function(f) {
          o <- tryCatch(readRDS(f), error = function(e) NULL)
          is.null(o) || is.null(o$basis_sig) ||
            !identical(as.character(o$basis_sig), sig) }, logical(1)))
      }
      if (stale > 0L)
        cat(sprintf("  %d of those were computed under a DIFFERENT contrast basis\n", stale),
            sprintf("    (no stamp = pre-bqq-0.5.1) and do NOT count. Current basis: %s\n",
                    .tv11_basis_sig()),
            "    They are left on disk; move or delete them if you want a clean folder.\n",
            sep = "")
    }
    cat(sprintf("  usable per config: min %d / median %g / max %d\n",
                min(have), stats::median(have), max(have)))
    cat(sprintf("  configs at >= target: %d/%d   TASKS STILL NEEDED: %d\n",
                sum(have >= target), n_configs, sum(need)))
    if (any(need > 0))
      cat("  deficit:", paste(sprintf("c%d:%d", which(need > 0), need[need > 0]),
                              collapse = " "), "\n")
    if (any(reps > have))
      cat("  (configs holding timed-out fits:",
          paste(sprintf("c%d:%d", which(reps > have), (reps - have)[reps > have]),
                collapse = " "), "-- topped up additively, nothing removed)\n")
  }
  invisible(data.frame(config = seq_len(n_configs), reps = reps, have = have, need = need))
}

# ---- preflight: validate a NEW workstation before committing to a run -------
arext12_preflight <- function(script = TV11_SCRIPT_DEFAULT, work_dir = TV11_WORK_DEFAULT,
                           box_dir = TV11_BOX_DEFAULT, cores = tv11_default_cores(),
                           target = 20L, sweep_stray = FALSE) {
  ok <- TRUE
  say <- function(lbl, pass, detail = "") {
    cat(sprintf("  [%s] %-34s %s\n", if (pass) "OK" else "!!", lbl, detail))
    if (!pass) ok <<- FALSE
  }
  cat("tv11 preflight\n")
  say("study script found", file.exists(script), script)
  if (file.exists(script))
    say("study script parses",
        !inherits(tryCatch(parse(script), error = function(e) e), "error"))
  if (file.exists(script))
    say("balanced-scheduling patch present",
        any(grepl("BALANCED MODE", readLines(script, warn = FALSE))))
  for (p in c("parallel", "cpm", "PearsonDS", "rstan", "bbq", "forecast")) {
    v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
    say(paste("package", p), !is.na(v), if (is.na(v)) "MISSING" else v)
  }
  bqqv <- tryCatch(utils::packageVersion("bbq"), error = function(e) NULL)
  # [ar_ext5 CHANGE 6] 0.6.0 REQUIRED: cv_copss_grid() gained loss = "score" (the
  # CV selection criterion of this study) and the both-loss return columns
  # (val_score / val_pinball). Under 0.5.x the child would die on its own version
  # guard AFTER cluster startup -- fail here instead, before any work is spent.
  # All earlier constraints (0.5.2: no iq_em_update; 0.5.1: full-interval lmom)
  # are inherited. Label and condition must state the SAME version; they drifted
  # once before.
  say("bbq >= 0.1.2", !is.null(bqqv) && bqqv >= "0.1.2",   # 0.6.7: one CV function, maxent retired, 376-row checkpoints
      if (is.null(bqqv)) "" else as.character(bqqv))
  sup_sig <- .tv11_basis_sig()
  say("basis fingerprint computable", !is.na(sup_sig), sup_sig)
  # The study script stamps checkpoints with its OWN copy of the formula. If the two
  # ever drift, every checkpoint this run writes would be rejected by the next one.
  sim_sig <- tryCatch({
    # evaluate the driver's own BASIS_SIG assignment with its helpers, pointing it at `script`
    e <- .tv11_driver_defs(script, c("TAUS", "SCHEMA_HASH_NAMES", ".tv11_schema_hash", ".tv11_basis_sig_from"))
    e$.tv11_this_script <- normalizePath(script, winslash = "/", mustWork = TRUE)
    ex <- parse(script, keep.source = FALSE)
    for (x in ex) if (is.call(x) && identical(as.character(x[[1]]), "<-") && is.name(x[[2]]) &&
                      identical(as.character(x[[2]]), "BASIS_SIG")) eval(x, e)
    e$BASIS_SIG
  }, error = function(err) NA_character_)
  say("study stamp matches supervisor", identical(sim_sig, sup_sig),
      if (identical(sim_sig, sup_sig)) "" else paste("study:", sim_sig))
  for (d in unique(c(work_dir, box_dir))) {
    if (!dir.exists(d)) try(dir.create(d, recursive = TRUE), silent = TRUE)
    t <- file.path(d, ".tv11_write_test")
    w <- isTRUE(tryCatch({ writeLines("x", t); file.remove(t) }, error = function(e) FALSE))
    say("writable", w, d)
  }
  say("work_dir is NOT the synced box_dir",
      normalizePath(work_dir, mustWork = FALSE) != normalizePath(box_dir, mustWork = FALSE),
      "hot I/O must stay on a local disk")
  nc <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)))
  say("cores", !is.na(nc) && cores >= 1 && cores <= nc,
      sprintf("using %d of %s detected (~%.0f GB peak at ~1 GB/worker)", cores,
              ifelse(is.na(nc), "?", nc), cores * 1))
  # The interactive requirement exists ONLY because the stray sweep is
  # `taskkill /F /IM Rscript.exe`, which under Rscript would kill the supervisor
  # itself. With sweep_stray = FALSE that sweep never fires, so headless is fine.
  if (isTRUE(sweep_stray))
    say("interactive session (NOT Rscript)", interactive(),
        "the supervisor sweeps Rscript.exe and would kill itself")
  else
    cat(sprintf("  [OK] %-34s %s\n", "headless allowed",
                "sweep_stray = FALSE -- the Rscript sweep can never fire"))
  st <- arext12_status(target = target, work_dir = work_dir, box_dir = box_dir, quiet = TRUE)
  cat(sprintf("  ..  %-34s %d tasks to reach %d reps/config\n", "work remaining",
              sum(st$need), target))
  cat(if (ok) "PREFLIGHT PASSED\n" else "PREFLIGHT FAILED -- fix the !! lines above\n")
  invisible(ok)
}

# ---- the supervisor --------------------------------------------------------
run_arext12_supervised <- function(
    script            = TV11_SCRIPT_DEFAULT,
    work_dir          = TV11_WORK_DEFAULT,
    box_dir           = TV11_BOX_DEFAULT,
    target            = NULL,     # balanced mode: level every config to this many reps
    n_configs         = TV11_N_CONFIGS,  # five-arm grid -> completion sentinel
    cores             = tv11_default_cores(),   # 80 on the workstation
    init_batch        = NULL,     # NULL -> ceiling(cores/4): exactly FOUR init waves
    init_delay        = 90L,      # seconds between the four waves
    sweep_stray       = FALSE,     # FALSE on a shared box: never taskkill Rscript.exe
    stall_busy_min    = 45,       # busy fleet all-silent this long -> frozen -> relaunch
    stall_tail_min    = 200,      # tail (few workers): allow one long fit (max legit ~3 h)
    busy_workers      = 8L,       # >= this many recently-beating workers = "busy" tier
    recycle_hours     = 8,        # planned restart to shed rstan worker memory growth
    startup_grace_min = 45,       # no stall verdicts until the fleet is up (~15 min ramp)
    poll_sec          = 60,       # how often to check
    sync_min          = 15,       # periodic backup of local checkpoints -> box_dir (minutes)
    launch_settle_min = 5,        # PID breadcrumb deadline: did R boot + load packages
    startup_deadline_min = 40,    # alive but no worker has beaten yet -> startup is wedged
    max_failed_starts = 3L,       # consecutive startups with no worker heartbeat -> STOP
    max_restarts      = 10000L,
    kill_existing     = TRUE      # clear any currently-frozen run on start
) {
  if (!file.exists(script))
    stop("study script not found: ", script,
         "\n  Pass script = \"...\" or set TV11_SCRIPT. Resolved dir: ", TV11_DIR)
  if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
  # Stan loads in exactly FOUR waves regardless of worker count.
  if (is.null(init_batch)) init_batch <- tv11_init_batch(cores)
  # Balanced mode: passed to the child via TV11_TARGET. NULL/0 -> original uniform run.
  if (!is.null(target) && !is.na(target) && target > 0L)
    Sys.setenv(TV11_TARGET = as.character(as.integer(target)))
  else
    Sys.unsetenv("TV11_TARGET")
  # the child hashes its own source for the fingerprint; tell it exactly which file it is
  if (exists("script", inherits = FALSE) && is.character(script) && file.exists(script))
    Sys.setenv(TV11_SCRIPT = normalizePath(script, winslash = "/", mustWork = TRUE))
  rscript <- file.path(R.home("bin"), if (.tv11_is_windows()) "Rscript.exe" else "Rscript")
  fn      <- file.path(work_dir, "arext12")  # MUST match the prefix simulation_arext12.R
                                            # writes; the completion sentinel is
                                            # file.exists(paste0(fn, ".rds")).
  pidfile <- file.path(work_dir, "arext12_run.pids")
  logfile <- file.path(work_dir, "arext_watchdog.log")

  # Every decision goes to the console AND a file, so the next post-mortem has data
  # even if this console is closed (defect 4: last outage left no supervisor trace).
  wlog <- function(fmt, ...) {
    line <- sprintf("[%s] %s", format(Sys.time(), "%m-%d %H:%M:%S"), sprintf(fmt, ...))
    message(line)
    try(cat(line, "\n", file = logfile, append = TRUE), silent = TRUE)
  }
  free_gb <- function() {   # free physical RAM; confirms/refutes the memory-decay theory
    if (!.tv11_is_windows()) return(NA_real_)
    kb <- suppressWarnings(try(system2(
      "powershell", c("-NoProfile", "-Command",
                      "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"),
      stdout = TRUE, stderr = FALSE), silent = TRUE))
    if (inherits(kb, "try-error") || !length(kb)) return(NA_real_)
    round(suppressWarnings(as.numeric(kb[1])) / 1024^2, 1)
  }

  n_task_files <- function() length(Sys.glob(file.path(work_dir, "arext12_task_config_*.rds")))
  n_cfg_files  <- function() length(Sys.glob(file.path(work_dir, "arext12_config_*.rds")))
  # Completion sentinel. In BALANCED mode the full-grid aggregate sentinel does
  # not apply, so completion is "every config has reached target".
  # Without this the supervisor would relaunch forever after the target was met.
  .balanced <- !is.null(target) && !is.na(target) && target > 0L
  is_done <- function() {
    if (.balanced)
      return(sum(arext12_status(target = target, work_dir = work_dir, box_dir = work_dir,
                             n_configs = n_configs, quiet = TRUE)$need) == 0L)
    file.exists(paste0(fn, ".rds")) && n_cfg_files() >= n_configs
  }

  # Liveness = newest mtime over heartbeats AND checkpoints. "busy" = how many
  # workers were still beating in the 30 min before the silence began; that picks
  # the stall tier (a busy fleet all-silent for 45 min is frozen; a lone tail
  # worker inside a 1-3 h fit is not).
  signals <- function() {
    hb <- Sys.glob(file.path(work_dir, "tv11_hb_*"))
    ck <- c(Sys.glob(file.path(work_dir, "arext12_task_config_*.rds")),
            Sys.glob(file.path(work_dir, "arext12_config_*.rds")))
    mt_hb  <- if (length(hb)) file.info(hb)$mtime else as.POSIXct(character(0))
    mt_all <- c(mt_hb, if (length(ck)) file.info(ck)$mtime else as.POSIXct(character(0)))
    if (!length(mt_all)) return(list(newest = as.POSIXct(NA), busy = 0L))
    newest <- max(mt_all, na.rm = TRUE)
    list(newest = newest,
         busy = sum(as.numeric(difftime(newest, mt_hb, units = "mins")) <= 30))
  }

  # Periodic INCREMENTAL backup of the resumable checkpoints from the LOCAL work_dir
  # to box_dir (the durable home = update_tv_shift_11). Only files changed since the
  # last sync are copied, so each sync stays cheap even though hot I/O remains local.
  # Runs in THIS supervisor process (never the workers), so a slow Box write only
  # delays a poll -- it can't freeze the fleet. Returns the timestamp for the next call.
  sync_to_box <- function(since) {
    mark <- Sys.time()
    n <- tv11_sync(work_dir = work_dir, box_dir = box_dir, since = since, quiet = TRUE)
    if (n > 0L) wlog("box-sync: %d changed checkpoint file(s) -> box_dir", n)
    mark
  }

  launch_n <- 0L
  launch <- function() {
    # Pre-launch hook: the child builds its task list ONCE at startup, so anything
    # imported afterwards is invisible to it until the next relaunch. Defining
    # tv11_pre_launch_hook() lets a second machine's results be pulled in right
    # before each (re)launch, so every recycle picks up the other box's progress.
    if (exists("tv11_pre_launch_hook", mode = "function"))
      try(get("tv11_pre_launch_hook")(), silent = TRUE)
    if (file.exists(pidfile)) try(unlink(pidfile), silent = TRUE)
    # Heartbeats from the PREVIOUS run outlive it, and signals() reads their mtimes
    # relative to each other -- so a dead fleet reports a constant, plausible-looking
    # "43 busy workers" forever (it did, for 4.5 h, on 2026-08-04). Clearing them makes
    # "a worker has beaten since THIS launch" a truthful startup signal.
    for (f in Sys.glob(file.path(work_dir, "tv11_hb_*"))) try(unlink(f), silent = TRUE)
    launch_n <<- launch_n + 1L
    console <- file.path(work_dir, sprintf("arext_console_%03d.log", launch_n))
    .tv11_state$console <- console
    Sys.setenv(TV11_WORK_DIR = work_dir, TV11_BOX_DIR = box_dir,
               TV11_CORES = as.character(cores),
               TV11_INIT_BATCH = as.character(init_batch),
               TV11_INIT_DELAY = as.character(init_delay),
               TV11_CONSOLE = console)
    # No shell redirection here: system(wait = FALSE) on Windows bypasses the shell,
    # so '>' would reach Rscript as a literal argument (defect 4). The child
    # redirects its own stdout/stderr to TV11_CONSOLE (see simulation_arext12.R).
    system(sprintf('"%s" "%s"', rscript, script), wait = FALSE)
    .tv11_state$launch_time <- Sys.time()
    wlog("launched study (#%d): %s  [console -> %s]", launch_n, script, basename(console))
  }

  master_pid <- function() {
    if (!file.exists(pidfile)) return(NA_integer_)
    suppressWarnings(as.integer(readLines(pidfile, warn = FALSE))[1])
  }

  # Has ANY worker beaten since this launch? The only unambiguous proof that the child
  # got past package loading, the Box import and makeCluster. Startup progress used to
  # be inferred from the PID manifest, which does not appear until the cluster is up --
  # so a child that was merely slow looked exactly like a child that never ran.
  beat_since_launch <- function() {
    hb <- Sys.glob(file.path(work_dir, "tv11_hb_*"))
    if (!length(hb)) return(FALSE)
    isTRUE(any(file.info(hb)$mtime > .tv11_state$launch_time, na.rm = TRUE))
  }

  # Free space on the work_dir volume. free_gb() reports RAM only, so a full local disk
  # -- which stops the child dead before it can log anything -- was invisible in the
  # watchdog log during the 2026-08-04 outage.
  free_disk_gb <- function() {
    if (!.tv11_is_windows()) return(NA_real_)
    drv <- substr(normalizePath(work_dir, winslash = "/", mustWork = FALSE), 1, 1)
    b <- suppressWarnings(try(system2(
      "powershell", c("-NoProfile", "-Command", sprintf("(Get-PSDrive %s).Free", drv)),
      stdout = TRUE, stderr = FALSE), silent = TRUE))
    if (inherits(b, "try-error") || !length(b)) return(NA_real_)
    round(suppressWarnings(as.numeric(b[1])) / 1024^3, 1)
  }

  # Copy the child's own console into THIS log. Without it a repeated startup failure
  # leaves nothing but identical one-line relaunch records -- 18 of them, in the case
  # this was written for -- and the cause has to be guessed at.
  dump_console <- function(n = 40L) {
    path <- .tv11_state$console
    if (is.null(path) || !file.exists(path)) {
      wlog("  child console not found (%s)", if (is.null(path)) "<unset>" else basename(path))
      return(invisible())
    }
    txt <- suppressWarnings(try(readLines(path, warn = FALSE), silent = TRUE))
    if (inherits(txt, "try-error")) { wlog("  child console %s unreadable.", basename(path)); return(invisible()) }
    if (!length(txt)) {
      wlog("  child console %s is EMPTY -- Rscript.exe died before writing a single line.",
           basename(path))
      return(invisible())
    }
    wlog("  ---- last %d line(s) of %s ----", min(n, length(txt)), basename(path))
    for (l in utils::tail(txt, n)) wlog("  | %s", l)
    wlog("  ---- end %s ----", basename(path))
  }

  kill_everything <- function() {
    .tv11_kill_run(pidfile)        # targeted: manifest PIDs (master first, /T cascades)
    # The stray sweep is `taskkill /F /IM Rscript.exe` -- it kills EVERY Rscript on the
    # machine, including unrelated jobs. Safe on the dedicated workstation, NOT on a
    # shared desktop, so it is opt-out via sweep_stray.
    if (isTRUE(sweep_stray)) .tv11_kill_stray_workers()
  }

  wlog("================ tv11 watchdog v2 ================")
  wlog("script   : %s", script)
  wlog("work_dir : %s   (local, hot I/O)", work_dir)
  wlog("box_dir  : %s   (periodic + final sync)", box_dir)
  wlog("cores=%d  stall busy/tail=%g/%g min  busy_workers=%d  recycle=%gh  grace=%g  poll=%gs  box-sync=%gmin",
       cores, stall_busy_min, stall_tail_min, busy_workers, recycle_hours,
       startup_grace_min, poll_sec, sync_min)
  wlog("==================================================")

  if (kill_existing) {
    wlog("clearing any existing/frozen run (%s)...",
         if (isTRUE(sweep_stray)) "manifest PIDs + Rscript.exe sweep"
         else "manifest PIDs ONLY -- unrelated Rscript jobs left running")
    kill_everything()
    Sys.sleep(3)
  }

  if (is_done()) { wlog("Study already complete -- nothing to do."); return(invisible(TRUE)) }

  launch()
  restarts <- 0L
  fail_starts <- 0L      # consecutive launches that never produced a worker heartbeat
  last_sync <- Sys.time()
  repeat {
    Sys.sleep(poll_sec)
    if (as.numeric(difftime(Sys.time(), last_sync, units = "mins")) >= sync_min)
      last_sync <- sync_to_box(last_sync)
    if (is_done()) {
      kill_everything()        # the child has nothing left to do; stop it cleanly
      sync_to_box(last_sync)   # final incremental backup of the latest checkpoints
      if (.balanced) wlog("COMPLETE: every config has reached target = %d.", target)
      else wlog("COMPLETE: %d/%d configs aggregated and synced to Box.", n_cfg_files(), n_configs)
      break
    }
    now          <- Sys.time()
    since_launch <- as.numeric(difftime(now, .tv11_state$launch_time, units = "mins"))
    has_pid      <- file.exists(pidfile)
    alive        <- has_pid && .tv11_pid_alive(master_pid())
    sg           <- signals()
    # silence is bounded by launch age: signals older than this launch don't count
    since_sig    <- if (is.na(sg$newest)) since_launch
                    else min(since_launch,
                             as.numeric(difftime(now, sg$newest, units = "mins")))
    stall_lim    <- if (sg$busy >= busy_workers) stall_busy_min else stall_tail_min

    # Startup is now judged on worker heartbeats, not on the PID manifest. The child
    # writes its PID immediately (breadcrumb) and its workers beat as soon as the fleet
    # is up, so the three startup outcomes are finally distinguishable:
    #   no breadcrumb        -> Rscript.exe never ran            (seconds to detect)
    #   breadcrumb, no PID   -> child died during startup        (one poll)
    #   alive, but no beats  -> startup wedged (Box mount, disk, makeCluster)
    started <- beat_since_launch()
    if (started) fail_starts <- 0L

    # One decision. Five restart reasons; NULL reason = healthy.
    reason <- NULL; startup_failure <- FALSE
    if (has_pid && !alive) {
      reason <- "child process died"; startup_failure <- !started
    } else if (!has_pid && since_launch > launch_settle_min) {
      # No breadcrumb has THREE causes, not two, and they need different responses.
      # The third one bit on 08-10: the child ran perfectly, found 0 tasks to do,
      # printed "TARGET MET", and exited -- so it never left a PID behind. Reporting
      # that as "never got past package load" sent the reader to look at rstan and
      # package loading, and relaunching it just repeated the no-op every 8 minutes.
      # Read the console before deciding.
      tail_txt <- tryCatch(paste(utils::tail(readLines(.tv11_state$console, warn = FALSE), 60),
                                 collapse = "\n"),
                           error = function(e) "")
      finished_clean <- grepl("TARGET MET for every config", tail_txt, fixed = TRUE) ||
                        grepl("Distributing 0 tasks",        tail_txt, fixed = TRUE)
      if (finished_clean) {
        wlog("child exited with NOTHING TO DO (it printed 'TARGET MET' / 'Distributing 0 tasks').")
        wlog("  This is a terminal condition, not a startup failure -- relaunching would")
        wlog("  just repeat the no-op. Stopping the watchdog.")
        wlog("  If you expected work: every checkpoint on disk is being counted as")
        wlog("  usable. Checkpoints from a different contrast basis (no stamp = pre-bqq")
        wlog("  0.5.1) must not count; current stamp is %s.", .tv11_basis_sig())
        kill_everything()
        break
      }
      # The breadcrumb is written just after the package block, so this means Rscript
      # never ran OR it died loading parallel/cpm/bqq/PearsonDS. The console dump below
      # separates the two: empty file = never ran, error text = package load failure.
      reason <- sprintf("no PID breadcrumb %.1f min after launch (never got past package load)",
                        since_launch)
      startup_failure <- TRUE
    } else if (!started && since_launch > startup_deadline_min) {
      reason <- sprintf("alive but no worker heartbeat %.1f min after launch (startup wedged)",
                        since_launch)
      startup_failure <- TRUE
    } else if (alive && since_launch > startup_grace_min && since_sig > stall_lim) {
      reason <- sprintf("no heartbeat/checkpoint for %.1f min (limit %g, %d busy workers)",
                        since_sig, stall_lim, sg$busy)
    } else if (alive && since_launch > recycle_hours * 60) {
      reason <- sprintf("planned recycle after %.1f h (shed worker memory growth)",
                        since_launch / 60)
    }

    if (is.null(reason)) {
      wlog("ok (%s): %d task + %d cfg files, %.1f min since last signal, %d busy workers, %s GB RAM / %s GB disk free",
           if (!started) "starting up" else "running", n_task_files(), n_cfg_files(),
           since_sig, sg$busy, format(free_gb()), format(free_disk_gb()))
      next
    }

    restarts <- restarts + 1L
    if (startup_failure) fail_starts <- fail_starts + 1L
    wlog("RELAUNCH #%d -- %s. (%d task, %d cfg files so far, %s GB RAM / %s GB disk free)",
         restarts, reason, n_task_files(), n_cfg_files(), format(free_gb()), format(free_disk_gb()))
    kill_everything()
    # The child's file handle is released by the kill above, so whatever it managed to
    # flush is readable now. Always dump on a startup failure -- that is the case where
    # the supervisor has no other information at all.
    if (startup_failure) dump_console()

    # Relaunching a startup that has already failed the same way three times in a row
    # cannot succeed: the cause is outside this process (Box mount, disk, R library,
    # a wedged makeCluster) and only a human can clear it. Looping instead of stopping
    # is what turned a one-off failure into a 4.5 h silent outage on 2026-08-04.
    if (fail_starts >= max_failed_starts) {
      wlog("STOPPING: %d consecutive launches failed to start. This is an environment",
           fail_starts)
      wlog("  problem, not a transient one -- relaunching again would not help. Check, in order:")
      wlog("  1. free disk on the work_dir volume (%s GB free now)", format(free_disk_gb()))
      wlog("  2. Box Drive is signed in and '%s' is browsable", box_dir)
      wlog("  3. library(bbq); library(rstan) load in a fresh R session")
      wlog("  4. the child console tail printed above")
      wlog("  Fix the cause, then re-source run_sim.R. Nothing on disk was lost:")
      wlog("  %d task checkpoints are intact and the run resumes from them.", n_task_files())
      break
    }
    if (restarts > max_restarts) { wlog("max_restarts reached; stopping watchdog."); break }
    Sys.sleep(10)
    launch()
  }
  invisible(is_done())
}
