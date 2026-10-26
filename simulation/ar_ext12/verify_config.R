# Static checks only. Does not run models or modify study checkpoints.
# Every ar_ext12 cell must equal the corresponding cell of ar_ext10 (ten settings) or
# ar_ext11 (HK), and the data-generating process must be the one those studies used.
here <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    p <- get0("ofile", envir = sys.frame(i), inherits = FALSE)
    if (is.character(p) && length(p) == 1L && basename(p) == "verify_config.R")
      return(dirname(normalizePath(p, winslash = "/")))
  }
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) return(dirname(normalizePath(sub("^--file=", "", arg[1]), winslash = "/")))
  getwd()
})
ok <- function(label, x) { stopifnot(isTRUE(x)); cat("PASS:", label, "\n") }
defs <- function(file, wanted) {
  e <- new.env(parent = globalenv())
  for (x in parse(file))
    if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]) &&
        as.character(x[[2]]) %in% wanted) eval(x, e)
  e
}
rhs <- function(path, wanted) {
  out <- list()
  for (x in parse(path, keep.source = FALSE))
    if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]) &&
        as.character(x[[2]]) %in% wanted) out[[as.character(x[[2]])]] <- x[[3]]
  out
}
driver <- file.path(here, "simulation_arext12.R")
ext10  <- file.path(here, "..", "ar_ext10", "simulation_arext10.R")
ext11  <- file.path(here, "..", "ar_ext11", "simulation_arext11.R")
ext11n <- file.path(here, "..", "ar_ext11", "simulation_arext11_null.R")
invisible(lapply(list.files(here, "\\.R$", full.names = TRUE), parse))
ok("all study scripts parse", TRUE)

cfg_names <- c("DEFAULT_PRIORS", "get_default_config", "create_config", "make_configs",
               "CONFIG_MAP", ".assert_config_map", "IQ_LAMBDA_GRID", "TAUS", "build_cv_grid",
               "BQQ_BASE_FIELDS", "LMOM_METHODS", "BQQ_METHODS", "CPM_TYPES")
cur <- defs(driver, cfg_names); old10 <- defs(ext10, cfg_names)
old11 <- defs(ext11, cfg_names); old11n <- defs(ext11n, cfg_names)
cfg <- cur$make_configs(n_sim = 120L)
ok("33 cells: 11 settings x {null, mean, variance}", length(cfg) == 33L &&
     identical(unique(vapply(cfg, `[[`, "", "shift_type")), c("null", "mean", "variance")))
key <- function(x) paste(x$setting, x$shift_type)
c10 <- old10$make_configs(n_sim = 120L); names(c10) <- vapply(c10, key, "")
c11 <- c(old11$make_configs(n_sim = 120L), old11n$make_configs(n_sim = 120L))
names(c11) <- vapply(c11, key, "")
for (x in cfg) {
  ref <- if (x$setting == "HK") c11[[key(x)]] else c10[[key(x)]]
  ok(paste("cell", key(x), "identical to", if (x$setting == "HK") "ar_ext11" else "ar_ext10"),
     !is.null(ref) && identical(x, ref))
}
ok("priors, lambda_iq grid and benchmark types unchanged",
   identical(cur$DEFAULT_PRIORS, old10$DEFAULT_PRIORS) &&
     identical(cur$IQ_LAMBDA_GRID, old10$IQ_LAMBDA_GRID) && identical(cur$CPM_TYPES, old10$CPM_TYPES))
for (prior in cur$DEFAULT_PRIORS)
  ok(paste("CV grid", prior), identical(cur$build_cv_grid(prior, cfg[[1]]), old10$build_cv_grid(prior, c10[[1]])))
ok("three fitted quantiles", identical(cur$TAUS, c(.25, .5, .75)))
ok("30 inference methods per prior", identical(cur$BQQ_METHODS, old10$BQQ_METHODS))
sup <- defs(file.path(here, "run_arext12_supervised.R"), "TV11_N_CONFIGS")
ok("supervisor covers all 33 cells", identical(sup$TV11_N_CONFIGS, 33L))

# DGP: functions unchanged from ar_ext10; parameters equal to ar_ext10 / ar_ext11.
fns <- c("simulate_data", "draw_baseline", "draw_post", ".scalemix_s2", ".draw_scalemix",
         ".scalemix_density", "ARIMA_ARGS", "PEARSON_SKEW", "KURT_MIX_NORMAL", "KURT_MIX_SKEWED",
         "POST_MEAN", "POST_SD")
ok("DGP functions, normal-arm constants and ARMA options identical to ar_ext10",
   identical(rhs(driver, fns), rhs(ext10, fns)))
par_names <- c(".pf", "POST_MEAN", "POST_SD", "HK_POST_MEAN", "HK_POST_SD", "PEARSON_SKEW",
               "KURT_MIX_NORMAL", "KURT_MIX_SKEWED", "BASE_PARAMS", "POST_PARAMS")
p12 <- defs(driver, par_names)
p10 <- defs(ext10, par_names)
p11 <- defs(ext11, c(".pf", "POST_MEAN", "POST_SD", "BASE_PARAMS", "POST_PARAMS"))
ok("normal and right-skewed parameters identical to ar_ext10",
   identical(p12$BASE_PARAMS[c("normal", "right_skewed")], p10$BASE_PARAMS[c("normal", "right_skewed")]) &&
     identical(p12$POST_PARAMS[c("normal", "right_skewed")], p10$POST_PARAMS[c("normal", "right_skewed")]))
ok("HK constants equal ar_ext11", identical(p12$HK_POST_MEAN, p11$POST_MEAN) &&
     identical(p12$HK_POST_SD, p11$POST_SD))
ok("HK parameters identical to ar_ext11",
   identical(p12$BASE_PARAMS$high_kurtosis, p11$BASE_PARAMS$high_kurtosis) &&
     identical(p12$POST_PARAMS$high_kurtosis, p11$POST_PARAMS$high_kurtosis))
cat("AR_EXT12 CONFIGURATION CHECKS PASSED\n")
