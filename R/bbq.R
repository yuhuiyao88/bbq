# BBQ: Bayesian bi-aspect quantile change-point detection (manuscript workflow)
#
# Thin wrappers that run the procedure of the BBQ manuscript (three quantiles,
# Q and QSS bases, unadjusted UI decisions) by calling the existing functions
# with the manuscript's settings. No estimation or inference code lives here:
#   design          -> getSustainedShift / getIsolatedShift          (Sec. 2.1)
#   tuning          -> cv_copss_grid, score loss                     (Sec. 2.5, Eq. 18)
#   final fit       -> getModel, MAP + Hessian Laplace draws         (Secs. 2.2-2.4)
#   detection       -> detectChangepoints_gamma, UI, raw, score      (Sec. 3, Eqs. 20-25)

.bbq_tau_levels <- c(0.25, 0.50, 0.75)

# The prior-specific shrinkage hyperparameter tuned by cross-validation (Sec. 2.5):
# the spike standard deviation for the spike-and-slab-type priors, the gamma
# hyperprior rate b_lambda (a_lambda = 1, Appendix B) for the LASSO-type priors.
.bbq_shrinkage_arg <- function(prior) {
  if (prior %in% c("spike_slab", "spike_slab_lasso")) "spike_sd" else "lambda_lasso2_b"
}

.bbq_shrinkage_grid <- function(prior) {
  if (prior %in% c("spike_slab", "spike_slab_lasso")) c(0.10, 0.15, 0.25, 0.40)
  else c(0.01, 0.05, 0.10, 0.50, 1.00)
}

#' Design shift matrix for BBQ
#'
#' The block design of manuscript Section 2.1: \eqn{r = \lceil (T - w)/l \rceil}
#' columns after a warm-up of \code{w} observations. Calls
#' \code{\link{getSustainedShift}} or \code{\link{getIsolatedShift}}.
#'
#' @param n Number of observations \eqn{T}.
#' @param l Block length.
#' @param w Warm-up (reference) period.
#' @param pattern \code{"sustained"} (default) or \code{"isolated"}.
#' @return The \code{n x r} design matrix \code{H}.
#' @examples
#' H <- bbqDesign(365, l = 30, w = 30)
#' dim(H)
#' @export
bbqDesign <- function(n, l, w, pattern = c("sustained", "isolated")) {
  pattern <- match.arg(pattern)
  if (pattern == "sustained") getSustainedShift(n, l, w) else getIsolatedShift(n, l, w)
}

#' Cross-validation for BBQ hyperparameters
#'
#' Order-preserved two-fold (odd/even) cross-validation of manuscript Section 2.5,
#' ranked by the negative unrestricted score likelihood of Eq. (18). The grid
#' crosses the interquantile weight \eqn{\lambda_{iq}} with the prior-specific
#' shrinkage hyperparameter (Section 4): the spike standard deviation for
#' \code{"spike_slab"} and \code{"spike_slab_lasso"}, the hyperprior rate
#' \eqn{b_\lambda} (with \eqn{a_\lambda = 1}) for the LASSO-type priors. The
#' non-crossing weight is fixed. Calls \code{\link{cv_copss_grid}} with
#' \code{loss = "score"} at \eqn{\tau = (0.25, 0.50, 0.75)}.
#'
#' @param y Numeric series.
#' @param l Block length.
#' @param w Warm-up period.
#' @param prior Prior on the shift coefficients (Section 2.3, Appendix B).
#' @param pattern \code{"sustained"} or \code{"isolated"} (see \code{\link{bbqDesign}}).
#' @param lambda_iq Grid of interquantile weights \eqn{\lambda_{iq}} of Eq. (15).
#'   Passed to \code{getModel} as \code{lambda_iq2 = lambda_iq^2}, since
#'   \code{getModel} applies the square root of \code{lambda_iq2}.
#' @param shrinkage Grid of the prior-specific shrinkage hyperparameter; the
#'   Section 4 grid for \code{prior} when \code{NULL}.
#' @param lambda_nc Non-crossing weight \eqn{\lambda_{nc}} of Eq. (14).
#' @param seed Random seed.
#' @param verbose Print progress.
#' @return The \code{cv_copss_grid} data frame sorted by validation loss (row 1
#'   is the selected configuration), with a \code{lambda_iq} column added and
#'   attributes \code{prior}, \code{shrinkage_arg}, \code{pattern}, \code{l}, \code{w}.
#' @export
bbqTune <- function(y, l, w, prior = c("spike_slab", "spike_slab_lasso", "lasso", "group_lasso",
                             "het_group_lasso", "adaptive_lasso"),
                    pattern = c("sustained", "isolated"),
                    lambda_iq = c(0, 1, 10, 100, 1000), shrinkage = NULL,
                    lambda_nc = 50, seed = 1, verbose = FALSE) {
  prior <- match.arg(prior)
  pattern <- match.arg(pattern)
  hp <- .bbq_shrinkage_arg(prior)
  if (is.null(shrinkage)) shrinkage <- .bbq_shrinkage_grid(prior)
  H <- bbqDesign(length(y), l, w, pattern)
  grid <- expand.grid(lambda_nc = lambda_nc, shrinkage = shrinkage,
                      lambda_iq2 = lambda_iq^2)
  names(grid)[names(grid) == "shrinkage"] <- hp
  cv <- cv_copss_grid(y, .bbq_tau_levels, H = H, w = w, grid = grid,
                      base_args = list(prior_gamma = prior, adaptive_iq = FALSE),
                      loss = "score", seed = seed, verbose = verbose)
  map_par <- attr(cv, "map_par")
  cv$lambda_iq <- sqrt(cv$lambda_iq2)
  attr(cv, "map_par") <- map_par
  attr(cv, "prior") <- prior
  attr(cv, "shrinkage_arg") <- hp
  attr(cv, "pattern") <- pattern
  attr(cv, "l") <- l
  attr(cv, "w") <- w
  cv
}

#' Fit the BBQ model
#'
#' MAP estimation of the penalized log-posterior of manuscript Eq. (17) at
#' \eqn{\tau = (0.25, 0.50, 0.75)}, initialized at the \code{rq.fit.lasso} pilot,
#' with posterior draws from the Hessian-based Laplace approximation. Calls
#' \code{\link{getModel}} with \code{fit_method = "map"} and
#' \code{map_hessian = TRUE}; \code{getModel} stops when the Hessian-based
#' Laplace approximation is unavailable. With the same \code{seed}, repeated
#' calls return identical fits and draws.
#'
#' @inheritParams bbqTune
#' @param lambda_iq Interquantile weight \eqn{\lambda_{iq}} (scalar).
#' @param shrinkage Prior-specific shrinkage hyperparameter (scalar): the spike
#'   standard deviation or \eqn{b_\lambda}, as in \code{\link{bbqTune}}.
#' @param n_draws Number of Laplace posterior draws \eqn{S}. The default 20000 is the
#'   number used in the simulation study: the decorrelation of Eq. (20) uses a
#'   covariance estimated from these draws, and fewer draws inflate the false alarm
#'   probability above its nominal level.
#' @return The \code{getModel} fit, with the design matrix in \code{$H}.
#' @export
bbqFit <- function(y, l, w, prior = c("spike_slab", "spike_slab_lasso", "lasso", "group_lasso",
                             "het_group_lasso", "adaptive_lasso"),
                   pattern = c("sustained", "isolated"),
                   lambda_iq, shrinkage, lambda_nc = 50,
                   n_draws = 20000, seed = 1) {
  prior <- match.arg(prior)
  pattern <- match.arg(pattern)
  H <- bbqDesign(length(y), l, w, pattern)
  args <- list(y = y, taus = .bbq_tau_levels, H = H, w = w, prior_gamma = prior,
               lambda_nc = lambda_nc, lambda_iq2 = lambda_iq^2, adaptive_iq = FALSE,
               fit_method = "map", map_hessian = TRUE, map_init = "pilot",
               laplace_n_samples = n_draws, seed = seed)
  args[[.bbq_shrinkage_arg(prior)]] <- shrinkage
  do.call(getModel, args)
}

#' BBQ change-point detection
#'
#' Section 3 of the manuscript: posterior decorrelation (Eq. 20), the block UI
#' statistic with its posterior probability and charting constant (Eq. 21) on
#' the quantile basis and on the QSS basis (median and IQR shifts, Eqs. 23-24),
#' unadjusted block decisions at \eqn{FAP_0}, cell signals within OOC blocks at
#' the block constant, and localization within each OOC block by the score
#' criterion of Eq. (25). Calls \code{\link{detectChangepoints_gamma}} with
#' \code{basis = c("quantile", "qss")}, \code{statistic = "ui"},
#' \code{adjust = "raw"} and \code{signal_position = "score"}.
#'
#' @param fit A fit from \code{\link{bbqFit}}.
#' @param l Block length.
#' @param w Warm-up period.
#' @param y The series (defaults to \code{fit$y}).
#' @param FAP0 Nominal false alarm probability.
#' @return The \code{detectChangepoints_gamma} result. The OOC blocks are
#'   \code{$tests$quantile$ui$raw} and \code{$tests$qss$ui$raw}, the cell
#'   signals \code{$tests[[basis]]$ui$cells$raw}, and the localized change
#'   points \code{$detected_blocks$signal_obs}.
#' @export
bbqDetect <- function(fit, l, w, y = fit$y, FAP0 = 0.05) {
  detectChangepoints_gamma(fit, .bbq_tau_levels, l = l, w = w, y = y,
                           signal_position = "score", alpha = FAP0,
                           basis = c("quantile", "qss"), statistic = "ui",
                           adjust = "raw")
}

#' BBQ procedure: tuning, final fit and detection
#'
#' Runs the manuscript workflow end to end: \code{\link{bbqTune}} over the
#' Section 4 grid, \code{\link{bbqFit}} on the full series at the selected
#' configuration, and \code{\link{bbqDetect}}.
#'
#' @inheritParams bbqTune
#' @inheritParams bbqFit
#' @param FAP0 Nominal false alarm probability.
#' @return A list with \code{cv} (the tuning table), \code{fit}, \code{detection},
#'   and \code{settings} (prior, pattern, l, w, selected \code{lambda_iq} and
#'   shrinkage, \code{lambda_nc}, \code{FAP0}).
#' @examples
#' \donttest{
#' # AR(1) series, phi = 0.5, location shift of 1 from observation 151
#' set.seed(1)
#' e <- rnorm(180); e[151:180] <- e[151:180] + 1
#' y <- as.numeric(stats::filter(e, 0.5, method = "recursive"))
#' fit <- bbqFit(y, l = 30, w = 30, prior = "spike_slab_lasso",
#'               lambda_iq = 10, shrinkage = 0.1)
#' det <- bbqDetect(fit, l = 30, w = 30)
#' det$tests$qss$ui$raw
#' }
#' @export
bbq <- function(y, l, w, prior = c("spike_slab", "spike_slab_lasso", "lasso", "group_lasso",
                             "het_group_lasso", "adaptive_lasso"),
                pattern = c("sustained", "isolated"),
                lambda_iq = c(0, 1, 10, 100, 1000), shrinkage = NULL,
                lambda_nc = 50, FAP0 = 0.05, n_draws = 20000, seed = 1,
                verbose = FALSE) {
  prior <- match.arg(prior)
  pattern <- match.arg(pattern)
  cv <- bbqTune(y, l, w, prior, pattern, lambda_iq, shrinkage, lambda_nc, seed, verbose)
  hp <- attr(cv, "shrinkage_arg")
  best <- cv[1, ]
  fit <- bbqFit(y, l, w, prior, pattern, lambda_iq = best$lambda_iq,
                shrinkage = best[[hp]], lambda_nc = lambda_nc,
                n_draws = n_draws, seed = seed)
  list(cv = cv, fit = fit, detection = bbqDetect(fit, l, w, y, FAP0),
       settings = list(prior = prior, pattern = pattern, l = l, w = w,
                       lambda_iq = best$lambda_iq, shrinkage = best[[hp]],
                       shrinkage_arg = hp, lambda_nc = lambda_nc, FAP0 = FAP0))
}
