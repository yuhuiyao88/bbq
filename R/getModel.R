# --- Session-level cache for compiled Stan model ---
# The Stan code is a fixed string; only stan_data changes between calls.
# Compiling once per R session avoids redundant 30-60s compilations.
.bbq_stan_cache <- new.env(parent = emptyenv())


# ====================================================================
# EM for the squared interquantile penalty weight, lambda_iq2
# ====================================================================
# Derivation: the EM fixed point for lambda_iq2 (see .bbq_iq_em_step below).
#
# The Stan model subtracts  lambda * P  from the target, where
# lambda = sqrt(lambda_iq2) and P = pen_iq_total is the NORMALIZED weighted
# penalty actually accumulated in the model block:
#
#   P = (1/n_comp) * [ (1/(r(m-1))) sum_{j,q} w^g_{q-1,j} |gamma[q,j]-gamma[q-1,j]|
#                    + (1/(p_s(m-1))) sum_{j,q} w^b_{q-1,j} |beta[q,j+1]-beta[q-1,j+1]| ]
#
# Writing the contribution as -lambda * sum_k c_k |d_k| over the N differences,
# the Laplace rate of difference k is a_k = lambda * c_k. Reading the fused factor
# as a Laplace prior on the differences and applying the Gaussian scale-mixture
# augmentation d_k | t_k ~ N(0, t_k), t_k | lambda ~ Exp(a_k^2 / 2) gives the M-step
#
#   lambda^2 = 2N / sum_k c_k^2 E[t_k],     E[t_k] = E|d_k| / a_k + 1 / a_k^2
#
# and, since sum_k c_k^2 E[t_k] = Sbar / lambda^(s) + N / (lambda^(s))^2 with
# Sbar := sum_k c_k E|d_k| = E[P], the closed-form recursion on the SQUARED scale:
#
#   lambda_iq2^(s+1) = 2 N lambda_iq2^(s) / ( sqrt(lambda_iq2^(s)) * Sbar + N )
#
# whose fixed point is lambda* = N / Sbar, i.e. lambda_iq2* = (N / Sbar)^2.
#
# Sbar is therefore just the posterior mean of the very quantity the Stan model
# computes, which is what .bbq_iq_Sbar() evaluates from the draws.
#
# IMPORTANT (see EM_FOR_LAMBDA_IQ_math.md sections 4 and 7): the scale-mixture
# augmentation is exact, but this E-step uses a Laplace approximation to
# p(gamma | y, lambda), not the exact conditional. The procedure is therefore an
# APPROXIMATE empirical-Bayes EM and carries no monotone-ascent guarantee.

# Number of penalized adjacent-quantile differences (the N above).
.bbq_iq_n_diff <- function(r, p_slope, m) {
  if (m < 2L) return(0L)
  as.integer((if (r > 0) r * (m - 1L) else 0L) +
             (if (p_slope > 0) p_slope * (m - 1L) else 0L))
}

# Posterior mean of pen_iq_total, i.e. Sbar = sum_k c_k E|d_k|.
# draws$gamma is [S, m, r]; draws$beta is [S, m, p_total] with column 1 the
# intercept, which is NOT penalized (mirrors the Stan model block exactly).
.bbq_iq_Sbar <- function(draws, w_iq_gamma, w_iq_beta, r, p_slope, m) {
  if (m < 2L) return(NA_real_)
  qs <- 2:m
  qm <- 1:(m - 1L)
  total <- 0
  n_comp <- 0L

  if (r > 0 && !is.null(draws$gamma)) {
    g <- draws$gamma
    dg <- abs(g[, qs, , drop = FALSE] - g[, qm, , drop = FALSE])   # [S, m-1, r]
    Ed <- apply(dg, c(2, 3), mean)                                  # [m-1, r]
    total <- total + sum(w_iq_gamma * Ed) / (r * (m - 1L))
    n_comp <- n_comp + 1L
  }

  if (p_slope > 0 && !is.null(draws$beta)) {
    b <- draws$beta[, , 2:(p_slope + 1L), drop = FALSE]             # drop intercept
    db <- abs(b[, qs, , drop = FALSE] - b[, qm, , drop = FALSE])
    Ed <- apply(db, c(2, 3), mean)
    total <- total + sum(w_iq_beta * Ed) / (p_slope * (m - 1L))
    n_comp <- n_comp + 1L
  }

  if (n_comp == 0L) return(NA_real_)
  total / n_comp
}

# Closed-form E-step (author, 2026-09-04): under the Laplace approximation each adjacent
# difference is d ~ N(dhat, s^2), so |d| = s * (noncentral chi_1 with noncentrality
# |dhat|/s) and E|d| is the folded-normal mean
#   E|d| = s sqrt(2/pi) exp(-dhat^2 / (2 s^2)) + |dhat| (1 - 2 Phi(-|dhat|/s)),
# with dhat the difference of the modes and s^2 = S_qq + S_{q-1,q-1} - 2 S_{q,q-1} from
# the Laplace covariance of the coefficient block. Exact under the same approximation
# the draws come from, and free of Monte Carlo noise. Same normalization as
# .bbq_iq_Sbar(). Returns NA when the fit carries no Laplace covariance (mcmc paths).
.bbq_folded_mean <- function(dhat, s) {
  s <- pmax(s, 1e-12); a <- abs(dhat)
  s * sqrt(2 / pi) * exp(-a^2 / (2 * s^2)) + a * (1 - 2 * stats::pnorm(-a / s))
}
.bbq_iq_Sbar_closed <- function(ls, w_iq_gamma, w_iq_beta, r, p_slope, m) {
  if (is.null(ls) || is.null(ls$eta_cov) || is.null(ls$eta_names) || m < 2L) return(NA_real_)
  nm <- ls$eta_names; mu <- ls$eta_mean; S <- ls$eta_cov
  pair_E <- function(n1, n2) {
    i1 <- match(n1, nm); i2 <- match(n2, nm)
    if (is.na(i1) || is.na(i2)) return(NA_real_)
    .bbq_folded_mean(mu[i1] - mu[i2], sqrt(max(S[i1, i1] + S[i2, i2] - 2 * S[i1, i2], 0)))
  }
  total <- 0; n_comp <- 0L
  if (r > 0) {
    Ed <- matrix(NA_real_, m - 1L, r)
    for (j in seq_len(r)) for (q in 2:m) Ed[q - 1L, j] <- pair_E(sprintf("gamma[%d,%d]", q, j), sprintf("gamma[%d,%d]", q - 1L, j))
    if (anyNA(Ed)) return(NA_real_)
    total <- total + sum(w_iq_gamma * Ed) / (r * (m - 1L)); n_comp <- n_comp + 1L
  }
  if (p_slope > 0) {
    Ed <- matrix(NA_real_, m - 1L, p_slope)
    for (j in seq_len(p_slope)) for (q in 2:m) Ed[q - 1L, j] <- pair_E(sprintf("betaX[%d,%d]", q, j), sprintf("betaX[%d,%d]", q - 1L, j))
    if (anyNA(Ed)) return(NA_real_)
    total <- total + sum(w_iq_beta * Ed) / (p_slope * (m - 1L)); n_comp <- n_comp + 1L
  }
  if (n_comp == 0L) return(NA_real_)
  total / n_comp
}

# One update of lambda_iq2 (manuscript Appendix C).
#   "recursion"  : the M-step recursion, Eq. (C.8), same fixed point, partial moves.
#   "fixedpoint" : DEFAULT since 0.6.6 -- Eq. (C.9), the maximizer of the expected
#                  complete-data log posterior over lambda. No other update is offered.
.bbq_iq_em_step <- function(lambda_iq2, Sbar, N, step = c("fixedpoint", "recursion")) {
  step <- match.arg(step)
  if (!is.finite(Sbar) || !is.finite(lambda_iq2) || N <= 0) return(NA_real_)
  if (Sbar <= 0) return(NA_real_)          # no fusion signal; caller keeps current value
  if (step == "recursion") {
    # Appendix C, Eq. (C.8): the M-step recursion (C.7) with (C.6) substituted,
    #   lambda2_{s+1} = 2 N lambda2_s / (lambda_s Sbar + N),
    # named the primary update in the manuscript; its fixed point is Eq. (C.9). It has the
    # same fixed point but moves only part of the way each iteration, which keeps the
    # iteration from leaping across a fused/unfused boundary and cycling (seen at
    # r = 5 and r = 10 blocks on a real daily series, 2026-09-02).
    return(2 * N * lambda_iq2 / (sqrt(lambda_iq2) * Sbar + N))
  }
  (N / Sbar)^2
}

# Posterior draws of beta/gamma from whichever machinery the fit_method used.
# Laplace draws (fit_method = "map") and MCMC draws share the [S, m, .] layout.
# ====================================================================
# Exact penalized log posterior of manuscript Eq. (17), recomputed in R
# ====================================================================
# Author's ruling (2026-09-04): along the EM chain, monitor exactly Eq. (17),
#   log pi(theta | Y, H) = l*(theta | Y, H) - P_NC(eta0, gamma) - P_IQ(gamma)
#                          + log pi(gamma) + log pi(eta0),
# recomputed here from the returned MAP coefficients and the Stan data. Stan's
# returned value is NOT used: the `~` statements drop normalizing constants
# and the experimental build smooths |d| in the IQ penalty, whereas (17) uses
# Eq. (15) with the exact absolute value.
#   l*      score likelihood, -0.5 ||B||^2 / n, as in the model block
#   P_NC    lambda_nc x mean positive crossing of adjacent quantile curves
#   P_IQ    sqrt(lambda_iq2) x normalized weighted sum of EXACT |d| (Eq. 15)
#   pi(gamma) prior of gamma, written through its hierarchy latents and their
#           hyperpriors exactly as the model block does (the mode is over them)
#   pi(eta0) prior of beta0 and betaX (with their latents)
# The data Jacobian of the log transform is added when log_flag = 1 (it is
# constant in theta unless jittering). All log densities carry their constants.
# `par` is the named vector returned by rstan::optimizing(as_vector = TRUE).
# Named parameter vector in Stan's as_vector naming (name, name[i], name[i,j]) from
# either the vector rstan::optimizing(as_vector = TRUE) returns or the list of
# as_vector = FALSE. Since 0.6.9 the only bare scalar of the Stan program is smooth_T.
.bbq_par_as_vector <- function(pl) {
  if (is.numeric(pl) && !is.null(names(pl))) return(pl)
  if (!is.list(pl)) return(NULL)
  scalars <- c("smooth_T")
  out <- numeric(0)
  for (nm in names(pl)) {
    v <- pl[[nm]]
    if (!is.numeric(v) || length(v) == 0L) next
    if (!is.null(dim(v)) && length(dim(v)) == 2L) {
      nms <- outer(seq_len(nrow(v)), seq_len(ncol(v)), function(i, j) sprintf("%s[%d,%d]", nm, i, j))
      out <- c(out, stats::setNames(as.vector(v), as.vector(nms)))
    } else if (nm %in% scalars && length(v) == 1L) {
      out <- c(out, stats::setNames(as.numeric(v), nm))
    } else {
      out <- c(out, stats::setNames(as.numeric(v), sprintf("%s[%d]", nm, seq_along(v))))
    }
  }
  out
}

.bbq_lp17 <- function(par, sd) {
  par <- .bbq_par_as_vector(par)
  if (!is.numeric(par) || is.null(names(par))) return(NULL)
  gv <- function(prefix, k) {
    if (k == 0L) return(numeric(0))
    unname(par[match(sprintf("%s[%d]", prefix, seq_len(k)), names(par))])
  }
  gm <- function(prefix, nr, nc) {
    if (nr == 0L || nc == 0L) return(matrix(0, nr, nc))
    out <- matrix(NA_real_, nr, nc)
    for (j in seq_len(nc))
      out[, j] <- unname(par[match(sprintf("%s[%d,%d]", prefix, seq_len(nr), j), names(par))])
    out
  }
  gs <- function(nm) unname(par[[nm]])
  n <- as.integer(sd$n); m <- as.integer(sd$m); r <- as.integer(sd$r); px <- as.integer(sd$px)
  p_slope <- as.integer(sd$p_slope)
  tau <- as.numeric(sd$tau_q)

  beta0 <- gv("beta0", m)
  betaX <- gm("betaX", m, px)
  gamma <- gm("gamma", m, r)
  beta  <- cbind(beta0, betaX)                       # m x (px + 1)
  X_design <- cbind(rep(1, n), if (px > 0) sd$X else NULL)
  H <- if (r > 0) sd$H else NULL

  # ---- y_eff (jitter, log) ----
  y <- as.numeric(sd$y)
  jit <- isTRUE(as.numeric(sd$jittering) == 1)
  logf <- isTRUE(as.numeric(sd$log_flag) == 1)
  u <- if (jit) gv("u", n) else NULL
  y_eff <- y
  if (jit) y_eff <- y_eff + u
  if (logf) y_eff <- log(y_eff)

  # ---- eta (n x m) ----
  ETA <- X_design %*% t(beta) + as.numeric(sd$offset)
  if (r > 0) ETA <- ETA + H %*% t(gamma)

  # ---- l*: score likelihood ----
  Z <- cbind(X_design, H)
  pr <- ncol(Z)
  R_i <- y_eff - ETA                                 # n x m residuals
  z <- pmin(pmax(-R_i / as.numeric(sd$base_scale), -20), 20)   # matrix first: pmin/pmax keep dims of arg 1
  PSI <- sweep(-plogis(z), 2, tau, "+")              # tau_q - inv_logit(z)
  S <- crossprod(Z, PSI)                             # pr x m
  Gs <- crossprod(Z) / n + diag(1e-8, pr)
  L_Gs <- t(chol(Gs))
  Qk <- outer(tau, tau, pmin) - outer(tau, tau)
  L_Q <- t(chol(Qk))
  A <- forwardsolve(L_Gs, S)                         # pr x m
  B <- forwardsolve(L_Q, t(A))                       # m x pr
  ll <- -0.5 * sum(B^2) / n

  # ---- P_NC ----
  dtau <- diff(tau)
  DF <- sweep(ETA[, 2:m, drop = FALSE] - ETA[, 1:(m - 1L), drop = FALSE], 2, dtau, "/")
  pen_nc <- as.numeric(sd$lambda_nc) * sum(pmax(0, -DF)) / (n * (m - 1L))

  # ---- P_IQ, Eq. (15): exact |d| ----
  lam_iq <- sqrt(as.numeric(sd$lambda_iq2))
  pen_iq <- 0
  if (lam_iq > 0) {
    tot <- 0; ncomp <- 0L
    if (r > 0) {
      dg <- abs(gamma[2:m, , drop = FALSE] - gamma[1:(m - 1L), , drop = FALSE])
      tot <- tot + sum(sd$w_iq_gamma * dg) / (r * (m - 1L)); ncomp <- ncomp + 1L
    }
    if (p_slope > 0) {
      db <- abs(betaX[2:m, 1:p_slope, drop = FALSE] - betaX[1:(m - 1L), 1:p_slope, drop = FALSE])
      tot <- tot + sum(sd$w_iq_beta * db) / (p_slope * (m - 1L)); ncomp <- ncomp + 1L
    }
    if (ncomp > 0L) pen_iq <- lam_iq * tot / ncomp
  }

  # ---- helpers for the hierarchies ----
  dinvgamma_log <- function(x, a, b) a * log(b) - lgamma(a) - (a + 1) * log(x) - b / x
  dlaplace_log  <- function(x, s) -log(2 * s) - abs(x) / s
  log_mix <- function(pi, l1, l2) { mx <- pmax(l1, l2); mx + log(pi * exp(l1 - mx) + (1 - pi) * exp(l2 - mx)) }
  # marginal group-LASSO prior of a column block G (m x k): lambda^m exp(-lambda ||g||) with its constant
  dgrouplap_log <- function(G, lam) {
    c_grp <- -0.5 * (m + 1) * log(2) - 0.5 * (m - 1) * log(2 * pi) - lgamma(0.5 * (m + 1))
    sum(m * log(lam) - lam * sqrt(colSums(G^2)) + c_grp)
  }

  # ---- log pi(gamma) ----
  lp_gamma <- 0
  if (r > 0) {
    pc <- as.integer(sd$prior_code)
    eff <- as.numeric(sd$lambda_lasso2_fixed)
    if (!(pc %in% c(3L, 5L, 6L)) && as.integer(sd$adaptive_gamma) == 1L) {
      ll2 <- gv("lambda_lasso2", 1L)
      lp_gamma <- lp_gamma + dgamma(ll2, shape = sd$lambda_lasso2_a, rate = sd$lambda_lasso2_b, log = TRUE)
      eff <- ll2
    }
    if (pc == 1L) {
      lp_gamma <- lp_gamma + dgrouplap_log(gamma, sqrt(eff))
    } else if (pc == 2L) {
      lp_gamma <- lp_gamma + sum(dlaplace_log(gamma, 1 / sqrt(eff)))
    } else if (pc == 5L) {
      l2l <- gm("lambda2_gamma_local", m, r)
      lp_gamma <- lp_gamma + sum(dgamma(l2l, shape = sd$lambda_lasso2_a, rate = sd$lambda_lasso2_b, log = TRUE)) +
        sum(dlaplace_log(gamma, 1 / sqrt(l2l)))
    } else if (pc == 3L) {
      pi_s <- gv("pi_slab", 1L)
      lp_gamma <- lp_gamma + dbeta(pi_s, sd$slab_pi_a, sd$slab_pi_b, log = TRUE) +
        sum(log_mix(pi_s, dnorm(gamma, 0, sd$slab_sd, log = TRUE), dnorm(gamma, 0, sd$spike_sd, log = TRUE)))
    } else if (pc == 4L) {
      om <- gv("omega_group", r)
      lp_gamma <- lp_gamma + sum(dinvgamma_log(om, 0.5, 0.5 * eff)) +
        sum(dlaplace_log(gamma, 1 / sqrt(matrix(om, m, r, byrow = TRUE))))
    } else if (pc == 6L) {
      pi_s <- gv("pi_slab", 1L)
      lp_gamma <- lp_gamma + dbeta(pi_s, sd$slab_pi_a, sd$slab_pi_b, log = TRUE) +
        sum(log_mix(pi_s, dlaplace_log(gamma, sd$slab_sd / sqrt(2)), dlaplace_log(gamma, sd$spike_sd / sqrt(2))))   # variance = sd^2, Appendix B
    }
  }

  # ---- log pi(eta0): beta0 and betaX ----
  lp_eta0 <- sum(dnorm(beta0, as.numeric(sd$beta0_loc), as.numeric(sd$beta0_scale), log = TRUE))
  if (px > 0) {
    pb <- as.integer(sd$prior_beta_code)
    effb <- as.numeric(sd$lambda_beta2_fixed)
    if (pb %in% c(2L, 4L, 5L) && as.integer(sd$adaptive_beta) == 1L) {
      lb2 <- gv("lambda_beta2", 1L)
      lp_eta0 <- lp_eta0 + dgamma(lb2, shape = sd$lambda_beta2_a, rate = sd$lambda_beta2_b, log = TRUE)
      effb <- lb2
    }
    if (pb == 1L) {
      lp_eta0 <- lp_eta0 + sum(dnorm(betaX, 0, as.numeric(sd$beta_sd), log = TRUE))
    } else if (pb == 2L) {
      lp_eta0 <- lp_eta0 + sum(dlaplace_log(betaX, 1 / sqrt(effb)))
    } else if (pb == 6L) {
      l2l <- gm("lambda2_beta_local", m, px)
      lp_eta0 <- lp_eta0 + sum(dgamma(l2l, shape = sd$lambda_beta2_a, rate = sd$lambda_beta2_b, log = TRUE)) +
        sum(dlaplace_log(betaX, 1 / sqrt(l2l)))
    } else if (pb == 3L) {
      pi_b <- gv("pi_slab_beta", 1L)
      lp_eta0 <- lp_eta0 + dbeta(pi_b, sd$beta_slab_pi_a, sd$beta_slab_pi_b, log = TRUE) +
        sum(log_mix(pi_b, dnorm(betaX, 0, sd$beta_slab_sd, log = TRUE), dnorm(betaX, 0, sd$beta_spike_sd, log = TRUE)))
    } else if (pb == 4L) {
      lp_eta0 <- lp_eta0 + dgrouplap_log(betaX, sqrt(effb))
    } else if (pb == 5L) {
      om <- gv("omega_beta_group", px)
      lp_eta0 <- lp_eta0 + sum(dinvgamma_log(om, 0.5, 0.5 * effb)) +
        sum(dlaplace_log(betaX, 1 / sqrt(matrix(om, m, px, byrow = TRUE))))
    }
  }
  # u ~ beta(1, 1) has log density 0 on (0, 1).

  # ---- data Jacobian of the log transform (as in the model block) ----
  jac <- 0
  if (logf) jac <- if (jit) -sum(log(y + u)) else -sum(log(y))

  total <- ll - pen_nc - pen_iq + lp_gamma + lp_eta0 + jac
  list(total = total, loglik = ll, pen_nc = pen_nc, pen_iq = pen_iq,
       lp_gamma = lp_gamma, lp_eta0 = lp_eta0, jacobian = jac)
}

.bbq_iq_draws <- function(res) {
  ls <- res$laplace_samples
  if (!is.null(ls) && (!is.null(ls$gamma) || !is.null(ls$beta))) {
    return(list(beta = ls$beta, gamma = ls$gamma))
  }
  if (!is.null(res$fit)) {
    d <- tryCatch(rstan::extract(res$fit, pars = c("beta", "gamma")),
                  error = function(e) NULL)
    if (!is.null(d)) return(list(beta = d$beta, gamma = d$gamma))
  }
  NULL
}

#' Smoothed Quantile Regression with Interquantile Shrinkage (Stan)
#'
#' Fits a multi-quantile (\eqn{m}) regression model where the conditional
#' quantile function is modeled with a per-quantile intercept (in time or index)
#' with optional fixed effects \eqn{X} and structured effects \eqn{H}.
#' The \eqn{H}-coefficients are shrunk via **Bayesian LASSO-type priors**
#' (group, elementwise, heterogeneous-group, or adaptive), and adjacent quantiles
#' are softly penalized to discourage crossings. **Interquantile shrinkage**
#' stabilizes outer quantiles by penalizing differences between adjacent quantile
#' coefficients.
#'
#' @section Model (high level):
#' \describe{
#'   \item{Data & design}{
#'     \itemize{
#'       \item \eqn{y_i} is optionally jittered (\code{u ~ Beta(1,1)}) and/or log-transformed.
#'       \item Combined linear predictor
#'       \eqn{\eta_{qi} = \beta_{0,q} + x_i^\top \beta_{X,q} + h_i^\top \gamma_q + \mathrm{offset}_i}.
#'       \item \eqn{\beta_{0,q}} is a per-quantile intercept whose default prior is
#'       centered at the warm-up-window empirical quantiles with unit-information
#'       scale (Kass & Wasserman 1995) — equivalently, a power prior on the warm-up
#'       window with discount \eqn{a_0 = 1/w} (Ibrahim & Chen 2000).
#'     }
#'   }
#'   \item{Interquantile shrinkage}{
#'     Penalizes differences between adjacent quantile coefficients for gamma and beta slopes
#'     using data-driven adaptive weights from pilot quantile regression estimates:
#'     \eqn{\text{pen}_{\text{IQ}} = \sum_{q=2}^m \sum_j w_{q,j} |\theta_{q,j} - \theta_{q-1,j}|}
#'     where \eqn{w_{q,j} = (|\tilde{\theta}_{q,j} - \tilde{\theta}_{q-1,j}| + \epsilon_w)^{-1}}
#'     and \eqn{\tilde{\theta}} are pilot estimates from separate quantile regressions
#'     (Jiang, Wang, & Bondell 2013).
#'     Weights are applied separately to gamma and beta slopes.
#'     Falls back to uniform weights (all 1) when \pkg{quantreg} is not available.
#'     Note: Intercept (beta0) is NOT penalized (per Jiang, Wang, & Bondell 2013).
#'   }
#'   \item{Non-crossing penalty}{
#'     Adds an L1 hinge on the finite-difference derivative in \eqn{\tau},
#'     scaled by \code{lambda_nc}.
#'   }
#' }
#'
#' @param y Numeric vector of responses of length \eqn{n}.
#' @param taus Numeric vector of target quantile levels in \eqn{(0,1)}, length \eqn{m}.
#' @param H Numeric matrix \eqn{n \times r} of structured predictors for group-lasso
#'   coefficients \eqn{\gamma}. If \eqn{r = 0}, pass a zero-column matrix.
#' @param w Integer \eqn{\ge 1}. Used for initial quantile estimation from first w observations.
#' @param X Optional numeric matrix \eqn{n \times p_x} of user-supplied covariates.
#'   An intercept column is added internally and assigned its own weakly informative prior.
#' @param offset Optional numeric vector of length \eqn{n} added to the linear predictor.
#' @param alpha Deprecated scalar retained for backward compatibility. It is no
#'   longer used in the \code{adaptive_lasso} or \code{het_group_lasso} prior
#'   construction.
#' @param eps_w Positive scalar FLOOR on the pilot adjacent-difference used in the
#'   adaptive IQ weight. Appendix C / manuscript Eq. (15) define
#'   \eqn{w_{q,j} = |\gamma^{pilot}_{q,j} - \gamma^{pilot}_{q-1,j}|^{-1}}; this floor reproduces
#'   that exactly whenever the pilot difference is at least eps_w, and only caps the
#'   weight at 1/eps_w for degenerate (near-zero) differences. Before 0.6.1 the code
#'   used \eqn{(diff + eps_w)^{-1}}, which perturbed EVERY weight rather than only the
#'   degenerate ones
#'   in the IQ shrinkage weights (default 1e-6).
#' @param c_sigma Positive scalar scaling factor for the base scale (default 1.0).
#' @param base_scale Optional positive scalar smoothing bandwidth. Defaults to the
#'   Fernandes, Guerre & Horta (2021) rule-of-thumb
#'   \code{1.06 * s * length(y)^(-1/5)} of manuscript Eq. (12), with
#'   \code{s = min(RMSE_0, IQR_0/1.38898)}, \code{RMSE_0 = sqrt(sum(e^2)/(length(y) - 1))} and
#'   \code{IQR_0} the interquartile range of the median residuals \code{e = y - median(y)}
#'   (the median (tau = 0.5) regression residuals when \code{X} is given), fitted on the ordinary
#'   predictors [intercept | X] only (the change-point design H is excluded). Used
#'   as the smoothing temperature \code{smooth_T = base_scale}; it does NOT set
#'   \code{beta0_scale}.
#' @param beta0_loc Prior location for the per-quantile intercept \code{beta0};
#'   \code{NULL} (default) or a length-\code{m} numeric vector on the modeling
#'   scale of \code{y}. When \code{NULL}, set to the empirical \code{taus}-quantiles
#'   of the warm-up period \code{y_model[1:w]} — the empirical-quantile anchoring
#'   of the tau-specific intercept of Yang and He (2012, Annals of Statistics).
#' @param beta0_scale Prior std dev for the per-quantile intercept \code{beta0};
#'   \code{NULL} (default), a positive scalar (recycled), or a length-\code{m}
#'   vector, on the modeling scale of \code{y} (log scale when \code{log_flag = 1}).
#'   When \code{NULL}, set to the unit-information prior of Kass and Wasserman
#'   (1995): \code{sqrt(taus * (1 - taus)) / f_hat}, with \code{f_hat} a kernel
#'   density estimate of the warm-up period evaluated at \code{beta0_loc}, so the
#'   prior carries the information of a single warm-up observation about each
#'   quantile. Together with the default \code{beta0_loc} this equals the power
#'   prior on the warm-up period with discount \code{a0 = 1/w} (Ibrahim and Chen
#'   2000; Bourazas, Kiagias and Tsiamyrtzis 2022). The prior is
#'   \code{beta0[q] ~ Normal(beta0_loc[q], beta0_scale[q])}.
#' @param beta_sd Positive scalar prior std dev for \code{betaX} coefficients under
#'   \code{prior_beta = "normal"} (default 1.0).
#' @param lambda_nc Positive scalar weight for the non-crossing penalty (larger is stricter).
#' @param lambda_iq2 \code{NULL} (default) or a non-negative scalar, the \strong{squared}
#'   interquantile (IQ) shrinkage weight \eqn{\lambda_{iq}^2}. The penalty applied to the
#'   target is \eqn{\sqrt{\lambda_{iq}^2}\sum_q |\gamma_q - \gamma_{q-1}|}, i.e. the
#'   effective L1 fusion rate is the \emph{square root} of this argument, matching the
#'   \code{lambda_lasso2} / \code{lambda_beta2} convention. Larger = stronger fusion;
#'   0 = none. \code{NULL} means 1. The value is chosen by cross-validation on a grid
#'   (see \code{\link{cv_copss_grid}}); the README workflow uses
#'   \code{c(0, 1, 10, 100, 1000)^2}. When \code{adaptive_iq = TRUE} it is instead the
#'   starting value of the EM and must be strictly positive.
#' @param adaptive_iq Logical; default \code{FALSE} (since 0.6.11): \code{lambda_iq2} is
#'   used as a fixed value and the model is fitted once. If \code{TRUE},
#'   \eqn{\lambda_{iq}^2} is learned by the EM of Appendix C, run between refits.
#'   The EM is no longer the documented workflow: on a change-point design its fixed
#'   point is total fusion (\eqn{\lambda_{iq}} of order \eqn{10^4}), at which scale and
#'   shape shifts are not detectable; the grid search on the held-out score selects
#'   moderate values instead. Each EM iteration is a full refit, so
#'   \code{adaptive_iq = TRUE} costs up to \code{iq_em_max_iter} times a single fit.
#' @param iq_em_max_iter Maximum number of EM iterations (refits) when
#'   \code{adaptive_iq = TRUE} (default 60). Each iteration recomputes \eqn{\bar S}
#'   from a fit at the updated \eqn{\lambda}. The fixed-point update usually needs a
#'   handful of iterations; after a switch to the recursion, allow several dozen.
#' @param iq_smooth Smoothing constant of the absolute value in the interquantile
#'   penalty inside the optimizer: |d| is replaced by sqrt(d^2 + iq_smooth^2) so that
#'   the objective is differentiable at fused solutions and warm restarts move
#'   (default 1e-4; 0 restores the exact absolute value). The EM monitor \code{lp17}
#'   always uses the exact absolute value of Eq. (15).
#' @param iq_em_inner_iter \code{NULL} (default): within every EM iteration the optimizer
#'   runs to its own stopping rule. An integer caps the L-BFGS iterations of EVERY
#'   EM step's optimizer call at that value; \code{1} gives the pattern one optimizer
#'   iteration, one EM update, repeat (each call continues from the previous solution
#'   but rstan restarts the L-BFGS memory). Tested in the ar_ext5 smoke test
#'   (2026-09-04); the E-step then evaluates the Laplace approximation at a point that
#'   is not the mode, so the trace should be read with that in mind.
#' @param iq_em_estep How the E-step expectation E|d| of Appendix C, Eq. (C.6), is
#'   evaluated. \code{"closed"} (default): the folded-normal mean of each adjacent
#'   difference under the Laplace approximation, E|d| = s sqrt(2/pi) exp(-dhat^2/(2 s^2))
#'   + |dhat| (1 - 2 Phi(-|dhat|/s)), with dhat the difference of the modes and s its
#'   Laplace standard deviation -- exact under the approximation and free of Monte
#'   Carlo noise (|d|/s is a noncentral chi with one degree of freedom).
#'   \code{"draws"}: the average of |d| over the \code{laplace_n_samples} draws.
#'   Fits without a Laplace covariance (\code{"mcmc"}, \code{"map_mcmc"}, Hessian
#'   failure) always use the draws. The trace records both values.
#' @param iq_em_lp_tol Stopping rule of the EM chain: the chain stops as soon as the
#'   relative gain in the complete-data log posterior -- manuscript Eq. (17) plus the
#'   normalizing term \eqn{r(m-1) \log \lambda_{iq}} of the fused prior (Appendix C,
#'   Eq. C.1), recomputed exactly at the step's solution (trace columns \code{lp_cd},
#'   \code{lp_cd_gain}; Eq. 17 alone in \code{lp17}) -- is below this fraction of the
#'   previous value's magnitude (floored at 1). Default 1e-2, the classic loose EM
#'   criterion. The term is constant within a
#'   fit, so the inner optimization is unchanged, and it is a constant when
#'   \code{adaptive_iq = FALSE}. The only other exits are a non-finite update and
#'   \code{iq_em_max_iter}.
#' @param iq_em_step How \eqn{\lambda_{iq}^2} is updated at each EM iteration:
#'   \code{"fixedpoint"} (default), the fixed point of the M-step, Appendix C Eq. (C.9),
#'   which the maximization of the expected complete-data log posterior gives directly;
#'   \code{"recursion"}, the augmented-model recursion Eq. (C.8), which has the same
#'   fixed point and moves part of the way per iteration.
#' @param iq_em_warm_jitter The estimation is one chain: the optimizer starts at the
#'   pilot initialization of Section 2.4 once, and every EM iteration after the first
#'   starts the MAP optimization at the previous iteration's full solution
#'   (coefficients and hierarchy latents). If the optimizer returns that start
#'   unchanged (a failed first line search on a kink of the fused penalty), it is
#'   restarted from the same solution with a small jitter on the unbounded
#'   coefficients; \code{iq_em_warm_jitter} gives the jitter standard deviations
#'   tried in order. The trace column \code{warm_status} records "start" (first
#'   iteration), "moved", "jitter<k>" or "stalled". The chain applies to
#'   \code{fit_method = "map"} and to the MAP stage of \code{"map_mcmc"}.

#' @param adaptive_beta Logical; if TRUE (default), the beta-side shrinkage level
#'   \eqn{\lambda_\beta^2} is learned from data for LASSO-type priors. If FALSE,
#'   \code{lambda_beta2_fixed} is used.
#' @param lambda_beta2_a,lambda_beta2_b Positive shape/rate hyperparameters for the
#'   beta-side LASSO-type shrinkage hierarchy.
#' @param lambda_beta2_fixed Positive scalar; fixed value for the beta-side shrinkage
#'   level \eqn{\lambda_\beta^2} when \code{adaptive_beta = FALSE} (default 1).
#' @param lambda_lasso2_a,lambda_lasso2_b Positive shape/rate hyperparameters of the
#'   conjugate gamma hyperprior on the rates of the LASSO-type priors (Appendix B).
#'   For \code{"lasso"}, \code{"group_lasso"} and \code{"het_group_lasso"} they are
#'   the prior of the global rate \eqn{\lambda^2 \sim Gamma(a_\lambda, b_\lambda)}
#'   (Park and Casella, 2008, Sec. 3.2). For \code{"adaptive_lasso"} they are the prior
#'   of every local rate directly, \eqn{\lambda^2_{q,j} \sim Gamma(a_\lambda, b_\lambda)}
#'   (Leng, Tran and Nott, 2014, Eq. 7); the adaptive LASSO has no global rate, so
#'   \code{adaptive_gamma} and \code{lambda_lasso2_fixed} do not affect it.
#' @param adaptive_gamma Logical; if TRUE (default), the global shrinkage level
#'   \eqn{\lambda^2} is learned from data via a Gamma prior. If FALSE, the fixed
#'   value \code{lambda_lasso2_fixed} is used.
#' @param lambda_lasso2_fixed Positive scalar; fixed value for the global shrinkage
#'   level \eqn{\lambda^2} when \code{adaptive_gamma = FALSE} (default 1).
#' @param log_flag Integer \code{0/1}. If 1, fit on \code{log(y)} (or
#'   \code{log(y + u)} with jittering). All data-driven inputs -- the \code{beta0}
#'   prior location, the \code{base_scale} bandwidth, and the pilot fits behind the
#'   adaptive IQ weights -- are then computed on the log scale as well (jitter
#'   approximated by its mean 0.5), so user-supplied \code{beta0_scale} /
#'   \code{base_scale} must also be given on the log scale.
#' @param jittering Integer \code{0/1}. If 1, add \eqn{u \sim \mathrm{Beta}(1,1)} to \eqn{y}.
#' @param chains Number of MCMC chains.
#' @param iter Total iterations per chain.
#' @param warmup Warmup iterations per chain.
#' @param control Optional list passed to \code{rstan::sampling()}.
#' @param seed RNG seed.
#' @param verbose show the log.
#' @param map_hessian Logical; if \code{TRUE} (default), compute the Hessian in the
#'   MAP step. The Hessian enables the standard Laplace approximation for posterior
#'   sampling (see \code{\link{getLaplaceSamples}}). When \code{FALSE}, the Hessian
#'   is not computed and no posterior draws are generated (\code{laplace_samples} is
#'   \code{NULL}); only the MAP point estimate is returned, as in the tuning fits of
#'   \code{cv_copss_grid}. There is no heuristic substitute for the Hessian-based draws.
#' @param map_tol_obj,map_tol_grad,map_tol_rel_grad,map_tol_rel_obj,map_tol_param
#'   MAP optimizer (L-BFGS) tolerances. Defaults are deliberately tight
#'   (\code{map_tol_rel_grad = 1e2}, \code{map_tol_rel_obj = 1e2}, with
#'   \code{map_iter = 10000}): the smoothed-score objective has near-flat plateaus
#'   (see \code{map_init}), and looser relative tolerances (e.g., rstan's default
#'   \code{tol_rel_obj = 1e4}) can terminate there prematurely while reporting
#'   convergence, leaving outer-quantile coefficients stranded far from the optimum.
#' @param map_history_size L-BFGS history size (default 25; rstan's default is 5).
#'   A larger history improves the curvature approximation on this ridge-shaped
#'   objective, typically reducing the iteration count without loosening any
#'   stopping tolerance.
#' @param map_iter Maximum iterations for MAP optimization.
#' @param map_init_values Optional named list of starting values for the MAP
#'   optimization, typically \code{cv_winner_init(cv)}: the fold-averaged
#'   \code{beta0}, \code{betaX} and \code{gamma} of the cross-validation winner.
#'   Components whose dimensions match the model replace the corresponding entries
#'   of the \code{map_init} start; everything else keeps that start. Ignored for
#'   \code{fit_method = "mcmc"}.
#' @param map_init Initialization for MAP optimization. \code{"pilot"} (default)
#'   fits a marginal LASSO quantile regression at each quantile level on
#'   \code{[1 | X | H]} (\code{quantreg::rq.fit.lasso}; intercept unpenalized,
#'   slope penalty at the \code{sqrt(tau(1-tau) n log d)} rate) and uses those
#'   estimates as the starting coefficients (clamped to a generous data range).
#'   The sparse, data-proximal start begins where the smoothed score has live
#'   gradients, matches the spike-slab prior's sparsity, and typically reduces
#'   the iteration count. It falls back to \code{"prior_center"} (with a
#'   warning) when the pilot fits are unavailable (no \code{X}/\code{H} columns,
#'   or \pkg{quantreg} missing).
#'   \code{"prior_center"} starts at the prior mode -- \code{beta0 = beta0_loc},
#'   \code{betaX = 0}, \code{gamma = 0}, hierarchy scales at 1 -- i.e., the
#'   in-control state. Mode-selection nuance: \code{"pilot"} approaches the
#'   optimum from the data-fitting side while \code{"prior_center"} grows the
#'   solution from the null; for multimodal spike-slab posteriors the two can
#'   select different modes, so compare the final log-posterior
#'   (\code{$map$value}) when in doubt. \code{"random"} restores rstan's default
#'   random initialization. Background: the logistic score saturates (clamped at
#'   |z| = 20) once a curve is ~20 bandwidths away from the data, leaving zero
#'   gradient; random inits can start (or wander) into such flat regions and
#'   produce runaway MAP modes.
#' @param fit_method One of "mcmc", "map_mcmc", or "map":
#'   \itemize{
#'     \item "mcmc": Estimators are posterior median from MCMC; posterior draws from MCMC.
#'     \item "map_mcmc": Estimators are MAP; posterior draws from MCMC (MAP used as init).
#'     \item "map": Estimators are MAP; posterior draws from Laplacian approximation.
#'   }
#' @param laplace_n_samples Number of samples for Laplacian approximation (when fit_method = "map").
#' @param prior_beta Prior type for \code{betaX}: \code{"normal"}, \code{"lasso"},
#'   \code{"spike_slab"}, \code{"group_lasso"}, \code{"het_group_lasso"}, or
#'   \code{"adaptive_lasso"}. The intercept \code{beta0} always retains its own
#'   normal prior \code{beta0[q] ~ Normal(quantile(y[1:w], tau_q), beta0_scale[q])}.
#' @param prior_gamma Prior type for gamma: \code{"spike_slab"} (default),
#'   \code{"group_lasso"}, \code{"lasso"}, \code{"het_group_lasso"}, \code{"adaptive_lasso"}, or
#'   \code{"spike_slab_lasso"} (a continuous Laplace spike + Laplace slab mixture,
#'   Rockova & George 2018; reuses \code{spike_sd}/\code{slab_sd} as Laplace scales).
#'   The \code{"adaptive_lasso"} option follows a Leng et al. (2014)-style hierarchy
#'   with coefficient-specific local shrinkage parameters. The
#'   \code{"het_group_lasso"} option combines group-level shrinkage with
#'   coefficient-level local scales, without pilot weights.
#'   Since 0.6.8 every LASSO-type prior is fitted with its normal-scale latents
#'   integrated out: \code{"lasso"} is the Laplace prior with rate
#'   \eqn{\lambda} (Park and Casella, 2008), \code{"group_lasso"} the multivariate
#'   Laplace \eqn{\lambda^m \exp(-\lambda \|\gamma_j\|_2)} of Kyung et al. (2010),
#'   \code{"adaptive_lasso"} the Laplace with local rate \eqn{\lambda_{q,j}} whose
#'   square keeps its gamma hyperprior, and \code{"het_group_lasso"} the Laplace
#'   with block rate \eqn{\sqrt{\omega_j}} whose Levy mixer stays a parameter. The
#'   MAP is therefore the mode of the marginal posterior; the joint mode over the
#'   scale latents does not exist (the density is unbounded as a scale tends to 0).
#' @param beta_spike_sd,beta_slab_sd,beta_slab_pi_a,beta_slab_pi_b Spike-and-slab
#'   hyperparameters for \code{betaX}.
#' @param spike_sd,slab_sd,slab_pi_a,slab_pi_b Spike-and-slab hyperparameters.
#'   NOTE: \code{spike_sd}/\code{slab_sd} are STANDARD DEVIATIONS; the
#'   manuscript's spike-and-slab specification is written in variances
#'   (variance = sd^2). For \code{prior_gamma = "spike_slab_lasso"} the Laplace
#'   components are parametrized so that their variances are likewise
#'   \code{spike_sd^2} and \code{slab_sd^2} (Laplace scale = sd / sqrt(2)), as in
#'   Appendix B.
#'
#'   \code{spike_sd} defaults to \strong{0.1} (raised from 0.05 in 0.5.2). A
#'   narrower spike makes the null mixture component close to a point mass, so
#'   any noise-driven \code{gamma} is pushed into the slab and flagged. In the
#'   lmom_3cfg simulation this dominated the false-alarm rate: on the null arm,
#'   \code{spike_sd <= 0.05} flagged 71 of 155 replications (0.458) against 1 of
#'   85 (0.012) for \code{spike_sd >= 0.10} -- Fisher exact p = 6.5e-16, odds
#'   ratio 70. Values below 0.1 are still accepted but will inflate the
#'   false-alarm rate of \code{prior_gamma = "spike_slab"} and
#'   \code{"spike_slab_lasso"}; the effect is not removable by any multiplicity
#'   correction, since those detections survive Bonferroni, Holm, BH and the
#'   calibrated constant identically.
#'
#'   \code{beta_spike_sd} was raised to 0.1 alongside it, for consistency of the
#'   spike width across the two coefficient blocks. Note the direct evidence
#'   above concerns \code{spike_sd} only: it governs the block coefficients that
#'   drive detection, whereas \code{beta_spike_sd} governs \code{betaX}, and no
#'   equivalent false-alarm study was run for the covariate side.
#'
#' @return A list with components:
#'   \itemize{
#'     \item fit: stanfit object (NULL if fit_method = "map")
#'     \item map: MAP estimates (contains $par with parameter values)
#'     \item y, H, X: Input data and design matrices
#'     \item hessian: Hessian at MAP (if computed)
#'     \item fit_method: The estimation method used
#'     \item laplace_samples: Pre-generated Laplacian samples (if fit_method = "map")
#'     \item stan_data: The data list passed to Stan. Its \code{lambda_iq2} entry is
#'       the value the returned fit was actually computed at.
#'     \item iq_em: Diagnostics for the \eqn{\lambda_{iq}^2} EM: \code{adaptive},
#'       \code{update}, \code{n_diff} (the number of penalized differences \eqn{N}),
#'       \code{lambda_iq2} (final, matching the returned fit), \code{lambda_iq2_next}
#'       (the update the last iteration proposed), \code{Sbar}, \code{converged},
#'       \code{n_iter}, \code{lp_tol} (the \code{iq_em_lp_tol} used), \code{note}, and
#'       \code{trace} -- a data frame with one row per
#'       EM iteration recording \code{lambda_iq2}, \code{lambda_iq}, \code{Sbar},
#'       \code{lambda_iq2_next}, \code{rel_change}, \code{warm_status}, \code{lp17}
#'       (Eq. 17 recomputed at the step's solution) and \code{lp17_gain}.
#'   }
#'
#' @references
#' Jiang, L., Wang, H. J., & Bondell, H. D. (2013). Interquantile Shrinkage in Regression Models.
#' Journal of Computational and Graphical Statistics, 22(4), 970-986.
#'
#' Leng, C., Tran, M. N., & Nott, D. (2014). Bayesian adaptive lasso.
#' Annals of the Institute of Statistical Mathematics, 66(2), 221-244.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 100
#' y <- rnorm(n)
#' taus <- c(0.25, 0.5, 0.75)
#' H <- getIsolatedShift(n, l = 20, w = 20)
#' fit <- getModel(y, taus, H = H, w = 20, fit_method = "map",
#'                 map_hessian = FALSE, map_iter = 500)
#' }
#'
#' @importFrom rstan stan_model sampling stan
#' @importFrom stats lm resid sd rnorm quantile mad median coef cov pchisq qchisq p.adjust
#' @importFrom MASS ginv
#' @export

getModel <- function(y, taus, H = NULL, X = NULL, offset = NULL, w = 0,
                        alpha = 0.75, eps_w = 1e-6, c_sigma = 1.0,
                        base_scale = NULL, beta0_loc = NULL, beta0_scale = NULL,
                        beta_sd = 1.0,
                        lambda_nc = 2,
                        lambda_iq2 = NULL,
                        adaptive_iq = FALSE,
                        iq_em_max_iter = 60,
                        iq_em_step = c("fixedpoint", "recursion"),
                        adaptive_beta = TRUE,
                        lambda_beta2_a = 1, lambda_beta2_b = 0.05,
                        lambda_beta2_fixed = 1,
                        adaptive_gamma = TRUE,
                        lambda_lasso2_a = 1, lambda_lasso2_b = 0.05,
                        lambda_lasso2_fixed = 1,
                        log_flag = 0, jittering = 0,
                        chains = 1, iter = 1500, warmup = 500,
                        control = list(adapt_delta = 0.99),
                        seed = 123, verbose = FALSE,
                        map_hessian = TRUE,
                        map_init = c("pilot", "prior_center", "random"),
                        map_init_values = NULL,
                        map_tol_obj = 1e-12, map_tol_grad = 1e-8,
                        map_tol_rel_grad = 1e2, map_tol_param = 1e-8,
                        map_tol_rel_obj = 1e2,
                        map_iter = 10000,
                        map_history_size = 25,
                        fit_method = c("mcmc", "map_mcmc", "map"),
                        laplace_n_samples = 1000,
                        prior_beta = c("normal", "lasso", "spike_slab",
                                       "group_lasso", "het_group_lasso", "adaptive_lasso"),
                        prior_gamma = c("spike_slab", "group_lasso", "lasso",
                                        "het_group_lasso", "adaptive_lasso",
                                        "spike_slab_lasso"),
                        beta_spike_sd = 0.1, beta_slab_sd = 2.0,
                        beta_slab_pi_a = 1, beta_slab_pi_b = 1,
                        spike_sd = 0.1, slab_sd = 2.0,
                        slab_pi_a = 1, slab_pi_b = 1,
                        # arguments new in the experimental build, kept last so that
                        # positional calls written for 0.6.5 bind unchanged
                        iq_smooth = 1e-4,
                        iq_em_lp_tol = 1e-2,
                        iq_em_warm_jitter = c(0.01, 0.02, 0.05),
                        iq_em_estep = c("closed", "draws"),
                        iq_em_inner_iter = NULL) {

  if (isTRUE(adaptive_iq) && identical(match.arg(fit_method), "map") && !isTRUE(map_hessian)) {
    stop("adaptive_iq = TRUE with fit_method = \"map\" needs the Laplace draws: set map_hessian = TRUE.")
  }
  prior_beta   <- match.arg(prior_beta)
  prior_gamma  <- match.arg(prior_gamma)
  iq_em_estep  <- match.arg(iq_em_estep)
  fit_method   <- match.arg(fit_method)
  map_init     <- match.arg(map_init)

  if (!is.null(lambda_iq2) && (!is.numeric(lambda_iq2) || length(lambda_iq2) != 1L ||
      is.na(lambda_iq2) || lambda_iq2 < 0)) {
    stop("lambda_iq2 must be NULL or a single non-negative scalar (it is lambda^2, not lambda).")
  }
  if (isTRUE(adaptive_iq) && !is.null(lambda_iq2) && lambda_iq2 <= 0) {
    stop("lambda_iq2 must be > 0 when adaptive_iq = TRUE: it is the EM starting value, ",
         "and the recursion has an absorbing fixed point at 0.")
  }
  if (!is.numeric(iq_smooth) || length(iq_smooth) != 1L || !is.finite(iq_smooth) || iq_smooth < 0) {
    stop("iq_smooth must be a single finite non-negative number (0 = exact absolute value).")
  }
  if (!is.numeric(iq_em_lp_tol) || length(iq_em_lp_tol) != 1L || !is.finite(iq_em_lp_tol) || iq_em_lp_tol < 0) {
    stop("iq_em_lp_tol must be a single finite non-negative number.")
  }
  if (!is.null(iq_em_inner_iter) && (!is.numeric(iq_em_inner_iter) || length(iq_em_inner_iter) != 1L ||
      !is.finite(iq_em_inner_iter) || iq_em_inner_iter < 1)) {
    stop("iq_em_inner_iter must be NULL (optimizer runs to its own stopping rule) or a single integer >= 1.")
  }
  if (!is.numeric(iq_em_warm_jitter) || any(!is.finite(iq_em_warm_jitter)) || any(iq_em_warm_jitter < 0)) {
    stop("iq_em_warm_jitter must be a numeric vector of finite non-negative standard deviations ",
         "(an empty vector disables the jittered retries).")
  }
  prior_beta_code <- switch(
    prior_beta,
    normal            = 1L,
    lasso             = 2L,
    spike_slab        = 3L,
    group_lasso       = 4L,
    het_group_lasso   = 5L,
    adaptive_lasso    = 6L
  )
  prior_code <- switch(
    prior_gamma,
    group_lasso        = 1L,
    lasso              = 2L,
    spike_slab         = 3L,
    het_group_lasso    = 4L,
    adaptive_lasso     = 5L,
    spike_slab_lasso   = 6L
  )

  safe_pilot_coefs <- function(y, Z, tau, eps_w = 1e-3, lambda_lasso = NULL) {
    d <- ncol(Z)
    if (d == 0) return(numeric(0))

    fit_q <- try(
      suppressWarnings(quantreg::rq(y ~ Z - 1, tau = tau, method = "fn")),
      silent = TRUE
    )

    if (inherits(fit_q, "try-error")) {
      n <- length(y)
      if (is.null(lambda_lasso)) {
        lambda_lasso <- sqrt(log(d + 1L) / n)
      }
      fit_q <- try(
        suppressWarnings(
          quantreg::rq(y ~ Z - 1, tau = tau,
                       method = "lasso", lambda = lambda_lasso)
        ),
        silent = TRUE
      )
      if (inherits(fit_q, "try-error")) {
        warning("Both rq(method = 'fn') and rq(method = 'lasso') failed at tau=",
                tau, "; using zero pilot coefficients.")
        return(rep(0, d))
      }
    }

    coef_hat <- as.numeric(stats::coef(fit_q))
    if (length(coef_hat) != d) {
      warning("Pilot fit returned length mismatch at tau=", tau,
              "; using zero pilot coefficients.")
      return(rep(0, d))
    }
    coef_hat
  }

  n <- length(y)

  stan_code <- "
  data {
      int<lower=1> n;                  // observations
      int<lower=0> px;                 // user predictors in eta (X), excluding intercept
      int<lower=2> m;                  // quantiles
      int<lower=0> r;                  // predictors in eta (H)

      matrix[n, px] X;                 // n x px
      matrix[n, r] H;                  // n x r
      vector[n] y;
      vector[n] offset;
      vector[m] tau_q;
      vector[m] beta0_loc;
      vector<lower=0>[m] beta0_scale;

      real<lower=1e-12> base_scale;
      real<lower=0>      c_sigma;
      real<lower=0>      beta_sd;

      real<lower=0> lambda_nc;         // non-crossing penalty weight

      // Interquantile (fused lasso) shrinkage: SQUARED penalty weight.
      // Convention matches lambda_lasso2 / lambda_beta2: the stored quantity is
      // lambda^2 and the Laplace rate applied to |gamma[q]-gamma[q-1]| is its
      // square root. Held fixed within a fit; when adaptive_iq = 1 the value is
      // updated between fits by the EM recursion in getModel() (see .bbq_iq_em_step).
      real<lower=0> lambda_iq2;        // SQUARED L1 fusion weight on adjacent-quantile differences
      real<lower=0> iq_smooth2;        // squared smoothing constant of the absolute value in the IQ penalty

      real<lower=0> lambda_lasso2_a;
      real<lower=0> lambda_lasso2_b;
      int<lower=0, upper=1> adaptive_gamma;      // 1 = data-adaptive, 0 = fixed
      real<lower=0> lambda_lasso2_fixed;         // fixed value when adaptive_gamma = 0

      real<lower=0> lambda_beta2_a;
      real<lower=0> lambda_beta2_b;
      int<lower=0, upper=1> adaptive_beta;       // 1 = data-adaptive, 0 = fixed
      real<lower=0> lambda_beta2_fixed;          // fixed value when adaptive_beta = 0

      int<lower=0, upper=1> jittering;
      int<lower=0, upper=1> log_flag;

      // prior selectors
      int<lower=1, upper=6> prior_beta_code;
      int<lower=1, upper=6> prior_code;

      // beta spike-and-slab hyperparameters
      real<lower=0> beta_spike_sd;
      real<lower=0> beta_slab_sd;
      real<lower=0> beta_slab_pi_a;
      real<lower=0> beta_slab_pi_b;

      // spike-and-slab hyperparameters
      real<lower=0> spike_sd;
      real<lower=0> slab_sd;
      real<lower=0> slab_pi_a;
      real<lower=0> slab_pi_b;

      // Data-driven IQ shrinkage weights (from pilot quantile regressions)
      int<lower=0> p_slope;
      matrix[m-1, r] w_iq_gamma;
      matrix[m-1, p_slope] w_iq_beta;
  }

  transformed data {
      // Quantile kernel Q[a,b] = min(tau_a, tau_b) - tau_a * tau_b
      matrix[m, m] Q;
      for (a in 1:m)
        for (b in 1:m)
          Q[a, b] = fmin(tau_q[a], tau_q[b]) - tau_q[a] * tau_q[b];

      // Combined X design with explicit intercept
      int p = px + 1;
      matrix[n, p] X_design;
      for (i in 1:n) {
        X_design[i, 1] = 1;
        if (px > 0) {
          for (j in 1:px)
            X_design[i, j + 1] = X[i, j];
        }
      }

      // Combined design Z = [X | H] (n x pr)
      int pr = p + r;
      matrix[n, pr] Z;
      {
        for (j in 1:p)
          for (i in 1:n)
            Z[i, j] = X_design[i, j];

        if (r > 0) {
          for (j in 1:r)
            for (i in 1:n)
              Z[i, p + j] = H[i, j];
        }
      }

      // Gram for score: Gs = Z'Z / n and its Cholesky
      matrix[pr, pr] Gs;
      matrix[pr, pr] L_Gs;
      if (pr > 0) {
        matrix[pr, n] Zt = Z';
        Gs = (Zt * Z) / n;
        for (k in 1:pr) Gs[k, k] = Gs[k, k] + 1e-8;  // tiny ridge
        L_Gs = cholesky_decompose(Gs);
      } else {
        Gs   = rep_matrix(0, 0, 0);
        L_Gs = rep_matrix(0, 0, 0);
      }

      // IQ shrinkage weights are now data-driven (passed as w_iq_gamma, w_iq_beta)
      // Computed from pilot quantile regression estimates in R
  }

  parameters {
      // Per-quantile intercept and user X-coefficients
      vector[m] beta0;
      matrix[m, px] betaX;

      // H-coefficients
      matrix[m, r] gamma;

      // Hierarchy latents exist only for the prior that uses them (size 0 otherwise),
      // so the Hessian is over the active parameters only.
      // 0.6.8: the normal-scale latents (sigma2_*) of the LASSO-type hierarchies are
      // INTEGRATED OUT (their joint mode does not exist: the density is unbounded as a
      // scale -> 0 with its coefficients); the MAP is taken on the marginal prior. Only
      // the local rates of the adaptive LASSO and the block mixers of the heterogeneous
      // group LASSO remain as parameters.
      // Local adaptive shrinkage rates for beta adaptive lasso
      matrix<lower=0>[m, px * (prior_beta_code == 6)] lambda2_beta_local;

      // Local adaptive shrinkage rates for adaptive lasso
      matrix<lower=0>[m, r * (prior_code == 5)] lambda2_gamma_local;

      // 0.6.9: the global rates and the mixing weights exist only for the priors that use
      // them (length 1) and are absent otherwise (length 0), so no unidentified parameter
      // is left in the Hessian. Stan names them lambda_lasso2[1], pi_slab[1], ...
      // Global beta shrinkage rate (learned when adaptive_beta = 1; betaX LASSO-type priors)
      vector<lower=0>[(px > 0 && adaptive_beta == 1 && (prior_beta_code == 2 || prior_beta_code == 4 || prior_beta_code == 5)) ? 1 : 0] lambda_beta2;

      // Global LASSO rate (learned when adaptive_gamma = 1; gamma LASSO-type priors)
      vector<lower=0>[(r > 0 && adaptive_gamma == 1 && (prior_code == 1 || prior_code == 2 || prior_code == 4)) ? 1 : 0] lambda_lasso2;

      // Spike-and-slab mixing weights
      vector<lower=0, upper=1>[(px > 0 && prior_beta_code == 3) ? 1 : 0] pi_slab_beta;
      vector<lower=0, upper=1>[(prior_code == 3 || prior_code == 6) ? 1 : 0] pi_slab;

      // Group-level mixer for beta hetero group lasso
      vector<lower=0>[px * (prior_beta_code == 5)] omega_beta_group;

      // Group-level mixer for hetero group lasso (Levy)
      // One per time block (consistent with group lasso grouping)
      vector<lower=0>[r * (prior_code == 4)] omega_group;

      // jitter variable (only when jittering = 1)
      vector<lower=1e-12, upper = 1>[n * jittering] u;
  }

  transformed parameters {
      matrix[m, px + 1] beta;
      for (q in 1:m) {
        beta[q, 1] = beta0[q];
        if (px > 0) {
          for (j in 1:px)
            beta[q, j + 1] = betaX[q, j];
        }
      }

      // Smoothing bandwidth (FGH rule-of-thumb; computed on the R side as base_scale)
      real<lower=1e-12> smooth_T = base_scale;

      vector[n] y_eff;
      y_eff = y;

      if (jittering == 1) {
        y_eff = y_eff + u;
      }
      if (log_flag == 1) {
        y_eff = log(y_eff);
      }

      vector[m-1] dtau;
      for (q in 1:(m-1)) dtau[q] = tau_q[q+1] - tau_q[q];
  }

  model {

      // jitter prior
      u ~ beta(1, 1);

      // beta0 prior (per-quantile intercept; informative warm-up-window prior)
      beta0 ~ normal(beta0_loc, beta0_scale);

      // betaX prior (user covariates only)
      if (px > 0) {
        if (prior_beta_code == 1) {
          to_vector(betaX) ~ normal(0, beta_sd);
        } else if (prior_beta_code == 2) {
          real lambda_beta2_eff;
          if (adaptive_beta == 1) {
            lambda_beta2[1] ~ gamma(lambda_beta2_a, lambda_beta2_b);
            lambda_beta2_eff = lambda_beta2[1];
          } else {
            lambda_beta2_eff = lambda_beta2_fixed;
          }
          // marginal of the Park & Casella hierarchy: Laplace with rate sqrt(lambda_beta2)
          for (j in 1:m) {
            for (i in 1:px) {
              target += double_exponential_lpdf(betaX[j, i] | 0, inv_sqrt(lambda_beta2_eff));
            }
          }
        } else if (prior_beta_code == 6) {
          // adaptive lasso on the slopes: local rates with the gamma hyperprior
          // directly (Leng, Tran & Nott 2014, Eq. 7); no global rate in this branch;
          // the normal scale is integrated out -> Laplace with rate sqrt(lambda2_local)
          for (j in 1:m) {
            for (i in 1:px) {
              lambda2_beta_local[j, i] ~ gamma(lambda_beta2_a, lambda_beta2_b);
              target += double_exponential_lpdf(betaX[j, i] | 0, inv_sqrt(lambda2_beta_local[j, i]));
            }
          }
        } else if (prior_beta_code == 3) {
          pi_slab_beta[1] ~ beta(beta_slab_pi_a, beta_slab_pi_b);
          for (j in 1:m) {
            for (i in 1:px) {
              target += log_mix(
                pi_slab_beta[1],
                normal_lpdf(betaX[j, i] | 0, beta_slab_sd),
                normal_lpdf(betaX[j, i] | 0, beta_spike_sd)
              );
            }
          }
        } else if (prior_beta_code == 4) {
          real lambda_beta2_eff;
          if (adaptive_beta == 1) {
            lambda_beta2[1] ~ gamma(lambda_beta2_a, lambda_beta2_b);
            lambda_beta2_eff = lambda_beta2[1];
          } else {
            lambda_beta2_eff = lambda_beta2_fixed;
          }
          // marginal of the Kyung et al. (2010) group hierarchy: lambda^m exp(-lambda ||betaX_i||)
          {
            real lam_b = sqrt(lambda_beta2_eff);
            real c_grp_b = -0.5 * (m + 1) * log2() - 0.5 * (m - 1) * log(2 * pi()) - lgamma(0.5 * (m + 1));
            for (i in 1:px)
              target += m * log(lam_b) - lam_b * sqrt(dot_self(col(betaX, i)) + iq_smooth2) + c_grp_b;
          }
        } else if (prior_beta_code == 5) {
          real c_levy_beta;
          real lambda_beta2_eff;
          if (adaptive_beta == 1) {
            lambda_beta2[1] ~ gamma(lambda_beta2_a, lambda_beta2_b);
            lambda_beta2_eff = lambda_beta2[1];
          } else {
            lambda_beta2_eff = lambda_beta2_fixed;
          }
          c_levy_beta = lambda_beta2_eff;
          for (i in 1:px) {
            omega_beta_group[i] ~ inv_gamma(0.5, 0.5 * c_levy_beta);
            // normal scale integrated out -> Laplace with rate sqrt(omega_i)
            for (j in 1:m) {
              target += double_exponential_lpdf(betaX[j, i] | 0, inv_sqrt(omega_beta_group[i]));
            }
          }
        }
      }

      // ----- Score-based likelihood using Z = [X | H] with logit smoothing -----
      {
        if ((p + r) > 0) {
          matrix[pr, m] S;

          for (q in 1:m) {
            vector[pr] s_q = rep_vector(0.0, pr);
            for (i in 1:n) {
              real xb = dot_product(to_vector(row(X_design, i)), to_vector(beta[q]));
              real hb = (r > 0) ? dot_product(to_vector(row(H, i)), to_vector(gamma[q])) : 0;

              real eta = xb + hb + offset[i];
              real r_i = y_eff[i] - eta;

              real z  = fmin(20, fmax(-20, -r_i / smooth_T));
              real Ilt = inv_logit(z);
              real psi = tau_q[q] - Ilt;

              s_q[1:p]      += to_vector(row(X_design, i)) * psi;
              if (r > 0) s_q[(p+1):pr] += to_vector(row(H, i)) * psi;
            }
            S[, q] = s_q;
          }

          matrix[m, m] L_Q = cholesky_decompose(Q);
          matrix[pr, m] A = mdivide_left_tri_low(L_Gs, S);
          matrix[m, pr] B = mdivide_left_tri_low(L_Q, A');

          target += -0.5 * dot_self(to_vector(B)) / n;
        }
      }

      // ----- Priors on gamma (H-coefficients) -----
      if (r > 0) {

        // Determine effective lambda_lasso2 value for lasso-type priors.
        // Global rate lambda^2 with its gamma hyperprior (Park & Casella 2008, Sec. 3.2):
        // used by the LASSO, group LASSO and heterogeneous group LASSO. The adaptive
        // LASSO (code 5) has no global rate: each local rate carries the gamma
        // hyperprior directly (Leng, Tran & Nott 2014, Eq. 7) -- Appendix B.
        real lambda_lasso2_eff;
        if (prior_code != 3 && prior_code != 5 && prior_code != 6) {
          if (adaptive_gamma == 1) {
            lambda_lasso2[1] ~ gamma(lambda_lasso2_a, lambda_lasso2_b);
            lambda_lasso2_eff = lambda_lasso2[1];
          } else {
            lambda_lasso2_eff = lambda_lasso2_fixed;
          }
        } else {
          lambda_lasso2_eff = lambda_lasso2_fixed;
        }

        // 1 = group lasso: marginal of the Kyung et al. (2010) hierarchy (6),
        //     sigma2_j ~ Gamma((m+1)/2, lambda^2/2), gamma_j | sigma2_j ~ N_m(0, sigma2_j I):
        //     p(gamma_j) = lambda^m 2^{-(m+1)/2} (2 pi)^{-(m-1)/2} / Gamma((m+1)/2) exp(-lambda ||gamma_j||_2).
        //     The norm is smoothed with the IQ constant (exact when iq_smooth2 = 0).
        if (prior_code == 1) {
          real lam = sqrt(lambda_lasso2_eff);
          real c_grp = -0.5 * (m + 1) * log2() - 0.5 * (m - 1) * log(2 * pi()) - lgamma(0.5 * (m + 1));
          for (i in 1:r) {
            target += m * log(lam) - lam * sqrt(dot_self(col(gamma, i)) + iq_smooth2) + c_grp;
          }

        // 2 = lasso: marginal of the Park & Casella (2008) hierarchy, Laplace with rate lambda
        } else if (prior_code == 2) {
          for (j in 1:m) {
            for (i in 1:r) {
              target += double_exponential_lpdf(gamma[j, i] | 0, inv_sqrt(lambda_lasso2_eff));
            }
          }

        // 5 = adaptive lasso with coefficient-specific local rates, each with the
        //     conjugate gamma hyperprior Gamma(a, b) of Leng, Tran & Nott (2014), Eq. (7)
        //     (their r = a, delta = b); no global rate in this branch (Appendix B).
        //     The normal scale is integrated out: gamma_qj | lambda2_qj ~ Laplace(rate sqrt(lambda2_qj)).
        } else if (prior_code == 5) {
          for (j in 1:m) {
            for (i in 1:r) {
              lambda2_gamma_local[j, i] ~ gamma(lambda_lasso2_a, lambda_lasso2_b);
              target += double_exponential_lpdf(gamma[j, i] | 0, inv_sqrt(lambda2_gamma_local[j, i]));
            }
          }

        // 3 = spike-and-slab
        } else if (prior_code == 3) {
          pi_slab[1] ~ beta(slab_pi_a, slab_pi_b);
          for (j in 1:m) {
            for (i in 1:r) {
              target += log_mix(
                pi_slab[1],
                normal_lpdf(gamma[j, i] | 0, slab_sd),
                normal_lpdf(gamma[j, i] | 0, spike_sd)
              );
            }
          }

        // 4 = heterogeneous group lasso with Levy mixing
        // Groups by time block (consistent with group lasso prior_code=1)
        // omega_group[i]: block-level Levy scale (one per time block), kept as a parameter
        // element scales integrated out: gamma_qj | omega_i ~ Laplace(rate sqrt(omega_i))
        } else if (prior_code == 4) {
          real c_levy = lambda_lasso2_eff;
          for (i in 1:r) {
            omega_group[i] ~ inv_gamma(0.5, 0.5 * c_levy);
            for (j in 1:m) {
              target += double_exponential_lpdf(gamma[j, i] | 0, inv_sqrt(omega_group[i]));
            }
          }

        // 6 = spike-and-slab LASSO (Laplace spike + Laplace slab; Rockova & George 2018).
        // The Laplace spike contributes zero curvature away from 0, so the Laplace-
        // approximation posterior variance is driven by the likelihood (honest) rather
        // than a tiny prior variance -- avoids the near-Dirac Gaussian-spike degeneracy.
        } else if (prior_code == 6) {
          pi_slab[1] ~ beta(slab_pi_a, slab_pi_b);
          for (j in 1:m) {
            for (i in 1:r) {
              // Appendix B writes the components as Laplace(0, sigma^2) with sigma^2 a
              // VARIANCE; Stan's double_exponential takes the scale b, Var = 2 b^2,
              // so b = sigma / sqrt(2) (author's alignment ruling, 2026-09-04).
              target += log_mix(
                pi_slab[1],
                double_exponential_lpdf(gamma[j, i] | 0, slab_sd / sqrt2()),
                double_exponential_lpdf(gamma[j, i] | 0, spike_sd / sqrt2())
              );
            }
          }
        }
      }

      // ---- Non-crossing penalty ----
      {
        real pen = 0;

        for (i in 1:n) {
          vector[m] eta_row;
          for (q in 1:m) {
            real xb = dot_product(to_vector(row(X_design, i)), to_vector(beta[q]));
            real hb = (r > 0) ? dot_product(to_vector(row(H, i)), to_vector(gamma[q])) : 0;
            eta_row[q] = xb + hb + offset[i];
          }

          // Non-crossing penalty: penalize negative derivatives
          for (q in 1:(m-1)) {
            real dfdtau = (eta_row[q + 1] - eta_row[q]) / dtau[q];
            pen += fmax(0, -dfdtau);
          }
        }

        pen /= (n * (m - 1));
        target += - lambda_nc * pen;
      }

      // ---- Interquantile (fused lasso) shrinkage penalty ----
      // Penalizes |coef[q] - coef[q-1]| with IQ weights (more shrinkage at outer quantiles)
      // Applied to gamma and betaX slopes (NOT the intercept beta0)
      {
        // Effective IQ penalty weight is the SQUARE ROOT of lambda_iq2, mirroring
        // the lambda_lasso2 convention where the stored quantity is lambda^2 and
        // the Laplace rate is lambda. Fixed within a fit; updated across fits by EM.
        real lambda_iq_eff = sqrt(lambda_iq2);

        if (lambda_iq_eff > 0) {
          real pen_iq_gamma = 0;
          real pen_iq_beta = 0;
          int n_components = 0;

          // Penalty on gamma (H-coefficients / change-point effects)
          if (r > 0) {
            for (j in 1:r) {
              for (q in 2:m) {
                pen_iq_gamma += w_iq_gamma[q-1, j] * sqrt(square(gamma[q, j] - gamma[q-1, j]) + iq_smooth2);
              }
            }
            pen_iq_gamma /= (r * (m - 1));
            n_components += 1;
          }

          // Penalty on beta (X-coefficients EXCLUDING intercept)
          // Note: Column 1 of X is the intercept - do NOT penalize it
          // (Jiang, Wang, & Bondell 2013 only penalize slope coefficients)
          if (p_slope > 0) {
            for (j in 1:p_slope) {
              for (q in 2:m) {
                pen_iq_beta += w_iq_beta[q-1, j] * sqrt(square(beta[q, j+1] - beta[q-1, j+1]) + iq_smooth2);
              }
            }
            pen_iq_beta /= (p_slope * (m - 1));
            n_components += 1;
          }

          // Average across components
          if (n_components > 0) {
            real pen_iq_total = (pen_iq_gamma + pen_iq_beta) / n_components;
            target += - lambda_iq_eff * pen_iq_total;
          }
        }
      }

      // Jacobian adjustment for log transform
      if (log_flag == 1) {
        if (jittering == 1) {
          target += -sum(log(y + u));
        } else {
          target += -sum(log(y));
        }
      }
  }
  "

  # ---------------- R-side pre-processing ----------------

  n <- length(y)
  m <- length(taus)

  if (is.null(offset)) {
    offset <- rep(0, n)
  }

  if (is.null(X)) {
    px <- 0L
    X <- matrix(0, nrow = n, ncol = 0)
  } else {
    X <- as.matrix(X)
    px <- ncol(X)
  }

  if (is.null(H)) {
    r <- 0
    H <- matrix(0, n, 0)
  } else {
    r <- ncol(H)
    if (r == 0) {
      H <- matrix(0, n, 0)
    }
  }

  # Modeling-scale response: the Stan likelihood operates on y_eff = log(y (+ u))
  # when log_flag = 1, so every data-driven quantity below -- the beta0 prior
  # location, the smoothing bandwidth, and the pilot fits for the adaptive IQ
  # weights -- must be computed on that same scale. The jitter u ~ Beta(1,1) is
  # approximated by its prior mean 0.5.
  if (log_flag == 1) {
    y_shift <- if (jittering == 1) 0.5 else 0
    if (any(y + y_shift <= 0)) {
      stop("log_flag = 1 requires y (plus the 0.5 jitter midpoint when ",
           "jittering = 1) to be strictly positive.")
    }
    y_model <- log(y + y_shift)
  } else {
    y_model <- y
  }

  # Smoothing bandwidth (manuscript Eq. 12; Fernandes, Guerre & Horta 2021):
  #   base_scale = 1.06 * s * n^(-1/5),  s = min(RMSE_0, IQR_0 / 1.38898),
  #   RMSE_0 = sqrt(sum((y_t - Q_2(y))^2) / (n - 1)),  IQR_0 = Q_3(y) - Q_1(y),
  # computed on the median (tau = 0.5) residuals y_t - Q_2(y); with X these are the
  # median regression residuals. The pilot median fit uses only the
  # ordinary predictors [intercept | X]; the change-point design H is excluded so the
  # shift structure being detected cannot absorb the residual scale. Overridable via
  # the base_scale argument.
  if (is.null(base_scale)) {
    if (!is.null(X) && ncol(as.matrix(X)) > 0 &&
        requireNamespace("quantreg", quietly = TRUE)) {
      Xd_med <- cbind(1, as.matrix(X))
      fit_med <- try(suppressWarnings(
        quantreg::rq(y_model ~ Xd_med - 1, tau = 0.5, method = "fn")), silent = TRUE)
      med_resid <- if (inherits(fit_med, "try-error")) {
        y_model - as.numeric(stats::quantile(y_model, 0.5))
      } else {
        y_model - as.numeric(Xd_med %*% as.numeric(stats::coef(fit_med)))
      }
    } else {
      med_resid <- y_model - as.numeric(stats::quantile(y_model, 0.5))
    }
    rmse0 <- sqrt(sum(med_resid^2) / (length(med_resid) - 1))
    s_disp <- min(rmse0, stats::IQR(med_resid) / 1.38898)
    base_scale <- max(1e-8, 1.06 * s_disp * length(y)^(-1 / 5))
  }

  # beta0 prior location: defaults to the per-quantile empirical level of the
  # warm-up period, on the modeling scale (log scale when log_flag = 1) --
  # empirical-quantile anchoring of the tau-specific intercept in the sense of
  # Yang & He (2012, Ann. Statist., p. 1107).
  if (is.null(beta0_loc)) {
    if (w > 0) {
      beta0_loc <- stats::quantile(y_model[1:w], probs = taus)
    } else {
      beta0_loc <- stats::quantile(y_model, probs = taus)
    }
  }
  beta0_loc <- as.vector(beta0_loc)
  if (length(beta0_loc) != m || anyNA(beta0_loc)) {
    stop("beta0_loc must be NULL or a length-m numeric vector without NAs.")
  }

  # beta0 prior scale: defaults to the unit-information prior (Kass & Wasserman
  # 1995) -- the prior carries the information of a single warm-up observation
  # about each quantile: sd = sqrt(tau (1 - tau)) / f_hat(beta0_loc), with f_hat
  # a kernel density estimate of the warm-up period. Equivalently, the power
  # prior on the warm-up period with discount a0 = 1/w (Ibrahim & Chen 2000;
  # Bourazas, Kiagias & Tsiamyrtzis 2022). Override with a scalar or length-m
  # vector on the modeling scale of y.
  if (is.null(beta0_scale)) {
    yw <- if (w > 0) y_model[1:w] else y_model
    if (length(unique(yw)) < 2L) {
      stop("The warm-up period is (nearly) constant; supply beta0_scale explicitly.")
    }
    dw <- stats::density(yw)
    f_hat <- stats::approx(dw$x, dw$y, xout = beta0_loc, rule = 2)$y
    f_hat <- pmax(f_hat, 0.05 / diff(range(yw)))   # numerical floor for sparse tails
    beta0_scale <- sqrt(taus * (1 - taus)) / f_hat
  } else if (length(beta0_scale) == 1L) {
    beta0_scale <- rep(beta0_scale, m)
  }
  if (length(beta0_scale) != m || any(beta0_scale <= 0)) {
    stop("beta0_scale must be a positive scalar or a length-m positive vector.")
  }

  # ---- IQ weight matrices construction (data-driven adaptive weights) ----
  p_total <- px + 1L
  p_slope <- px
  w_iq_gamma <- matrix(1, nrow = m - 1L, ncol = max(r, 0L))
  w_iq_beta  <- matrix(1, nrow = m - 1L, ncol = max(p_slope, 0L))

  if ((r > 0 || p_slope > 0) && requireNamespace("quantreg", quietly = TRUE)) {
    X_design <- cbind(Intercept = 1, X)
    Z_pilot <- if (r > 0) cbind(X_design, H) else X_design
    d_pilot <- ncol(Z_pilot)

    pilot_coefs <- matrix(NA, nrow = m, ncol = d_pilot)
    for (q in seq_len(m)) {
      pilot_coefs[q, ] <- safe_pilot_coefs(
        y = y_model, Z = Z_pilot, tau = taus[q], eps_w = eps_w
      )
    }

    if (r > 0) {
      gamma_pilot <- pilot_coefs[, (p_total + 1):(p_total + r), drop = FALSE]
      for (q in 2:m) {
        for (j in seq_len(r)) {
          diff_val <- abs(gamma_pilot[q, j] - gamma_pilot[q - 1, j])
          w_iq_gamma[q - 1, j] <- 1 / max(diff_val, eps_w)
        }
      }
    }

    if (p_slope > 0) {
      beta_pilot <- pilot_coefs[, 2:p_total, drop = FALSE]
      for (q in 2:m) {
        for (j in seq_len(p_slope)) {
          diff_val <- abs(beta_pilot[q, j] - beta_pilot[q - 1, j])
          w_iq_beta[q - 1, j] <- 1 / max(diff_val, eps_w)
        }
      }
    }
  } else if ((r > 0 || p_slope > 0) && !requireNamespace("quantreg", quietly = TRUE)) {
    warning("Package 'quantreg' not available; using uniform IQ shrinkage weights.")
  }


  # Pilot-scale start for the EM (default since 0.6.3). Appendix C's fixed point
  # (C.9), lambda = r(m-1)^2 / sum_k w_k E|d_k|, evaluated at the pilot fit: with the
  # adaptive weights w_k = 1/|d_k^pilot| the sum equals the number of penalized
  # differences N = r(m-1) (+ p_slope(m-1)), so the start is lambda = N, i.e.
  # lambda_iq2 = N^2. This is the analogue of Park & Casella (2008, Sec. 3.1), who
  # start their marginal-likelihood EM at the fixed point evaluated with least-squares
  # estimates. Without the EM, NULL falls back to the pre-0.6.3 default of 1.
  if (is.null(lambda_iq2)) {
    N0 <- .bbq_iq_n_diff(r, p_slope, m)
    lambda_iq2 <- if (isTRUE(adaptive_iq) && N0 > 0L) as.numeric(N0)^2 else 1
  }

  stan_data <- list(
    n = n, px = px, m = m, r = r,
    X = X, H = H,
    y = y, offset = offset, tau_q = taus,
    beta0_loc = beta0_loc, beta0_scale = beta0_scale,
    base_scale = base_scale, c_sigma = c_sigma, beta_sd = beta_sd,
    lambda_nc = lambda_nc,
    lambda_iq2 = lambda_iq2,                      # SQUARED IQ fusion weight (rate = sqrt)
    iq_smooth2 = iq_smooth^2,
    lambda_beta2_a = lambda_beta2_a, lambda_beta2_b = lambda_beta2_b,
    adaptive_beta = as.integer(adaptive_beta),
    lambda_beta2_fixed = lambda_beta2_fixed,
    lambda_lasso2_a = lambda_lasso2_a, lambda_lasso2_b = lambda_lasso2_b,
    adaptive_gamma = as.integer(adaptive_gamma),
    lambda_lasso2_fixed = lambda_lasso2_fixed,
    log_flag = as.integer(log_flag), jittering = as.integer(jittering),
    prior_beta_code = prior_beta_code,
    prior_code = prior_code,
    beta_spike_sd = beta_spike_sd,
    beta_slab_sd = beta_slab_sd,
    beta_slab_pi_a = beta_slab_pi_a,
    beta_slab_pi_b = beta_slab_pi_b,
    spike_sd = spike_sd,
    slab_sd  = slab_sd,
    slab_pi_a = slab_pi_a,
    slab_pi_b = slab_pi_b,
    p_slope = as.integer(p_slope),
    w_iq_gamma = if (r > 0) w_iq_gamma else matrix(0, m - 1L, 0L),
    w_iq_beta  = if (p_slope > 0) w_iq_beta else matrix(0, m - 1L, 0L)
  )

  # Prior-center initialization for MAP optimization: start beta0 at its prior
  # location, all shift/slope coefficients at 0, and the shrinkage hierarchy at
  # neutral values. The smoothed score saturates (|z| clamped at 20) once a curve
  # sits ~20 bandwidths from the data, leaving zero gradient there; random inits
  # can start in (or reach) these flat regions and yield runaway MAP modes.
  # Sizes of the hierarchy latents follow the Stan parameter block: a latent
  # exists only for the prior that uses it (size 0 otherwise); u only if jittering.
  n_bl <- if (prior_beta_code == 6L) px else 0L
  n_gl <- if (prior_code == 5L) r else 0L
  n_lb2 <- if (px > 0L && as.integer(adaptive_beta) == 1L && prior_beta_code %in% c(2L, 4L, 5L)) 1L else 0L
  n_ll2 <- if (r > 0L && as.integer(adaptive_gamma) == 1L && prior_code %in% c(1L, 2L, 4L)) 1L else 0L
  n_pib <- if (px > 0L && prior_beta_code == 3L) 1L else 0L
  n_pi  <- if (prior_code %in% c(3L, 6L)) 1L else 0L
  n_ob <- if (prior_beta_code == 5L) px else 0L
  n_og <- if (prior_code == 4L) r else 0L
  n_u  <- if (as.integer(jittering) == 1L) n else 0L
  init_prior_center <- list(
    beta0 = as.array(beta0_loc),
    betaX = matrix(0, m, px),
    gamma = matrix(0, m, r),
    lambda2_beta_local = matrix(1, m, n_bl),
    lambda2_gamma_local = matrix(1, m, n_gl),
    lambda_beta2 = as.array(rep(1, n_lb2)), lambda_lasso2 = as.array(rep(1, n_ll2)),
    pi_slab_beta = as.array(rep(0.5, n_pib)), pi_slab = as.array(rep(0.5, n_pi)),
    omega_beta_group = as.array(rep(1, n_ob)),
    omega_group = as.array(rep(1, n_og)),
    u = as.array(rep(0.5, n_u))
  )

  # Pilot-based initialization (map_init = "pilot"): marginal LASSO quantile
  # regression per quantile level on [1 | X | H] via quantreg::rq.fit.lasso,
  # with the intercept unpenalized and the slope penalty at the
  # sqrt(tau (1 - tau) n log d) rate (cf. Belloni & Chernozhukov 2011). The
  # sparse, data-proximal start begins inside the region where the smoothed
  # score has live gradients and matches the spike-slab prior's sparsity.
  # (The adaptive IQ weights keep their own pilot fits; this block is
  # initialization only.) Mode-selection nuance: "pilot" approaches the optimum
  # from the data-fitting side, whereas "prior_center" grows the solution from
  # the in-control state; for multimodal spike-slab posteriors compare final lp.
  init_pilot <- init_prior_center
  if (map_init == "pilot") {
    pilot_init <- NULL
    if ((px + r) > 0 && requireNamespace("quantreg", quietly = TRUE)) {
      Z_init <- cbind(1, X, H)
      d_init <- ncol(Z_init)
      pilot_init <- matrix(NA_real_, m, d_init)
      for (q in seq_len(m)) {
        lam_q <- sqrt(taus[q] * (1 - taus[q]) * n * log(d_init))
        lam_vec <- c(0, rep(lam_q, d_init - 1L))  # do not penalize the intercept
        fit_l <- try(suppressWarnings(
          quantreg::rq.fit.lasso(Z_init, y_model, tau = taus[q], lambda = lam_vec)
        ), silent = TRUE)
        cf <- if (inherits(fit_l, "try-error")) NULL else as.numeric(fit_l$coefficients)
        if (!is.null(cf) && length(cf) == d_init && all(is.finite(cf))) {
          pilot_init[q, ] <- cf
        }
      }
      if (anyNA(pilot_init)) pilot_init <- NULL
    }
    if (!is.null(pilot_init)) {
      cap <- 5 * diff(range(y_model))
      clamp <- function(x, ref = 0) ref + pmin(pmax(x - ref, -cap), cap)
      init_pilot$beta0 <- as.array(clamp(pilot_init[, 1], beta0_loc))
      if (px > 0) init_pilot$betaX <-
        matrix(clamp(pilot_init[, 2:(px + 1), drop = FALSE]), m, px)
      if (r > 0) init_pilot$gamma <-
        matrix(clamp(pilot_init[, (px + 2):(px + 1 + r), drop = FALSE]), m, r)
    } else {
      warning("map_init = 'pilot' requested but rq.fit.lasso pilot fits are ",
              "unavailable; falling back to 'prior_center'.")
      map_init <- "prior_center"
    }
  }

  # Starting values actually used: the map_init choice, overridden component-wise
  # by map_init_values (e.g. the cross-validation winner's coefficients) where the
  # dimensions match. Applied to both the direct and the EM MAP paths below.
  init_used <- if (map_init == "prior_center") init_prior_center else init_pilot
  map_init_used <- map_init
  if (!is.null(map_init_values) && fit_method == "map") {
    stopifnot(is.list(map_init_values))
    taken <- character(0)
    for (nm in names(map_init_values)) {
      v <- map_init_values[[nm]]
      if (nm %in% names(init_used) && is.numeric(v) && all(is.finite(v)) &&
          identical(as.integer(dim(as.array(v))), as.integer(dim(as.array(init_used[[nm]]))))) {
        init_used[[nm]] <- if (is.matrix(init_used[[nm]])) matrix(v, nrow(init_used[[nm]]), ncol(init_used[[nm]])) else as.array(as.numeric(v))
        taken <- c(taken, nm)
      }
    }
    if (length(taken)) map_init_used <- paste0(map_init, "+values(", paste(taken, collapse = ","), ")")
    else warning("map_init_values supplied but no component matched the model dimensions; using map_init = '", map_init, "'.")
  }

  # Run rstan::optimizing with its generic "non-zero return code" warning
  # replaced by an informative termination record. Under the tight default
  # tolerances, L-BFGS commonly exits with code 70 (line search cannot further
  # improve the objective) -- the expected exit at a converged optimum for this
  # objective. Exit code 0 only means a stopping tolerance fired first; it is
  # not by itself evidence of a better fit (see map_init/tolerance docs).
  run_optimizing <- function(opt_args) {
    res <- withCallingHandlers(
      do.call(rstan::optimizing, opt_args),
      warning = function(w) {
        if (grepl("non-zero return code", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
    rc <- res$return_code
    res$termination <- if (identical(rc, 0L)) {
      "converged: a stopping tolerance was met"
    } else if (identical(rc, 70L)) {
      paste0("line search exhausted (exit code 70): no direction improves the ",
             "objective further; expected at a converged optimum")
    } else {
      sprintf("optimizer exit code %s", rc)
    }
    if (!identical(rc, 0L) && !identical(rc, 70L)) {
      warning("MAP optimization ended with unusual exit code ", rc,
              "; check fit quality via fit$map$coverage and fit$map$value.")
    }
    res
  }

  # Fit-quality diagnostic: empirical coverage of the fitted quantile curves,
  # colMeans(y < eta_tau), which should sit near tau. Guards against silently
  # bad MAP modes (curves stranded away from the data) -- loud where the
  # optimizer's exit code is uninformative.
  map_coverage_check <- function(par_vec) {
    tryCatch({
      b0 <- par_vec[sprintf("beta0[%d]", seq_len(m))]
      eta <- matrix(rep(b0, each = n), n, m)
      if (px > 0) {
        bX <- matrix(par_vec[grep("^betaX\\[", names(par_vec))], nrow = m)
        eta <- eta + X %*% t(bX)
      }
      if (r > 0) {
        gm <- matrix(par_vec[grep("^gamma\\[", names(par_vec))], nrow = m)
        eta <- eta + H %*% t(gm)
      }
      eta <- eta + offset
      cov_hat <- colMeans(y_model < eta)
      names(cov_hat) <- paste0("tau=", taus)
      tol_cov <- pmax(0.10, 4 * sqrt(taus * (1 - taus) / n))
      if (any(abs(cov_hat - taus) > tol_cov)) {
        warning("MAP fit-quality check: empirical coverage of the fitted ",
                "quantile curves deviates from the target levels (",
                paste(sprintf("%.3f vs %.3f", cov_hat, taus), collapse = "; "),
                "). The optimizer may have stopped in a poor mode; consider ",
                "map_init = 'prior_center', a larger map_iter, or tighter ",
                "tolerances.", call. = FALSE)
      }
      cov_hat
    }, error = function(e) NULL)
  }

  # Compile Stan model once per session (cached)
  if (is.null(.bbq_stan_cache$sm)) {
    if (verbose) message("Compiling BBQ Stan model (one-time per session)...")
    .bbq_stan_cache$sm <- rstan::stan_model(model_code = stan_code)
  }
  sm <- .bbq_stan_cache$sm

  # Initialize outputs
  fit <- NULL
  map_fit <- NULL
  hessian <- NULL
  laplace_samples <- NULL

  # Helper: parse 2D parameter indices from names like "gamma[1,2]"
  parse_2d_idx <- function(par_names, idx, prefix) {
    dims_str <- gsub(paste0(prefix, "\\[|\\]"), "", par_names[idx])
    dims_split <- strsplit(dims_str, ",")
    list(
      row = as.integer(sapply(dims_split, `[`, 1)),
      col = as.integer(sapply(dims_split, `[`, 2))
    )
  }

  # Helper: scatter flat sample vector into 3D array
  scatter_to_array <- function(samples_mat, row_idx, col_idx, n_row, n_col, n_samples) {
    arr <- array(NA, dim = c(n_samples, n_row, n_col))
    for (i in seq_along(row_idx)) {
      arr[, row_idx[i], col_idx[i]] <- samples_mat[, i]
    }
    arr
  }

  # Helper: generate Laplace-approximation samples from the MAP fit. The Hessian
  # from rstan::optimizing() is on the UNCONSTRAINED scale (k x k, raw parameters
  # only). The eta-relevant parameters are the intercept beta0, user slopes betaX,
  # and the change-point effects gamma; the baseline lives entirely in beta0.
  generate_laplace_samples <- function(par_map, hessian, n_samples, seed_val) {
    if (!is.null(seed_val)) set.seed(seed_val)

    par_names <- names(par_map)
    n_par <- length(par_map)

    # --- Identify raw parameter indices (first k entries of par_map) ---
    # The Hessian has k rows/cols corresponding to the raw (unconstrained) parameters.
    # par_map[1:k] = raw parameters on CONSTRAINED scale
    # par_map[(k+1):n_par] = transformed parameters
    k <- if (!is.null(hessian)) nrow(hessian) else 0
    raw_par_names <- if (k > 0) names(par_map)[1:k] else character(0)

    # Indices within raw parameters for the eta-relevant components
    beta0_idx  <- grep("^beta0\\[",  raw_par_names)
    betaX_idx  <- grep("^betaX\\[",  raw_par_names)
    gamma_idx  <- grep("^gamma\\[",  raw_par_names)
    eta_param_idx <- c(beta0_idx, betaX_idx, gamma_idx)

    # Parse dimensions from parameter names
    beta0_parsed  <- if (length(beta0_idx) > 0) {
      list(idx = as.integer(gsub("beta0\\[|\\]", "", raw_par_names[beta0_idx])))
    } else NULL
    betaX_parsed  <- if (length(betaX_idx) > 0) parse_2d_idx(raw_par_names, betaX_idx, "betaX") else NULL
    gamma_parsed  <- if (length(gamma_idx) > 0) parse_2d_idx(raw_par_names, gamma_idx, "gamma") else NULL

    # --- Try proper Hessian-based Laplace approximation ---
    laplace_ok <- FALSE
    eta_mean <- NULL; eta_cov <- NULL; eta_names <- NULL
    beta_array <- NULL
    gamma_array <- NULL

    if (!is.null(hessian) && k > 0 && length(eta_param_idx) > 0) {
      laplace_result <- tryCatch({

        # Build unconstrained mean vector. beta0/betaX/gamma are unbounded;
        # lower-bounded scale params are log-transformed below for correct inversion.
        theta_unc_full <- as.numeric(par_map[1:k])

        # lambda_lasso2, lambda_beta2, local rates, omega_* (<lower=0>): log
        # (0.6.8: no sigma2_* latents remain; the LASSO-type priors are marginal)
        for (pat in c("^lambda_lasso2\\[", "^lambda_beta2\\[", "^lambda2_beta_local\\[", "^lambda2_gamma_local\\[", "^omega_group", "^omega_beta_group")) {
          idx_tmp <- grep(pat, raw_par_names)
          if (length(idx_tmp) > 0) theta_unc_full[idx_tmp] <- log(pmax(par_map[idx_tmp], 1e-10))
        }
        # pi_slab_beta (<lower=0, upper=1>): logit
        pi_beta_idx <- grep("^pi_slab_beta\\[", raw_par_names)
        if (length(pi_beta_idx) > 0) {
          pv <- pmin(pmax(par_map[pi_beta_idx], 1e-10), 1 - 1e-10)
          theta_unc_full[pi_beta_idx] <- log(pv / (1 - pv))
        }
        # pi_slab (<lower=0, upper=1>): logit
        pi_idx <- grep("^pi_slab\\[", raw_par_names)
        if (length(pi_idx) > 0) {
          pv <- pmin(pmax(par_map[pi_idx], 1e-10), 1 - 1e-10)
          theta_unc_full[pi_idx] <- log(pv / (1 - pv))
        }
        # u (<lower=1e-12, upper=1>): logit (approx)
        u_idx <- grep("^u\\[", raw_par_names)
        if (length(u_idx) > 0) {
          uv <- pmin(pmax(par_map[u_idx], 1e-10), 1 - 1e-10)
          theta_unc_full[u_idx] <- log(uv / (1 - uv))
        }

        # Invert full Hessian to get posterior covariance on unconstrained scale.
        # Moore-Penrose pseudo-inverse is robust to the zero-curvature directions
        # of unused-prior scale latents (their flat Hessian rows otherwise make
        # solve() report a computationally singular system). eta is decoupled from
        # those latents, so its submatrix below is the exact marginal covariance.
        H_neg <- -(hessian + t(hessian)) / 2  # ensure symmetry
        H_neg_reg <- H_neg + diag(1e-6, k)
        Sigma_full <- MASS::ginv(H_neg_reg)

        # Extract marginal covariance for the eta-related parameters
        Sigma_sub <- Sigma_full[eta_param_idx, eta_param_idx]

        # Ensure positive definiteness
        eig <- eigen(Sigma_sub, symmetric = TRUE)
        eig$values <- pmax(eig$values, 1e-8)
        n_eta <- length(eig$values)
        L <- t(eig$vectors %*% diag(sqrt(eig$values), nrow = n_eta, ncol = n_eta))
        # the covariance the draws are actually generated from (floored), kept for the
        # closed-form E-step (.bbq_iq_Sbar_closed)
        Sigma_pd <- eig$vectors %*% diag(eig$values, nrow = n_eta, ncol = n_eta) %*% t(eig$vectors)

        # Sample from MVN on unconstrained scale
        theta_unc_sub <- theta_unc_full[eta_param_idx]
        z_mat <- matrix(rnorm(n_samples * length(eta_param_idx)), n_samples, length(eta_param_idx))
        samples_unc <- sweep(z_mat %*% L, 2, theta_unc_sub, "+")

        # --- Map columns back to parameter blocks (order of eta_param_idx) ---
        beta0_cols  <- seq_along(beta0_idx)
        betaX_cols  <- length(beta0_idx) + seq_along(betaX_idx)
        gamma_cols  <- length(beta0_idx) + length(betaX_idx) + seq_along(gamma_idx)

        # --- Assemble beta and gamma arrays ---
        beta_arr <- if (length(beta0_idx) > 0 || length(betaX_idx) > 0) {
          m_beta <- if (!is.null(beta0_parsed)) {
            max(beta0_parsed$idx)
          } else {
            max(betaX_parsed$row)
          }
          p_beta <- 1 + if (!is.null(betaX_parsed)) max(betaX_parsed$col) else 0
          arr <- array(0, dim = c(n_samples, m_beta, p_beta))
          if (!is.null(beta0_parsed)) {
            for (i in seq_along(beta0_parsed$idx)) {
              arr[, beta0_parsed$idx[i], 1] <- samples_unc[, beta0_cols[i]]
            }
          }
          if (!is.null(betaX_parsed)) {
            for (i in seq_along(betaX_parsed$row)) {
              arr[, betaX_parsed$row[i], betaX_parsed$col[i] + 1] <- samples_unc[, betaX_cols[i]]
            }
          }
          arr
        } else NULL

        gamma_arr <- if (length(gamma_idx) > 0) {
          scatter_to_array(samples_unc[, gamma_cols, drop = FALSE],
                           gamma_parsed$row, gamma_parsed$col,
                           max(gamma_parsed$row), max(gamma_parsed$col), n_samples)
        } else NULL

        list(beta = beta_arr, gamma = gamma_arr,
             eta_mean = theta_unc_full[eta_param_idx], eta_cov = Sigma_pd,
             eta_names = raw_par_names[eta_param_idx])

      }, error = function(e) {
        stop("Hessian-based Laplace approximation failed: ", conditionMessage(e),
             call. = FALSE)
      })

      if (!is.null(laplace_result)) {
        beta_array <- laplace_result$beta
        gamma_array <- laplace_result$gamma
        eta_mean <- laplace_result$eta_mean; eta_cov <- laplace_result$eta_cov; eta_names <- laplace_result$eta_names
        laplace_ok <- TRUE
      }
    }

    if (!laplace_ok)
      stop("Hessian-based Laplace approximation produced no draws (no Hessian or no ",
           "eta parameters); fit with map_hessian = TRUE.", call. = FALSE)
    list(mu = NULL, beta = beta_array, gamma = gamma_array,
         eta_mean = eta_mean, eta_cov = eta_cov, eta_names = eta_names)
  }

  # ------------------------------------------------------------------
  # One complete fit at a given stan_data -- in particular, at a given
  # lambda_iq2. Factored out so the EM recursion below can call it
  # repeatedly. EM sits OUTSIDE the fit: it updates lambda_iq2 between
  # refits and touches no Stan code.
  # ------------------------------------------------------------------
  run_one_fit <- function(stan_data, warm = NULL, inner_iter = NULL) {
    fit <- NULL
    map_fit <- NULL
    hessian <- NULL
    laplace_samples <- NULL

  # ------------------------------------------------------------------
  # fit_method = "mcmc": MCMC only, estimators are posterior median
  # ------------------------------------------------------------------
  if (fit_method == "mcmc") {
    fit <- rstan::sampling(
      sm, data = stan_data,
      chains = chains, iter = iter, warmup = warmup,
      control = control, seed = seed, verbose = verbose
    )

    # Extract posterior median as point estimates
    draws <- rstan::extract(fit, pars = c("beta", "gamma"))
    map_fit <- list(par = list())
    if (!is.null(draws$beta)) map_fit$par$beta <- apply(draws$beta, c(2, 3), median)
    if (!is.null(draws$gamma)) map_fit$par$gamma <- apply(draws$gamma, c(2, 3), median)
    map_fit$estimator <- "posterior_median"
    # Posterior medians of EVERY parameter, Stan-named, so that the EM monitor
    # (Eq. 17 at the point estimate) is available for this method as well.
    map_fit$par_vec <- tryCatch({
      sm_ <- rstan::summary(fit)$summary
      v <- stats::setNames(sm_[, "50%"], rownames(sm_)); v[names(v) != "lp__"]
    }, error = function(e) NULL)

  # ------------------------------------------------------------------
  # fit_method = "map_mcmc": MAP estimators, MCMC posterior draws
  # ------------------------------------------------------------------
  } else if (fit_method == "map_mcmc") {
    # First get MAP estimates
    opt_args <- list(
      object = sm,
      data = stan_data,
      hessian = map_hessian,
      as_vector = FALSE,
      seed = seed,
      verbose = verbose
    )
    if (!is.null(warm)) opt_args$init <- warm
    else if (map_init != "random") opt_args$init <- init_used
    if (!is.null(map_tol_obj))       opt_args$tol_obj       <- map_tol_obj
    if (!is.null(map_tol_grad))      opt_args$tol_grad      <- map_tol_grad
    if (!is.null(map_tol_rel_grad))  opt_args$tol_rel_grad  <- map_tol_rel_grad
    if (!is.null(map_tol_param))     opt_args$tol_param     <- map_tol_param
    if (!is.null(map_tol_rel_obj))   opt_args$tol_rel_obj   <- map_tol_rel_obj
    if (!is.null(map_history_size)) opt_args$history_size  <- map_history_size
    if (!is.null(map_iter))          opt_args$iter          <- map_iter
    if (!is.null(inner_iter))        opt_args$iter          <- as.integer(inner_iter)
    map_fit <- run_optimizing(opt_args)
    map_fit$estimator <- "map"
    map_fit$par_vec <- .bbq_par_as_vector(map_fit$par)

    hessian <- if (map_hessian && !is.null(map_fit$hessian)) map_fit$hessian else NULL

    # Then run MCMC with MAP as initialization
    init_theta <- map_fit$par
    fit <- rstan::sampling(
      sm, data = stan_data,
      chains = chains, iter = iter, warmup = warmup,
      init = function() init_theta, init_r = 0.01,
      control = control, seed = seed, verbose = verbose
    )

  # ------------------------------------------------------------------
  # fit_method = "map": MAP estimators, Laplacian approximation draws
  # ------------------------------------------------------------------
  } else if (fit_method == "map") {
    opt_args <- list(
      object = sm,
      data = stan_data,
      hessian = map_hessian,
      as_vector = TRUE,
      seed = seed,
      verbose = verbose
    )
    if (!is.null(warm)) opt_args$init <- warm
    else if (map_init != "random") opt_args$init <- init_used
    if (!is.null(map_tol_obj))       opt_args$tol_obj       <- map_tol_obj
    if (!is.null(map_tol_grad))      opt_args$tol_grad      <- map_tol_grad
    if (!is.null(map_tol_rel_grad))  opt_args$tol_rel_grad  <- map_tol_rel_grad
    if (!is.null(map_tol_param))     opt_args$tol_param     <- map_tol_param
    if (!is.null(map_tol_rel_obj))   opt_args$tol_rel_obj   <- map_tol_rel_obj
    if (!is.null(map_history_size)) opt_args$history_size  <- map_history_size
    if (!is.null(map_iter))          opt_args$iter          <- map_iter
    if (!is.null(inner_iter))        opt_args$iter          <- as.integer(inner_iter)
    map_fit <- run_optimizing(opt_args)
    map_fit$estimator <- "map"
    map_fit$par_vec <- map_fit$par
    map_fit$coverage <- map_coverage_check(map_fit$par)

    if (map_hessian && is.null(map_fit$hessian))
      stop("map_hessian = TRUE but the optimizer returned no Hessian.", call. = FALSE)
    hessian <- if (map_hessian) map_fit$hessian else NULL

    # Hessian-based Laplace draws; without the Hessian the fit is the MAP point only
    laplace_samples <- if (!is.null(hessian)) generate_laplace_samples(
      par_map = map_fit$par,
      hessian = hessian,
      n_samples = laplace_n_samples,
      seed_val = seed
    ) else NULL
  }

    list(fit = fit, map = map_fit, hessian = hessian,
         laplace_samples = laplace_samples)
  }

  # ------------------------------------------------------------------
  # lambda_iq2: EM across refits, or a single fit at the supplied value
  # ------------------------------------------------------------------
  n_iq_diff <- .bbq_iq_n_diff(r, p_slope, m)
  iq_em_step <- match.arg(iq_em_step)
  iq_em_mode <- iq_em_step
  iq_em <- list(
    adaptive        = isTRUE(adaptive_iq),
    update          = iq_em_step,
    n_diff          = n_iq_diff,
    lp_tol          = iq_em_lp_tol,
    lambda_iq2      = lambda_iq2,
    lambda_iq2_next = NA_real_,
    Sbar            = NA_real_,
    converged       = NA,
    stop_reason     = NA_character_,
    n_iter          = 0L,
    trace           = NULL,
    note            = NULL
  )

  if (isTRUE(adaptive_iq) && n_iq_diff > 0L) {
    cur <- lambda_iq2
    lam_used <- cur
    nxt <- NA_real_
    Sbar <- NA_real_
    converged <- FALSE
    tr <- NULL
    res <- NULL
    stop_reason <- "max_iter"

    for (s in seq_len(iq_em_max_iter)) {
      stan_data$lambda_iq2 <- cur
      warm <- NULL
      # One chain (author's ruling 2026-09-04): the pilot initialization of Sec. 2.4
      # starts the iterations once; every step after the first starts the optimizer at
      # the previous step's full MAP (coefficients AND every hierarchy latent; since 0.6.9
      # the global rates and mixing weights are length-0/1 vectors named with brackets).
      # Applies to "map" (named vector, as_vector = TRUE) and to
      # the MAP stage of "map_mcmc" (list, as_vector = FALSE). Zero-length blocks
      # (px = 0 or r = 0) are skipped.
      if (s > 1L && fit_method %in% c("map", "map_mcmc") && !is.null(res)) {
        pv <- res$map$par
        if (is.list(pv)) {
          keep <- intersect(names(pv), names(init_used))
          keep <- keep[vapply(keep, function(nm) length(init_used[[nm]]) > 0L, logical(1))]
          if (length(keep)) warm <- pv[keep]
        } else if (is.numeric(pv) && !is.null(names(pv))) {
          warm <- list()
          for (nm in names(init_used)) {
            tmpl <- init_used[[nm]]
            if (length(tmpl) == 0L) next
            if (nm %in% names(pv)) { warm[[nm]] <- unname(pv[[nm]]); next }   # true scalar
            idx <- which(startsWith(names(pv), paste0(nm, "[")))
            if (length(idx) == length(tmpl)) {
              v <- unname(pv[idx])
              warm[[nm]] <- if (is.matrix(tmpl)) matrix(v, nrow(tmpl), ncol(tmpl)) else as.array(v)
            }
          }
          if (!length(warm)) warm <- NULL
        }
      }
      # No per-step optimizer cap: a capped fit is not at its optimum, so Sbar from
      # it is meaningless and the chain oscillates (lp_test 2026-09-04). Every step
      # runs the optimizer to its own stopping rule.
      res <- run_one_fit(stan_data, warm = warm, inner_iter = iq_em_inner_iter)
      warm_status <- if (is.null(warm)) "start" else "moved"
      if (!is.null(warm)) {
        # Stall guard: a warm-started L-BFGS that fails its first line search returns
        # the start unchanged (the solution sits on kinks of the fused penalty).
        # Detect it by comparing the returned unbounded coefficients with the start
        # and retry from the same start with a small jitter.
        same_as <- function(r0, w0) {
          pv <- .bbq_par_as_vector(r0$map$par); if (!is.numeric(pv)) return(FALSE)
          dmax <- 0
          for (nm in intersect(names(w0), c("beta0", "betaX", "gamma"))) {
            if (length(w0[[nm]]) == 0L) next
            idx <- which(startsWith(names(pv), paste0(nm, "[")))
            if (length(idx) == length(w0[[nm]])) dmax <- max(dmax, max(abs(unname(pv[idx]) - as.vector(w0[[nm]]))))
          }
          dmax < 1e-8
        }
        if (same_as(res, warm)) {
          warm_status <- "stalled"
          for (k in seq_along(iq_em_warm_jitter)) {
            wj <- warm
            for (nm in intersect(names(wj), c("beta0", "betaX", "gamma"))) {
              if (length(wj[[nm]]) == 0L) next
              wj[[nm]] <- wj[[nm]] + stats::rnorm(length(wj[[nm]]), 0, iq_em_warm_jitter[k])
            }
            res_j <- run_one_fit(stan_data, warm = wj, inner_iter = iq_em_inner_iter)
            if (!same_as(res_j, wj)) { res <- res_j; warm_status <- paste0("jitter", k); break }
          }
          if (verbose && warm_status == "stalled")
            message(sprintf("[iq-EM %02d] warm start stalled after %d jittered retries; keeping the previous solution", s, length(iq_em_warm_jitter)))
        }
      }
      lam_used <- cur

      draws <- .bbq_iq_draws(res)
      Sbar_draws  <- if (is.null(draws)) NA_real_ else .bbq_iq_Sbar(draws, w_iq_gamma, w_iq_beta, r, p_slope, m)
      Sbar_closed <- .bbq_iq_Sbar_closed(res$laplace_samples, w_iq_gamma, w_iq_beta, r, p_slope, m)
      Sbar <- if (iq_em_estep == "closed" && is.finite(Sbar_closed)) Sbar_closed else Sbar_draws
      nxt <- .bbq_iq_em_step(cur, Sbar, n_iq_diff, step = iq_em_mode)
      rel <- if (is.finite(nxt)) abs(nxt - cur) / max(cur, .Machine$double.eps) else NA_real_

      # Monitor (author's ruling 2026-09-04): exactly Eq. (17), recomputed in R
      # from the returned coefficients at this step's lambda; and its gain over
      # the previous step of the chain. Nothing else is compared.
      lp17_obj <- tryCatch(.bbq_lp17(res$map$par_vec, stan_data), error = function(e) NULL)
      lp17 <- if (!is.null(lp17_obj) && is.finite(lp17_obj$total)) lp17_obj$total else NA_real_
      lp17_gain <- if (!is.null(tr) && is.finite(lp17) && is.finite(tr$lp17[nrow(tr)])) lp17 - tr$lp17[nrow(tr)] else NA_real_
      # Complete-data log posterior (author's ruling 2026-09-04, late): the EM's objective
      # is Eq. (17) plus the normalizing term N log lambda_iq of the fused prior
      # (Appendix C, Eq. C.1), constant within a fit, decisive across EM steps. Its gain
      # vanishes at the EM fixed point (-P + N/lambda = 0 with P = Sbar); this is the
      # stopping rule. Eq. (17) alone is kept in the trace for reference.
      lp_cd <- if (is.finite(lp17)) lp17 + n_iq_diff * log(sqrt(lam_used)) else NA_real_
      lp_cd_gain <- if (!is.null(tr) && is.finite(lp_cd) && is.finite(tr$lp_cd[nrow(tr)])) lp_cd - tr$lp_cd[nrow(tr)] else NA_real_
      tr <- rbind(tr, data.frame(
        iter = s, step = iq_em_mode, lambda_iq2 = lam_used, lambda_iq = sqrt(lam_used),
        Sbar = Sbar, Sbar_draws = Sbar_draws, Sbar_closed = Sbar_closed,
        lambda_iq2_next = nxt, rel_change = rel, warm_status = warm_status,
        lp17 = lp17, lp17_gain = lp17_gain, lp_cd = lp_cd, lp_cd_gain = lp_cd_gain
      ))
      if (verbose) {
        message(sprintf("[iq-EM %02d] lambda_iq2 = %.6g  (lambda = %.6g)  Sbar = %.6g (closed %.6g, draws %.6g)  -> %.6g  | lp_cd = %.6g  gain = %s  (lp17 = %.6g)",
                        s, lam_used, sqrt(lam_used), Sbar, Sbar_closed, Sbar_draws, nxt, lp_cd,
                        if (is.finite(lp_cd_gain)) sprintf("%+.6g", lp_cd_gain) else "-", lp17))
      }

      if (!is.finite(nxt)) {
        stop_reason <- "no_usable_update"
        iq_em$note <- sprintf(
          "EM halted at iteration %d: Sbar = %s produced no usable update; keeping lambda_iq2 = %.6g.",
          s, format(Sbar), lam_used)
        break
      }
      # The only convergence rule (author's ruling 2026-09-04/05): the complete-data log
      # posterior lp_cd has stopped gaining along the chain -- relative change below
      # iq_em_lp_tol (default 1e-2), guarded by max(1, |lp_cd|). The other exits are a
      # non-finite update (above) and iq_em_max_iter.
      if (is.finite(lp_cd_gain) && abs(lp_cd_gain) < iq_em_lp_tol * max(1, abs(tr$lp_cd[nrow(tr)] - lp_cd_gain))) {
        converged <- TRUE
        stop_reason <- "lp_gain"
        break
      }
      cur <- nxt
    }

    # The returned fit is the one actually computed at lam_used, so the model
    # and the reported lambda_iq2 are self-consistent; no extra refit is done.
    stan_data$lambda_iq2 <- lam_used
    fit             <- res$fit
    map_fit         <- res$map
    hessian         <- res$hessian
    laplace_samples <- res$laplace_samples

    iq_em$lambda_iq2      <- lam_used
    iq_em$lambda_iq2_next <- nxt
    iq_em$Sbar            <- Sbar
    iq_em$converged       <- converged
    iq_em$stop_reason     <- stop_reason
    iq_em$n_iter          <- if (is.null(tr)) 0L else nrow(tr)
    iq_em$trace           <- tr

    if (!converged && is.null(iq_em$note)) {
      iq_em$note <- sprintf(
        paste0("EM stopped at iq_em_max_iter = %d before the relative gain in the complete-data ",
               "log posterior fell below iq_em_lp_tol = %g (last gain %s). Reported lambda_iq2 is ",
               "the last fitted value."),
        iq_em_max_iter, iq_em_lp_tol,
        if (is.finite(tr$lp_cd_gain[nrow(tr)])) sprintf("%.3g", tr$lp_cd_gain[nrow(tr)]) else "NA")
      warning(iq_em$note, call. = FALSE)
    }

  } else {
    if (isTRUE(adaptive_iq) && n_iq_diff == 0L) {
      iq_em$adaptive <- FALSE
      iq_em$note <- "adaptive_iq = TRUE ignored: there are no penalized adjacent-quantile differences (need m >= 2 and r > 0 or p_slope > 0)."
    }
    res <- run_one_fit(stan_data)
    fit             <- res$fit
    map_fit         <- res$map
    hessian         <- res$hessian
    laplace_samples <- res$laplace_samples
  }

  list(
    fit = fit,
    map = map_fit,
    y = y, H = H, X = X, taus = taus,
    hessian = hessian,
    fit_method = fit_method,
    map_init_used = if (fit_method == "map") map_init_used else NA_character_,
    laplace_samples = laplace_samples,
    stan_data = stan_data,
    iq_em = iq_em
  )
}
