#' Cross-Validation for Interquantile Shrinkage Model (getModel)
#'
#' COPSS-style order-preserved 2-fold CV for tuning hyperparameters
#' in the interquantile shrinkage quantile regression model. Two validation
#' losses are available: the exact unconstrained score-likelihood criterion
#' (default) and the average pinball (check) loss.
#'
#' @name cv_copss
NULL

#' Pinball (check) loss for quantile regression
#'
#' Computes the average pinball loss across observations and quantile levels.
#'
#' @param y_val Numeric vector of validation responses.
#' @param qhat Numeric matrix of predicted quantiles (n x m).
#' @param taus Numeric vector of quantile levels.
#' @return Scalar average pinball loss.
#' @keywords internal
pinball_loss <- function(y_val, qhat, taus) {
  n <- length(y_val)
  m <- length(taus)
  losses <- matrix(0, n, m)
  for (j in seq_len(m)) {
    u <- y_val - qhat[, j]
    losses[, j] <- u * (taus[j] - as.numeric(u < 0))
  }
  mean(losses)
}


#' Exact unconstrained score-likelihood loss for quantile regression
#'
#' Held-out version of the unrestricted score likelihood (manuscript Eq. (6)):
#' the quadratic form \eqn{vec(S)' (Q \otimes G)^{-1} vec(S) / (2 n)}, where S stacks the
#' EXACT (unsmoothed) quantile scores psi_qt = tau_q - 1\{y_t < qhat_qt\} along
#' the design directions, Q[a,b] = min(tau_a, tau_b) - tau_a tau_b is the
#' quantile kernel (the covariance of the exact score, which is why no smoothing
#' is applied here), and G = Z'Z/n is the Gram of the evaluation fold's design.
#' Equals the negative held-out score log-likelihood up to the model's own 1/(2n)
#' scaling, so lower is better. The non-crossing indicator and all penalties are
#' deliberately excluded: validation evaluates the data-fit term only.
#'
#' Mirrors the Stan likelihood block of \code{getModel} (same Cholesky solves,
#' same 1e-8 Gram ridge), with the smoothed sigmoid replaced by the exact
#' indicator because no gradients are needed at validation.
#'
#' @param y_val Numeric vector of validation responses.
#' @param qhat Numeric matrix of predicted quantiles (n x m).
#' @param taus Numeric vector of quantile levels.
#' @param Z_val Design matrix of the evaluation fold (n x (p + r)), the same
#'   \code{[1 | X | H]} layout used in the Stan score.
#' @return Scalar loss (negative held-out score log-likelihood), or NA if the
#'   whitening factorization fails.
#' @keywords internal
score_loss <- function(y_val, qhat, taus, Z_val) {
  n <- length(y_val)
  m <- length(taus)
  psi <- matrix(0, n, m)
  for (j in seq_len(m)) {
    psi[, j] <- taus[j] - as.numeric(y_val < qhat[, j])
  }
  tryCatch({
    S <- crossprod(Z_val, psi)                       # (p+r) x m score matrix
    G <- crossprod(Z_val) / n
    diag(G) <- diag(G) + 1e-8                        # same tiny ridge as Stan
    Qk <- outer(taus, taus, pmin) - tcrossprod(taus)
    A <- forwardsolve(t(chol(G)), S)                 # L_G^{-1} S
    B <- forwardsolve(t(chol(Qk)), t(A))             # L_Q^{-1} (L_G^{-1} S)'
    0.5 * sum(B * B) / n                             # tr(S' G^-1 S Q^-1) / (2n)
  }, error = function(e) {
    warning("score_loss failed: ", e$message)
    NA_real_
  })
}


# Both losses are computed for every fit (the fits dominate the cost); `loss`
# only decides which one fills train_loss/val_loss and ranks the grid.
.bbq_cv_losses <- function(y_tr, eta_tr, y_val, eta_val, taus, Z_tr, Z_val) {
  list(
    train_pinball = pinball_loss(y_tr, eta_tr, taus),
    val_pinball   = pinball_loss(y_val, eta_val, taus),
    train_score   = score_loss(y_tr, eta_tr, taus, Z_tr),
    val_score     = score_loss(y_val, eta_val, taus, Z_val)
  )
}

.bbq_cv_losses_na <- function() {
  list(train_pinball = NA_real_, val_pinball = NA_real_,
       train_score = NA_real_, val_score = NA_real_)
}


# Tuning grids and base_args are matched to getModel() formals BY NAME, and any
# name that is not a formal is silently dropped when the argument list is built.
# That silence is dangerous after a rename: a grid still carrying the old
# `lambda_iq` column would tune nothing while appearing to run normally. Fail
# loudly instead, and point at the rename when that is what happened.
.bbq_check_tuning_names <- function(grid_names, base_names, formal_names) {
  unknown <- setdiff(c(grid_names, base_names), formal_names)
  unknown <- unknown[nzchar(unknown)]
  if (!length(unknown)) return(invisible(TRUE))
  hint <- ""
  if ("lambda_iq" %in% unknown) {
    hint <- paste0(
      "\n  `lambda_iq` was renamed to `lambda_iq2` and is now the SQUARED weight: ",
      "the effective fusion rate is sqrt(lambda_iq2).\n  Replace lambda_iq = v with ",
      "lambda_iq2 = v^2 to reproduce previous results.")
  }
  stop("Unknown tuning argument(s) for getModel(): ",
       paste(sQuote(unknown), collapse = ", "), hint, call. = FALSE)
}



#' General grid search CV for getModel
#'
#' More flexible version that accepts a data.frame grid of hyperparameters.
#'
#' @param y Numeric vector of responses.
#' @param taus Numeric vector of quantile levels.
#' @param H,X Design matrices.
#' @param w Integer warm-up period.
#' @param grid data.frame with columns for hyperparameters (lambda_nc, etc.)
#' @param base_args Named list of additional arguments passed to getModel.
#' @param loss Validation loss used to rank the grid: \code{"score"} (default)
#'   is the exact unconstrained score-likelihood criterion of \code{score_loss};
#'   \code{"pinball"} restores the previous average pinball (check) loss. Both
#'   are always computed and returned; \code{loss} decides which one fills
#'   \code{train_loss}/\code{val_loss} and sorts the result.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return data.frame with grid and CV losses (\code{train_loss}, \code{val_loss},
#'   \code{train_score}, \code{val_score}, \code{train_pinball}, \code{val_pinball}), plus \code{lambda_iq2_fit}, the mean over folds of the
#'   \eqn{\lambda_{iq}^2} each row's EM settled on (NA when the EM is off), which a
#'   caller can pass as the start of the final fit.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 200
#' y <- rnorm(n)
#' taus <- c(0.25, 0.5, 0.75)
#' H <- getSustainedShift(n, l = 20, w = 30)
#' grid <- data.frame(lambda_nc = c(2, 5, 10))
#' cv_result <- cv_copss_grid(y = y, taus = taus, H = H, w = 30, grid = grid)
#' }
#'
#' @export
cv_copss_grid <- function(y, taus, H, X = NULL, w, grid,
                              base_args = list(),
                              loss = c("score", "pinball"),
                              seed = 123,
                              verbose = TRUE) {

  loss <- match.arg(loss)

  # Get getModel formals for default filling
  gm_formals <- as.list(formals(getModel))
  .bbq_check_tuning_names(names(grid), names(base_args), names(gm_formals))

  fit_and_score <- function(idx_train, idx_val, row_args, init_values = NULL) {
    # Convert factors to character
    row_args <- lapply(row_args, function(x) if (is.factor(x)) as.character(x) else x)
    base_args_l <- lapply(base_args, function(x) if (is.factor(x)) as.character(x) else x)

    H_tr <- H[idx_train, , drop = FALSE]
    H_val <- H[idx_val, , drop = FALSE]

    # Build full argument list
    full_args <- gm_formals
    full_args["y"] <- list(y[idx_train])
    full_args["taus"] <- list(taus)
    full_args["H"] <- list(H_tr)
    X_tr <- if (is.null(X)) NULL else as.matrix(X[idx_train, , drop = FALSE])
    full_args["X"] <- list(X_tr)
    full_args["w"] <- list(w)
    full_args["fit_method"] <- list("map")
    full_args["seed"] <- list(seed)
    full_args["verbose"] <- list(FALSE)
    # The tuning fits use the SAME estimation procedure as the final fit: adaptive_iq
    # follows getModel()'s default (FALSE since 0.6.11), so lambda_iq2 is a grid column
    # and each tuning fit is one MAP fit at its row's values. Passing adaptive_iq = TRUE
    # in base_args runs the EM inside every fold instead.
    full_args["laplace_n_samples"] <- list(2000L)   # E-step draws for tuning; base_args may override
    # The tuning fits call getModel with getModel's own defaults for everything not in
    # base_args or the grid: the same chain, to the same stopping rule, as the final fit.

    # Override with base_args then row_args
    for (nm in names(base_args_l)) full_args[nm] <- list(base_args_l[[nm]])
    for (nm in names(row_args)) full_args[nm] <- list(row_args[[nm]])
    if (!is.null(init_values) && is.null(full_args[["map_init_values"]]))
      full_args["map_init_values"] <- list(init_values)
    # The EM's E-step averages |d| over Laplace draws, so the Hessian is required
    # whenever the EM is on; without the EM the tuning fit needs no draws.
    full_args["map_hessian"] <- list(isTRUE(full_args[["adaptive_iq"]]))

    # Keep only valid arguments
    full_args <- full_args[intersect(names(full_args), names(gm_formals))]

    fit <- tryCatch({
      do.call(getModel, full_args)
    }, error = function(e) {
      # A wall-clock / CPU time-limit hit (from setTimeLimit in a caller) must
      # abort the whole CV rather than be silently downgraded to a skipped grid
      # point: the caller treats a timeout as an invalid run and redoes the
      # simulation from fresh data. Genuine numerical fit failures are still
      # tolerated (skip this grid point and continue).
      if (grepl("reached elapsed time limit|reached CPU time limit",
                conditionMessage(e))) stop(e)
      warning("Model fitting failed: ", e$message)
      return(NULL)
    })

    if (is.null(fit) || is.null(fit$map)) {
      return(.bbq_cv_losses_na())
    }

    par <- fit$map$par

    # Reconstruct design
    X_tr_design <- {
      X0 <- matrix(1, nrow = length(idx_train), ncol = 1)
      X_raw <- if (is.null(X)) NULL else as.matrix(X[idx_train, , drop = FALSE])
      out <- if (is.null(X_raw) || ncol(X_raw) == 0) X0 else cbind(X0, X_raw)
      storage.mode(out) <- "double"
      out
    }

    m <- length(taus)
    p <- ncol(X_tr_design)
    r <- ncol(H_tr)

    beta_vec <- par[grep("^beta\\[", names(par))]
    gamma_vec <- par[grep("^gamma\\[", names(par))]

    if (length(beta_vec) == 0) beta_vec <- rep(0, m * p)
    if (length(gamma_vec) == 0) gamma_vec <- rep(0, m * r)

    beta <- matrix(beta_vec, m, p, byrow = FALSE)
    gamma <- if (r > 0) matrix(gamma_vec, m, r, byrow = FALSE) else matrix(0, m, 0)

    # Training predictions: X*beta[q,] (incl. intercept) + H*gamma[q,]
    eta_tr <- matrix(0, nrow = length(idx_train), ncol = m)
    for (j in seq_len(m)) {
      eta_tr[, j] <- as.numeric(X_tr_design %*% beta[j, ]) +
        if (r > 0) as.numeric(H_tr %*% gamma[j, ]) else 0
    }

    # Validation predictions: X*beta[q,] (incl. intercept) + H*gamma[q,]
    X_val_raw <- if (is.null(X)) NULL else as.matrix(X[idx_val, , drop = FALSE])
    X_val <- {
      X0 <- matrix(1, nrow = length(idx_val), ncol = 1)
      out <- if (is.null(X_val_raw) || ncol(X_val_raw) == 0) X0 else cbind(X0, X_val_raw)
      storage.mode(out) <- "double"
      out
    }

    n_val <- length(idx_val)

    eta_val <- matrix(0, n_val, m)
    for (j in seq_len(m)) {
      eta_val[, j] <- as.numeric(X_val %*% beta[j, ]) +
        if (r > 0) as.numeric(H_val %*% gamma[j, ]) else 0
    }

    res <- .bbq_cv_losses(y[idx_train], eta_tr, y[idx_val], eta_val, taus,
                          Z_tr = cbind(X_tr_design, H_tr),
                          Z_val = cbind(X_val, H_val))
    res$lambda_iq2 <- if (!is.null(fit$iq_em) && isTRUE(fit$iq_em$adaptive))
      fit$iq_em$lambda_iq2 else NA_real_
    res$map_par <- par                      # fold MAP parameters, for cv_winner_init()
    res
  }

  idx_odd  <- seq(1, length(y), by = 2)
  idx_even <- seq(2, length(y), by = 2)

  tr_sel <- paste0("train_", loss)
  val_sel <- paste0("val_", loss)

  train_pinballs <- val_pinballs <- numeric(nrow(grid))
  train_scores <- val_scores <- numeric(nrow(grid))
  lambda_iq2_fit <- rep(NA_real_, nrow(grid))
  map_par <- vector("list", nrow(grid))

  report <- function(k, row_args, a, b) {
    if (!verbose) return(invisible(NULL))
    hp_str <- paste(
      sprintf("%s=%s", names(row_args), vapply(row_args, function(x) format(x, digits = 3), "")),
      collapse = ", "
    )
    message(sprintf("[cv_copss_grid] iter %d/%d | %s | %s: train=%.4f | val=%.4f",
                    k, nrow(grid), hp_str, loss,
                    (a[[tr_sel]] + b[[tr_sel]]) / 2, (a[[val_sel]] + b[[val_sel]]) / 2))
  }
  for (k in seq_len(nrow(grid))) {
    row_args <- as.list(grid[k, , drop = FALSE])
    # every grid row starts from the pilot and runs its own chain
    a <- fit_and_score(idx_odd, idx_even, row_args)
    b <- fit_and_score(idx_even, idx_odd, row_args)
    train_pinballs[k] <- (a$train_pinball + b$train_pinball) / 2
    val_pinballs[k]   <- (a$val_pinball + b$val_pinball) / 2
    train_scores[k]   <- (a$train_score + b$train_score) / 2
    val_scores[k]     <- (a$val_score + b$val_score) / 2
    lambda_iq2_fit[k] <- mean(c(a$lambda_iq2, b$lambda_iq2), na.rm = TRUE)
    map_par[[k]] <- list(a$map_par, b$map_par)
    report(k, row_args, a, b)
  }

  out <- cbind(grid,
               train_loss = if (loss == "score") train_scores else train_pinballs,
               val_loss   = if (loss == "score") val_scores else val_pinballs,
               train_pinball = train_pinballs, val_pinball = val_pinballs,
               train_score = train_scores, val_score = val_scores,
               lambda_iq2_fit = lambda_iq2_fit,
               cv_loss = loss)
  ord <- order(out$val_loss)
  out <- out[ord, ]
  rownames(out) <- NULL
  attr(out, "map_par") <- map_par[ord]   # fold MAP parameters per row, in the sorted order

  out
}

#' Starting values for the final fit from a cross-validation winner
#'
#' Averages the two fold MAP fits of grid row \code{k} of a \code{cv_copss_grid}
#' result (row 1 is the winner) and returns them as a list \code{beta0},
#' \code{betaX}, \code{gamma} for \code{getModel(map_init_values = )}. The folds
#' use the same design as the full series, so the averaged coefficients are on the
#' scale of the final fit and start the optimizer close to its solution.
#'
#' @param cv Result of \code{cv_copss_grid}.
#' @param k Grid row (in the sorted result); default 1, the winner.
#' @return A named list of starting values, or \code{NULL} when the fold fits were
#'   not recorded.
#' @noRd
cv_winner_init <- function(cv, k = 1L) {
  mp <- attr(cv, "map_par")
  if (is.null(mp) || length(mp) < k || is.null(mp[[k]])) return(NULL)
  pars <- Filter(Negate(is.null), mp[[k]])
  if (!length(pars)) return(NULL)
  nm <- Reduce(intersect, lapply(pars, names))
  avg <- Reduce(`+`, lapply(pars, function(p) p[nm])) / length(pars)
  .bbq_par_to_init(avg)
}

# Named MAP parameter vector ("beta0[1]", "gamma[2,3]", ...) -> list(beta0, betaX, gamma)
.bbq_par_to_init <- function(avg) {
  nm <- names(avg)
  pick <- function(prefix) {
    v <- avg[grep(paste0("^", prefix, "\\["), nm)]
    if (!length(v)) return(NULL)
    idx <- do.call(rbind, lapply(regmatches(names(v), regexpr("\\[.*\\]", names(v))),
                                 function(z) as.integer(strsplit(gsub("[][]", "", z), ",")[[1]])))
    if (ncol(idx) == 1L) return(as.array(as.numeric(v[order(idx[, 1])])))
    M <- matrix(0, max(idx[, 1]), max(idx[, 2]))
    M[idx] <- v
    M
  }
  out <- list(beta0 = pick("beta0"), betaX = pick("betaX"), gamma = pick("gamma"))
  Filter(Negate(is.null), out)
}

