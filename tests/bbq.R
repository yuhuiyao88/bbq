# Checks that the BBQ wrappers are the existing functions at the manuscript settings.
# Run from the package root. Uses the installed package selected by R_LIBS.
suppressPackageStartupMessages(library(bbq))
ok <- function(label, cond) {
  if (!isTRUE(cond)) stop("FAIL: ", label)
  cat("PASS:", label, "\n")
}

ok("sustained design is getSustainedShift",
   identical(bbqDesign(365, 30, 30), getSustainedShift(365, 30, 30)))
ok("isolated design is getIsolatedShift",
   identical(bbqDesign(365, 30, 30, "isolated"), getIsolatedShift(365, 30, 30)))
ok("r = ceiling((T - w) / l)", ncol(bbqDesign(1826, 30, 90)) == 58L)

# Synthetic posterior draws with a location shift in block 3 and a scale shift in block 6.
set.seed(7)
taus <- c(0.25, 0.5, 0.75); n <- 365; l <- 30; w <- 30
H <- bbqDesign(n, l, w); r <- ncol(H); draws <- 1000L
gamma <- array(rnorm(draws * 3 * r, sd = 0.1), c(draws, 3, r))
gamma[, , 3] <- gamma[, , 3] + 1
gamma[, 1, 6] <- gamma[, 1, 6] - 0.8; gamma[, 3, 6] <- gamma[, 3, 6] + 0.8
beta <- matrix(qnorm(taus), 3, 1)
fit <- list(fit_method = "map", taus = taus, y = rnorm(n), H = H, X = NULL,
            map = list(par = list(beta = beta, gamma = apply(gamma, c(2, 3), mean))),
            laplace_samples = list(gamma = gamma))

det <- bbqDetect(fit, l, w)
ref <- detectChangepoints_gamma(fit, taus, l = l, w = w, y = fit$y,
                                signal_position = "score", alpha = 0.05,
                                basis = c("quantile", "qss"), statistic = "ui",
                                adjust = "raw")
ok("bbqDetect is detectChangepoints_gamma at the manuscript settings", identical(det, ref))
ok("decision rule is the unadjusted block rule", identical(det$adjust, "raw"))
for (b in c("quantile", "qss")) {
  ui <- det$tests[[b]]$ui; m <- nrow(det$tests[[b]]$z)
  ok(paste(b, "block constant is Eq. (21)"),
     isTRUE(all.equal(ui$c_block, qnorm((1 + 0.95^(1 / m)) / 2))))
  ok(paste(b, "block probability is Eq. (21)"),
     isTRUE(all.equal(ui$pvalue, 1 - (2 * pnorm(ui$stat) - 1)^m)))
  ok(paste(b, "OOC blocks are p_j < FAP0"), identical(ui$raw, which(ui$pvalue < 0.05)))
}
ok("QSS basis has location and scale cells", identical(rownames(det$tests$qss$z), c("LS", "ScS")))
ok("location shift block flagged with its location cell",
   3L %in% det$tests$qss$ui$raw && det$tests$qss$ui$cells$raw["LS", 3])
ok("scale shift block flagged with its scale cell",
   6L %in% det$tests$qss$ui$raw && det$tests$qss$ui$cells$raw["ScS", 6])

# Table 2: QSS profiles are the median and the IQR of the fitted quantiles.
eta <- array(NA_real_, c(draws, 3, n))
for (s in seq_len(draws)) eta[s, , ] <- beta[, 1] + gamma[s, , ] %*% t(H)
qss <- getQSS(eta)
ok("QSS location is the fitted median", isTRUE(all.equal(qss[, 1, ], eta[, 2, ])))
ok("QSS scale is the fitted IQR", isTRUE(all.equal(qss[, 2, ], eta[, 3, ] - eta[, 1, ])))

# Hessian-based Laplace draws only: no heuristic substitute.
ok("getLaplaceSamples requires a Hessian",
   inherits(tryCatch(getLaplaceSamples(list(par = c(a = 1)), NULL), error = identity), "error"))
ok("no noise-scale argument remains",
   !any(c("noise_scale", "laplace_noise_scale") %in%
          c(names(formals(getLaplaceSamples)), names(formals(getModel)))))
cat("BBQ WRAPPER CHECKS PASSED\n")
