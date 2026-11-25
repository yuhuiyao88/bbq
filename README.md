# bbq

Bayesian bi-aspect quantile (BBQ) change-point detection for retrospective
(Phase I) monitoring of individual autocorrelated observations, as described in
the manuscript *A Bayesian Bi-aspect Quantile Change-Point Detection Method for
Retrospective Location and Scale in Individual Autocorrelated Data*.

The three quantile levels (0.25, 0.50, 0.75) are modeled jointly by Bayesian
multiple quantile regression with a block shift design, a smoothed score
likelihood, non-crossing and interquantile penalties, and a sparsity-inducing
prior on the shift coefficients. The posterior is summarized by its mode and a
Hessian-based Laplace approximation. Decorrelated union-intersection statistics
then flag shifted blocks on the quantile (Q) basis and on the location (median)
and scale (IQR) basis (QSS), each change is localized within its block, and the
cells responsible are identified.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yuhuiyao88/bbq")
```

The package compiles its Stan model on first use in each R session (a few minutes).

## The manuscript workflow

| Manuscript | Function |
|---|---|
| Sec. 2.1, design shift matrix | `bbqDesign()` |
| Sec. 2.5, order-preserved cross-validation (Eq. 18) over the Sec. 4 grid | `bbqTune()` |
| Secs. 2.2-2.4, MAP fit with Laplace draws (Eqs. 5-17) | `bbqFit()` |
| Sec. 3, decorrelated UI detection and localization (Eqs. 20-25) | `bbqDetect()` |
| Steps 1-6 end to end | `bbq()` |
| Figures 6-8 | `plotQuantileProcess()`, `plotQSSProcess()`, `plotGammaHeatmap()` |

```r
library(bbq)   # the Stan model compiles on the first fit of each session

# Synthetic AR(1) series as in the reference setting P05 of Section 4.1:
# phi = 0.5, T = 365, l = 30, w = 30, and a location shift of 1 in the
# innovations from observation 292 (0.8 T) onward.
set.seed(2026)
n <- 365; e <- rnorm(n); e[292:n] <- e[292:n] + 1
y <- as.numeric(stats::filter(e, 0.5, method = "recursive"))

# End to end (Steps 1-6): cross-validation over the Section 4 grid, refit of
# the selected configuration on all T observations, and detection.
out <- bbq(y, l = 30, w = 30, prior = "spike_slab_lasso")
out$settings[c("lambda_iq", "shrinkage")]   # selected by Eq. (18): lambda_iq = 0, spike sd = 0.1
det <- out$detection
det$tests$qss$ui$raw                        # OOC blocks on the QSS basis (Eq. 21): 9 10
det$tests$qss$ui$cells$raw[, det$tests$qss$ui$raw, drop = FALSE]   # location / scale cells
det$detected_blocks$signal_obs[det$tests$qss$ui$raw]               # localized change-points (Eq. 25): 277 302
det$tests$qss$overall_ui_p                  # overall UI p-value (Eq. 22)

# Figures 6-8 for the QSS basis
plotQuantileProcess(out$fit, detection = det, basis = "qss")
plotQSSProcess(out$fit, detection = det)
plotGammaHeatmap(out$fit, detection = det, basis = "qss")
```

The same steps one at a time; `bbqTune` sorts the grid by validation loss, so
row 1 is the selected configuration, and the results are identical to `bbq()`:

```r
cv  <- bbqTune(y, l = 30, w = 30, prior = "spike_slab_lasso")
fit <- bbqFit(y, l = 30, w = 30, prior = "spike_slab_lasso",
              lambda_iq = cv$lambda_iq[1], shrinkage = cv[[attr(cv, "shrinkage_arg")]][1])
det <- bbqDetect(fit, l = 30, w = 30)
```

Blocks 9 and 10 cover observations 271-300 and 301-330, so the shift at
observation 292 is found in the block that contains it and in the next one.
The end-to-end run takes about five minutes on a laptop, most of it the 20
cross-validation configurations.

The shrinkage hyperparameter is the spike standard deviation for
`"spike_slab"` and `"spike_slab_lasso"`, and the gamma hyperprior rate
`b_lambda` for `"lasso"`, `"group_lasso"`, `"het_group_lasso"` and
`"adaptive_lasso"`. `lambda_iq` is the interquantile weight of Eq. (15). All
functions take a `seed` (default 1), so repeated runs give identical results.
