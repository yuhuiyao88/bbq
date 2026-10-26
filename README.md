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
library(bbq)

# synthetic AR(1) series (phi = 0.5), one year of daily data, with a location
# shift of 1 standard deviation from day 292 (0.8 T), as in Section 4
set.seed(2026)
n <- 365; e <- rnorm(n); e[292:n] <- e[292:n] + 1
y <- as.numeric(stats::filter(e, 0.5, method = "recursive"))

# end to end: cross-validation over the Section 4 grid, refit, detection
out <- bbq(y, l = 30, w = 30, prior = "spike_slab_lasso")

# or step by step at a chosen configuration
fit <- bbqFit(y, l = 30, w = 30, prior = "spike_slab_lasso", lambda_iq = 100, shrinkage = 0.1)
det <- bbqDetect(fit, l = 30, w = 30)
det$tests$qss$ui$raw                         # blocks flagged on the QSS basis
det$tests$qss$ui$cells$raw                   # location / scale cells within them
det$detected_blocks$signal_obs[det$tests$qss$ui$raw]   # localized change-points
```

The shrinkage hyperparameter is the spike standard deviation for
`"spike_slab"` and `"spike_slab_lasso"`, and the gamma hyperprior rate
`b_lambda` for `"lasso"`, `"group_lasso"`, `"het_group_lasso"` and
`"adaptive_lasso"`. `lambda_iq` is the interquantile weight of Eq. (15). All
functions take a `seed` (default 1), so repeated runs give identical results.

## Simulation study

`simulation/ar_ext12` holds the Section 4 study (11 settings × no shift,
location and scale; 33 configurations) with its own copy of the package. It is
not part of the package build; see `simulation/ar_ext12/README.md`.
