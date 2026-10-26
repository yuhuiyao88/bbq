# ar_ext12 — Section 4 simulation for the BBQ manuscript, on the bbq package

One study with every setting in Section 4: 11 settings × {no shift, location,
scale} = 33 configurations, 100 usable replications each (120 simulation IDs
available per configuration for replacements). It is the `ar_ext10` code with

- the package swapped from bqq 0.6.12 to the bundled **bbq 0.1.2**
  (`package/bbq_0.1.2.tar.gz`, no ARCOS data inside);
- the HK setting added from `ar_ext11` (Pearson Type VII innovations with
  moments (0, 1, 0, 4), calibrated location and scale arms at TV 0.383);
- the skewness and kurtosis arms dropped.

Everything else is reused unchanged: data generation, seeds, six priors, CV
grids, both streams (raw series and ARMA residuals), cpm benchmarks, the
supervisor with resume from checkpoints, batch initialization, retries,
timeouts, completion ladder and aggregation. `verify_config.R` checks every one
of the 33 cells against `../ar_ext10` or `../ar_ext11`.

| IDs | Setting (label in output) | Manuscript |
|---|---|---|
| 1–3 | P00 | P00, independent |
| 4–6 | D0 | P05, reference (AR 0.5) |
| 7–9 | P02 | P02 |
| 10–12 | P08 | P08 |
| 13–15 | K | RS, right-skewed |
| 16–18 | T18 | T18 |
| 19–21 | T09 | T09 |
| 22–24 | L14 | L14 |
| 25–27 | L07 | L07 |
| 28–30 | ISO | ISO |
| 31–33 | HK | HK |

Within each setting the IDs are no shift, location, scale. The manuscript methods
are `q_ui_raw_w` and `qss_ui_raw_w` (unadjusted UI on the Q and QSS bases); the
other methods are computed as before. What changes results relative to
ar_ext10/ar_ext11 comes from the package: the bandwidth rule
s = min(RMSE₀, IQR₀/1.38898) of Eq. (12), and no heuristic posterior draws.

## Run on the workstation (fresh R session)

```r
study <- "<path>/simulation_study/ar_ext12"
source(file.path(study, "install_arext12_package.R"))   # private library ~/arext12_library
source(file.path(study, "verify_config.R"))              # static checks; needs ../ar_ext10, ../ar_ext11
source(file.path(study, "smoke_driver_arext12.R"))      # one replication, one worker (T09 null)
Sys.setenv(AREXT12_CORES = "80")                         # match the workstation
source(file.path(study, "run_arext12.R"))
```

Overrides: `AREXT12_R_LIBRARY` (package library), `AREXT12_WORK_DIR` (hot output,
default `~/arext12_run`), `AREXT12_CORES`; `SMOKE_CONFIG=<id>` picks another smoke
cell. After loading the supervisor, `arext12_status(target = 100)` reports progress
and `run_arext12_kill()` stops the study's own processes. Aggregate with
`build_app_csv_from_tasks.R`.

Workstation packages: rstan (with a C++ toolchain), quantreg, MASS, cpm,
PearsonDS, forecast; ggplot2 and patchwork for figures.
