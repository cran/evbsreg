# evbsreg 1.2.3

* `envelope_qq()` did not restore the graphical parameter `pty` on exit: it
  saved `par("pty")` without its name, so `par(oldpar)` restored nothing and
  warned '"m" is not a graphical parameter'. It now uses
  `oldpar <- par(pty = "s"); on.exit(par(oldpar))`. A test checks that `pty`
  is restored and that no warning is emitted.
* `envelope_qq()` failed for `nrep` below 40 (an order-statistic index of 0).
  The indices are now kept valid; the bands are unchanged for `nrep >= 40`.

# evbsreg 1.2.2

* Removed the dependency on 'SpatialExtremes', which was archived on CRAN on
  2026-10-21 (and with it evbsreg). The package used only
  `SpatialExtremes::rgev()`, inside `revbs()`; it is replaced by an internal
  generator built on the same inversion of an Exp(1) draw, so `revbs()` and
  `evbsreg.fit.mc()` return identical values for a given seed.
* Examples and vignette are fully reproducible: the only random steps
  (`envelope_qq()` envelope, `evbs_block_boot()` resampling) now run after
  `set.seed()`, and numbers quoted in the vignette text are computed inline.
* The vignette now also covers the upper endpoint, return levels and expected
  shortfall, block bootstrap standard errors, `evbs_monitor()`, and the
  comparison with GEV regression using the common statistic `Bbar`.
* `DESCRIPTION`, `README.md` and `inst/CITATION` updated to the current title
  of the accompanying paper, with the arXiv preprint; the version history was
  removed from the `Description` field.
* pkgdown site: favicons regenerated in the format pkgdown >= 2.1 links to
  (`favicon.svg`, `favicon-96x96.png`, `site.webmanifest`), which were missing.

# evbsreg 1.2.1

## Bug fixes and documentation

* `devbs()` documentation stated that the extreme-value index of the response
  is always 2*gamma. That holds only for gamma > 0; for gamma <= 0 the index
  is gamma itself (Ferreira, Gomes and Leiva, 2012, Theorems 4.1-4.2).
* `evbs_return_level()` documentation now states that the level is
  conditional on the covariate value, that `period` is in units of one
  observation (months, for monthly maxima), and what reading it as an annual
  return level assumes.
* `evbs_monitor()` documentation now states that the endpoint is evaluated at
  the covariate value of the running maximum, as the code does.
* `inst/scripts/`: the simulation scripts used gamma in {-0.2, -0.1, 0.1, 0.2}
  whereas the paper's tables use {-0.2, 0, 0.2}; Scenario 3 used alpha = 0.5,
  which with the contamination alpha - 0.5 gave a degenerate alpha = 0. Both
  now match the paper (Scenario 3: beta = (1, 2), alpha = 1.2).

## New features

* `cnc_diagnostics_gev()` also returns `Bbar`, the mean aggregate contribution
  computed exactly as `colMeans(cnc_diagnostics(fit)$Bj)` is for the EVBS
  model, so that influence rankings can be compared across the two families.
  The existing `Bj` (normalised diagonal of the influence matrix) is unchanged.
* `gev_monitor()` computes the prospective identifiability index of
  `evbs_monitor()` under the GEV regression.

# evbsreg 1.2.0

## Bug fixes

* `plot_aggregate_contributions()` labelled flagged observations with a fixed
  offset above each point, so labels for observations close together on the
  index axis overlapped and became unreadable. Labels are now staggered
  vertically for points that fall within a closeness threshold of each other,
  leaving isolated points unaffected.
* `envelope_qq()` saved and restored the caller's *entire* `par()` state
  (`par(no.readonly = TRUE)`), which includes the multi-panel layout position
  (`mfg`). When called as one step of a multi-panel figure set up by the
  caller (`par(mfrow = c(1, 2))` followed by `envelope_qq()` then a second
  `plot()`), restoring the saved state on exit rewound the layout position,
  so the second plot overwrote the first panel instead of advancing to the
  second, leaving one panel blank. The function now saves and restores only
  `pty`, the one parameter it actually changes.

## New features

* `evbs_monitor()` computes the prospective endpoint-identifiability index
  $F_t$, combining the conformal normal curvature with the finite upper endpoint,
  and `plot.evbs_monitor()` draws the corresponding control chart. Values below
  one indicate that the endpoint is determined by a single observation and should
  not be quoted as a design value.

# evbsreg 1.1.0

## New features

* `logEVBS()` provides a `gamlss.family` implementation of the log-EVBS
  distribution, allowing every parameter, including the tail-shape parameter,
  to depend on covariates or smooth terms. With a constant predictor it
  reproduces the fixed-parameter fit of `evbsreg.fit()`. Companion functions
  `dlogEVBS()`, `plogEVBS()` and `qlogEVBS()` are also provided.
* `devbs()`, `pevbs()` and `qevbs()` complete the distribution family. The
  package previously exported only `revbs()`, so users had no way to compute
  quantiles or return levels.
* `evbs_endpoint()` returns the finite upper endpoint of the fitted model when
  `gama < 0` (Weibull max-domain of attraction).
* `evbs_return_level()` computes return levels and expected shortfall from the
  exact EVBS quantile function.
* `gevreg.fit()`, `gev_scores()` and `cnc_diagnostics_gev()` extend the local
  influence framework to the generalized extreme-value regression model.
* `evbs_block_boot()` provides moving-block bootstrap standard errors for
  series with residual serial dependence.

## Important notes

* The extreme-value index of the EVBS *response* is `2 * gama`, not `gama`. The
  transformation preserves the max-domain of attraction but doubles the tail
  index. Documentation updated accordingly.
* The GEV return-level formula does **not** apply to the EVBS response. Use
  `qevbs()` or `evbs_return_level()`.
* `gevreg.fit()` centres are strongly recommended: the GEV likelihood is poorly
  conditioned for uncentred covariates and may silently fail to converge.

# evbsreg 1.0.0

## New Features

* **Initial CRAN Release** - Complete implementation of local influence diagnostics for Extreme-Value Birnbaum-Saunders (EVBS) regression models

* **Estimation** - `evbsreg.fit()` function for joint maximum likelihood estimation of EVBS regression models with flexible parameter specification

* **Diagnostics** - Conformal normal curvature-based local influence diagnostics under three perturbation schemes:
  - Case-weight perturbation
  - Response variable perturbation
  - Explanatory variable perturbation

* **Residuals** - Randomized quantile residuals (`rcoxsnell()`, `rqrandomized()`) with simulation envelopes for model validation

* **Visualization** - Publication-quality diagnostic and density plots:
  - `plot_cnc()` for local influence plots
  - `envelope_qq()` for quantile-quantile plots with envelopes
  - `plot_evbs_alpha()` and `plot_evbs_gama()` for parameter density visualization
  - `plot_aggregate_contributions()` for influence aggregation
  - `plot_normalized_eigenvalues()` for eigenvalue analysis

* **Monte Carlo Utilities** - `generate_evbs_data()` and `generate_logevbs_data()` for simulation studies

* **Random Number Generation** - `revbs()` for generating random variates from EVBS distributions with flexible GEV parent distributions

## Methodology

The methods implemented in this package are described in:
- Ospina, Lima, Barros, and Macedo (2026, submitted)

Application to real-world data:
- Monthly maximum wind gust data from Itajai, Brazil (included in `itajai` dataset)

## Documentation

* Complete function reference with examples
* Comprehensive vignette demonstrating workflow on real data
* CITATION file with proper attribution

## Dependencies

* **Imports**: stats, graphics, SpatialExtremes, ggplot2
* **Suggests**: gamlss, grDevices, knitr, rmarkdown, testthat

---

For more information, visit: https://raydonal.github.io/evbsreg/
