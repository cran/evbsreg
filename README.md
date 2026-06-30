# evbsreg <img src="man/figures/logo.png" align="right" height="139" alt="evbsreg website" />

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/Raydonal/evbsreg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Raydonal/evbsreg/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Local influence diagnostics for the Extreme-Value Birnbaum–Saunders
(EVBS) regression model.**

This package implements the methodology of:

> Ospina, R., Lima, J. I. C., Barros, M., and Macêdo, A. M. S. (2026).
> *Local influence diagnostics for the extreme-value Birnbaum–Saunders
> regression model: methodology, validation, and application to anomalous
> wind gusts.* Submitted.

It provides joint maximum likelihood estimation, conformal normal curvature
(CNC) diagnostics under three perturbation schemes, randomized quantile
residuals with simulation envelope, Monte Carlo utilities, and
publication-quality density and diagnostic plots.

## Installation

From a local source tarball:

```r
install.packages("evbsreg_1.0.0.tar.gz", repos = NULL, type = "source")
```

Or from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Raydonal/evbsreg")
```

Dependencies: `SpatialExtremes` (GEV random number generation) and
`ggplot2` (density plots).

## Quick start

```r
library(evbsreg)
data(itajai)

# 1. Fit the EVBS regression model
X   <- cbind(1, itajai$pressure)
fit <- evbsreg.fit(X, itajai$wind)
round(fit$coeff, 4)

# 2. Local influence diagnostics
diag <- cnc_diagnostics(fit)
plot_cnc(diag, q = 7)

# 3. Most influential observation
which(diag$Bj[7, ] > diag$bq[7])

# 4. Refit without it and measure the impact
fit82 <- evbsreg.fit(X[-82, ], itajai$wind[-82])
round(100 * (fit82$coeff - fit$coeff) / abs(fit$coeff), 2)
```

The tail-shape parameter changes by about **−73.67%** when the catastrophic
event of 26 April 2017 (observation 82) is removed, while the regression
structure remains stable.

## Main functions

| Function | Purpose |
|---|---|
| `evbsreg.fit()` | Joint maximum likelihood fit of the EVBS regression model |
| `cnc_diagnostics()` | Conformal normal curvature diagnostics |
| `plot_cnc()` | Two-panel diagnostic figure (eigenvalues + contributions) |
| `rqrandomized()`, `rcoxsnell()` | Quantile and Cox–Snell residuals |
| `envelope_qq()` | Normal probability plot with simulation envelope |
| `revbs()` | EVBS random number generation |
| `evbsreg.fit.mc()` | Monte Carlo simulation study |
| `plot_evbs_alpha()` … | Density plots (Figures 1–2 of the paper) |

See `vignette("evbsreg")` for the full worked example.

## Reproducing the paper

Five standalone scripts reproduce every figure, table, and simulation:

```r
source(system.file("scripts/script_01_density_figures.R",    package = "evbsreg"))
source(system.file("scripts/script_02_itajai_application.R", package = "evbsreg"))
source(system.file("scripts/script_03_simulation_scenario1.R", package = "evbsreg"))
source(system.file("scripts/script_04_simulation_scenario2.R", package = "evbsreg"))
source(system.file("scripts/script_05_simulation_scenario3.R", package = "evbsreg"))
```

Each simulation script defaults to `m = 5000` replicates (matching the
paper). Set `m <- 500` at the top of a script for a quick check.

## Documentation

Full documentation, including the reference index and the "Get started"
vignette, is available at the package website:
<https://raydonal.github.io/evbsreg/>.

## License

MIT © Raydonal Ospina