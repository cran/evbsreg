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
* **Suggests**: grDevices, knitr, rmarkdown, testthat

---

For more information, visit: https://raydonal.github.io/evbsreg/
