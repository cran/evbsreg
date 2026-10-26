## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(evbsreg)

## ----data---------------------------------------------------------------------
data(itajai)
str(itajai)
summary(itajai$wind)

## ----fit----------------------------------------------------------------------
X <- cbind(1, itajai$pressure)
fit <- evbsreg.fit(X, itajai$wind)

data.frame(
  Parameter = c("beta0", "beta1", "alpha", "gamma"),
  Estimate  = round(fit$coeff, 4),
  SE        = round(c(fit$stderrors, fit$stderroralpha, fit$stderrorgama), 4)
)

## ----pvalues------------------------------------------------------------------
round(fit$pvalues, 4)

## ----diag---------------------------------------------------------------------
diag <- cnc_diagnostics(fit)

## Top four normalized eigenvalues
round(head(diag$eigenvalues_norm, 4), 5)

## Observations flagged at q = 7
which(diag$Bj[7, ] > diag$bq[7])

## ----plot-cnc, fig.width = 10, fig.height = 4.5-------------------------------
plot_cnc(diag, q = 7)

## ----deletion-----------------------------------------------------------------
fit82 <- evbsreg.fit(X[-82, ], itajai$wind[-82])
rc <- 100 * (fit82$coeff - fit$coeff) / abs(fit$coeff)
names(rc) <- c("beta0", "beta1", "alpha", "gamma")
round(rc, 2)

## ----residuals----------------------------------------------------------------
r <- rqrandomized(X, itajai$wind)
shapiro.test(r)

## ----envelope, fig.width = 6, fig.height = 6----------------------------------
set.seed(2023)
envelope_qq(X, itajai$wind, nrep = 100)

## ----endpoint-----------------------------------------------------------------
eta82 <- sum(X[82, ] * fit$betahat)
c(endpoint = unname(evbs_endpoint(eta82, fit$alphahat, fit$gamahat)),
  observed = itajai$wind[82])

## ----return-level-------------------------------------------------------------
evbs_return_level(fit, x = c(1, mean(itajai$pressure)),
                  period = 12 * 50, es = TRUE)
## 10-, 50- and 100-year levels at the same covariate value
qevbs(1 - 1 / (12 * c(10, 50, 100)),
      eta = sum(c(1, mean(itajai$pressure)) * fit$betahat),
      alpha = fit$alphahat, gama = fit$gamahat)

## ----blockboot----------------------------------------------------------------
set.seed(2023)
bb <- evbs_block_boot(X, itajai$wind, L = 4, B = 100)
round(bb$se, 4)

## ----monitor, fig.width = 7, fig.height = 4.5---------------------------------
ct <- evbs_monitor(X, itajai$wind, burn = 40)
ct[ct$F < 1, c("t", "flagged", "gama", "endpoint", "F")][1:4, ]
plot(ct)

## ----gev----------------------------------------------------------------------
Xc <- cbind(1, itajai$pressure - mean(itajai$pressure))
gfit <- gevreg.fit(Xc, itajai$wind)
gfit$convergence
gdiag <- cnc_diagnostics_gev(gfit)

## ----gev-compare--------------------------------------------------------------
Bbar_evbs <- colMeans(diag$Bj)
head(order(Bbar_evbs, decreasing = TRUE), 3)
head(order(gdiag$Bbar, decreasing = TRUE), 3)
cor(Bbar_evbs, gdiag$Bbar, method = "spearman")

## ----densities, fig.width = 7, fig.height = 5---------------------------------
plot_evbs_alpha()

## ----scripts, eval = FALSE----------------------------------------------------
# # Density figures
# source(system.file("scripts/script_01_density_figures.R", package = "evbsreg"))
# 
# # Itajai application
# source(system.file("scripts/script_02_itajai_application.R", package = "evbsreg"))
# 
# # Monte Carlo study; set m <- 500 inside for a quick check
# source(system.file("scripts/script_03_simulation_scenario1.R", package = "evbsreg"))
# source(system.file("scripts/script_04_simulation_scenario2.R", package = "evbsreg"))
# source(system.file("scripts/script_05_simulation_scenario3.R", package = "evbsreg"))

