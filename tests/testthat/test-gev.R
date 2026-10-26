test_that("GEV analytic score matches numerical differentiation", {
  data(itajai)
  X <- cbind(1, itajai$pressure - mean(itajai$pressure))
  y <- itajai$wind
  f <- gevreg.fit(X, y)
  S <- gev_scores(y, X, f$beta, f$sigma, f$xi)
  lli <- function(par, i) {
    mu <- sum(X[i, ] * par[1:2]); s <- par[3]; x <- par[4]
    z <- (y[i] - mu) / s; u <- 1 + x * z
    -log(s) - (1 + 1 / x) * log(u) - u^(-1 / x)
  }
  p0 <- c(f$beta, f$sigma, f$xi); h <- 1e-6
  num <- vapply(1:4, function(k) {
    pp <- p0; pm <- p0; pp[k] <- pp[k] + h; pm[k] <- pm[k] - h
    (lli(pp, 1) - lli(pm, 1)) / (2 * h)
  }, numeric(1))
  expect_equal(unname(S[1, ]), num, tolerance = 1e-4)
})
test_that("GEV fit converges: scores sum to zero at the MLE", {
  data(itajai)
  X <- cbind(1, itajai$pressure - mean(itajai$pressure))
  f <- gevreg.fit(X, itajai$wind)
  S <- gev_scores(itajai$wind, X, f$beta, f$sigma, f$xi)
  expect_lt(max(abs(colSums(S, na.rm = TRUE))), 1e-4)
})
test_that("GEV diagnostics flag observation 82, as EVBS does", {
  data(itajai)
  X <- cbind(1, itajai$pressure - mean(itajai$pressure))
  d <- cnc_diagnostics_gev(gevreg.fit(X, itajai$wind))
  expect_identical(which.max(d$Bj), 82L)
  expect_true(82L %in% d$flagged)
})

test_that("Bbar for the GEV matches a finite-difference curvature computation", {
  data(itajai)
  X <- cbind(1, itajai$pressure - mean(itajai$pressure)); y <- itajai$wind
  f <- gevreg.fit(X, y)
  lli <- function(par) {
    mu <- as.vector(X %*% par[1:2]); s <- exp(par[3]); x <- par[4]
    u <- 1 + x * (y - mu) / s
    -log(s) - (1 + 1 / x) * log(u) - u^(-1 / x)
  }
  par <- c(f$beta, log(f$sigma), f$xi); h <- 1e-5; k <- 4
  S <- sapply(1:k, function(a) { e <- replace(0 * par, a, h)
    (lli(par + e) - lli(par - e)) / (2 * h) })
  L <- matrix(0, k, k)
  for (a in 1:k) for (b in 1:k) {
    ea <- replace(0 * par, a, h); eb <- replace(0 * par, b, h)
    L[a, b] <- (sum(lli(par + ea + eb)) - sum(lli(par + ea - eb)) -
                sum(lli(par - ea + eb)) + sum(lli(par - ea - eb))) / (4 * h^2)
  }
  Bnum <- evbsreg:::.mean_aggregate(S %*% solve(-L) %*% t(S))
  d <- cnc_diagnostics_gev(f)
  expect_equal(d$Bbar, Bnum, tolerance = 1e-4)
  expect_identical(which.max(d$Bbar), 82L)
})

test_that(".mean_aggregate reproduces colMeans(cnc_diagnostics()$Bj)", {
  data(itajai)
  fit <- evbsreg.fit(cbind(1, itajai$pressure), itajai$wind)
  expect_equal(evbsreg:::.mean_aggregate(fit$B),
               colMeans(cnc_diagnostics(fit)$Bj), tolerance = 1e-8)
})

test_that("gev_monitor: the Itajai chart crosses the limit at April 2017 (t = 82)", {
  skip_on_cran()
  data(itajai)
  X <- cbind(1, itajai$pressure - mean(itajai$pressure))
  gt <- gev_monitor(X, itajai$wind, burn = 40)
  expect_equal(nrow(gt), 85L)
  expect_true(all(gt$F[gt$t %in% 82:95] < 1))
  expect_true(all(gt$F[gt$t %in% 44:81] > 1))
  expect_true(all(gt$flagged[gt$t %in% 82:95] == 82L))
})
