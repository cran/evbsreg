test_that("evbsreg.fit reproduces the published Itajai estimates", {
  data(itajai)
  X <- cbind(1, itajai$pressure)
  fit <- evbsreg.fit(X, itajai$wind)

  # Published estimates (paper Table 2)
  expect_equal(unname(fit$coeff[1]), 25.5148, tolerance = 1e-3)
  expect_equal(unname(fit$coeff[2]), -0.0227, tolerance = 1e-3)
  expect_equal(unname(fit$alphahat),  0.1857, tolerance = 1e-3)
  expect_equal(unname(fit$gamahat),  -0.1551, tolerance = 1e-3)
})

test_that("cnc_diagnostics flags observation 82", {
  data(itajai)
  X <- cbind(1, itajai$pressure)
  fit <- evbsreg.fit(X, itajai$wind)
  diag <- cnc_diagnostics(fit)

  # Observation 82 has the largest aggregate contribution at q=7
  expect_equal(which.max(diag$Bj[7, ]), 82L)

  # Top normalized eigenvalue matches the paper
  expect_equal(diag$eigenvalues_norm[1], 0.69678, tolerance = 1e-4)
})

test_that("deletion of obs 82 changes gamma by about -73.67%", {
  data(itajai)
  X <- cbind(1, itajai$pressure)
  fit <- evbsreg.fit(X, itajai$wind)
  fit82 <- evbsreg.fit(X[-82, ], itajai$wind[-82])
  rc <- 100 * (fit82$gamahat - fit$gamahat) / abs(fit$gamahat)
  expect_equal(unname(rc), -73.67, tolerance = 0.1)
})

test_that("revbs returns the requested number of positive variates", {
  set.seed(2023)
  x <- revbs(50, alpha = 0.5, beta = 1, gama = 0.2)
  expect_length(x, 50)
  expect_true(all(x > 0))
})

test_that("revbs draws match the former SpatialExtremes::rgev construction", {
  for (g in c(-0.2, 0, 0.2)) {
    set.seed(1)
    x <- revbs(20, alpha = 0.5, beta = 2, gama = g)
    set.seed(1)
    e <- rexp(20)
    z <- if (g == 0) -log(e) else (e^(-g) - 1) / g
    expect_identical(x, 2 * (1 + 0.25 * z^2 / 2 + 0.5 * z * sqrt(1 + 0.25 * z^2 / 4)))
  }
})

test_that("envelope_qq restores par('pty') without warnings and accepts small nrep", {
  data(itajai)
  X <- cbind(1, itajai$pressure)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  graphics::par(pty = "m")
  expect_silent(env <- envelope_qq(X, itajai$wind, nrep = 100))
  expect_identical(graphics::par("pty"), "m")
  expect_true(all(env$e1 <= env$e2))
  expect_silent(envelope_qq(X, itajai$wind, nrep = 10))
})
