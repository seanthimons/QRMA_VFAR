test_that("exact beta-Poisson response is a valid dose-response curve", {
  d <- c(0, 1, 10, 100, 1000)
  P <- exact_beta_poisson_response(d, alpha = 0.265, beta = 5)
  expect_equal(P[[1L]], 0)                 # P(0) = 0
  expect_true(all(P >= 0 & P <= 1))        # probabilities
  expect_true(all(diff(P) > 0))            # strictly increasing in dose
})

test_that("exact beta-Poisson survival matches a high-precision oracle", {
  # Slow (~minutes): the oracle sums at least `dose` terms at 240-bit precision.
  # Gated off by default so routine test()/check() stays fast; opt in on CI major
  # builds via QRMAVFAR_SLOW_TESTS=1.
  if (!nzchar(Sys.getenv("QRMAVFAR_SLOW_TESTS"))) {
    skip("High-precision oracle is slow; set QRMAVFAR_SLOW_TESTS=1 to run.")
  }
  skip_if_not_installed("Rmpfr")

  # Ground-truth survival S = M(a, a+b, -d) via the all-positive Kummer
  # transform summed at high precision (no cancellation, no asymptotic):
  #   M(a, a+b, -d) = exp(-d) * M(beta, a+b, d)
  oracle_survival <- function(a, b, d, prec = 240L) {
    a <- Rmpfr::mpfr(a, prec); b <- Rmpfr::mpfr(b, prec); d <- Rmpfr::mpfr(d, prec)
    A <- b; B <- a + b
    term <- Rmpfr::mpfr(1, prec); s <- term
    tol <- Rmpfr::mpfr(1e-60, prec)
    n <- 0L
    repeat {
      n <- n + 1L
      term <- term * ((A + (n - 1L)) / (B + (n - 1L))) * (d / n)
      s <- s + term
      if (term < s * tol && n > as.numeric(d)) break
    }
    as.numeric(exp(-d) * s)
  }

  # dose capped at 1e3 to keep the high-precision sum tractable; the extreme-dose
  # / extreme-beta regime is certified separately in dev/validation (mpmath).
  grid <- expand.grid(
    a = c(0.05, 0.265, 1, 2),
    b = c(0.1, 1, 5, 1000),
    d = c(1e-3, 10, 100, 1e3)
  )
  err <- mapply(
    function(a, b, d) abs(exact_beta_poisson_survival(d, a, b) - oracle_survival(a, b, d)) /
      abs(oracle_survival(a, b, d)),
    grid$a, grid$b, grid$d
  )
  expect_lt(max(err), 1e-7)
})

test_that("exact beta-Poisson fits and effective_dose inverts the curve", {
  ward <- ward_fixture()
  fit <- fit_dose_response(ward, "exact_beta_poisson")

  expect_identical(fit$model, "exact_beta_poisson")
  expect_named(fit$coefficients, c("alpha", "beta"))
  expect_true(fit$convergence == 0L)

  ed50 <- effective_dose(fit, 0.5)
  expect_equal(
    exact_beta_poisson_response(ed50, fit$coefficients[["alpha"]], fit$coefficients[["beta"]]),
    0.5,
    tolerance = 1e-6
  )
})

test_that("exact and approximate beta-Poisson agree when beta is large", {
  # Large beta is the regime where the approximation is valid; the two forms
  # should coincide. (Small beta is where they diverge and the exact model earns
  # its keep.)
  d <- c(1, 10, 100, 1000)
  alpha <- 0.5
  beta <- 1e5
  n50 <- beta * (2^(1 / alpha) - 1)              # approximate-model N50 for these shapes
  P_exact <- exact_beta_poisson_response(d, alpha, beta)
  P_approx <- beta_poisson_response(d, alpha, n50)
  expect_equal(P_exact, P_approx, tolerance = 1e-4)
})
