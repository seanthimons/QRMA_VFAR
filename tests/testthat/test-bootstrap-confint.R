# Deterministic qdr_bootstrap object so the percentile math is exactly checkable
# without depending on random refits (which test-bootstrap-plot.R already covers).
fake_bootstrap <- function(estimate, converged = rep(TRUE, length(estimate))) {
  tbl <- tibble::tibble(
    model = "exponential",
    replicate = seq_along(estimate),
    converged = converged,
    k = estimate
  )
  structure(tbl, class = c("qdr_bootstrap", class(tbl)), parameter_columns = "k")
}

test_that("bootstrap_confint returns the requested percentile bounds", {
  boot <- fake_bootstrap(1:100)
  ci <- bootstrap_confint(boot, levels = 0.95)

  expect_named(ci, c("model", "term", "level", "lower", "upper"))
  expect_equal(nrow(ci), 1L)
  expect_equal(ci$lower, stats::quantile(1:100, 0.025, names = FALSE))
  expect_equal(ci$upper, stats::quantile(1:100, 0.975, names = FALSE))
})

test_that("wider confidence levels give wider intervals", {
  boot <- fake_bootstrap(1:100)
  ci <- bootstrap_confint(boot, levels = c(0.95, 0.99))

  expect_equal(nrow(ci), 2L)
  lo <- ci$lower[ci$level == 0.99]
  hi <- ci$upper[ci$level == 0.99]
  expect_lt(lo, ci$lower[ci$level == 0.95])
  expect_gt(hi, ci$upper[ci$level == 0.95])
})

test_that("non-converged replicates are excluded from the interval", {
  # 100 good replicates plus two non-converged extreme outliers
  boot <- fake_bootstrap(
    estimate = c(1:100, 1e9, 1e9),
    converged = c(rep(TRUE, 100), FALSE, FALSE)
  )
  ci <- bootstrap_confint(boot, levels = 0.95)

  # bounds match the converged-only sample; the 1e9 outliers are ignored
  expect_equal(ci$upper, stats::quantile(1:100, 0.975, names = FALSE))
  expect_lt(ci$upper, 1e6)
})

test_that("bootstrap_confint validates its inputs", {
  boot <- fake_bootstrap(1:100)
  expect_error(bootstrap_confint(boot, levels = 1.5), "between zero and one")
  expect_error(bootstrap_confint(boot, levels = numeric(0)), "between zero and one")
  expect_error(bootstrap_confint(tibble::tibble(k = 1:3)), "qdr_bootstrap")
})
