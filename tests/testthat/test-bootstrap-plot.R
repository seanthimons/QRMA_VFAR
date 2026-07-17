test_that("bootstrap results are reproducible and tidy", {
  fit <- fit_dose_response(ward_fixture(), "beta_poisson")
  first <- bootstrap_dose_response(fit, times = 20, seed = 42)
  second <- bootstrap_dose_response(fit, times = 20, seed = 42)
  intervals <- bootstrap_confint(first)

  expect_s3_class(first, "qdr_bootstrap")
  expect_identical(attr(first, "backend"), "sequential")
  expect_equal(tibble::as_tibble(first), tibble::as_tibble(second))
  expect_equal(nrow(first), 20L)
  expect_true(all(first$converged))
  expect_setequal(intervals$term, c("alpha", "n50", "ed10", "ed50"))
  expect_true(all(intervals$lower <= intervals$upper))
})

test_that("mirai and sequential backends use identical bootstrap samples", {
  skip_if_not_installed("mirai", minimum_version = "2.5.0")
  skip_if(Sys.getenv("_R_CHECK_PACKAGE_NAME_") == "", "requires an installed package namespace")
  fit <- fit_dose_response(ward_fixture(), "beta_poisson")
  sequential <- bootstrap_dose_response(fit, times = 8, seed = 42)

  mirai::daemons(1, dispatcher = FALSE)
  on.exit(
    {
      mirai::daemons(0)
      Sys.sleep(1)
    },
    add = TRUE
  )
  parallel <- bootstrap_dose_response(fit, times = 8, seed = 42, backend = "mirai")

  expect_identical(attr(parallel, "backend"), "mirai")
  expect_equal(lapply(parallel, identity), lapply(sequential, identity))
})

test_that("prediction curves and plots include bootstrap uncertainty", {
  fit <- fit_dose_response(ward_fixture(), "beta_poisson")
  bootstrap <- bootstrap_dose_response(fit, times = 20, seed = 42)
  curve <- prediction_curve(fit, bootstrap, levels = c(0.95, 0.99), points = 30)

  expect_equal(nrow(curve), 60L)
  expect_setequal(curve$level, c(0.95, 0.99))
  expect_true(all(curve$lower <= curve$upper))
  model_plot <- plot_dose_response(fit, bootstrap)
  bootstrap_plot <- ggplot2::autoplot(bootstrap)
  expect_s3_class(model_plot, "ggplot")
  expect_s3_class(bootstrap_plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(model_plot))
  expect_no_error(ggplot2::ggplot_build(bootstrap_plot))
})

test_that("the complete analysis workflow returns plot-ready results", {
  analysis <- analyze_dose_response(ward_fixture(), bootstrap_times = 10, seed = 2026)

  expect_s3_class(analysis, "qdr_analysis")
  expect_named(analysis, c("data", "trend", "fits", "goodness_of_fit", "comparison", "assessment", "bootstraps"))
  expect_named(analysis$fits, c("exponential", "beta_poisson"))
  expect_named(analysis$bootstraps, c("exponential", "beta_poisson"))
  expect_equal(analysis$assessment, assess_dose_response_models(analysis))
  analysis_plot <- ggplot2::autoplot(analysis, points = 30)
  expect_s3_class(analysis_plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(analysis_plot))
})
