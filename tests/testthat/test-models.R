test_that("response functions are vectorized and honor their parameterization", {
  dose <- c(0, 1, 10)

  expect_equal(exponential_response(dose, k = 0.2), 1 - exp(-0.2 * dose))
  expect_equal(beta_poisson_response(10, alpha = 0.5, n50 = 10), 0.5)
  expect_equal(beta_poisson_response(dose, alpha = 0.5, n50 = 10)[[1L]], 0)
})

test_that("Ward model estimates remain stable", {
  ward <- ward_fixture()
  fits <- fit_dose_response_models(ward)

  expect_equal(unname(coef(fits$exponential)[["k"]]), 0.00108819, tolerance = 1e-6)
  expect_equal(unname(coef(fits$beta_poisson)[["alpha"]]), 0.2650067, tolerance = 1e-5)
  expect_equal(unname(coef(fits$beta_poisson)[["n50"]]), 5.597221, tolerance = 1e-5)
  expect_equal(effective_dose(fits$beta_poisson, 0.5), coef(fits$beta_poisson)[["n50"]])
  expect_true(all(vapply(fits, function(fit) fit$convergence == 0L, logical(1))))
})

test_that("Ward diagnostics identify the beta-Poisson model", {
  ward <- ward_fixture()
  trend <- dose_trend_test(ward)
  fits <- fit_dose_response_models(ward)
  goodness <- goodness_of_fit(fits)
  comparison <- compare_dose_response_models(fits)
  assessment <- assess_dose_response_models(fits)

  expect_equal(trend$statistic, 5.035736, tolerance = 1e-6)
  expect_true(trend$passes)
  expect_false(goodness$good_fit[goodness$model == "exponential"])
  expect_true(goodness$good_fit[goodness$model == "beta_poisson"])
  expect_identical(goodness$assessment[goodness$model == "exponential"], "inadequate")
  expect_identical(goodness$assessment[goodness$model == "beta_poisson"], "adequate")
  expect_match(goodness$conclusion[goodness$model == "beta_poisson"], "adequate fit")
  expect_equal(comparison$model[[1L]], "beta_poisson")
  expect_true(comparison$preferred[[1L]])
  expect_true(comparison$preferred_by_AIC[[1L]])
  expect_true(comparison$preferred_by_BIC[[1L]])
  expect_equal(comparison$chi_square_df[[1L]], 1)
  expect_equal(comparison$chi_square_critical[[1L]], stats::qchisq(0.95, 1))
  expect_lt(comparison$chi_square_p_value[[1L]], 0.05)
  expect_identical(comparison$selection[[1L]], "preferred")
  expect_identical(assessment$recommendation[assessment$model == "beta_poisson"], "recommended")
  expect_identical(assessment$recommendation[assessment$model == "exponential"], "not_recommended")
  expect_true(assessment$appropriate[assessment$model == "beta_poisson"])
  expect_false(assessment$appropriate[assessment$model == "exponential"])
})

test_that("chi-squared and information-criterion preferences remain distinct", {
  fits <- fit_dose_response_models(ward_fixture())
  fits$exponential$deviance <- 10
  fits$exponential$log_lik <- -10
  fits$beta_poisson$deviance <- 7
  fits$beta_poisson$log_lik <- -8.5

  comparison <- compare_dose_response_models(fits)
  assessment <- assess_dose_response_models(fits)

  expect_identical(comparison$model[comparison$preferred], "exponential")
  expect_identical(comparison$model[comparison$preferred_by_AIC], "beta_poisson")
  expect_false(unique(comparison$significant_improvement))
  expect_identical(assessment$recommendation[assessment$model == "exponential"], "recommended")
  expect_identical(assessment$recommendation[assessment$model == "beta_poisson"], "acceptable_alternative")
})

test_that("relative preference does not imply adequate fit", {
  fits <- fit_dose_response_models(ward_fixture())
  fits$exponential$deviance <- 30
  fits$exponential$log_lik <- -20
  fits$beta_poisson$deviance <- 20
  fits$beta_poisson$log_lik <- -15

  assessment <- assess_dose_response_models(fits)
  beta_poisson <- assessment[assessment$model == "beta_poisson", ]

  expect_true(beta_poisson$preferred)
  expect_false(beta_poisson$appropriate)
  expect_identical(beta_poisson$recommendation, "preferred_but_inadequate")
  expect_match(beta_poisson$conclusion, "does not show an adequate fit")
})

test_that("Ward model-selection criteria give a unanimous consensus", {
  analysis <- analyze_dose_response(ward_fixture(), bootstrap_times = 0)
  consensus <- consensus_model_decision(analysis)
  beta_poisson <- consensus[consensus$model == "beta_poisson", ]

  expect_identical(beta_poisson$votes, 3L)
  expect_true(beta_poisson$consensus_selected)
  expect_identical(beta_poisson$decision, "selected")
  expect_identical(beta_poisson$agreement, "unanimous")
  expect_match(beta_poisson$conclusion, "selected unanimously")
})

test_that("consensus identifies a two-of-three majority", {
  fits <- fit_dose_response_models(ward_fixture())
  fits$exponential$deviance <- 10
  fits$exponential$log_lik <- -10
  fits$beta_poisson$deviance <- 7
  fits$beta_poisson$log_lik <- -8.5

  consensus <- consensus_model_decision(fits)
  exponential <- consensus[consensus$model == "exponential", ]
  beta_poisson <- consensus[consensus$model == "beta_poisson", ]

  expect_identical(exponential$votes, 2L)
  expect_true(exponential$consensus_selected)
  expect_identical(exponential$agreement, "majority")
  expect_identical(beta_poisson$votes, 1L)
  expect_false(beta_poisson$consensus_selected)
})

test_that("consensus abstains when fewer than two criteria agree", {
  comparison <- compare_dose_response_models(fit_dose_response_models(ward_fixture()))
  comparison$preferred <- FALSE
  comparison$preferred_by_AIC <- comparison$model == "beta_poisson"
  comparison$preferred_by_BIC <- comparison$model == "exponential"

  consensus <- consensus_model_decision(comparison)

  expect_false(any(consensus$consensus_selected))
  expect_true(all(consensus$decision == "no_consensus"))
  expect_true(all(consensus$agreement == "no_consensus"))
  expect_true(all(consensus$criteria_available == 2L))
})

test_that("fit objects provide tidy modeling interfaces", {
  fit <- fit_dose_response(ward_fixture(), "beta_poisson")
  tidied <- generics::tidy(fit, conf.int = TRUE)
  glanced <- generics::glance(fit)
  augmented <- generics::augment(fit)

  expect_named(tidied, c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"))
  expect_equal(glanced$n_groups, 8L)
  expect_equal(glanced$nobs, 59)
  expect_named(augmented, c("dose", "positive", "negative", "total", "response", ".fitted", ".resid"))
  expect_true(all(augmented$.fitted >= 0 & augmented$.fitted <= 1))
  expect_equal(as.numeric(logLik(fit)), fit$log_lik)
})
