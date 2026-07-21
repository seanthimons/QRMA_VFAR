test_that("response functions are vectorized and honor their parameterization", {
  dose <- c(0, 1, 10)

  expect_equal(exponential_response(dose, k = 0.2), 1 - exp(-0.2 * dose))
  expect_equal(beta_poisson_response(10, alpha = 0.5, n50 = 10), 0.5)
  expect_equal(beta_poisson_response(dose, alpha = 0.5, n50 = 10)[[1L]], 0)
})

test_that("Ward model estimates remain stable", {
  ward <- ward_fixture()
  fits <- fit_dose_response_models(ward)

  # tolerances allow for optimizer-path variation at the same MLE (multi-start
  # keeps the best-likelihood seed, which can settle ~1e-5 from a single-start fit)
  expect_equal(unname(coef(fits$exponential)[["k"]]), 0.00108819, tolerance = 1e-5)
  expect_equal(unname(coef(fits$beta_poisson)[["alpha"]]), 0.2650067, tolerance = 1e-4)
  expect_equal(unname(coef(fits$beta_poisson)[["n50"]]), 5.597221, tolerance = 1e-4)
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

test_that("low-infectivity data fits without collapsing to a degenerate optimum", {
  # responses appear only at the high-dose end (small k); the old geometric-mean
  # dose start sent optim() to k ~ 0 with a deviance in the hundreds
  low_inf <- as_dose_response(data.frame(
    dose = 10^(0:7),
    pos = c(0, 0, 0, 0, 0, 1, 3, 9),
    neg = c(10, 10, 10, 10, 10, 9, 7, 1)
  ))
  fit <- fit_dose_response(low_inf, "exponential")
  expect_true(fit$convergence == 0L)
  expect_lt(fit$deviance, 15)
  ed50 <- effective_dose(fit, 0.5)
  expect_gt(ed50, 1e5)
  expect_lt(ed50, 1e8)
})

test_that("trend pre-screen is one-sided: only an increasing trend passes", {
  # increasing response with dose -> passes (matches Weir Zca > 1.645 gate)
  increasing <- as_dose_response(data.frame(
    dose = c(1, 10, 100, 1000),
    pos = c(0, 2, 5, 8),
    neg = c(8, 6, 3, 0)
  ))
  up <- dose_trend_test(increasing)
  expect_gt(up$statistic, 1.645)
  expect_equal(up$p_value, stats::pnorm(up$statistic, lower.tail = FALSE))
  expect_true(up$passes)

  # the same data reversed (decreasing response) must NOT pass the pre-screen,
  # even though its trend is equally strong in magnitude
  decreasing <- as_dose_response(data.frame(
    dose = c(1, 10, 100, 1000),
    pos = c(8, 5, 2, 0),
    neg = c(0, 3, 6, 8)
  ))
  down <- dose_trend_test(decreasing)
  expect_lt(down$statistic, 0)
  expect_false(down$passes)
})

test_that("comparison generalizes to three models with the exponential baseline", {
  ward <- ward_fixture()
  fits <- list(
    exponential = fit_dose_response(ward, "exponential"),
    beta_poisson = fit_dose_response(ward, "beta_poisson"),
    exact_beta_poisson = fit_dose_response(ward, "exact_beta_poisson")
  )
  cmp <- compare_dose_response_models(fits)

  expect_equal(nrow(cmp), 3L)

  # the exponential is the nested baseline: no self-test
  exp_row <- cmp[cmp$model == "exponential", ]
  expect_true(is.na(exp_row$chi_square_df))
  expect_false(exp_row$significant_improvement)

  # both 2-parameter beta-Poisson forms are tested against the exponential on 1 df
  bp_rows <- cmp[cmp$model != "exponential", ]
  expect_true(all(bp_rows$chi_square_df == 1))
  expect_true(all(bp_rows$significant_improvement))

  # exactly one preferred model; it is a beta-Poisson form and, among the
  # non-nested 2-parameter models, the lowest-AIC one (the global AIC minimum)
  preferred <- cmp$model[cmp$preferred]
  expect_length(preferred, 1L)
  expect_true(preferred %in% c("beta_poisson", "exact_beta_poisson"))
  expect_identical(preferred, cmp$model[cmp$AIC == min(cmp$AIC)])
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
