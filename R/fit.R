#' Pre-screen for an increasing dose-response trend
#'
#' Performs a one-sided Cochran-Armitage-type trend test using log dose as the
#' score. This is a pre-screening step that confirms a monotonic increasing
#' dose-response relationship (response rising with dose) before models are fit,
#' mirroring the `Zca > 1.645` gate in the Weir CAMRA workflow. The test is
#' directional by design: only an increasing trend passes.
#'
#' @param data Grouped dose-response data accepted by [as_dose_response()].
#' @param alpha Significance level used for the `passes` indicator.
#'
#' @return A one-row tibble containing the test statistic, one-sided p-value,
#'   significance level, and pass indicator.
#' @export
dose_trend_test <- function(data, alpha = 0.05) {
  data <- as_dose_response(data)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one number between zero and one.", call. = FALSE)
  }

  score <- log(data$dose)
  score_mean <- sum(data$total * score) / sum(data$total)
  response_mean <- sum(data$positive) / sum(data$total)
  variance <- response_mean * (1 - response_mean) * sum(data$total * (score - score_mean)^2)
  statistic <- if (variance > 0) {
    sum((score - score_mean) * data$positive) / sqrt(variance)
  } else {
    NA_real_
  }
  p_value <- if (is.na(statistic)) NA_real_ else stats::pnorm(statistic, lower.tail = FALSE)

  tibble::tibble(
    statistic = statistic,
    p_value = p_value,
    alpha = alpha,
    passes = !is.na(p_value) & p_value < alpha,
    alternative = "increasing"
  )
}

#' Fit a mechanistic dose-response model
#'
#' Fits grouped binomial data by maximum likelihood. Parameters are optimized
#' on the log scale to enforce positivity.
#'
#' @param data Grouped dose-response data accepted by [as_dose_response()].
#' @param model Either `"exponential"` or `"beta_poisson"`.
#' @param start Optional named vector of positive parameter starting values.
#' @param control Control list passed to [stats::optim()].
#'
#' @return A `qdr_fit` object.
#' @export
fit_dose_response <- function(
  data,
  model = c("exponential", "beta_poisson", "exact_beta_poisson"),
  start = NULL,
  control = list(maxit = 2000, reltol = 1e-10)
) {
  data <- as_dose_response(data)
  model <- match.arg(model)
  parameter_names <- model_parameters(model)
  starts <- if (is.null(start)) {
    candidate_starts(data, model)
  } else {
    list(validate_start(start, parameter_names))
  }

  objective <- function(log_parameters) {
    coefficients <- stats::setNames(exp(log_parameters), parameter_names)
    probability <- clamp_probability(model_probability(model, data$dose, coefficients))
    -sum(stats::dbinom(data$positive, size = data$total, prob = probability, log = TRUE))
  }

  # Fit from each candidate start and keep the highest-likelihood result: a
  # misspecified model's likelihood can have start-dependent local optima
  # (e.g. an exponential fit to beta-Poisson-shaped data). A supplied `start` is
  # used on its own, so warm-started bootstrap refits remain single-shot.
  optimization <- NULL
  for (seed in starts) {
    candidate <- tryCatch(
      stats::optim(par = log(seed), fn = objective, method = "BFGS", hessian = TRUE, control = control),
      error = function(cnd) NULL
    )
    if (
      !is.null(candidate) &&
        is.finite(candidate$value) &&
        (is.null(optimization) || candidate$value < optimization$value)
    ) {
      optimization <- candidate
    }
  }
  if (is.null(optimization)) {
    stop(
      sprintf("The %s model could not be fit from any starting value.", model_label(model)),
      call. = FALSE
    )
  }
  log_coefficients <- stats::setNames(optimization$par, parameter_names)
  coefficients <- exp(log_coefficients)
  fitted_probability <- clamp_probability(model_probability(model, data$dose, coefficients))
  log_lik <- sum(stats::dbinom(data$positive, size = data$total, prob = fitted_probability, log = TRUE))
  saturated_log_lik <- sum(stats::dbinom(
    data$positive,
    size = data$total,
    prob = data$response,
    log = TRUE
  ))
  vcov_log <- tryCatch(
    solve(optimization$hessian),
    error = function(cnd) matrix(NA_real_, length(parameter_names), length(parameter_names))
  )
  dimnames(vcov_log) <- list(parameter_names, parameter_names)

  fit <- structure(
    list(
      call = match.call(),
      model = model,
      data = data,
      coefficients = coefficients,
      log_coefficients = log_coefficients,
      fitted = fitted_probability,
      log_lik = log_lik,
      deviance = 2 * (saturated_log_lik - log_lik),
      df_residual = nrow(data) - length(coefficients),
      convergence = optimization$convergence,
      message = optimization$message %||% "",
      hessian = optimization$hessian,
      vcov_log = vcov_log,
      counts = optimization$counts
    ),
    class = "qdr_fit"
  )

  if (optimization$convergence != 0L) {
    warning(
      sprintf("The %s model optimizer did not converge (code %d).", model_label(model), optimization$convergence),
      call. = FALSE
    )
  }
  fit
}

#' Fit a set of dose-response models
#'
#' @param models Character vector of models to fit, a subset of
#'   `"exponential"`, `"beta_poisson"`, and `"exact_beta_poisson"`. Defaults to
#'   the exponential and approximate beta-Poisson models; add
#'   `"exact_beta_poisson"` to include the exact model.
#' @param check_trend If `TRUE`, warn when the increasing trend test does not
#'   pass at `alpha`.
#' @inheritParams dose_trend_test
#'
#' @return A named `qdr_model_set` containing one fit per requested model, with
#'   the trend result stored in the `trend` attribute.
#' @export
fit_dose_response_models <- function(
  data,
  models = c("exponential", "beta_poisson"),
  check_trend = TRUE,
  alpha = 0.05
) {
  data <- as_dose_response(data)
  models <- validate_model_set(models)
  trend <- dose_trend_test(data, alpha = alpha)
  if (check_trend && !trend$passes) {
    warning("No statistically significant increasing dose-response trend was detected.", call. = FALSE)
  }

  fits <- stats::setNames(lapply(models, function(model) fit_dose_response(data, model)), models)
  structure(fits, class = c("qdr_model_set", "list"), trend = trend, data = data)
}

validate_model_set <- function(models) {
  valid <- c("exponential", "beta_poisson", "exact_beta_poisson")
  models <- unique(models)
  if (!is.character(models) || length(models) < 1L || !all(models %in% valid)) {
    stop(
      sprintf("`models` must be a non-empty subset of: %s.", paste(valid, collapse = ", ")),
      call. = FALSE
    )
  }
  models
}

model_parameters <- function(model) {
  switch(
    model,
    exponential = "k",
    beta_poisson = c("alpha", "n50"),
    exact_beta_poisson = c("alpha", "beta")
  )
}

# Candidate starting values tried in turn when no `start` is supplied; the fit
# keeps whichever converges to the highest likelihood. Multiple seeds are needed
# because a single data-driven seed cannot avoid every start-dependent local
# optimum across the QMRA-wiki dose-response set:
#   - capped ID50 seed: anchors low-infectivity fits (responses only at the top
#     doses) near their true tiny rate; the cap stops it saturating the max dose.
#   - geometric-mean-dose seed: the historical default, robust for shallow data.
#   - fixed k = exp(-5): the CAMRA reference's fixed exponential seed.
candidate_starts <- function(data, model) {
  parameter_names <- model_parameters(model)
  middle_dose <- exp(stats::median(log(data$dose)))
  median_dose <- estimate_median_dose(data)
  seeds <- switch(
    model,
    exponential = list(
      c(k = min(log(2) / median_dose, 10 / max(data$dose))),
      c(k = log(2) / middle_dose),
      c(k = exp(-5))
    ),
    beta_poisson = list(
      c(alpha = 1, n50 = middle_dose),
      c(alpha = 0.3, n50 = median_dose),
      c(alpha = 0.15, n50 = median_dose)
    ),
    exact_beta_poisson = list(
      c(alpha = 0.5, beta = middle_dose / 3)
    )
  )
  lapply(seeds, validate_start, expected_names = parameter_names)
}

# Seed the optimizer with a data-driven ID50 (dose at ~50% response), found by
# interpolating the observed response over log dose. A geometric-mean-dose seed
# is many orders of magnitude off when responses concentrate at the high-dose
# end (low-infectivity pathogens), which sends optim() into a degenerate flat
# optimum; anchoring the start near the response transition avoids that.
estimate_median_dose <- function(data) {
  ordering <- order(data$dose)
  dose <- data$dose[ordering]
  response <- data$response[ordering]
  if (all(response < 0.5)) {
    return(max(dose))
  }
  crossing <- which(response >= 0.5)[1L]
  if (crossing == 1L) {
    return(dose[1L])
  }
  window <- c(crossing - 1L, crossing)
  log_dose <- stats::approx(response[window], log(dose[window]), xout = 0.5)$y
  exp(log_dose)
}

validate_start <- function(start, expected_names) {
  if (!is.numeric(start) || !identical(names(start), expected_names) || any(!is.finite(start)) || any(start <= 0)) {
    stop(
      sprintf("`start` must be a positive named numeric vector: %s.", paste(expected_names, collapse = ", ")),
      call. = FALSE
    )
  }
  start
}

model_label <- function(model) {
  switch(
    model,
    exponential = "Exponential",
    beta_poisson = "Beta-Poisson",
    exact_beta_poisson = "Exact Beta-Poisson"
  )
}

#' Calculate an effective dose
#'
#' @param object A fitted `qdr_fit` object.
#' @param probability Response probabilities strictly between zero and one.
#'
#' @return A numeric vector of doses.
#' @export
effective_dose <- function(object, probability = 0.5) {
  check_qdr_fit(object)
  if (!is.numeric(probability) || any(!is.finite(probability)) || any(probability <= 0 | probability >= 1)) {
    stop("`probability` must contain values strictly between zero and one.", call. = FALSE)
  }

  if (object$model == "exponential") {
    return(-log1p(-probability) / object$coefficients[["k"]])
  }

  if (object$model == "exact_beta_poisson") {
    alpha <- object$coefficients[["alpha"]]
    beta <- object$coefficients[["beta"]]
    doses <- object$data$dose
    lo <- log(min(doses[doses > 0])) - 5
    hi <- log(max(doses)) + 15
    solve_one <- function(p) {
      target <- 1 - p
      f <- function(log_d) exact_beta_poisson_survival(exp(log_d), alpha, beta) - target
      exp(stats::uniroot(f, interval = c(lo, hi), extendInt = "downX")$root)
    }
    return(vapply(probability, solve_one, numeric(1)))
  }

  alpha <- object$coefficients[["alpha"]]
  n50 <- object$coefficients[["n50"]]
  log_numerator <- log_expm1(-log1p(-probability) / alpha)
  log_denominator <- log_expm1(log(2) / alpha)
  exp(log(n50) + log_numerator - log_denominator)
}

#' Goodness-of-fit statistics
#'
#' @param object A `qdr_fit` or `qdr_model_set` object.
#' @param alpha Significance level for the chi-squared reference test.
#'
#' @return A tibble with one row per model. `good_fit` is the logical decision,
#'   `assessment` is a stable machine-readable code, and `conclusion` is a
#'   human-readable interpretation.
#' @export
goodness_of_fit <- function(object, alpha = 0.05) {
  fits <- as_fit_list(object)
  purrr::map_dfr(fits, function(fit) {
    critical <- stats::qchisq(1 - alpha, df = fit$df_residual)
    result <- tibble::tibble(
      model = fit$model,
      deviance = fit$deviance,
      df = fit$df_residual,
      critical_value = critical,
      p_value = stats::pchisq(fit$deviance, df = fit$df_residual, lower.tail = FALSE),
      good_fit = fit$deviance < critical
    )
    result |>
      dplyr::mutate(
        assessment = ifelse(.data$good_fit, "adequate", "inadequate"),
        conclusion = ifelse(
          .data$good_fit,
          paste(model_label(fit$model), "shows an adequate fit to the data."),
          paste(model_label(fit$model), "does not show an adequate fit to the data.")
        )
      )
  })
}

#' Compare fitted dose-response models
#'
#' Compares two or more fitted models. The single-hit exponential is the common
#' nested sub-model, so each higher-parameter model is tested against it by the
#' chi-squared deviance difference on the extra degrees of freedom. `preferred`
#' is the lowest-AIC model among those that significantly improve on the
#' exponential, or the exponential itself when none do. Models with equal
#' parameter counts (approximate vs exact beta-Poisson) are not nested in each
#' other and are separated only by AIC/BIC. Because the exponential is a limiting
#' case rather than a regular nested model, the chi-squared comparison is an
#' approximation. With exactly two models this reduces to the original
#' simpler-vs-fuller deviance test.
#'
#' @param object A `qdr_model_set` or list of two or more `qdr_fit` objects.
#' @param alpha Significance level for the chi-squared deviance comparison.
#'
#' @return A tibble with one row per model: fit metrics, a `converged` flag,
#'   AIC/BIC and their deltas, the nested chi-squared test of each model against
#'   the exponential baseline (`NA` on the baseline row), logical preference
#'   flags, a stable `selection` code, and a human-readable `conclusion`. A
#'   warning is issued when any compared model did not converge, since a poorly
#'   converged richer model can make a simpler model appear preferred by default.
#' @export
compare_dose_response_models <- function(object, alpha = 0.05) {
  fits <- as_fit_list(object)
  if (length(fits) < 2L) {
    stop("At least two fitted models are required for comparison.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one number between zero and one.", call. = FALSE)
  }
  result <- purrr::imap_dfr(fits, function(fit, name) {
    parameters <- length(fit$coefficients)
    observations <- sum(fit$data$total)
    tibble::tibble(
      model = fit$model,
      parameters = parameters,
      converged = fit$convergence == 0L,
      log_lik = fit$log_lik,
      deviance = fit$deviance,
      AIC = -2 * fit$log_lik + 2 * parameters,
      BIC = -2 * fit$log_lik + log(observations) * parameters
    )
  })
  if (anyDuplicated(result$model)) {
    stop("Each model in a comparison must be a distinct model type.", call. = FALSE)
  }
  # A model that only reached a poor local optimum can make a simpler model look
  # preferred by default. Surface it (the `converged` column carries the flag per
  # model) so an "exponential preferred" result is not silently trusted when a
  # richer model failed to converge rather than genuinely failing to improve.
  if (!all(result$converged)) {
    warning(
      sprintf(
        "These models did not converge; the comparison may be unreliable: %s.",
        paste(vapply(result$model[!result$converged], model_label, character(1)), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # The exponential is the common nested sub-model (unique fewest-parameter model
  # and the limiting case of both beta-Poisson families). Each higher-parameter
  # model is tested against it by the deviance difference on `parameters - 1` df.
  # Equal-parameter models (approximate vs exact beta-Poisson) are not nested and
  # are separated only by AIC/BIC. Without a unique baseline no nested test runs.
  min_parameters <- min(result$parameters)
  has_baseline <- sum(result$parameters == min_parameters) == 1L
  baseline_deviance <- if (has_baseline) {
    result$deviance[result$parameters == min_parameters]
  } else {
    NA_real_
  }

  nested <- has_baseline & result$parameters > min_parameters
  parameter_difference <- ifelse(nested, result$parameters - min_parameters, NA_integer_)
  deviance_difference <- ifelse(nested, pmax(baseline_deviance - result$deviance, 0), NA_real_)
  chi_square_critical <- ifelse(nested, stats::qchisq(1 - alpha, df = parameter_difference), NA_real_)
  chi_square_p_value <- ifelse(
    nested,
    stats::pchisq(deviance_difference, df = parameter_difference, lower.tail = FALSE),
    NA_real_
  )
  significant_improvement <- nested & !is.na(deviance_difference) & deviance_difference > chi_square_critical

  # Preferred: among models that significantly improve on the exponential
  # baseline, the lowest-AIC one; otherwise the baseline (or, with no baseline,
  # the lowest-AIC model). With two models this is exactly "the fuller model if
  # its extra parameter significantly cuts deviance, else the simpler model".
  candidates <- which(significant_improvement)
  selected_model <- if (length(candidates) > 0L) {
    result$model[candidates][which.min(result$AIC[candidates])]
  } else if (has_baseline) {
    result$model[result$parameters == min_parameters]
  } else {
    result$model[which.min(result$AIC)]
  }

  result |>
    dplyr::mutate(
      delta_AIC = .data$AIC - min(.data$AIC),
      delta_BIC = .data$BIC - min(.data$BIC),
      deviance_difference = deviance_difference,
      chi_square_df = parameter_difference,
      chi_square_critical = chi_square_critical,
      chi_square_p_value = chi_square_p_value,
      significant_improvement = significant_improvement,
      preferred = .data$model == selected_model,
      preferred_by_AIC = .data$AIC == min(.data$AIC),
      preferred_by_BIC = .data$BIC == min(.data$BIC),
      selection = ifelse(.data$preferred, "preferred", "not_preferred"),
      conclusion = ifelse(
        .data$preferred,
        paste(vapply(.data$model, model_label, character(1)), "is the preferred model."),
        paste(vapply(.data$model, model_label, character(1)), "is not the preferred model.")
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$preferred), .data$AIC)
}

as_fit_list <- function(object) {
  if (inherits(object, "qdr_fit")) {
    return(list(object))
  }
  if (inherits(object, "qdr_model_set") || (is.list(object) && all(vapply(object, inherits, logical(1), "qdr_fit")))) {
    return(unclass(object))
  }
  stop("`object` must be a qdr_fit or a collection of qdr_fit objects.", call. = FALSE)
}

check_qdr_fit <- function(object) {
  if (!inherits(object, "qdr_fit")) {
    stop("`object` must be a qdr_fit object.", call. = FALSE)
  }
  invisible(object)
}
