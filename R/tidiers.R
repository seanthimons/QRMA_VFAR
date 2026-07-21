#' @export
print.qdr_fit <- function(x, ...) {
  cat("<qdr_fit>", model_label(x$model), "dose-response model\n")
  print(tidy.qdr_fit(x), row.names = FALSE)
  cat(sprintf("Deviance: %.4f on %d df\n", x$deviance, x$df_residual))
  invisible(x)
}

#' @export
print.qdr_model_set <- function(x, ...) {
  cat("<qdr_model_set>\n")
  print(compare_dose_response_models(x), n = Inf)
  invisible(x)
}

#' Tidy a dose-response fit
#'
#' @param x A `qdr_fit` object.
#' @param conf.int Include Wald confidence intervals.
#' @param conf.level Confidence level for Wald intervals.
#' @param ... Unused.
#'
#' @return A tibble with one row per model parameter.
#' @exportS3Method generics::tidy
tidy.qdr_fit <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  check_qdr_fit(x)
  standard_error <- x$coefficients * sqrt(diag(x$vcov_log))
  result <- tibble::tibble(
    term = names(x$coefficients),
    estimate = unname(x$coefficients),
    std.error = unname(standard_error)
  ) |>
    dplyr::mutate(
      statistic = .data$estimate / .data$std.error,
      p.value = 2 * stats::pnorm(abs(.data$statistic), lower.tail = FALSE)
    )
  if (conf.int) {
    quantile <- stats::qnorm(1 - (1 - conf.level) / 2)
    result <- result |>
      dplyr::mutate(
        conf.low = exp(x$log_coefficients - quantile * sqrt(diag(x$vcov_log))),
        conf.high = exp(x$log_coefficients + quantile * sqrt(diag(x$vcov_log)))
      )
  }
  result
}

#' Glance at a dose-response fit
#'
#' @param x A `qdr_fit` object.
#' @param ... Unused.
#'
#' @return A one-row model summary tibble.
#' @exportS3Method generics::glance
glance.qdr_fit <- function(x, ...) {
  check_qdr_fit(x)
  parameters <- length(x$coefficients)
  observations <- sum(x$data$total)
  tibble::tibble(
    model = x$model,
    nobs = observations,
    n_groups = nrow(x$data),
    logLik = x$log_lik,
    AIC = -2 * x$log_lik + 2 * parameters,
    BIC = -2 * x$log_lik + log(observations) * parameters,
    deviance = x$deviance,
    df.residual = x$df_residual,
    convergence = x$convergence
  )
}

#' Augment dose-response data with fitted values
#'
#' @param x A `qdr_fit` object.
#' @param data Optional grouped dose-response data. Defaults to the fitting
#'   data.
#' @param ... Unused.
#'
#' @return The standardized data with `.fitted` and `.resid` columns.
#' @exportS3Method generics::augment
augment.qdr_fit <- function(x, data = NULL, ...) {
  check_qdr_fit(x)
  data <- if (is.null(data)) x$data else as_dose_response(data)
  fitted <- stats::predict(x, newdata = data)
  data |>
    dplyr::mutate(
      .fitted = fitted,
      .resid = grouped_deviance_residual(.data$positive, .data$total, .data$.fitted)
    )
}

grouped_deviance_residual <- function(positive, total, probability) {
  expected <- total * probability
  negative <- total - positive
  expected_negative <- total - expected
  contribution <- 2 * (xlog_ratio(positive, expected) + xlog_ratio(negative, expected_negative))
  sign(positive - expected) * sqrt(pmax(contribution, 0))
}

xlog_ratio <- function(observed, expected) {
  ifelse(observed == 0, 0, observed * log(observed / expected))
}

#' @export
predict.qdr_fit <- function(object, newdata = NULL, ...) {
  check_qdr_fit(object)
  dose <- if (is.null(newdata)) {
    object$data$dose
  } else if (is.numeric(newdata)) {
    newdata
  } else if (is.data.frame(newdata) && "dose" %in% names(newdata)) {
    newdata$dose
  } else {
    stop("`newdata` must be a numeric dose vector or a data frame with a `dose` column.", call. = FALSE)
  }
  model_probability(object$model, dose, object$coefficients)
}

#' @export
coef.qdr_fit <- function(object, ...) {
  object$coefficients
}

#' @export
logLik.qdr_fit <- function(object, ...) {
  structure(
    object$log_lik,
    class = "logLik",
    df = length(object$coefficients),
    nobs = sum(object$data$total)
  )
}

#' @export
deviance.qdr_fit <- function(object, ...) {
  object$deviance
}
