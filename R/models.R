#' Evaluate an exponential dose-response model
#'
#' @param dose A non-negative numeric vector of doses.
#' @param k A positive exponential rate parameter.
#'
#' @return A numeric vector of response probabilities.
#' @export
exponential_response <- function(dose, k) {
  validate_prediction_arguments(dose, k = k)
  -expm1(-k * dose)
}

#' Evaluate an approximate beta-Poisson dose-response model
#'
#' Uses the median-dose parameterization where `n50` is the dose associated
#' with a response probability of 0.5.
#'
#' @param dose A non-negative numeric vector of doses.
#' @param alpha A positive shape parameter.
#' @param n50 A positive median-dose parameter.
#'
#' @return A numeric vector of response probabilities.
#' @export
beta_poisson_response <- function(dose, alpha, n50) {
  validate_prediction_arguments(dose, alpha = alpha, n50 = n50)
  log_ratio <- log(dose) - log(n50) + log_expm1(log(2) / alpha)
  -expm1(-alpha * log1pexp(log_ratio))
}

validate_prediction_arguments <- function(dose, ...) {
  if (!is.numeric(dose) || any(!is.finite(dose)) || any(dose < 0)) {
    stop("`dose` must contain only finite, non-negative values.", call. = FALSE)
  }
  parameters <- list(...)
  for (parameter in names(parameters)) {
    value <- parameters[[parameter]]
    if (!is.numeric(value) || any(!is.finite(value)) || any(value <= 0)) {
      stop(sprintf("`%s` must contain only finite, positive values.", parameter), call. = FALSE)
    }
  }
  invisible(NULL)
}

log_expm1 <- function(x) {
  result <- numeric(length(x))
  large <- x > 50
  result[large] <- x[large] + log1p(-exp(-x[large]))
  result[!large] <- log(expm1(x[!large]))
  attributes(result) <- attributes(x)
  result
}

log1pexp <- function(x) {
  result <- numeric(length(x))
  positive <- x > 0
  result[positive] <- x[positive] + log1p(exp(-x[positive]))
  result[!positive] <- log1p(exp(x[!positive]))
  attributes(result) <- attributes(x)
  result
}

clamp_probability <- function(probability) {
  pmin(pmax(probability, .Machine$double.eps), 1 - .Machine$double.eps)
}

model_probability <- function(model, dose, coefficients) {
  switch(
    model,
    exponential = exponential_response(dose, coefficients[["k"]]),
    beta_poisson = beta_poisson_response(dose, coefficients[["alpha"]], coefficients[["n50"]]),
    exact_beta_poisson = exact_beta_poisson_response(dose, coefficients[["alpha"]], coefficients[["beta"]])
  )
}
