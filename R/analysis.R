#' Run the complete microbial dose-response analysis
#'
#' Standardizes the input, evaluates its trend, fits the requested models,
#' compares their fit, and optionally bootstraps uncertainty.
#'
#' @param data Grouped dose-response data accepted by [as_dose_response()].
#' @param models Character vector of models to fit and compare, a subset of
#'   `"exponential"`, `"beta_poisson"`, and `"exact_beta_poisson"`. Defaults to
#'   the exponential and approximate beta-Poisson models; add
#'   `"exact_beta_poisson"` to include the exact model.
#' @param bootstrap_times Number of bootstrap replicates for the exponential and
#'   approximate beta-Poisson models. Use zero to skip their bootstrapping.
#' @param exact_bootstrap_times Number of bootstrap replicates for the exact
#'   beta-Poisson model, which is roughly an order of magnitude slower to fit.
#'   Only used when `"exact_beta_poisson"` is in `models`; use zero to fit and
#'   compare the exact model without bootstrapping it.
#' @param resample Bootstrap method passed to [bootstrap_dose_response()].
#' @param seed Optional integer random seed.
#' @param backend Bootstrap execution backend passed to
#'   [bootstrap_dose_response()].
#' @param compute Optional mirai compute profile name.
#' @param check_trend Warn when an increasing trend is not detected.
#' @param alpha Significance level for trend and goodness-of-fit tests.
#'
#' @return A `qdr_analysis` object containing standardized data, trend test,
#'   model fits, diagnostics, parsed assessment, comparison, and bootstrap
#'   results.
#' @export
analyze_dose_response <- function(
  data,
  models = c("exponential", "beta_poisson"),
  bootstrap_times = 1000L,
  exact_bootstrap_times = 10000L,
  resample = c("observed", "fitted"),
  seed = NULL,
  backend = c("sequential", "mirai"),
  compute = NULL,
  check_trend = TRUE,
  alpha = 0.05
) {
  data <- as_dose_response(data)
  resample <- match.arg(resample)
  backend <- match.arg(backend)
  validate_bootstrap_times(bootstrap_times, "bootstrap_times")
  validate_bootstrap_times(exact_bootstrap_times, "exact_bootstrap_times")
  fits <- fit_dose_response_models(data, models = models, check_trend = check_trend, alpha = alpha)

  # Bootstrap each fitted model with its own replicate count: the exact model is
  # far slower to fit, so it gets a separate (opt-out via zero) budget. Models
  # with a zero budget are fit and compared but not bootstrapped.
  bootstraps <- purrr::imap(fits, function(fit, name) {
    times <- if (fit$model == "exact_beta_poisson") exact_bootstrap_times else bootstrap_times
    if (times <= 0L) {
      return(NULL)
    }
    model_seed <- if (is.null(seed)) NULL else seed + match(name, names(fits)) - 1L
    bootstrap_dose_response(
      fit,
      times = times,
      resample = resample,
      seed = model_seed,
      backend = backend,
      compute = compute
    )
  })
  bootstraps <- purrr::compact(bootstraps)

  goodness <- goodness_of_fit(fits, alpha = alpha)
  comparison <- compare_dose_response_models(fits, alpha = alpha)
  assessment <- build_model_assessment(goodness, comparison)

  structure(
    list(
      data = data,
      trend = attr(fits, "trend"),
      fits = fits,
      goodness_of_fit = goodness,
      comparison = comparison,
      assessment = assessment,
      bootstraps = bootstraps
    ),
    class = "qdr_analysis"
  )
}

validate_bootstrap_times <- function(times, name) {
  if (length(times) != 1L || !is.numeric(times) || !is.finite(times) || times < 0 || times != floor(times)) {
    stop(sprintf("`%s` must be a non-negative whole number.", name), call. = FALSE)
  }
  invisible(times)
}

#' @export
print.qdr_analysis <- function(x, ...) {
  cat("<qdr_analysis> microbial dose-response analysis\n")
  cat(sprintf(
    "Trend: Z = %.3f, one-sided p = %.4g\n\n",
    x$trend$statistic,
    x$trend$p_value
  ))
  cat("Model assessment:\n")
  print(
    x$assessment[, c("model", "recommendation", "conclusion")],
    n = Inf
  )
  cat("\nModel comparison:\n")
  print(x$comparison, n = Inf)
  if (length(x$bootstraps) > 0L) {
    cat(sprintf("\nBootstrap replicates per model: %d\n", nrow(x$bootstraps[[1L]])))
  }
  invisible(x)
}
