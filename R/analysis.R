#' Run the complete microbial dose-response analysis
#'
#' Standardizes the input, evaluates its trend, fits both supported models,
#' compares their fit, and optionally bootstraps uncertainty.
#'
#' @param data Grouped dose-response data accepted by [as_dose_response()].
#' @param bootstrap_times Number of bootstrap replicates per model. Use zero to
#'   skip bootstrapping.
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
  bootstrap_times = 1000L,
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
  fits <- fit_dose_response_models(data, check_trend = check_trend, alpha = alpha)
  if (
    length(bootstrap_times) != 1L ||
      !is.numeric(bootstrap_times) ||
      !is.finite(bootstrap_times) ||
      bootstrap_times < 0 ||
      bootstrap_times != floor(bootstrap_times)
  ) {
    stop("`bootstrap_times` must be a non-negative whole number.", call. = FALSE)
  }
  bootstraps <- if (bootstrap_times > 0L) {
    bootstrap_dose_response_models(
      fits,
      times = bootstrap_times,
      resample = resample,
      seed = seed,
      backend = backend,
      compute = compute
    )
  } else {
    list()
  }

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
