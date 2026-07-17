#' Bootstrap a fitted dose-response model
#'
#' Generates grouped binomial bootstrap samples and refits the model. The
#' `"observed"` resampling method preserves the legacy CAMRA behavior by using
#' each dose group's observed response probability. The `"fitted"` method is a
#' model-based parametric bootstrap.
#'
#' @param object A `qdr_fit` object.
#' @param times Number of bootstrap replicates.
#' @param resample Either `"observed"` or `"fitted"`.
#' @param seed Optional integer random seed.
#' @param backend Execution backend. `"sequential"` runs refits in the current
#'   R process. `"mirai"` distributes refits to daemons previously configured
#'   with [mirai::daemons()].
#' @param compute Optional mirai compute profile name. Ignored by the sequential
#'   backend.
#'
#' @return A `qdr_bootstrap` tibble with one row per replicate.
#' @export
bootstrap_dose_response <- function(
  object,
  times = 1000L,
  resample = c("observed", "fitted"),
  seed = NULL,
  backend = c("sequential", "mirai"),
  compute = NULL
) {
  check_qdr_fit(object)
  resample <- match.arg(resample)
  backend <- match.arg(backend)
  if (length(times) != 1L || !is.numeric(times) || !is.finite(times) || times < 1 || times != floor(times)) {
    stop("`times` must be a positive whole number.", call. = FALSE)
  }
  if (!is.null(compute) && (!is.character(compute) || length(compute) != 1L || is.na(compute) || !nzchar(compute))) {
    stop("`compute` must be NULL or one non-empty mirai compute profile name.", call. = FALSE)
  }
  times <- as.integer(times)
  probability <- if (resample == "observed") object$data$response else object$fitted
  tasks <- generate_bootstrap_tasks(object, probability, times, seed)
  estimates <- if (backend == "sequential") {
    run_sequential_bootstrap(tasks, object)
  } else {
    run_mirai_bootstrap(tasks, object, compute)
  }
  parameter_columns <- c(model_parameters(object$model), "ed10", "ed50")
  structure(
    estimates,
    class = c("qdr_bootstrap", class(estimates)),
    fit = object,
    resample = resample,
    backend = backend,
    parameter_columns = parameter_columns
  )
}

generate_bootstrap_tasks <- function(object, probability, times, seed) {
  generate <- function() {
    lapply(seq_len(times), function(replicate) {
      positive <- stats::rbinom(nrow(object$data), size = object$data$total, prob = probability)
      list(
        replicate = replicate,
        sample_data = tibble::tibble(
          dose = object$data$dose,
          positive = positive,
          negative = object$data$total - positive
        )
      )
    })
  }

  if (is.null(seed)) generate() else withr::with_seed(seed, generate())
}

run_sequential_bootstrap <- function(tasks, object) {
  purrr::map_dfr(tasks, function(task) {
    bootstrap_refit(object, task$sample_data, task$replicate)
  })
}

run_mirai_bootstrap <- function(tasks, object, compute) {
  if (!requireNamespace("mirai", quietly = TRUE)) {
    stop("Package 'mirai' is required for `backend = \"mirai\"`.", call. = FALSE)
  }
  mirai::require_daemons(.compute = compute, call = environment())
  mapped <- mirai::mirai_map(
    tasks,
    function(task, original_fit, refit) {
      refit(original_fit, task$sample_data, task$replicate)
    },
    .args = list(original_fit = object, refit = bootstrap_refit),
    .compute = compute
  )
  results <- mapped[]
  failed <- vapply(results, function(result) inherits(result, c("miraiError", "errorValue")), logical(1))
  if (any(failed)) {
    first_failure <- results[[which(failed)[[1L]]]]
    stop(sprintf("A mirai bootstrap task failed: %s", as.character(first_failure)), call. = FALSE)
  }
  dplyr::bind_rows(results)
}

bootstrap_refit <- function(original_fit, sample_data, replicate) {
  tryCatch(
    {
      fit <- suppressWarnings(fit_dose_response(
        sample_data,
        model = original_fit$model,
        start = original_fit$coefficients
      ))
      values <- c(
        list(
          model = original_fit$model,
          replicate = replicate,
          converged = fit$convergence == 0L,
          log_lik = fit$log_lik,
          deviance = fit$deviance
        ),
        as.list(fit$coefficients),
        list(
          ed10 = effective_dose(fit, 0.1),
          ed50 = effective_dose(fit, 0.5),
          error = NA_character_
        )
      )
      tibble::as_tibble_row(values)
    },
    error = function(cnd) {
      parameters <- stats::setNames(
        as.list(rep(NA_real_, length(model_parameters(original_fit$model)))),
        model_parameters(original_fit$model)
      )
      values <- c(
        list(
          model = original_fit$model,
          replicate = replicate,
          converged = FALSE,
          log_lik = NA_real_,
          deviance = NA_real_
        ),
        parameters,
        list(ed10 = NA_real_, ed50 = NA_real_, error = conditionMessage(cnd))
      )
      tibble::as_tibble_row(values)
    }
  )
}

#' Bootstrap both dose-response models
#'
#' @param object A `qdr_model_set` or list of `qdr_fit` objects.
#' @inheritParams bootstrap_dose_response
#'
#' @return A named list of `qdr_bootstrap` objects.
#' @export
bootstrap_dose_response_models <- function(
  object,
  times = 1000L,
  resample = c("observed", "fitted"),
  seed = NULL,
  backend = c("sequential", "mirai"),
  compute = NULL
) {
  fits <- as_fit_list(object)
  resample <- match.arg(resample)
  backend <- match.arg(backend)
  purrr::imap(fits, function(fit, index) {
    model_seed <- if (is.null(seed)) NULL else seed + match(index, names(fits)) - 1L
    bootstrap_dose_response(
      fit,
      times = times,
      resample = resample,
      seed = model_seed,
      backend = backend,
      compute = compute
    )
  })
}

#' Bootstrap percentile confidence intervals
#'
#' @param object A `qdr_bootstrap` object.
#' @param levels Confidence levels strictly between zero and one.
#'
#' @return A tibble with percentile interval bounds for each parameter and
#'   effective dose.
#' @export
bootstrap_confint <- function(object, levels = c(0.95, 0.99)) {
  check_qdr_bootstrap(object)
  validate_confidence_levels(levels)
  parameter_columns <- attr(object, "parameter_columns")
  estimates <- object |>
    tibble::as_tibble() |>
    dplyr::filter(.data$converged) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(parameter_columns),
      names_to = "term",
      values_to = "estimate"
    )

  purrr::map_dfr(sort(levels), function(level) {
    tail_probability <- (1 - level) / 2
    estimates |>
      dplyr::group_by(.data$model, .data$term) |>
      dplyr::summarise(
        level = level,
        lower = stats::quantile(.data$estimate, tail_probability, na.rm = TRUE, names = FALSE),
        upper = stats::quantile(.data$estimate, 1 - tail_probability, na.rm = TRUE, names = FALSE),
        .groups = "drop"
      )
  })
}

validate_confidence_levels <- function(levels) {
  if (!is.numeric(levels) || length(levels) == 0L || any(!is.finite(levels)) || any(levels <= 0 | levels >= 1)) {
    stop("`levels` must contain values strictly between zero and one.", call. = FALSE)
  }
  invisible(levels)
}

check_qdr_bootstrap <- function(object) {
  if (!inherits(object, "qdr_bootstrap")) {
    stop("`object` must be a qdr_bootstrap object.", call. = FALSE)
  }
  invisible(object)
}

#' @export
print.qdr_bootstrap <- function(x, ...) {
  model <- unique(x$model)
  failures <- sum(!x$converged)
  cat(sprintf(
    "<qdr_bootstrap> %s model: %d replicates [%s]",
    model_label(model),
    nrow(x),
    attr(x, "backend") %||% "sequential"
  ))
  if (failures > 0L) {
    cat(sprintf(" (%d did not converge)", failures))
  }
  cat("\n")
  print(bootstrap_confint(x, levels = 0.95), n = Inf)
  invisible(x)
}
