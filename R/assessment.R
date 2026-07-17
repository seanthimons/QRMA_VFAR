#' Assess fitted dose-response models
#'
#' Combines absolute chi-squared goodness-of-fit with the relative model
#' comparison. A preferred model is only marked `recommended` when it also has
#' an adequate absolute fit.
#'
#' @param object A `qdr_analysis`, `qdr_model_set`, or list of two `qdr_fit`
#'   objects.
#' @param alpha Significance level used when diagnostics must be calculated.
#'
#' @return A tibble with one row per model. `recommendation` is one of
#'   `"recommended"`, `"acceptable_alternative"`,
#'   `"preferred_but_inadequate"`, or `"not_recommended"`.
#' @export
assess_dose_response_models <- function(object, alpha = 0.05) {
  if (inherits(object, "qdr_analysis")) {
    return(build_model_assessment(object$goodness_of_fit, object$comparison))
  }

  fits <- as_fit_list(object)
  build_model_assessment(
    goodness_of_fit(fits, alpha = alpha),
    compare_dose_response_models(fits, alpha = alpha)
  )
}

build_model_assessment <- function(goodness, comparison) {
  goodness |>
    dplyr::transmute(
      model = .data$model,
      deviance = .data$deviance,
      gof_df = .data$df,
      gof_critical_value = .data$critical_value,
      gof_p_value = .data$p_value,
      appropriate = .data$good_fit,
      fit_assessment = .data$assessment
    ) |>
    dplyr::left_join(
      comparison |>
        dplyr::transmute(
          model = .data$model,
          preferred = .data$preferred,
          preferred_by_AIC = .data$preferred_by_AIC,
          preferred_by_BIC = .data$preferred_by_BIC,
          deviance_difference = .data$deviance_difference,
          chi_square_df = .data$chi_square_df,
          chi_square_critical = .data$chi_square_critical,
          chi_square_p_value = .data$chi_square_p_value,
          significant_improvement = .data$significant_improvement,
          selection_assessment = .data$selection
        ),
      by = "model"
    ) |>
    dplyr::mutate(
      recommendation = dplyr::case_when(
        .data$appropriate & .data$preferred ~ "recommended",
        .data$appropriate ~ "acceptable_alternative",
        .data$preferred ~ "preferred_but_inadequate",
        TRUE ~ "not_recommended"
      ),
      conclusion = assessment_conclusion(
        .data$model,
        .data$recommendation
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$preferred), dplyr::desc(.data$appropriate))
}

assessment_conclusion <- function(model, recommendation) {
  labels <- vapply(model, model_label, character(1))
  dplyr::case_when(
    recommendation == "recommended" ~ paste(
      labels,
      "shows an adequate fit and is preferred by the chi-squared deviance comparison."
    ),
    recommendation == "acceptable_alternative" ~ paste(
      labels,
      "shows an adequate fit but is not preferred by the chi-squared deviance comparison."
    ),
    recommendation == "preferred_but_inadequate" ~ paste(
      labels,
      "is preferred relative to the alternative but does not show an adequate fit to the data."
    ),
    TRUE ~ paste(
      labels,
      "does not show an adequate fit and is not preferred by the chi-squared deviance comparison."
    )
  )
}

#' Experimental consensus model decision
#'
#' Treats the chi-squared deviance comparison, AIC, and BIC selections as three
#' votes. This is an experimental relative model-selection summary; it does not
#' replace the absolute goodness-of-fit assessment from
#' [assess_dose_response_models()].
#'
#' A model selected by all three criteria receives `"unanimous"` agreement. A
#' model selected by two criteria receives `"majority"` agreement. Tied or
#' unavailable criteria abstain, and fewer than two votes produces
#' `"no_consensus"`.
#'
#' @param object A comparison tibble returned by
#'   [compare_dose_response_models()], a `qdr_analysis`, a `qdr_model_set`, or a
#'   list of two `qdr_fit` objects.
#' @param alpha Significance level used when a comparison must be calculated.
#'
#' @return A tibble with one row per model and columns for each criterion's
#'   vote, total votes, consensus selection, agreement level, and conclusion.
#' @export
consensus_model_decision <- function(object, alpha = 0.05) {
  comparison <- as_model_comparison(object, alpha)
  required <- c("model", "preferred", "preferred_by_AIC", "preferred_by_BIC")
  missing_columns <- setdiff(required, names(comparison))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf("Comparison output is missing required columns: %s.", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }
  if (anyDuplicated(comparison$model) || nrow(comparison) < 2L) {
    stop("Comparison output must contain at least two uniquely named models.", call. = FALSE)
  }
  flag_columns <- required[-1L]
  valid_flags <- vapply(comparison[flag_columns], function(value) is.logical(value) && !anyNA(value), logical(1))
  if (!all(valid_flags)) {
    stop("Model preference columns must contain non-missing logical values.", call. = FALSE)
  }

  models <- comparison$model
  choices <- c(
    chi_squared = criterion_choice(models, comparison$preferred),
    AIC = criterion_choice(models, comparison$preferred_by_AIC),
    BIC = criterion_choice(models, comparison$preferred_by_BIC)
  )
  valid_choices <- stats::setNames(choices %in% models, names(choices))
  votes <- unname(vapply(models, function(model) sum(choices[valid_choices] == model), integer(1)))
  highest_vote <- max(votes)
  winners <- models[votes == highest_vote & highest_vote >= 2L]
  selected_model <- if (length(winners) == 1L) winners[[1L]] else NA_character_
  agreement <- if (is.na(selected_model)) {
    "no_consensus"
  } else if (highest_vote == 3L) {
    "unanimous"
  } else {
    "majority"
  }
  chi_squared_vote <- valid_choices[["chi_squared"]] & models == choices[["chi_squared"]]
  aic_vote <- valid_choices[["AIC"]] & models == choices[["AIC"]]
  bic_vote <- valid_choices[["BIC"]] & models == choices[["BIC"]]
  consensus_selected <- !is.na(selected_model) & models == selected_model
  decision <- if (is.na(selected_model)) {
    rep("no_consensus", length(models))
  } else {
    ifelse(consensus_selected, "selected", "not_selected")
  }

  result <- tibble::tibble(
    model = models,
    chi_squared_vote = chi_squared_vote,
    AIC_vote = aic_vote,
    BIC_vote = bic_vote,
    votes = votes,
    consensus_selected = consensus_selected,
    decision = decision,
    agreement = agreement,
    criteria_available = sum(valid_choices),
    chi_squared_choice = choices[["chi_squared"]],
    AIC_choice = choices[["AIC"]],
    BIC_choice = choices[["BIC"]]
  )
  result |>
    dplyr::mutate(
      conclusion = consensus_conclusion(
        .data$model,
        .data$votes,
        .data$consensus_selected,
        agreement
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$consensus_selected), dplyr::desc(.data$votes), .data$model)
}

as_model_comparison <- function(object, alpha) {
  if (inherits(object, "qdr_analysis")) {
    return(object$comparison)
  }
  if (is.data.frame(object)) {
    return(tibble::as_tibble(object))
  }
  compare_dose_response_models(object, alpha = alpha)
}

criterion_choice <- function(models, selected) {
  choices <- models[selected]
  if (length(choices) == 1L) {
    choices[[1L]]
  } else if (length(choices) > 1L) {
    "tie"
  } else {
    "unavailable"
  }
}

consensus_conclusion <- function(model, votes, selected, agreement) {
  labels <- vapply(model, model_label, character(1))
  dplyr::case_when(
    selected & agreement == "unanimous" ~ paste(
      labels,
      "is selected unanimously by AIC, BIC, and chi-squared comparison."
    ),
    selected & agreement == "majority" ~ paste(labels, "is selected by a two-of-three consensus majority."),
    agreement == "no_consensus" ~ "No model receives the two votes required for consensus selection.",
    TRUE ~ paste(labels, "receives", votes, "of three votes and is not selected by consensus.")
  )
}
