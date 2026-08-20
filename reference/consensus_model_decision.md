# Experimental consensus model decision

Treats the chi-squared deviance comparison, AIC, and BIC selections as
three votes. This is an experimental relative model-selection summary;
it does not replace the absolute goodness-of-fit assessment from
[`assess_dose_response_models()`](https://seanthimons.github.io/singlehit/reference/assess_dose_response_models.md).

## Usage

``` r
consensus_model_decision(object, alpha = 0.05)
```

## Arguments

- object:

  A comparison tibble returned by
  [`compare_dose_response_models()`](https://seanthimons.github.io/singlehit/reference/compare_dose_response_models.md),
  a `qdr_analysis`, a `qdr_model_set`, or a list of two `qdr_fit`
  objects.

- alpha:

  Significance level used when a comparison must be calculated.

## Value

A tibble with one row per model and columns for each criterion's vote,
total votes, consensus selection, agreement level, and conclusion.

## Details

A model selected by all three criteria receives `"unanimous"` agreement.
A model selected by two criteria receives `"majority"` agreement. Tied
or unavailable criteria abstain, and fewer than two votes produces
`"no_consensus"`.
