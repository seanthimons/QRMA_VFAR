# Assess fitted dose-response models

Combines absolute chi-squared goodness-of-fit with the relative model
comparison. A preferred model is only marked `recommended` when it also
has an adequate absolute fit.

## Usage

``` r
assess_dose_response_models(object, alpha = 0.05)
```

## Arguments

- object:

  A `qdr_analysis`, `qdr_model_set`, or list of two `qdr_fit` objects.

- alpha:

  Significance level used when diagnostics must be calculated.

## Value

A tibble with one row per model. `recommendation` is one of
`"recommended"`, `"acceptable_alternative"`,
`"preferred_but_inadequate"`, or `"not_recommended"`.
