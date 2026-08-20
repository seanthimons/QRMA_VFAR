# Fit a set of dose-response models

Fit a set of dose-response models

## Usage

``` r
fit_dose_response_models(
  data,
  models = c("exponential", "beta_poisson"),
  check_trend = TRUE,
  alpha = 0.05
)
```

## Arguments

- data:

  Grouped dose-response data accepted by
  [`as_dose_response()`](https://seanthimons.github.io/singlehit/reference/as_dose_response.md).

- models:

  Character vector of models to fit, a subset of `"exponential"`,
  `"beta_poisson"`, and `"exact_beta_poisson"`. Defaults to the
  exponential and approximate beta-Poisson models; add
  `"exact_beta_poisson"` to include the exact model.

- check_trend:

  If `TRUE`, warn when the increasing trend test does not pass at
  `alpha`.

- alpha:

  Significance level used for the `passes` indicator.

## Value

A named `qdr_model_set` containing one fit per requested model, with the
trend result stored in the `trend` attribute.
