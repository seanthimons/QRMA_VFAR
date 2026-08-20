# Test whether dose-response datasets are statistically poolable

Runs the poolability likelihood-ratio test (Haas, Rose & Gerba): for
each model, fit every dataset separately (unpooled) and the stacked
combination (pooled), then compare the deviance difference to a
chi-squared distribution. Datasets are poolable when the pooled model
does not fit significantly worse than the separate models, i.e. one
shared parameter set is adequate.

## Usage

``` r
poolability_test(
  datasets,
  models = c("exponential", "beta_poisson"),
  alpha = 0.05
)
```

## Arguments

- datasets:

  A named list of dose-response data frames (accepted by
  [`as_dose_response()`](https://seanthimons.github.io/singlehit/reference/as_dose_response.md))
  or already-standardized tibbles. Unnamed elements are named
  `dataset_1`, `dataset_2`, and so on.

- models:

  Character vector of models, a subset of `"exponential"`,
  `"beta_poisson"`, and `"exact_beta_poisson"`.

- alpha:

  Significance level for the chi-squared test.

## Value

A tibble with one row per model: `model`, `n_datasets`,
`deviance_pooled`, `deviance_unpooled`, `lrt_statistic`, `df`,
`chi_square_critical`, `p_value`, `poolable`, `converged`, and a
human-readable `conclusion`.

## Details

Datasets are combined by stacking their rows (each dataset's dose groups
are kept as separate binomial observations); they are not aggregated, so
repeated doses across datasets are preserved.
