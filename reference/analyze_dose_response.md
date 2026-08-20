# Run the complete microbial dose-response analysis

Standardizes the input, evaluates its trend, fits the requested models,
compares their fit, and optionally bootstraps uncertainty.

## Usage

``` r
analyze_dose_response(
  data,
  models = c("exponential", "beta_poisson"),
  bootstrap_times = 10000L,
  exact_bootstrap_times = 10000L,
  resample = c("observed", "fitted"),
  seed = NULL,
  backend = c("auto", "sequential", "mirai"),
  compute = NULL,
  workers = NULL,
  check_trend = TRUE,
  alpha = 0.05
)
```

## Arguments

- data:

  Grouped dose-response data accepted by
  [`as_dose_response()`](https://seanthimons.github.io/singlehit/reference/as_dose_response.md).

- models:

  Character vector of models to fit and compare, a subset of
  `"exponential"`, `"beta_poisson"`, and `"exact_beta_poisson"`.
  Defaults to the exponential and approximate beta-Poisson models; add
  `"exact_beta_poisson"` to include the exact model.

- bootstrap_times:

  Number of bootstrap replicates for the exponential and approximate
  beta-Poisson models. Use zero to skip their bootstrapping.

- exact_bootstrap_times:

  Number of bootstrap replicates for the exact beta-Poisson model, which
  is roughly an order of magnitude slower to fit. Only used when
  `"exact_beta_poisson"` is in `models`; use zero to fit and compare the
  exact model without bootstrapping it.

- resample:

  Bootstrap method passed to
  [`bootstrap_dose_response()`](https://seanthimons.github.io/singlehit/reference/bootstrap_dose_response.md).

- seed:

  Optional integer random seed.

- backend:

  Bootstrap execution backend passed to
  [`bootstrap_dose_response()`](https://seanthimons.github.io/singlehit/reference/bootstrap_dose_response.md).

- compute:

  Optional mirai compute profile name.

- workers:

  Number of temporary local mirai daemons, passed to
  [`bootstrap_dose_response()`](https://seanthimons.github.io/singlehit/reference/bootstrap_dose_response.md).

- check_trend:

  Warn when an increasing trend is not detected.

- alpha:

  Significance level for trend and goodness-of-fit tests.

## Value

A `qdr_analysis` object containing standardized data, trend test, model
fits, diagnostics, parsed assessment, comparison, and bootstrap results.
