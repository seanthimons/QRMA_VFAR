# Fit a mechanistic dose-response model

Fits grouped binomial data by maximum likelihood. Parameters are
optimized on the log scale to enforce positivity.

## Usage

``` r
fit_dose_response(
  data,
  model = c("exponential", "beta_poisson", "exact_beta_poisson"),
  start = NULL,
  control = list(maxit = 2000, reltol = 1e-10)
)
```

## Arguments

- data:

  Grouped dose-response data accepted by
  [`as_dose_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/as_dose_response.md).

- model:

  Either `"exponential"` or `"beta_poisson"`.

- start:

  Optional named vector of positive parameter starting values.

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A `qdr_fit` object.
