# Plot a fitted dose-response model

Plot a fitted dose-response model

## Usage

``` r
plot_dose_response(
  object,
  bootstrap = NULL,
  levels = c(0.95, 0.99),
  points = 200L
)
```

## Arguments

- object:

  A `qdr_fit` object.

- bootstrap:

  Optional matching `qdr_bootstrap` object.

- levels:

  Confidence levels for bootstrap pointwise intervals.

- points:

  Number of log-spaced doses in the curve.

## Value

A ggplot object.
