# Create prediction curves and bootstrap intervals

Create prediction curves and bootstrap intervals

## Usage

``` r
prediction_curve(
  object,
  bootstrap = NULL,
  levels = c(0.95, 0.99),
  points = 200L,
  dose_range = NULL
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

- dose_range:

  Optional positive length-two dose range.

## Value

A tibble containing dose, fitted response, confidence level, and
pointwise bounds. Without a bootstrap object, level and bounds are `NA`.
