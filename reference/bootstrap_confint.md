# Bootstrap percentile confidence intervals

Bootstrap percentile confidence intervals

## Usage

``` r
bootstrap_confint(object, levels = c(0.95, 0.99))
```

## Arguments

- object:

  A `qdr_bootstrap` object.

- levels:

  Confidence levels strictly between zero and one.

## Value

A tibble with percentile interval bounds for each parameter and
effective dose.
