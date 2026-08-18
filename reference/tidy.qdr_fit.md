# Tidy a dose-response fit

Tidy a dose-response fit

## Usage

``` r
# S3 method for class 'qdr_fit'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)
```

## Arguments

- x:

  A `qdr_fit` object.

- conf.int:

  Include Wald confidence intervals.

- conf.level:

  Confidence level for Wald intervals.

- ...:

  Unused.

## Value

A tibble with one row per model parameter.
