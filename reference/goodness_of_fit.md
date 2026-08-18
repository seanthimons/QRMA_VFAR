# Goodness-of-fit statistics

Goodness-of-fit statistics

## Usage

``` r
goodness_of_fit(object, alpha = 0.05)
```

## Arguments

- object:

  A `qdr_fit` or `qdr_model_set` object.

- alpha:

  Significance level for the chi-squared reference test.

## Value

A tibble with one row per model. `good_fit` is the logical decision,
`assessment` is a stable machine-readable code, and `conclusion` is a
human-readable interpretation.
