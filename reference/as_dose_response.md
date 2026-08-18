# Standardize grouped dose-response data

Validates grouped binomial counts, combines rows with the same dose, and
calculates total subjects and observed response probabilities.

## Usage

``` r
as_dose_response(data, dose = NULL, positive = NULL, negative = NULL)
```

## Arguments

- data:

  A data frame containing dose, positive-response count, and
  negative-response count columns.

- dose, positive, negative:

  Optional source column names.

## Value

A tibble with one row per dose and columns `dose`, `positive`,
`negative`, `total`, and `response`.
