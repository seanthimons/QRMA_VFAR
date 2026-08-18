# Augment dose-response data with fitted values

Augment dose-response data with fitted values

## Usage

``` r
# S3 method for class 'qdr_fit'
augment(x, data = NULL, ...)
```

## Arguments

- x:

  A `qdr_fit` object.

- data:

  Optional grouped dose-response data. Defaults to the fitting data.

- ...:

  Unused.

## Value

The standardized data with `.fitted` and `.resid` columns.
