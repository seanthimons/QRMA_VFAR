# Evaluate an approximate beta-Poisson dose-response model

Uses the median-dose parameterization where `n50` is the dose associated
with a response probability of 0.5.

## Usage

``` r
beta_poisson_response(dose, alpha, n50)
```

## Arguments

- dose:

  A non-negative numeric vector of doses.

- alpha:

  A positive shape parameter.

- n50:

  A positive median-dose parameter.

## Value

A numeric vector of response probabilities.
