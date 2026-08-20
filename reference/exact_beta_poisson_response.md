# Exact beta-Poisson dose-response function

Computes the response probability `P(dose) = 1 - S(dose)` of the exact
beta-Poisson model. Registered in `model_probability()` alongside
[`exponential_response()`](https://seanthimons.github.io/singlehit/reference/exponential_response.md)
and
[`beta_poisson_response()`](https://seanthimons.github.io/singlehit/reference/beta_poisson_response.md).

## Usage

``` r
exact_beta_poisson_response(dose, alpha, beta)
```

## Arguments

- dose:

  A non-negative numeric vector of doses.

- alpha:

  A positive shape parameter of the Beta(alpha, beta) infectivity
  distribution.

- beta:

  A positive shape parameter of the Beta(alpha, beta) infectivity
  distribution.

## Value

A numeric vector of response probabilities.
