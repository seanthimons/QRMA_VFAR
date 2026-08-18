# Exact beta-Poisson survival function

Computes the survival
`S(dose) = 1 - P(dose) = M(alpha, alpha + beta, -dose)` of the exact
beta-Poisson dose-response model (Kummer confluent hypergeometric form),
vectorized over `dose`.

## Usage

``` r
exact_beta_poisson_survival(dose, alpha, beta)
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

A numeric vector of survival probabilities `1 - P(dose)`.
