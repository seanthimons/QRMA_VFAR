# Compare fitted dose-response models

Compares two or more fitted models. The single-hit exponential is the
common nested sub-model, so each higher-parameter model is tested
against it by the chi-squared deviance difference on the extra degrees
of freedom. `preferred` is the lowest-AIC model among those that
significantly improve on the exponential, or the exponential itself when
none do. Models with equal parameter counts (approximate vs exact
beta-Poisson) are not nested in each other and are separated only by
AIC/BIC. Because the exponential is a limiting case rather than a
regular nested model, the chi-squared comparison is an approximation.
With exactly two models this reduces to the original simpler-vs-fuller
deviance test.

## Usage

``` r
compare_dose_response_models(object, alpha = 0.05)
```

## Arguments

- object:

  A `qdr_model_set` or list of two or more `qdr_fit` objects.

- alpha:

  Significance level for the chi-squared deviance comparison.

## Value

A tibble with one row per model: fit metrics, a `converged` flag,
AIC/BIC and their deltas, the nested chi-squared test of each model
against the exponential baseline (`NA` on the baseline row), logical
preference flags, a stable `selection` code, and a human-readable
`conclusion`. A warning is issued when any compared model did not
converge, since a poorly converged richer model can make a simpler model
appear preferred by default.
