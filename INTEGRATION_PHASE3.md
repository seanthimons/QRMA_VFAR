# Phase 3 integration — register the exact beta-Poisson model

Makes `exact_beta_poisson` a fittable model with `effective_dose` support, plus
an `Rmpfr` oracle test. Does **not** touch `fit_dose_response_models`,
`compare_dose_response_models`, `assess_*`, or `consensus_*` — the three-model
comparison (Pinch Point 7) is a separate phase.

Apply on Windows against the clean baseline you just committed.

## 1. New files (drop in as-is)

| From (my outputs) | To (repo) |
|---|---|
| `R_exact_beta_poisson.R` | `R/exact_beta_poisson.R` |
| `test-exact-beta-poisson.R` | `tests/testthat/test-exact-beta-poisson.R` |

## 2. `R/models.R` — register in `model_probability()`

Add the third `switch` branch:

```r
model_probability <- function(model, dose, coefficients) {
  switch(
    model,
    exponential = exponential_response(dose, coefficients[["k"]]),
    beta_poisson = beta_poisson_response(dose, coefficients[["alpha"]], coefficients[["n50"]]),
    exact_beta_poisson = exact_beta_poisson_response(dose, coefficients[["alpha"]], coefficients[["beta"]])
  )
}
```

(Remember to add a comma after the `beta_poisson = ...` line.)

## 3. `R/fit.R` — four small edits

**(a)** `fit_dose_response()` signature — allow the new model in `match.arg`:

```r
# was:  model = c("exponential", "beta_poisson"),
model = c("exponential", "beta_poisson", "exact_beta_poisson"),
```

**(b)** `model_parameters()`:

```r
model_parameters <- function(model) {
  switch(model,
    exponential = "k",
    beta_poisson = c("alpha", "n50"),
    exact_beta_poisson = c("alpha", "beta")
  )
}
```

**(c)** `default_start()` — add a start for `(alpha, beta)`. `beta` is chosen so
the implied median dose sits near the data (`N50 ≈ beta·(2^(1/alpha) − 1)`):

```r
default_start <- function(data, model) {
  middle_dose <- exp(stats::median(log(data$dose)))
  switch(
    model,
    exponential = c(k = log(2) / middle_dose),
    beta_poisson = c(alpha = 1, n50 = middle_dose),
    exact_beta_poisson = c(alpha = 0.5, beta = middle_dose / 3)
  )
}
```

**(d)** `model_label()`:

```r
model_label <- function(model) {
  switch(model,
    exponential = "Exponential",
    beta_poisson = "Beta-Poisson",
    exact_beta_poisson = "Exact Beta-Poisson"
  )
}
```

**(e)** `effective_dose()` — add a root-solving branch. Insert it **before** the
existing `alpha <- object$coefficients[["alpha"]]` / `n50` block (that block
assumes the approximate model). The survival is monotone decreasing in dose, so
`uniroot` with `extendInt = "downX"` is robust:

```r
  if (object$model == "exact_beta_poisson") {
    alpha <- object$coefficients[["alpha"]]
    beta  <- object$coefficients[["beta"]]
    doses <- object$data$dose
    lo <- log(min(doses[doses > 0])) - 5
    hi <- log(max(doses)) + 15
    solve_one <- function(p) {
      target <- 1 - p                       # survival at the requested response p
      f <- function(log_d) exact_beta_poisson_survival(exp(log_d), alpha, beta) - target
      exp(stats::uniroot(f, interval = c(lo, hi), extendInt = "downX")$root)
    }
    return(vapply(probability, solve_one, numeric(1)))
  }
```

## 4. `DESCRIPTION` — add the test oracle to `Suggests`

```
Suggests:
    mirai (>= 2.5.0),
    Rmpfr,
    testthat (>= 3.0.0)
```

## 5. Regenerate docs, then check

```r
devtools::document()   # writes man/exact_beta_poisson_response.Rd, updates NAMESPACE
devtools::test()       # test-exact-beta-poisson.R (Rmpfr test skips if not installed)
devtools::check()      # full R CMD check
```

`devtools::document()` adds `export(exact_beta_poisson_response)` and
`export(exact_beta_poisson_survival)` to `NAMESPACE` from the roxygen tags.

## 6. Smoke test

```r
library(qrmavfar)
ward <- read_dose_response(system.file("extdata", "Ward_rotavirus.txt", package = "qrmavfar"))
fit  <- fit_dose_response(ward, "exact_beta_poisson")
coef(fit)                 # alpha, beta
effective_dose(fit, 0.5)  # ID50
# compare to the approximate model on the same data:
coef(fit_dose_response(ward, "beta_poisson"))
```

## What is intentionally deferred

- `fit_dose_response_models()` still fits only the two legacy models. Fitting the
  exact model is opt-in via `fit_dose_response(data, "exact_beta_poisson")` until
  Phase 5 generalizes the model set and comparison.
- `compare_dose_response_models()` / `assess_*` / `consensus_*` remain two-model
  (Pinch Point 7 decision pending).
- `bootstrap_dose_response()` will work on an exact-model fit as-is (it keys off
  `model_parameters()` and `effective_dose()`), but the 10k-replicate performance
  pass and the `bootstrap_prediction_matrix()` / `autoplot` exact branches are
  Phase 4.

## Notes for review

- Parameterization is native `(alpha, beta)` (Pinch Point 1) — matches the
  collaborators' `Probinfbpa(a, b, d)` and the literature.
- Evaluator is pure base R (Pinch Point 2); `Rmpfr` is test-only.
- `effective_dose` uses monotone `uniroot` in log-dose (Pinch Point 6).
- The Kummer/asymptotic routing continuity was checked (~1e-6 across the handoff);
  if BFGS shows any instability on real species, the Pinch-Point-3 single
  continuous evaluator for the fitting path is the fallback.
