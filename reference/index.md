# Package index

## End-to-end workflow

- [`analyze_dose_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/analyze_dose_response.md)
  : Run the complete microbial dose-response analysis

## Data import

- [`as_dose_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/as_dose_response.md)
  : Standardize grouped dose-response data
- [`read_dose_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/read_dose_response.md)
  : Read grouped dose-response data

## Trend test

- [`dose_trend_test()`](https://seanthimons.github.io/QRMA_VFAR/reference/dose_trend_test.md)
  : Pre-screen for an increasing dose-response trend

## Model fitting

- [`fit_dose_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/fit_dose_response.md)
  : Fit a mechanistic dose-response model
- [`fit_dose_response_models()`](https://seanthimons.github.io/QRMA_VFAR/reference/fit_dose_response_models.md)
  : Fit a set of dose-response models

## Model response functions

- [`exponential_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/exponential_response.md)
  : Evaluate an exponential dose-response model
- [`beta_poisson_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/beta_poisson_response.md)
  : Evaluate an approximate beta-Poisson dose-response model
- [`exact_beta_poisson_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/exact_beta_poisson_response.md)
  : Exact beta-Poisson dose-response function
- [`exact_beta_poisson_survival()`](https://seanthimons.github.io/QRMA_VFAR/reference/exact_beta_poisson_survival.md)
  : Exact beta-Poisson survival function
- [`effective_dose()`](https://seanthimons.github.io/QRMA_VFAR/reference/effective_dose.md)
  : Calculate an effective dose
- [`prediction_curve()`](https://seanthimons.github.io/QRMA_VFAR/reference/prediction_curve.md)
  : Create prediction curves and bootstrap intervals

## Diagnostics, comparison, and assessment

- [`goodness_of_fit()`](https://seanthimons.github.io/QRMA_VFAR/reference/goodness_of_fit.md)
  : Goodness-of-fit statistics
- [`compare_dose_response_models()`](https://seanthimons.github.io/QRMA_VFAR/reference/compare_dose_response_models.md)
  : Compare fitted dose-response models
- [`assess_dose_response_models()`](https://seanthimons.github.io/QRMA_VFAR/reference/assess_dose_response_models.md)
  : Assess fitted dose-response models
- [`consensus_model_decision()`](https://seanthimons.github.io/QRMA_VFAR/reference/consensus_model_decision.md)
  : Experimental consensus model decision

## Bootstrap uncertainty

- [`bootstrap_dose_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/bootstrap_dose_response.md)
  : Bootstrap a fitted dose-response model
- [`bootstrap_dose_response_models()`](https://seanthimons.github.io/QRMA_VFAR/reference/bootstrap_dose_response_models.md)
  : Bootstrap both dose-response models
- [`bootstrap_confint()`](https://seanthimons.github.io/QRMA_VFAR/reference/bootstrap_confint.md)
  : Bootstrap percentile confidence intervals

## broom / tidymodels methods

- [`tidy(`*`<qdr_fit>`*`)`](https://seanthimons.github.io/QRMA_VFAR/reference/tidy.qdr_fit.md)
  : Tidy a dose-response fit
- [`glance(`*`<qdr_fit>`*`)`](https://seanthimons.github.io/QRMA_VFAR/reference/glance.qdr_fit.md)
  : Glance at a dose-response fit
- [`augment(`*`<qdr_fit>`*`)`](https://seanthimons.github.io/QRMA_VFAR/reference/augment.qdr_fit.md)
  : Augment dose-response data with fitted values

## Pooling multiple datasets

- [`poolability_test()`](https://seanthimons.github.io/QRMA_VFAR/reference/poolability_test.md)
  : Test whether dose-response datasets are statistically poolable
- [`group_datasets()`](https://seanthimons.github.io/QRMA_VFAR/reference/group_datasets.md)
  : Group dose-response datasets into mutually poolable sets

## Plotting

- [`plot_dose_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/plot_dose_response.md)
  : Plot a fitted dose-response model

## Data

- [`ward_rotavirus`](https://seanthimons.github.io/QRMA_VFAR/reference/ward_rotavirus.md)
  : Ward rotavirus human-challenge dose-response data
