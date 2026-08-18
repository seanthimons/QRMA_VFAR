# Bootstrap both dose-response models

Bootstrap both dose-response models

## Usage

``` r
bootstrap_dose_response_models(
  object,
  times = 10000L,
  resample = c("observed", "fitted"),
  seed = NULL,
  backend = c("sequential", "mirai"),
  compute = NULL
)
```

## Arguments

- object:

  A `qdr_model_set` or list of `qdr_fit` objects.

- times:

  Number of bootstrap replicates.

- resample:

  Either `"observed"` or `"fitted"`.

- seed:

  Optional integer random seed.

- backend:

  Execution backend. `"sequential"` runs refits in the current R
  process. `"mirai"` distributes refits to daemons previously configured
  with
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html).

- compute:

  Optional mirai compute profile name. Ignored by the sequential
  backend.

## Value

A named list of `qdr_bootstrap` objects.
