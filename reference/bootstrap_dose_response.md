# Bootstrap a fitted dose-response model

Generates grouped binomial bootstrap samples and refits the model. The
`"observed"` resampling method preserves the legacy CAMRA behavior by
using each dose group's observed response probability. The `"fitted"`
method is a model-based parametric bootstrap.

## Usage

``` r
bootstrap_dose_response(
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

  A `qdr_fit` object.

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

A `qdr_bootstrap` tibble with one row per replicate.
