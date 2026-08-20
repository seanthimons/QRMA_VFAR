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
  backend = c("auto", "sequential", "mirai"),
  compute = NULL,
  workers = NULL
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

  Execution backend. `"auto"` uses mirai above 1,000 replicates when it
  is installed and otherwise runs sequentially. `"sequential"` runs
  refits in the current R process. `"mirai"` distributes refits to
  existing daemons or a temporary local pool.

- compute:

  Optional mirai compute profile name. Ignored by the sequential
  backend.

- workers:

  Number of temporary local mirai daemons. `NULL` uses 75% of detected
  physical cores, with a minimum of one. Ignored when using existing
  daemons or the sequential backend.

## Value

A `qdr_bootstrap` tibble with one row per replicate.
