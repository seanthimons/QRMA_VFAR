# Start a non-blocking dose-response bootstrap

Starts mirai bootstrap refits and immediately returns a job. Call
[`collect_bootstrap()`](https://seanthimons.github.io/singlehit/reference/collect_bootstrap.md)
to wait for and retrieve the completed result.

## Usage

``` r
bootstrap_dose_response_async(
  object,
  times = 10000L,
  resample = c("observed", "fitted"),
  seed = NULL,
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

- compute:

  Optional mirai compute profile name. Ignored by the sequential
  backend.

- workers:

  Number of temporary local mirai daemons. `NULL` uses 75% of detected
  physical cores, with a minimum of one. Ignored when using existing
  daemons or the sequential backend.

## Value

A `qdr_bootstrap_job` for
[`collect_bootstrap()`](https://seanthimons.github.io/singlehit/reference/collect_bootstrap.md).
