# qrmavfar

`qrmavfar` converts the CAMRA dose-response modeling script into a reusable R
package. It fits exponential and approximate beta-Poisson models to grouped
binomial microbial response data and returns tidy results instead of writing
files into the working directory.

The package provides:

- validated import of dose, positive-response, and negative-response counts;
- a log-dose trend test;
- binomial maximum-likelihood fits for both mechanistic models;
- broom/tidymodels-compatible `tidy()`, `glance()`, and `augment()` methods;
- chi-squared goodness-of-fit and model-deviance comparisons, with AIC/BIC;
- observed-proportion or fitted-model binomial bootstraps;
- percentile parameter intervals and pointwise confidence curves; and
- `ggplot2` model and bootstrap plots.

Development notes and owner-requested workflow follow-ups live in [TODO.md](TODO.md).

Boosterpak remains optional repository-development tooling. It is not a package
dependency and is not needed to install or use `qrmavfar`.

## Example

```r
library(qrmavfar)

ward_path <- system.file("extdata", "Ward_rotavirus.txt", package = "qrmavfar")
ward <- read_dose_response(ward_path)

analysis <- analyze_dose_response(
  ward,
  bootstrap_times = 1000,
  resample = "observed",
  seed = 2026
)

analysis$trend
analysis$goodness_of_fit
analysis$comparison
analysis$assessment
bootstrap_confint(analysis$bootstraps$beta_poisson)

ggplot2::autoplot(analysis)
ggplot2::autoplot(analysis$bootstraps$beta_poisson)
```

For long bootstrap runs, configure mirai workers once and select the parallel
backend. The package leaves daemon lifecycle management to the caller:

```r
mirai::daemons(4)

analysis <- analyze_dose_response(
  ward,
  bootstrap_times = 10000,
  seed = 2026,
  backend = "mirai"
)

mirai::daemons(0)
```

Bootstrap samples are generated before parallel dispatch, so the same `seed`
produces identical results with `backend = "sequential"` and
`backend = "mirai"`. Only the model refits are sent to workers.

For development, the original Ward source file remains at
`data/raw/Ward_rotavirus.txt`. The legacy script is retained under `dev/` as a
methodological reference; package functions do not source it or depend on its
global variables.
