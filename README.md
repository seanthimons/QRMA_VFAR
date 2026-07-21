# qrmavfar

`qrmavfar` converts the CAMRA dose-response modeling script into a reusable R
package. It fits exponential and approximate beta-Poisson models to grouped
binomial microbial response data and returns tidy results instead of writing
files into the working directory.

The package provides:

- validated import of dose, positive-response, and negative-response counts;
- a log-dose trend test;
- binomial maximum-likelihood fits for both mechanistic models, using a
  multi-start optimizer for reliable convergence (see below);
- broom/tidymodels-compatible `tidy()`, `glance()`, and `augment()` methods;
- chi-squared goodness-of-fit and model-deviance comparisons, with AIC/BIC;
- observed-proportion or fitted-model binomial bootstraps;
- percentile parameter intervals and pointwise confidence curves; and
- `ggplot2` model and bootstrap plots.

### Fitting robustness (multi-start)

When no starting values are supplied, `fit_dose_response()` fits the model from
several candidate starting values and keeps the highest-likelihood result. A
single starting value cannot avoid every start-dependent local optimum across
the range of real dose-response data — for example an exponential fit to
beta-Poisson-shaped data, or a low-infectivity pathogen whose responses appear
only at the highest doses. Validated against the QMRA-wiki reference dataset,
this multi-start strategy reproduces the CAMRA reference fits far more reliably
than a single fixed start. Supplying `start` explicitly (as the bootstrap
warm-start does) uses that value alone, so bootstrap performance is unaffected.

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
