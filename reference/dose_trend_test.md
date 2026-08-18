# Pre-screen for an increasing dose-response trend

Performs a one-sided Cochran-Armitage-type trend test using log dose as
the score. This is a pre-screening step that confirms a monotonic
increasing dose-response relationship (response rising with dose) before
models are fit, mirroring the `Zca > 1.645` gate in the Weir CAMRA
workflow. The test is directional by design: only an increasing trend
passes.

## Usage

``` r
dose_trend_test(data, alpha = 0.05)
```

## Arguments

- data:

  Grouped dose-response data accepted by
  [`as_dose_response()`](https://seanthimons.github.io/QRMA_VFAR/reference/as_dose_response.md).

- alpha:

  Significance level used for the `passes` indicator.

## Value

A one-row tibble containing the test statistic, one-sided p-value,
significance level, and pass indicator.
