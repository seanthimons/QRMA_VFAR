# Group dose-response datasets into mutually poolable sets

Assigns datasets to groups such that the datasets within a group are
statistically poolable (via the
[`poolability_test()`](https://seanthimons.github.io/singlehit/reference/poolability_test.md)
likelihood-ratio test) and datasets in different groups are not.
Grouping is performed per model, because poolability can differ between
the exponential and beta-Poisson models.

## Usage

``` r
group_datasets(
  datasets,
  models = c("exponential", "beta_poisson"),
  alpha = 0.05,
  method = c("agglomerative", "exhaustive"),
  max_exhaustive = 6L
)
```

## Arguments

- datasets:

  A named list of dose-response data frames (accepted by
  [`as_dose_response()`](https://seanthimons.github.io/singlehit/reference/as_dose_response.md))
  or already-standardized tibbles. Unnamed elements are named
  `dataset_1`, `dataset_2`, and so on.

- models:

  Character vector of models, a subset of `"exponential"`,
  `"beta_poisson"`, and `"exact_beta_poisson"`.

- alpha:

  Significance level for the chi-squared test.

- method:

  Grouping strategy. `"agglomerative"` (default) starts from singletons
  and repeatedly merges the most compatible poolable pair until no pair
  is poolable; it scales to many datasets. `"exhaustive"` enumerates all
  set partitions and picks the coarsest one whose every group is
  poolable; it is limited to `max_exhaustive` datasets.

- max_exhaustive:

  Maximum number of datasets allowed for the exhaustive method (the
  number of set partitions grows super-exponentially).

## Value

A long tibble with columns `dataset`, `model`, and integer `group`. The
per-model matrix of pairwise pooling p-values is attached as the
`"pairwise"` attribute.
