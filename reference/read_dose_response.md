# Read grouped dose-response data

Reads a delimited text file and standardizes its columns for modeling.
Common positive-response names (`positive`, `pos`, and
`positive_response`) and negative-response names (`negative`, `neg`, and
`negative_response`) are detected automatically.

## Usage

``` r
read_dose_response(
  path,
  delim = NULL,
  dose = NULL,
  positive = NULL,
  negative = NULL
)
```

## Arguments

- path:

  Path to a delimited text file.

- delim:

  Delimiter passed to
  [`readr::read_delim()`](https://readr.tidyverse.org/reference/read_delim.html).
  `NULL` asks readr to detect it.

- dose, positive, negative:

  Optional source column names.

## Value

A tibble with columns `dose`, `positive`, `negative`, `total`, and
`response`.
