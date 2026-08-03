## Prepare the shipped `ward_rotavirus` dataset from the raw source file.
## Run with R 4.5.1: Rscript data-raw/ward_rotavirus.R

ward_rotavirus <- readr::read_delim(
  "inst/extdata/Ward_rotavirus.txt",
  delim = "\t",
  show_col_types = FALSE
)

# Ship the raw 3-column form under the canonical column names so that
# `as_dose_response(ward_rotavirus)` remains a meaningful teaching step. Drop
# readr's `spec`/`problems` attributes so the shipped object is a plain tibble.
ward_rotavirus <- ward_rotavirus |>
  dplyr::rename(positive = pos, negative = neg) |>
  as.data.frame() |>
  tibble::as_tibble()

usethis::use_data(ward_rotavirus, overwrite = TRUE)
