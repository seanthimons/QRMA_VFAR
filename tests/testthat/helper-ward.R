ward_fixture <- function() {
  path <- system.file("extdata", "Ward_rotavirus.txt", package = "singlehit")
  read_dose_response(path)
}
