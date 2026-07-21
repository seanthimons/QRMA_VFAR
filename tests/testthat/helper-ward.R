ward_fixture <- function() {
  path <- system.file("extdata", "Ward_rotavirus.txt", package = "qrmavfar")
  read_dose_response(path)
}
