test_that("Ward data is standardized from legacy column names", {
  ward <- ward_fixture()

  expect_s3_class(ward, "tbl_df")
  expect_named(ward, c("dose", "positive", "negative", "total", "response"))
  expect_equal(nrow(ward), 8L)
  expect_equal(sum(ward$total), 59)
  expect_true(all(diff(ward$dose) > 0))
  expect_equal(ward$response, ward$positive / ward$total)
})

test_that("duplicate doses are combined", {
  input <- data.frame(
    dose = c(1, 1, 10),
    pos = c(1, 2, 4),
    neg = c(3, 2, 0)
  )
  result <- as_dose_response(input)

  expect_equal(nrow(result), 2L)
  expect_equal(result$positive[[1L]], 3)
  expect_equal(result$negative[[1L]], 5)
})

test_that("invalid grouped counts fail before fitting", {
  expect_error(
    as_dose_response(data.frame(dose = c(0, 1, 2), pos = 1, neg = 1)),
    "greater than zero"
  )
  expect_error(
    as_dose_response(data.frame(dose = 1:3, pos = c(1, 1.5, 1), neg = 1)),
    "whole numbers"
  )
  expect_error(
    as_dose_response(data.frame(dose = 1:3, cases = 1, controls = 1)),
    "detect the positive"
  )
})
