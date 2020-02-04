context("Loading data")

test_that("Data dimensions are correct", {
  data('tdMethods', package = 'TreeDistData')
  nMetrics <- length(tdMethods)

  LengthWithout <- function(...) nMetrics - length(list(...))

  data("sevenTipDistances", package='TreeDistData')
  expect_equal(LengthWithout('mafi'), length(sevenTipDistances))
  expect_equal(c(945L, 945L), dim(sevenTipDistances[[1]]))

  data("distanceDistribution25", package='TreeDistData')
  expect_equal(c(LengthWithout('mafi'), 10000L),
               dim(distanceDistribution25))
  expect_equal(dim(distanceDistribution50), dim(distanceDistribution25))

  #TODO update other datasets for all 20 metrics?
  data("randomTreeDistances", package='TreeDistData')
  nLeafMeasurements <- 197L
  expect_equal(c(LengthWithout('mafi', 'nni_t'), 13L, nLeafMeasurements),
               dim(randomTreeDistances))

})
