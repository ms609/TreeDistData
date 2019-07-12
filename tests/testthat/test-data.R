context("Loading data")

test_that("Data dimensions are correct", {
  nMetrics <- 9L

  data("sevenTipDistances", package='TreeDistData')
  expect_equal(nMetrics + 1L, length(sevenTipDistances))
  expect_equal(945L, length(sevenTipDistances$shapes))
  vapply(sevenTipDistances[seq_len(nMetrics)], function (item)
    expect_equal(c(945L, 945L), dim(item)), integer(2))

  data("randomTreeDistances", package='TreeDistData')
  nLeafMeasurements <- 197L
  expect_equal(c(nMetrics, 4L, nLeafMeasurements), dim(randomTreeDistances))

  data("distanceDistribution25", package='TreeDistData')
  expect_equal(c(nMetrics, 10000L), dim(distanceDistribution25))
})
