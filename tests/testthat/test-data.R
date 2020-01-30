context("Loading data")

test_that("Data dimensions are correct", {
  nMetrics <- 20L

  data("sevenTipDistances", package='TreeDistData')
  expect_equal(nMetrics, length(sevenTipDistances))
  expect_equal(c(945L, 945L), dim(sevenTipDistances[[1]]))

  data("distanceDistribution25", package='TreeDistData')
  expect_equal(c(nMetrics - 1L, 10000L), dim(distanceDistribution25))

  #TODO update other datasets for all 20 metrics?
  allDistMethods <- c('dpi', 'msid', 'cid', 'qd', 'nts', 'ja2', 'ja4', 'jna2',
                      'jna4', 'msd', 'mast', 'masti', 'nni_l', 'nni_t', 'nni_u',
                      'spr', 'tbr_l', 'tbr_u', 'rf', 'path')
  metrics <- allDistMethods[!allDistMethods %in% 'nni_t']
  nMetrics <- length(metrics)
  data("randomTreeDistances", package='TreeDistData')
  nLeafMeasurements <- 197L
  expect_equal(c(nMetrics, 13L, nLeafMeasurements),
               dim(randomTreeDistances))
  expect_equal(sort(metrics), sort(dimnames(randomTreeDistances)[[1]]))

})
