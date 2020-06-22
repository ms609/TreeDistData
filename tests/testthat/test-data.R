context("Loading data")

test_that("Data dimensions are correct", {
  data('tdMethods', package = 'TreeDistData')
  nMetrics <- length(tdMethods)

  LengthWithout <- function(...) nMetrics - length(list(...))

  data("sevenTipDistances", package = 'TreeDistData')
  expect_equal(LengthWithout('mafi'), length(sevenTipDistances))
  expect_equal(c(945L, 945L), dim(sevenTipDistances[[1]]))

  data("distanceDistribution25", package = 'TreeDistData')
  expect_equal(c(LengthWithout('mafi'), 10000L),
               dim(distanceDistribution25))
  expect_equal(dim(distanceDistribution50), dim(distanceDistribution25))

  noMafiNnit <- tdMethods[!tdMethods %in% c('mafi', 'nni_t')]

  AllDistsThere <- function (x) expect_true(all(noMafiNnit %in% x))

  data("randomTreeDistances", package = 'TreeDistData')
  nLeafMeasurements <- 197L
  AllDistsThere(dimnames(randomTreeDistances)[[1]])
  expect_equal(c(LengthWithout('mafi', 'nni_t'), 13L, nLeafMeasurements),
               dim(randomTreeDistances))

  lapply(bullseyeDistances, function (x) {
    AllDistsThere(names(x))
  })
  lapply(bullseyeMorphScores, function (x) {
    AllDistsThere(dimnames(x)[[2]])
  })
  lapply(bullMoDiScores, function (x) {
    AllDistsThere(dimnames(x)[[2]])
  })
  AllDistsThere(dimnames(distanceDistribution25)[[1]])
  AllDistsThere(dimnames(distanceDistribution50)[[1]])
  AllDistsThere(dimnames(linTestOneResults)[[2]])
  AllDistsThere(dimnames(linTestTwoResults)[[2]])
  AllDistsThere(dimnames(linTestSPRResults)[[2]])


  #TODO update other datasets for all 20 metrics?
})
