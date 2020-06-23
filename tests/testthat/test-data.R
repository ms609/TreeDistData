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

  AllDistsThere <- function (x, nni_t = TRUE, mast = TRUE, mafi = FALSE) {
    exclude <- character(0)
    if (!mafi) exclude <- c(exclude, 'mafi')
    if (!nni_t) exclude <- c(exclude, 'nni_t')
    if (!mast) exclude <- c(exclude, 'mast', 'masti')
    methods <- tdMethods[!tdMethods %in% exclude]
    expect_true(all(methods %in% x))
    expect_equal(methods, x[!x %in% exclude])
  }

  data("randomTreeDistances", package = 'TreeDistData')
  nLeafMeasurements <- 197L
  AllDistsThere(dimnames(randomTreeDistances)[[1]], nni_t = FALSE)
  expect_equal(c(LengthWithout('mafi', 'nni_t'), 13L, nLeafMeasurements),
               dim(randomTreeDistances))

  lapply(bullseyeDistances, function (x) {
    AllDistsThere(names(x), nni_t = TRUE)
  })
  lapply(bullseyeMorphScores, function (x) {
    AllDistsThere(dimnames(x)[[2]], nni_t = TRUE)
  })
  lapply(bullMoDiScores, function (x) {
    AllDistsThere(dimnames(x)[[2]], nni_t = TRUE)
  })
  AllDistsThere(names(sevenTipDistances))
  AllDistsThere(dimnames(pectinateDistances11)[[1]], mafi = TRUE)
  AllDistsThere(dimnames(distanceDistribution25)[[1]])
  AllDistsThere(dimnames(distanceDistribution50)[[1]])
  AllDistsThere(dimnames(linTestOneResults)[[2]], FALSE, FALSE)
  AllDistsThere(dimnames(linTestTwoResults)[[2]], FALSE, FALSE)
  AllDistsThere(dimnames(linTestSPRResults)[[2]], FALSE, FALSE)
  AllDistsThere(names(shapeEffect))
  AllDistsThere(names(sprDistances))

})
