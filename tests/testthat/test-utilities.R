context('utilities.R')

test_that('Pairwise distances calculated correctly', {
  set.seed(0)
  trees <- lapply(rep(16, 6), ape::rtree, br=NULL)
  lapply(CompareAllTrees(trees), function (dist) {
    expect_equal(c(6, 6), dim(dist))
  })
})
