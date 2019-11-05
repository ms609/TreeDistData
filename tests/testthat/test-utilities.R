context('utilities.R')

test_that('Pairwise distances calculated correctly', {
  trees <- lapply(rep(16, 6), ape::rtree, br=NULL)
  dists <- PairwiseDistances(trees, phangorn::RF.dist)
  expect_equal(dists, t(dists))
  expect_equal(as.integer(TreeDist::RobinsonFoulds(trees)), as.integer(dists))
})
