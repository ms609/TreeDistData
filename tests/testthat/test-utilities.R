context('utilities.R')

test_that("AllDists doesn't crash", {
  library("TreeTools")
  AllDists(as.phylo(0:5, 6), BalancedTree(6), verbose = TRUE)
})

test_that('Pairwise distances calculated correctly', {
  nTrees <- 6L
  nTip <- 16L

  set.seed(0)
  trees <- lapply(rep(nTip, nTrees), ape::rtree, br=NULL)
  trees[[1]] <- TreeTools::BalancedTree(nTip)
  trees[[nTrees - 1L]] <- TreeTools::PectinateTree(nTip)
  class(trees) <- 'multiPhylo'

  dists <- PairwiseDistances(trees, phangorn::RF.dist)
  expect_equivalent(phangorn::RF.dist(trees), dists)

  lapply(CompareAllTrees(trees), function (dist) {
    expect_equal(c(nTrees, nTrees), dim(as.matrix(dist)))
  })


  sprWalk <- vector('list', nTrees)
  sprWalk[[1]] <- lastTree <- TreeTools::PectinateTree(nTip)

  for (i in seq.int(2, nTrees)) {
    sprWalk[[i]] <- lastTree <- TreeSearch::SPR(lastTree)
  }

  lapply(CompareAllTrees(sprWalk), function(dist) {
    expect_equal(c(nTrees, nTrees), dim(as.matrix(dist)))
  })
})