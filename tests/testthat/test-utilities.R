context('utilities.R')

test_that('Pairwise distances calculated correctly', {
  nTrees <- 6L
  nTip <- 16L

  set.seed(0)
  trees <- lapply(rep(nTip, nTrees), ape::rtree, br=NULL)
  lapply(CompareAllTrees(trees), function (dist) {
    expect_equal(c(nTrees, nTrees), dim(dist))
  })


  sprWalk <- vector('list', nTrees)
  sprWalk[[1]] <- lastTree <- TreeTools::PectinateTree(nTip)

  for (i in seq.int(2, nTrees)) {
    sprWalk[[i]] <- lastTree <- TreeSearch::SPR(lastTree)
  }

  lapply(CompareAllTrees(sprWalk), function(dist) {
    expect_equal(c(nTrees, nTrees), dim(dist))
  })
})
