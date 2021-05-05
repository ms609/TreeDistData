context('utilities.R')

test_that("AllDists doesn't crash", {
  library("TreeTools")
  AllDists(as.phylo(0:5, 6), BalancedTree(6), verbose = FALSE)
  expect_null(NULL)
})

test_that('Pairwise distances calculated correctly', {
  nTrees <- 6L
  nTip <- 16L

  set.seed(0)
  trees <- lapply(rep(nTip, nTrees), TreeTools::RandomTree, root = TRUE)
  trees[[1]] <- TreeTools::BalancedTree(nTip)
  trees[[nTrees - 1L]] <- TreeTools::PectinateTree(nTip)
  class(trees) <- 'multiPhylo'

  lapply(CompareAllTrees(trees), function (dist) {
    expect_equal(c(nTrees, nTrees), dim(as.matrix(dist)))
  })

  sprWalk <- vector('list', nTrees)
  sprWalk[[1]] <- lastTree <- TreeTools::PectinateTree(nTip)

  for (i in seq.int(2L, nTrees)) {
    sprWalk[[i]] <- lastTree <- TreeSearch::SPR(lastTree)
  }

  lapply(CompareAllTrees(sprWalk), function(dist) {
    expect_equal(c(nTrees, nTrees), dim(as.matrix(dist)))
  })
})

test_that("Colours retrieved", {
  expect_equal('#1965B0', TreeDistCol('pid'))
  expect_equal('#90C98765', TreeDistCol('mafi', 65))
  expect_warning(expect_equal(c('#F7F056', 'NA'),
                              TreeDistCol(c('spr', 'whoops'))))
})
