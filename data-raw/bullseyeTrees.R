RNGversion('3.6.0')
set.seed(1)
nTrees <- 1000L
nTips <- c(5L, 10L, 20L, 50L)
treesNames <- paste(nTips, 'leaves')

# Generate trees:
bullseyeTrees <- lapply(nTips, function (nTip)
  lapply(rep(nTip, nTrees), ape::rtree)
)
names(bullseyeTrees) <- treesNames
usethis::use_data(bullseyeTrees, overwrite = TRUE)
