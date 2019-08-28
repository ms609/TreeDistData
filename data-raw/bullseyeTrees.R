set.seed(1)
nTrees <- 1000L
nTips <- c(5L, 10L, 20L, 50L)
treesNames <- paste(nTips, 'tips')
subsamples <- 10:1 * 200

# Generate trees:
bullseyeTrees <- lapply(nTips, function (nTip)
  lapply(seq_len(nTrees), function (XX) ape::rtree(nTip))
)
names(bullseyeTrees) <- treesNames
usethis::use_data(bullseyeTrees, overwrite=TRUE)
