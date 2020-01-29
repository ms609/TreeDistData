library('TreeTools')
set.seed(0L)
nTrees <- 75L
nTip <- 50L

sprWalk <- vector('list', nTrees)
sprWalk[[1]] <- lastTree <- PectinateTree(nTip)

message('Walking')
for (i in seq_len(nTrees)[-1]) {
  sprWalk[[i]] <- lastTree <- TreeSearch::SPR(lastTree)
}

message('Calculating distances')
sprDistances <- TreeDistData::CompareAllTrees(sprWalk, verbose = TRUE)

usethis::use_data(sprDistances, compress='xz', overwrite = TRUE)
