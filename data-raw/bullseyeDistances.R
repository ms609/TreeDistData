library('TreeDist')
data('bullseyeTrees', package = 'TreeDistData')
bullseyeDistances <- lapply(bullseyeTrees, function(x) NULL)

for (tips in names(bullseyeTrees)) {
  message(tips)
  trees <- bullseyeTrees[[tips]]
  bullseyeDistances[[tips]] <- TreeDistData::CompareAllTrees(trees, verbose = TRUE)
}

usethis::use_data(bullseyeDistances, compress='xz', overwrite = TRUE)
