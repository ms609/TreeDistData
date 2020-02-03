stop("No longer required. Object deleted.")

library('TreeDist')
suppressWarnings(RNGversion("3.5.0")) # Stopgap until we can require R 3.6.0
set.seed(0)

repls <-  1000000L

RandomTree <- function(nTip) ape::rtree(nTip, br=NULL)
randomTreePairs <- lapply(seq_len(repls), function (i) list(RandomTree(11), RandomTree(11)))
cat("Generated random trees.")

distanceDistribution11 <- vapply(seq_len(repls), function (i) {
  treePair <- randomTreePairs[[i]]
  TreeDistData:::AllDists(treePair[[1]], treePair[[2]])
}, c(vpi = 0, vmsi = 0, vci = 0, qd = 0, nts = 0, msd = 0, rf = 0, path = 0,
     spr = 0)
)

usethis::use_data(distanceDistribution11, compress='xz', overwrite=TRUE)
