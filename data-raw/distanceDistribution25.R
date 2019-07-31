library('TreeDist')
suppressWarnings(RNGversion("3.5.0")) # Stopgap until we can require R 3.6.0
set.seed(0)

repls <-  10000L

RandomTree <- function(nTip) ape::rtree(nTip, br=NULL)
randomTreePairs25 <- lapply(seq_len(repls), function (i) list(RandomTree(25), RandomTree(25)))
cat("Generated random trees.")

AllDists <- function (tr1, tr2) {
  cat('.')
  c(VariationOfPhylogeneticInfo(tr1, tr2, normalize=TRUE),
    VariationOfMatchingSplitInfo(tr1, tr2, normalize=TRUE),
    VariationOfClusteringInfo(tr1, tr2, normalize=TRUE),
    Quartet::QuartetDivergence(Quartet::QuartetStatus(tr1, tr2), similarity = FALSE),
    1 - NyeTreeSimilarity(tr1, tr2, normalize=TRUE),
    MatchingSplitDistance(tr1, tr2),
    phangorn::treedist(tr1, tr2),
    phangorn::SPR.dist(tr1, tr2)
  )
}

distanceDistribution25 <- vapply(seq_len(repls), function (i) {
    treePair <- randomTreePairs25[[i]]
    AllDists(treePair[[1]], treePair[[2]])
  }, c(vpi = 0, vmsi = 0, vci = 0, qd = 0, nts = 0, msd = 0, rf = 0, path = 0,
       spr = 0)
)

usethis::use_data(randomTreePairs25, compress='bzip2', overwrite=TRUE)
usethis::use_data(distanceDistribution25, compress='xz', overwrite=TRUE)
