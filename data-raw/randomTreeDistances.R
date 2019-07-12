library('TreeDist')
suppressWarnings(RNGversion("3.5.0")) # Stopgap until we can require R 3.6.0
set.seed(0)

RandomDistances <- function (nLeaves, repls) {
  RandomTree <- function(nTip) ape::rtree(nTip, br=NULL)
  ret <- vapply(nLeaves, function (n) {
    cat('\n', n, 'Leaves ')
    distances <- vapply(seq_len(repls),
                        function (XX) {
                          tr1 <- RandomTree(n)
                          tr2 <- RandomTree(n)
                          cat('.')

                          c(VariationOfArborealInfo(tr1, tr2, normalize=TRUE),
                            VariationOfPartitionInfo(tr1, tr2, normalize=TRUE),
                            VariationOfClusteringInfo(tr1, tr2, normalize=TRUE),
                            Quartet::QuartetDivergence(Quartet::QuartetStatus(tr1, tr2), similarity = FALSE),
                            1 - NyeTreeSimilarity(tr1, tr2, normalize=TRUE),
                            MatchingSplitDistance(tr1, tr2, normalize=FALSE),
                            phangorn::treedist(tr1, tr2) / c(n + n - 6L, 1), # No norm for path
                            phangorn::SPR.dist(tr1, tr2) / (n / 2) # crude normalization!
                            )
                        },
                        double(9L))
    cbind(mean = rowMeans(distances), sd = apply(distances, 1, sd),
          min = apply(distances, 1, min),
          max = apply(distances, 1, max))
    }, matrix(0, ncol=4L, nrow=9L)
  )
  dimnames(ret) <- list(c('vai', 'vpi', 'vci', 'qd', 'nts', 'msd', 'rf',
                          'path', 'spr'),
                        c('mean', 'sd', 'min', 'max'),
                        nLeaves)
  ret
}

nLeaves <- 4:200
randomTreeDistances <- RandomDistances(nLeaves, repls=1000L)
usethis::use_data(randomTreeDistances, compress='gzip', overwrite=TRUE)
