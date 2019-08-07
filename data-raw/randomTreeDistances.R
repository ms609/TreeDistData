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
                          sprDist <- phangorn::SPR.dist(tr1, tr2)
                          dists <- TreeDistData:::AllDists(tr1, tr2)
                          dists[c(1:9, 9)] /
                            c(rep(1L, 6L), rf = n + n - 6L,
                              path = 1L,
                              sprUpper = n - 3L, sprLower = (n - 2L) / 2L
                              # Upper and lower bound for SPR diameter
                              # (Allen & Steel 2001)
                            )
                        },
                        double(10L))
    t(rbind(apply(distances, 1L, summary),
          apply(distances, 1L, quantile, probs=c(0.01, 0.05, 0.95, 0.99)),
          sd = apply(distances, 1L, sd)
    )[c(1, 7, 8, 2, 3, 4, 9, 10, 6, 5, 11), ])
    }, matrix(0, ncol=11L, nrow=10L)
  )
  dimnames(ret) <- list(c('vpi', 'vmsi', 'vci', 'qd', 'nts', 'msd', 'rf',
                          'path', 'spr', 'sprLB'),
                        c('min', '1%', '5%', '25%', '50%', '75%', '95%', '99%',
                          'max', 'mean', 'sd'),
                        nLeaves)
  ret
}

nLeaves <- 4:200
randomTreeDistances <- RandomDistances(nLeaves, repls=1000L)
usethis::use_data(randomTreeDistances, compress='gzip', overwrite=TRUE)
