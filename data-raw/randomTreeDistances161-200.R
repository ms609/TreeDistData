RNGversion("3.6.0")
library('TreeDist')

RandomDistances <- function (nLeaves, repls) {
  set.seed(0)
  RandomTree <- function(nTip) ape::rtree(nTip, br=NULL)
  ret <- vapply(nLeaves, function (n) {
    cat('\n', n, 'Leaves ')
    distances <- vapply(seq_len(repls),
                        function (XX) {
                          tr1 <- RandomTree(n)
                          tr2 <- RandomTree(n)
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
          apply(distances, 1L, quantile, probs=c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)),
          sd = apply(distances, 1L, sd)
    )[c(1, 7:9, 2, 3, 5, 10:12, 6, 4, 13), ])
    }, matrix(0, ncol=13L, nrow=10L)
  )
  dimnames(ret) <- list(c('vpi', 'vmsi', 'vci', 'qd', 'nts', 'msd', 'rf',
                          'path', 'spr', 'sprLB'),
                        c('min', '1%', '5%', '10%', '25%', '50%', '75%',
                          '90%', '95%', '99%', 'max', 'mean', 'sd'),
                        nLeaves)
  ret
}

nLeaves <- 161:200
randomTreeDistances161 <- RandomDistances(nLeaves, repls=1000L)
usethis::use_data(randomTreeDistances161, compress='gzip', overwrite=TRUE)
