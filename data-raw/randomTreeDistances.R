library('TreeDistData')
library('usethis')
suppressWarnings(RNGversion("3.5.0")) # Stopgap until we can require R 3.6.0
repls <- 1000

# Look for existing data object
use_directory('data')
paths <- fs::path('data', 'randomTreeDistances', ext='rda')
if (file.exists(proj_path(paths))) {
  load(proj_path(paths))
} else {
  methodAc <- c('vpi', 'vmsi', 'vci', 'qd', 'nts', 'ja2', 'ja4', 'jna2', 'jna4',
                'msd', 'mast', 'masti',
                'nni_l', 'nni_t', 'nni_u', 'spr', 'tbr_l', 'tbr_u',
                'rf', 'path')
  randomTreeDistances <- array(NA, #dim=c(197, 10, 13),
                               dim = c(length(methodAc), 13, 197),
                               dimnames =
                                 list(
                                   methodAc,
                                   c('min', '1%', '5%', '10%', '25%', '50%', '75%',
                                     '90%', '95%', '99%', 'max', 'mean', 'sd'),
                                   4:200))
  usethis::use_data(randomTreeDistances, compress='gzip', overwrite=TRUE)
}

RandomDistances <- function (nLeaves, repls) {
  set.seed(0)
  RandomTree <- function (nTip) ape::rtree(nTip, br = NULL)
  distances <- vapply(seq_len(repls),
                      function (XX) {
                        if (XX %% 72 == 0) cat("\n           ")
                        tr1 <- RandomTree(nLeaves)
                        tr2 <- RandomTree(nLeaves)
                        dists <- TreeDistData:::AllDists(tr1, tr2)
                        dists[c(1:9, 9)] /
                          c(rep(1L, 6L),
                            rf = nLeaves + nLeaves - 6L,
                            path = 1L,
                            sprUpper = nLeaves - 3L,
                            sprLower = (nLeaves - 2L) / 2L
                            # Upper and lower bound for SPR diameter
                            # (Allen & Steel 2001)
                          )
                      },
                      double(10L))
  t(rbind(apply(distances, 1L, summary),
          apply(distances, 1L, quantile, probs=c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)),
          sd = apply(distances, 1L, sd)
  )[c(1, 7:9, 2, 3, 5, 10:12, 6, 4, 13), ])
}

# Build steadily so that partial dataset is available,
# and so that progress is not lost if script interrupted.
while (any(empty <- is.na(randomTreeDistances[1, 1, ]))) {
  cat(sum(empty), 'to go...\n')
  doNext <- sample(names(empty)[empty], 1L)
  cat('\n', doNext, 'Leaves ')
  dists <- RandomDistances(as.integer(doNext), repls)
  load(proj_path(paths))
  cat('\n', ifelse(empty, '-', 'X'))
  randomTreeDistances[, , doNext] <- dists
  # Compress=xz was better, but encoding errors kept wiping the file |-:
  usethis::use_data(randomTreeDistances, compress='gzip', overwrite=TRUE)
}
