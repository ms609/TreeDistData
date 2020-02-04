library('TreeDistData')
library('usethis')
suppressWarnings(RNGversion("3.5.0")) # Stopgap until we can require R 3.6.0
repls <- 1000
ourMethods <- tdMethods[!tdMethods %in% 'nni_t']

# Look for existing data object
use_directory('data')
paths <- fs::path('data', 'randomTreeDistances', ext='rda')
if (file.exists(proj_path(paths))) {
  load(proj_path(paths))
} else {
  randomTreeDistances <- array(NA,
                               dim = c(length(ourMethods), 13, 197),
                               dimnames =
                                 list(
                                   ourMethods,
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
                        cat(".")
                        if (XX %% 72 == 0) cat(' ', XX, "\n           ")
                        tr1 <- RandomTree(nLeaves)
                        tr2 <- RandomTree(nLeaves)
                        TreeDistData:::AllDists(tr1, tr2)
                      },
                      double(length(tdMethods)))
  distances <- distances[tdMethods %in% ourMethods, ]
  t(rbind(apply(distances, 1L, summary),
          apply(distances, 1L, quantile,
                probs = c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)),
          sd = apply(distances, 1L, sd, na.rm = TRUE)
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

cat("\n # # # COMPLETE # # # \n")

# Upper and lower bound for SPR diameter:
# sprUpper = nLeaves - 3L
# sprLower = (nLeaves - 2L) / 2L
# (Allen & Steel 2001)
