library('TreeTools', warn.conflicts = FALSE)
library('TreeDist')
RNGversion("3.6.0")
set.seed(0)

repls <- 10000L
randomTreePairs25 <- lapply(rep(25L, repls), function (nTip)
  list(RandomTree(nTip), RandomTree(nTip)))
cat("Generated", repls, "random tree pairs.")

distanceDistribution25 <- vapply(randomTreePairs25, function (treePair) {
    TreeDistData:::AllDists(treePair[[1]], treePair[[2]])
  }, c(pid = 0, msid = 0, cid = 0, qd = 0, nye = 0, jnc2 = 0, jnc4 = 0,
       jco2 = 0, jco4 = 0, ms = 0, mast = 0, masti = 0, nni_l = 0,
       nni_t = 0, nni_u = 0, spr = 0, tbr_l = 0, tbr_u = 0, rf = 0,
       icrf = 0, path = 0)
)

usethis::use_data(randomTreePairs25, compress = 'bzip2', overwrite = TRUE)
usethis::use_data(distanceDistribution25, compress = 'xz', overwrite = TRUE)
message("Complete.")
