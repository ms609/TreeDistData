library('ape')
library('TreeTools')
library('TreeDist')
library('TreeDistData')

# Lin use nTrees = 100L, nTip = 100L, replicates=1000L
# k1 = 40, 50, 60, 70
# k2 = 10, 20, 30, 40

nTrees = 50 # Quadratic effect on runtime
nTip = 20 # Hyperexponential effect on runtime
replicates = 10 # Linear effect on runtime
message("Running tests on ", nTrees, ' ', nTip, "-leaf trees; ",
        replicates, " replicates.")

LinTestOneSet <- function (nTip, k, nTrees) {
  skeleton <- RandomTree(seq_len(k))
  structure(lapply(seq_len(nTrees), function (XX) {
    tr <- skeleton
    for (i in k + seq_len(nTip - k))
      tr <- AddTip(tr, label = i)
    tr
  }), class='multiPhylo')
}

LinTestTwoSet <- function (nTip, k, nTrees) {
  startTree <- ape::rtree(nTip, br=NULL)

  SwapTwo <- function (x, length.x = length(x)) {
    swapsies <- sample.int(length.x, 2)
    x[swapsies] <- x[rev(swapsies)]
    x
  }

  RepeatLLI <- function (tr, k) {
    labels <- tr$tip.label
    for (i in seq_len(k)) labels <- SwapTwo(labels)
    tr$tip.label <- labels
    tr
  }

  structure(lapply(seq_len(nTrees), function (XX) RepeatLLI(startTree, k)),
            class='multiPhylo')
}

LinTestSPRSet <- function (nTip, k, nTrees) {
  startTree <- ape::rtree(nTip, br=NULL)

  RepeatSPR <- function (tr, k) {
    for (i in seq_len(k)) tr <- TreeSearch::SPR(tr)
    tr
  }

  structure(lapply(seq_len(nTrees), function (XX) RepeatSPR(startTree, k)),
            class='multiPhylo')
}

SpectralClustering <- function (dat, nClusters) {
  # More efficient version of anocva::spectralClustering
  n <- ncol(dat)
  L <- diag(rowSums(dat)) - dat
  eigenVectors <- eigen(L, symmetric = TRUE)$vectors
  cluster::pam(x = eigenVectors[, n - seq_len(nClusters) + 1L],
               k = nClusters, cluster.only=TRUE)
}


LinTest <- function(k, TestSet = LinTestOneSet, nTip = 100L, nTrees = 100L) {
  cat (".")
  trees <- c(TestSet(nTip, k, nTrees), TestSet(nTip, k, nTrees))
  comparison <- CompareAllTrees(trees, slow = FALSE, verbose = FALSE)
  # Too slow to compute
  comparison$mast <- NULL
  comparison$masti <- NULL
  comparison$qd <- NULL

  # NAs not supported
  comparison$nni_t <- NULL

  ClusterOK <- function (Func, ...) apply(
    vapply(comparison, Func, FUN.VALUE = integer(nTrees + nTrees), ...), 2L,
    identical, y = rep(1:2, each=nTrees))

  HClusters <- function (dat, method) {
    clusters <- hclust(as.dist(dat), method)
    cutree(clusters, k = 2L)
  }

  SClusters <- function (dat) {
    if (!is.null(dat)) {
      if (max(dat) <= 1L) dat <- 1 - dat else dat <- max(dat) - dat
      SpectralClustering(as.matrix(dat), 2L)
    } else {
      rep(0L, nTrees + nTrees)
    }
  }

  cbind(spc = ClusterOK(SClusters),
        pam = ClusterOK(cluster::pam, k=2L, diss=TRUE, cluster.only=TRUE),
        h.cmp = ClusterOK(HClusters, method='complete'),
        h.sng = ClusterOK(HClusters, method='single'),
        h.avg = ClusterOK(HClusters, method='average')
        )
}
linTestReturn <- matrix(FALSE, nrow=16L, ncol=5L,
                        dimnames = list(c('vpi', 'vmsi', 'vci', 'nts',
                                          'ja2', 'ja4', 'jna2', 'jna4',
                                          'msd', 'nni_l', 'nni_u', 'spr',
                                          'tbr_l', 'tbr_u', 'rf', 'path'),
                                        c('spc', 'pam', 'h.cmp', 'h.sng', 'h.avg')))
runLinTestReturn <- t(0 * linTestReturn)

RunLinTest <- function (percent, TestSet = LinTestOneSet,
                        nTip = 100L, nTrees = 100L, replicates= 1000L) {
  message("\n  k = ", percent, "% ")
  colSums(aperm(vapply(seq_len(replicates), function (XX)
    LinTest(percent * nTip / 100, TestSet, nTip, nTrees), linTestReturn)), c(3, 1, 2))
}

message("Lin et al. (2012) test one")
linTestOneResults <-
vapply(seq(30L, 70L, by = 10L), RunLinTest, TestSet = LinTestOneSet,
       nTip = nTip, nTrees = nTrees, replicates = replicates,
       FUN.VALUE = runLinTestReturn)
usethis::use_data(linTestOneResults, compress = 'xz', overwrite = TRUE)

message("Lin et al. (2012) test two")
linTestTwoResults <-
vapply(seq(10L, 40L, by = 10L), RunLinTest, TestSet = LinTestTwoSet,
       nTip = nTip, nTrees = nTrees, replicates = replicates,
       FUN.VALUE = runLinTestReturn)
usethis::use_data(linTestTwoResults, compress = 'xz', overwrite = TRUE)

message("SPR cluster recovery test")
linTestSPRResults <-
vapply(seq(10L, 40L, by = 10L), RunLinTest, TestSet = LinTestTwoSet,
       nTip = nTip, nTrees = nTrees, replicates = replicates,
       FUN.VALUE = runLinTestReturn)
usethis::use_data(linTestSPRResults, compress = 'xz', overwrite = TRUE)
