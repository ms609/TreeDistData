library('TreeSearch')
devtools::load_all()

# Lin use nTrees = 100L, nTip = 100L, replicates=1000L
# k1 = 40, 50, 60, 70
# k2 = 10, 20, 30, 40

LinTestOneSet <- function (nTip, k, nTrees) {
  skeleton <- ape::rtree(k, br=NULL, tip.label = seq_len(k))
  structure(lapply(seq_len(nTrees), function (XX) {
    tr <- skeleton
    for (i in k + seq_len(nTip - k))
      tr <- AddTip(tr, sample(tr$edge[, 2], 1L), label = i)
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

  structure(lapply(seq_len(nTrees), function (XX) RepeatLLI(startTree, k)), class='multiPhylo')
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
  comparison <- CompareAllTrees(trees)
  maxDist <- c(
        vpi = 1,
        vci = 1,
        qd = 1,
        nts = 1,
        msd = nTip * nTip, # crude
        rf = nTip + nTip - 6,
        path = max(comparison[['path']]), # crude
        spr = nTip + nTip - 2 # crude
      )
  ClusterOK <- function (Func, ...) apply(
    vapply(comparison, Func, FUN.VALUE = integer(nTrees + nTrees), ...), 2L,
    identical, y = rep(1:2, each=nTrees))
  HClusters <- function (dat, method) {
    clusters <- hclust(as.dist(dat), method)
    cutree(clusters, k = 2L)
  }
  SClusters <- function(dat) {
    if (max(dat) <= 1L) dat <- 1 - dat else dat <- max(dat) - dat
    SpectralClustering(dat, 2L)
  }

  cbind(spc = ClusterOK(SClusters),
        pam = ClusterOK(cluster::pam, k=2L, diss=TRUE, cluster.only=TRUE),
        h.cmp = ClusterOK(HClusters, method='complete'),
        h.sng = ClusterOK(HClusters, method='single'),
        h.avg = ClusterOK(HClusters, method='average')
        )
}
linTestReturn <- matrix(FALSE, nrow=9L, ncol=5L)

RunLinTest <- function (k, TestSet = LinTestOneSet, nTip = 100L, nTrees = 100L, replicates= 1000L) {
  cat("\n", k, " ")
  colSums(aperm(vapply(seq_len(replicates), function (XX) LinTest(k, TestSet, nTip, nTrees),
                linTestReturn)), c(3, 1, 2))
}

linTestOneResults <-
vapply(seq(40L, 70L, by=10L), RunLinTest, TestSet=LinTestOneSet,
       nTip=100L, nTrees=100L, replicates=25L, FUN.VALUE = linTestReturn)
usethis::use_data(linTestOneResults, compress='xz', overwrite=TRUE)

linTestTwoResults <-
vapply(seq(10L, 40L, by=10L), RunLinTest, TestSet=LinTestTwoSet,
       nTip=100L, nTrees=100L, replicates=25L, FUN.VALUE = linTestReturn)
usethis::use_data(linTestTwoResults, compress='xz', overwrite=TRUE)
