library('phangorn')
set.seed(0)
nTrees <- 3L
nTip <- 15L

bullseyeTrees <- lapply(seq_len(nTrees), function (XX) ape::rtree(20L))
usethis::use_data(bullseyeTrees, compress='xz', overwrite=TRUE)

seqs <- lapply(bullseyeTrees, simSeq, l = 2000)

SubSample <- function (dat, n) {
  at <- attributes(dat)
  chosenInd <- at$index[seq_len(n)]
  ret <- lapply(dat, function (x) x[seq_len(max(chosenInd))])
  attributes(ret) <- at
  attr(ret, 'index') <- chosenInd
  attr(ret, 'nr') <- n
  ret
}

bullseyeInferred <- lapply(seq_len(nTrees), function (i) {
  tr <- bullseyeTrees[[i]]
  sq <- seqs[[i]]
  subsamples <- c(2000, 1800, 1600, 1400, 1200, 1000, 800)
  lapply(subsamples, function (n) {
    fitStart <- pml(tr, SubSample(seqs[[i]], n))
    optim.pml(fitStart, model="GTR", optNni=TRUE, rearrangement="stochastic",
              trace=-1L)$tree
  })
})

usethis::use_data(bullseyeInferred, compress='xz', overwrite=TRUE)

vapply(seq_along(bullseyeInferred), function (i) {
  tr <- root(bullseyeTrees[[i]], 't1', resolve.root=TRUE)
  tr$edge.length  <- NULL
  trs <- lapply(bullseyeInferred[[i]], root, 't1', resolve.root=TRUE)

  normInfo <- PartitionInfo(tr)
  cbind(
    vpi = VariationOfPhylogeneticInfo(tr, trs, normalize=normInfo),
    vmsi = VariationOfMatchingSplitInfo(tr, trs, normalize=normInfo),
    vci = VariationOfClusteringInfo(tr, trs, normalize=normInfo),
    qd = Quartet::QuartetDivergence(Quartet::QuartetStatus(trs, cf=tr), similarity = FALSE),
    nts = 1 - NyeTreeSimilarity(tr, trs, normalize=TRUE),
    msd = MatchingSplitDistance(tr, trs),
    t(vapply(trs, phangorn::treedist, tree2=tr, double(2))),
    spr = vapply(trs, phangorn::SPR.dist, tree2=tr, double(1))
  )
}, matrix(0, nrow = length(bullseyeInferred), ncol=9L))



vapply(1:3, function (x) dists, double(9))
