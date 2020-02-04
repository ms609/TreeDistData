library('phangorn')
library('phylosim')
library('TreeDist')
library('TreeSearch')


data("bullseyeTrees", package='TreeDistData') # Generated in bullseyeTrees.R
data('bullseyeSeqs', package = 'TreeDistData') # Generated in bullseyeSeqs.R
tipsNames <- names(bullseyeTrees)
nTrees <- length(bullseyeTrees[[1]])

subsamples <- 10:1 * 200
RNGversion('3.6.0')
set.seed(1)


bullseyeInferred <- vector('list', length(nTips))
names(bullseyeInferred) <- tipsNames
for (trees in tipsNames[seq_len(useML)]) {
  seqs <- bullseyeSeqs[[trees]]
  theseTrees <- bullseyeTrees[[trees]]
  bullseyeInferred[[trees]] <- lapply(seq_len(nTrees), function (i) {
    tr <- theseTrees[[i]]
    sq <- seqs[[i]]

    lapply(subsamples, function (n) {
      fitStart <- pml(unroot(tr), MatrixToPhyDat(sq[, seq_len(n)]))
      optim.pml(fitStart, model="K80", optNni=TRUE, rearrangement="stochastic",
                control=pml.control(epsilon=1e-06, trace = 0))$tree
    })
  })
}

for (trees in tipsNames[-seq_len(useML)]) {
  seqs <- bullseyeSeqs[[trees]]
  bullseyeInferred[[trees]] <- lapply(seq_len(nTrees), function (i) {
    tr <- theseTrees[[i]]
    sq <- seqs[[i]]

    lapply(subsamples, function (n) {
      dists <- dist.dna(as.DNAbin(sq[, seq_len(n)]))
      if (any(is.nan(dists)) || any(is.infinite(dists)) || any(is.na(dists))) {
        dists <- dist.ml(MatrixToPhyDat(sq[, seq_len(n)]))
      }
      upgma(dists)
    })
  })
}
usethis::use_data(bullseyeInferred, compress='xz', overwrite=TRUE)

bullseyeScores <- vector('list', length(nTips))
names(bullseyeScores) <- tipsNames
for (trees in tipsNames) {
  inferred <- bullseyeInferred[[trees]]
  theseTrees <- bullseyeTrees[[trees]]
  theseScores <- vapply(seq_along(inferred), function (i) {
    trueTree <- theseTrees[[i]]
    rootTip <- trueTree$tip.label[1]
    tr <- root(trueTree, rootTip, resolve.root=TRUE)
    tr$edge.length  <- NULL
    trs <- lapply(inferred[[i]], root, rootTip, resolve.root=TRUE)

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
  }, matrix(0, nrow = 10L, ncol=9L, dimnames=list(
    10:1 * 200,
    c('vpi', 'vmsi', 'vci', 'qd', 'nts', 'msd', 'rf', 'path', 'spr')
  )))
  bullseyeScores[[trees]] <- theseScores
}
usethis::use_data(bullseyeScores, compress='xz', overwrite=TRUE)
