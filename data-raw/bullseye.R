library('phangorn')
library('phylosim')
library('TreeDist')
library('TreeSearch')
set.seed(1)
nTrees <- 10L
nTips <- c(5L, 10L, 20L, 50L)
treesNames <- paste(nTips, 'tips')
subsamples <- 10:1 * 20

# Generate trees:
bullseyeTrees <- lapply(nTips, function (nTip)
  lapply(seq_len(nTrees), function (XX) ape::rtree(nTip))
)
names(bullseyeTrees) <- treesNames
usethis::use_data(bullseyeTrees, overwrite=TRUE)

# Infer trees:
bullseyeInferred <- vector('list', length(nTips))
names(bullseyeInferred) <- treesNames
kimura <- list(list(K80(Alpha = 1, Beta = 2)))
for (trees in treesNames) {
  theseTrees <- bullseyeTrees[[trees]]
  seqs <- lapply(theseTrees, function (tree) {
    nTip <- length(tree$tip.label)
    states <- Simulate(PhyloSim(root.seq=sampleStates(
      NucleotideSequence(len=2000, proc=kimura)),
      phylo=tree
    ))$alignment
    states[paste0('t', seq_len(nTip)), ]
  })

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
names(bullseyeScores) <- treesNames
for (trees in treesNames) {
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
