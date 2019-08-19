library('phangorn')
library('TreeSearch')
set.seed(0)
nTrees <- 100L
nTips <- c(5L, 10L, 20L, 50L)
treesNames <- paste(nTips, 'tips')
subsamples <- 10:1 * 20

# Generate trees:
bullseyeTrees <- lapply(nTips, function (nTip)
  lapply(seq_len(nTrees), function (XX) ape::rtree(nTip))
)
names(bullseyeTrees) <- treesNames
usethis::use_data(bullseyeTrees, overwrite=TRUE)


# Define functions:
SubSampleSic <- function (dat, n) {
  at <- attributes(dat)
  chosenInd <- at$index[seq_len(n)]
  indexedChosen <- seq_len(max(chosenInd))
  ret <- lapply(dat, function (x) x[indexedChosen])
  attributes(ret) <- at
  attr(ret, 'weight') <- as.integer(table(chosenInd))
  attr(ret, 'index') <- indexedChosen
  attr(ret, 'nr') <- n
  ret
}

SubSample <- function (dat, n) {
  at <- attributes(dat)
  tokens <- t(vapply(seq_along(dat), function(taxon) {
    dat[[taxon]][at$index[seq_len(n)]]
  }, double(n)))
  rownames(tokens) <- names(dat)
  ret <- MatrixToPhyDat(tokens)
  attr(ret, 'levels') <- at$levels
  ret
}


# Infer trees:
bullseyeInferred <- vector('list', length(nTips))
names(bullseyeInferred) <- treesNames
for (trees in treesNames) {
  theseTrees <- bullseyeTrees[[trees]]
  seqs <- lapply(theseTrees, simSeq, l = 2000)

  bullseyeInferred[[trees]] <- lapply(seq_len(nTrees), function (i) {
    tr <- theseTrees[[i]]
    sq <- seqs[[i]]

    lapply(subsamples, function (n) {
      dm <- dist.ml(SubSample(sq, n))
      upgma(dm)
    })
  })
}
usethis::use_data(bullseyeInferred, compress='xz', overwrite=TRUE)

for (trees in treesNames) {
  inferred <- bullseyeInferred[[trees]]
  theseTrees <- bullseyeTrees[[trees]]
  allScores <- vapply(seq_along(inferred), function (i) {
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
  }, matrix(0, nrow = 10L, ncol=9L))




  rhos <- apply(allScores, 3L, function (scores) apply(scores, 2, function (x)
    cor.test(x,  y = seq_len(10L), method='spearman', exact=FALSE)$estimate))

  rowMeans(rhos, na.rm=TRUE)


  pairwise <- apply(allScores, 2L, function (scores) {
    rowSums(scores[10:4, ] > scores[7:1, ])
  })

  apply(pairwise, 1, sum)

  treeMeans <- apply(allScores, 1L, rowMeans)



  par(mar=c(1,1, 1, 1))
  plot(type='n', xlim=c(0, 10), ylim=0:1, axes=F, x=0, y=0)
  lapply(1:9, function (i) {
    iScores <- treeMeans[i, ] / max(treeMeans[i, ])
    lines(iScores, col=i)
  }
  ) -> XX
  legend('topright', legend=rownames(treeMeans), col=1:9, lty=1, bty='n')


}
