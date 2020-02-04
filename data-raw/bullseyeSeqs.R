library('phangorn')
library('phylosim')

data("bullseyeTrees", package='TreeDistData') # Generated in bullseyeTrees.R
tipsNames <- names(bullseyeTrees)
nTrees <- length(bullseyeTrees[[1]])

useML <- c(1, 2)
subsamples <- 10:1 * 200
RNGversion('3.6.0')
set.seed(1)

# Infer trees:
bullseyeSeqs <- vector('list', length(tipsNames))
names(bullseyeSeqs) <- tipsNames
kimura <- list(list(K80(Alpha = 1, Beta = 2)))
for (trees in tipsNames) {
  message("== ", trees, " ==\n")
  theseTrees <- bullseyeTrees[[trees]]
  seqs <- lapply(theseTrees, function (tree) {
    nTip <- length(tree$tip.label)
    states <- Simulate(PhyloSim(root.seq=sampleStates(
      NucleotideSequence(len=2000, proc=kimura)),
      phylo = tree
    ))$alignment
    states[paste0('t', seq_len(nTip)), ]
  })
  bullseyeSeqs[[trees]] <- seqs
}
usethis::use_data(bullseyeSeqs, compress = 'xz', overwrite = TRUE)
message("== COMPLETE == ")
