library('TreeDist')
suppressWarnings(RNGversion("3.5.0")) # Stopgap until we can require R 3.6.0
set.seed(0)

pectinateTree <- structure(list(edge = structure(c(12L, 13L, 14L, 15L, 16L, 17L,
                                                    17L, 16L, 15L, 14L, 13L, 12L, 18L, 18L, 19L, 19L, 20L, 20L, 21L,
                                                    21L, 13L, 14L, 15L, 16L, 17L, 1L, 2L, 3L, 4L, 5L, 6L, 18L, 7L,
                                                    19L, 8L, 20L, 9L, 21L, 10L, 11L), .Dim = c(20L, 2L)), Nnode = 10L,
                                 tip.label = c("1", "2", "3", "4", "5", "6", "7", "8", "9",
                                               "10", "11")), class = "phylo", order = "cladewise")

repls <-  100000L

tipLabel <- pectinateTree$tip.label
randomTrees <- lapply(seq_len(repls), function (i) ape::rtree(11L, br=NULL, tip.label=tipLabel))
cat("Generated comparison trees.")

pectinateDistances11 <- rbind(
  VariationOfPhylogeneticInfo(pectinateTree, randomTrees, normalize=TRUE),
  VariationOfMatchingSplitInfo(pectinateTree, randomTrees, normalize=TRUE),
  VariationOfClusteringInfo(pectinateTree, randomTrees, normalize=TRUE),
  Quartet::QuartetDivergence(Quartet::QuartetStatus(randomTrees, cf=pectinateTree), similarity = FALSE),
  1 - NyeTreeSimilarity(pectinateTree, randomTrees, normalize=TRUE),
  MatchingSplitDistance(pectinateTree, randomTrees),
  vapply(randomTrees, phangorn::treedist, tree2=pectinateTree, double(2)),
  vapply(randomTrees, phangorn::SPR.dist, tree2=pectinateTree, double(1))
)
cat("Calculated and bound.")

rownames(pectinateDistances11) <- names(c(vpi = 0, vmsi = 0, vci = 0, qd = 0, nts = 0, msd = 0, rf = 0, path = 0,
     spr = 0))

usethis::use_data(pectinateDistances11, compress='xz', overwrite=TRUE)
cat("Used data. Complete.")
