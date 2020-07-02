library('TreeTools', quietly = TRUE, warn.conflict = FALSE)
library('TreeDist')

nTip <- 11L
tipLabel <- seq_len(nTip)
pectinateTree <- PectinateTree(tipLabel)

RNGversion("3.6.0")
set.seed(0)
repls <-  100000L
message("Generating ", repls, " random trees:")
randomTreeIds <- unique(runif(repls * 1.5, 0, NUnrooted(nTip) - 1L))[seq_len(repls)]
message("> Found ", repls, ' unique integers. Arboricating...')
randomTrees <- as.phylo(randomTreeIds, nTip, tipLabel)
message("> Converted integers to trees. Sorting...")
randomTrees <- structure(lapply(randomTrees, Postorder), class = 'multiPhylo')
message("> Sorted trees into postorder.\n")

message('Calculating info... ')
pid <- DifferentPhylogeneticInfo(pectinateTree, randomTrees, normalize = TRUE)
msid <- MatchingSplitInfoDistance(pectinateTree, randomTrees, normalize = TRUE)
cid <- ClusteringInfoDistance(pectinateTree, randomTrees, normalize = TRUE)
message('QD... ')
qd <- Quartet::QuartetDivergence(
  Quartet::QuartetStatus(randomTrees, cf = pectinateTree),
  similarity = FALSE)
message('NTS... ')
nye <- 1 - NyeSimilarity(pectinateTree, randomTrees, normalize = TRUE)
message('JRF... ')
jnc2 <- JaccardRobinsonFoulds(pectinateTree, randomTrees, normalize = TRUE,
                            k = 2, allowConflict = FALSE)
jco2 <- JaccardRobinsonFoulds(pectinateTree, randomTrees, normalize = TRUE,
                            k = 2, allowConflict = TRUE)
jnc4 <- JaccardRobinsonFoulds(pectinateTree, randomTrees, normalize = TRUE,
                            k = 4, allowConflict = FALSE)
jco4 <- JaccardRobinsonFoulds(pectinateTree, randomTrees, normalize = TRUE,
                            k = 4, allowConflict = TRUE)

message('MSD... ')
ms <- MatchingSplitDistance(pectinateTree, randomTrees)
message('MAST... ')
mast <- MASTSize(pectinateTree, randomTrees, rooted = FALSE)
message('NNI... ')
nni <- NNIDist(pectinateTree, randomTrees)[c('lower', 'best_lower',
                                             'tight_upper', 'best_upper',
                                             'loose_upper'), ]
rownames(nni) <- c('nni_l', 'nni_L', 'nni_t', 'nni_U', 'nni_u')
message('SPR... ')
spr <- phangorn::SPR.dist(pectinateTree, randomTrees)
message('TBR... ')
tbr <- TBRDist::TBRDist(pectinateTree, randomTrees)
message('MAFi... ')
mafi <- TBRDist::MAFInfo(pectinateTree, randomTrees)
message('RF... ')
rf <- RobinsonFoulds(pectinateTree, randomTrees)
message('ICRF... ')
icrf <- InfoRobinsonFoulds(pectinateTree, randomTrees)
message('Path... ')
path <- phangorn::path.dist(pectinateTree, randomTrees)

pectinateDistances11 <- rbind(pid = pid, msid = msid, cid = cid, qd = qd,
                              nye = nye,
                              jnc2 = jnc2, jnc4 = jnc4,
                              jco2 = jco2, jco4 = jco4,
                              ms = ms,
                              mast = mast, masti = LnUnrooted(mast) / log(2),
                              nni, spr = spr,
                              tbr_l = tbr$tbr_min, tbr_u = tbr$tbr_max,
                              rf = rf, icrf = icrf, path = path, mafi = mafi)


message("Calculated and bound.")

normalizers <- c(pid = 1, msid = 1, cid = 1, qd = 1, nye = 1, ms = 1,
                 mast = nTip, masti = LnUnrooted.int(nTip) / log(2),
                 nni_l = 18, nni_t = 18, nni_u = 18,
                 spr = 18, tbr_l = 18, tbr_u = 18,
                 rf = 1, icrf = NA, path = 1, mafi = 1)

usethis::use_data(pectinateDistances11, compress = 'xz', overwrite = TRUE)
message("Complete.")
