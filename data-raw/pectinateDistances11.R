library('TreeTools')
library('TreeDist')
suppressWarnings(RNGversion("3.5.0")) # Stopgap until we can require R 3.6.0

nTip <- 11L
tipLabel <- seq_len(nTip)
pectinateTree <- PectinateTree(tipLabel)

repls <-  100000L
set.seed(0)
randomTreeIds <- unique(floor(runif(repls * 2) * NUnrooted(nTip)))[seq_len(repls)]
randomTrees <- as.phylo(randomTreeIds, nTip, tipLabel)
randomTrees <- structure(lapply(randomTrees, Postorder), class = 'multiPhylo')
message("Generated ", repls, " comparison trees.")

message('Calculating VoI... ')
vpi <- VariationOfPhylogeneticInfo(pectinateTree, randomTrees, normalize=TRUE)
vmsi <- VariationOfMatchingSplitInfo(pectinateTree, randomTrees, normalize=TRUE)
vci <- VariationOfClusteringInfo(pectinateTree, randomTrees, normalize=TRUE)
message('QD... ')
qd <- Quartet::QuartetDivergence(Quartet::QuartetStatus(randomTrees, cf=pectinateTree), similarity = FALSE)
message('NTS... ')
nts <- 1 - NyeTreeSimilarity(pectinateTree, randomTrees, normalize=TRUE)
message('MSD... ')
msd <- MatchingSplitDistance(pectinateTree, randomTrees)
message('MAST... ')
mast <- MASTSize(pectinateTree, randomTrees)
message('NNI... ')
nni <- NNIDist(pectinateTree, randomTrees)
rownames(nni) <- c('nni_l', 'nni_t', 'nni_u')
message('SPR... ')
spr <- phangorn::SPR.dist(pectinateTree, randomTrees)
message('TBR... ')
tbr <- TBRDist::TBRDist(pectinateTree, randomTrees)
message('MAFi... ')
mafi <- TBRDist::MAFInfo(pectinateTree, randomTrees)
message('RF... ')
rf <- RobinsonFoulds(pectinateTree, randomTrees)
message('Path... ')
path <- phangorn::path.dist(pectinateTree, randomTrees)

pectinateDistances11 <- rbind(vpi = vpi, vmsi = vmsi, vci = vci, qd = qd,
                              nts = nts, mst = msd,
                              mast = mast, masti = LnUnrooted(mast) / log(2),
                              nni, spr = spr,
                              tbr_l = tbr$tbr_min, tbr_u = tbr$tbr_max,
                              mafi = mafi, rf = rf, path = path)


message("Calculated and bound.")

normalizers <- c(vpi = 1, vmsi = 1, vci = 1, qd = 1, nts = 1, msd = 1,
                 mast=nTip, masti=LnUnrooted.int(nTip) / log(2),
                 nni_l=18, nni_t = 18, nni_u= 18,
                 spr = 18, tbr_l =18, tbr_u = 18,
                 mafi = 1,
                 rf = 1, path = 1)

usethis::use_data(pectinateDistances11, compress='xz', overwrite=TRUE)
message("Used data. Complete.")
