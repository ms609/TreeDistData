library('TreeTools')
library('TreeDist')

nTip <- 11L
tipLabel <- seq_len(nTip)
pectinateTree <- PectinateTree(tipLabel)

RNGversion("3.5.0")
set.seed(0)
repls <-  100000L
randomTreeIds <- unique(floor(runif(repls * 2) * NUnrooted(nTip)))[seq_len(repls)]
randomTrees <- as.phylo(randomTreeIds, nTip, tipLabel)
randomTrees <- structure(lapply(randomTrees, Postorder), class = 'multiPhylo')
message("Generated ", repls, " comparison trees.")

message('Calculating info... ')
dpi <- DifferentPhylogeneticInfo(pectinateTree, randomTrees, normalize = TRUE)
msid <- MatchingSplitInfoDistance(pectinateTree, randomTrees, normalize = TRUE)
cid <- ClusteringInfoDistance(pectinateTree, randomTrees, normalize = TRUE)
message('QD... ')
qd <- Quartet::QuartetDivergence(
  Quartet::QuartetStatus(randomTrees, cf=pectinateTree),
  similarity = FALSE)
message('NTS... ')
nts <- 1 - NyeTreeSimilarity(pectinateTree, randomTrees, normalize=TRUE)
message('MSD... ')
msd <- MatchingSplitDistance(pectinateTree, randomTrees)
message('MAST... ')
mast <- MASTSize(pectinateTree, randomTrees, rooted = FALSE)
message('NNI... ')
nni <- matrix(unlist(NNIDist(pectinateTree, randomTrees)), nrow = 3,
              dimnames = list(c('nni_l', 'nni_t', 'nni_u'), NULL))
message('SPR... ')
spr <- phangorn::SPR.dist(pectinateTree, randomTrees)
message('TBR... ')
tbr <- TBRDist::TBRDist(pectinateTree, randomTrees)
message('MAFi... ')
mafi <- TBRDist::MAFInfo(pectinateTree, randomTrees)
message('RF... ')
rf <- RobinsonFoulds(pectinateTree, randomTrees)
message('RFI... ')
rfi <- RobinsonFouldsInfo(pectinateTree, randomTrees)
message('Path... ')
path <- phangorn::path.dist(pectinateTree, randomTrees)

pectinateDistances11 <- rbind(dpi = dpi, msid = msid, cid = cid, qd = qd,
                              nts = nts, msd = msd,
                              mast = mast, masti = LnUnrooted(mast) / log(2),
                              nni, spr = spr,
                              tbr_l = tbr$tbr_min, tbr_u = tbr$tbr_max,
                              mafi = mafi, rf = rf, rfi = rfi, path = path)


message("Calculated and bound.")

normalizers <- c(dpi = 1, msid = 1, cid = 1, qd = 1, nts = 1, msd = 1,
                 mast=nTip, masti=LnUnrooted.int(nTip) / log(2),
                 nni_l=18, nni_t = 18, nni_u= 18,
                 spr = 18, tbr_l = 18, tbr_u = 18,
                 mafi = 1,
                 rf = 1, rfi = NA, path = 1)

usethis::use_data(pectinateDistances11, compress='xz', overwrite=TRUE)
message("Used data. Complete.")
