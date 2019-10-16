library('TreeDist')
data('bullseyeTrees', package='TreeDistData')
bullseyeDistances <- lapply(bullseyeTrees, function(x) NULL)

for (tips in names(bullseyeTrees)) {
  message(tips)
  trees <- bullseyeTrees[[tips]][1:100]
  elementStatus <- Quartet::ManyToManyQuartetAgreement(trees)
  qd <- elementStatus[, , 'd'] / elementStatus[1, 1, 1]
  message("Calculated Quartet Distances")

  treeDists <- vapply(trees, function (tr1) vapply(trees, function (tr2) {
    c(phangorn::treedist(tr1, tr2)[c('symmetric.difference', 'path.difference')],
      phangorn::SPR.dist(tr1, tr2))
  }, double(3)), matrix(0, nrow=3, ncol=length(trees)))

  bullseyeDistances[[tips]] <- list(
    vpi = VariationOfPhylogeneticInfo(trees, trees, normalize = TRUE),
    vmsi = VariationOfMatchingSplitInfo(trees, trees, normalize = TRUE),
    vci = VariationOfClusteringInfo(trees, trees, normalize = TRUE),
    qd = qd,
    nts = 1 - NyeTreeSimilarity(trees, trees, normalize = TRUE),
    msd = MatchingSplitDistance(trees, trees),
    rf = treeDists['symmetric.difference', , ],
    path = treeDists['path.difference', , ],
    spr = treeDists['spr', , ]
  )
}

usethis::use_data(bullseyeDistances, compress='xz', overwrite=TRUE)
