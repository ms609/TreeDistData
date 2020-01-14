set.seed(1)
for (i in 1:500) {
  cat('.')
  trees <- lapply(rep(40, 50), ape::rtree, br = NULL)
  # Safest to re-order, as postordering avoids crash in path.dist
  trees <- structure(lapply(trees, TreeTools::Postorder), class='multiPhylo')

  TreeDistData::PairwiseDistances(trees, TreeDist::NNIDist, 3L)

  #phangorn::SPR.dist(trees)

  #tbr <- TBRDist::TBRDist(trees, exact = exact)
}


