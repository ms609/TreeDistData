#' All distances between a pair of trees
#' @param tr1,tr2 Phylogenetic trees of class `phylo`.
#' @author Martin R. Smith
#' @importFrom TreeDist VariationOfPhylogeneticInfo VariationOfMatchingSplitInfo
#' NyeTreeSimilarity MatchingSplitDistance
#' VariationOfClusteringInfo RobinsonFoulds
#' @importFrom Quartet QuartetDivergence QuartetStatus
#' @importFrom phangorn treedist SPR.dist
#' @family pairwise tree distances
#' @export
AllDists <- function (tr1, tr2) {
  cat('.')
  c(VariationOfPhylogeneticInfo(tr1, tr2, normalize=TRUE),
    VariationOfMatchingSplitInfo(tr1, tr2, normalize=TRUE),
    VariationOfClusteringInfo(tr1, tr2, normalize=TRUE),
    QuartetDivergence(QuartetStatus(tr1, tr2), similarity = FALSE),
    1 - NyeTreeSimilarity(tr1, tr2, normalize=TRUE),
    MatchingSplitDistance(tr1, tr2),
    RobinsonFoulds(tr1, tr2),
    path.dist(tr1, tr2),
    SPR.dist(tr1, tr2)
  )
}

#' All distances between each pair of trees
#' @param trees List of bifurcating trees of class `phylo`.
#' @author Martin R. Smith
#' @importFrom TreeTools as.Splits
#' @importFrom TreeDist VariationOfPhylogeneticInfo VariationOfMatchingSplitInfo
#' NyeTreeSimilarity MatchingSplitDistance
#' VariationOfClusteringInfo RobinsonFoulds
#' @importFrom Quartet ManyToManyQuartetAgreement
#' @importFrom phangorn path.dist SPR.dist
#' @family pairwise tree distances
#' @export
CompareAllTrees <- function (trees) {
  elementStatus <- ManyToManyQuartetAgreement(trees)
  qd <- elementStatus[, , 'd'] / elementStatus[1, 1, 's']

  splits <- as.Splits(trees)
  if(!inherits(trees, 'multiPhylo')) {
    trees <- structure(trees, class='multiPhylo')
  }

  list(
    vpi = VariationOfPhylogeneticInfo(splits, normalize=TRUE),
    vmsi = VariationOfMatchingSplitInfo(splits, normalize=TRUE),
    vci = VariationOfClusteringInfo(splits, normalize=TRUE),
    qd = qd,
    nts = 1 - NyeTreeSimilarity(splits, normalize=TRUE),
    msd = MatchingSplitDistance(splits),
    rf = RobinsonFoulds(splits),
    path = as.matrix(path.dist(trees)),
    spr = as.matrix(SPR.dist(trees))
  )
}
