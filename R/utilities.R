#' All distances between a pair of trees
#' @param tr1,tr2 Phylogenetic trees of class `phylo`.
#' @author Martin R. Smith
#' @importFrom TreeDist VariationOfPhylogeneticInfo VariationOfMatchingSplitInfo
#' NyeTreeSimilarity MatchingSplitDistance
#' VariationOfClusteringInfo
#' @importFrom Quartet QuartetDivergence QuartetStatus
#' @importFrom phangorn treedist SPR.dist
#' @keywords internal
#' @export
AllDists <- function (tr1, tr2) {
  cat('.')
  c(VariationOfPhylogeneticInfo(tr1, tr2, normalize=TRUE),
    VariationOfMatchingSplitInfo(tr1, tr2, normalize=TRUE),
    VariationOfClusteringInfo(tr1, tr2, normalize=TRUE),
    QuartetDivergence(QuartetStatus(tr1, tr2), similarity = FALSE),
    1 - NyeTreeSimilarity(tr1, tr2, normalize=TRUE),
    MatchingSplitDistance(tr1, tr2),
    treedist(tr1, tr2),
    SPR.dist(tr1, tr2)
  )
}

#' All distances between each pair of trees
#' @param trees List of bifurcating trees of class `phylo`.
#' @author Martin R. Smith
#' @importFrom TreeDist VariationOfPhylogeneticInfo VariationOfMatchingSplitInfo
#' NyeTreeSimilarity MatchingSplitDistance
#' VariationOfClusteringInfo
#' @importFrom Quartet ManyToManyQuartetAgreement
#' @importFrom phangorn treedist SPR.dist
#' @keywords internal
#' @export
CompareAllTrees <- function (trees) {
  elementStatus <- ManyToManyQuartetAgreement(trees)
  qd <- elementStatus[, , 'd'] / elementStatus[1, 1, 's']

  treeDists <- vapply(trees, function (tr1) vapply(trees, function (tr2) {
    c(treedist(tr1, tr2)[c('symmetric.difference', 'path.difference')],
      SPR.dist(tr1, tr2))
  }, double(3)), matrix(0, nrow=3, ncol=length(trees)))

  list(
    vpi = VariationOfPhylogeneticInfo(trees, trees, normalize=TRUE),
    vmsi = VariationOfMatchingSplitInfo(trees, trees, normalize=TRUE),
    vci = VariationOfClusteringInfo(trees, trees, normalize=TRUE),
    qd = qd,
    nts = 1 - NyeTreeSimilarity(trees, trees, normalize=TRUE),
    msd = MatchingSplitDistance(trees, trees),
    rf = treeDists['symmetric.difference', , ],
    path = treeDists['path.difference', , ],
    spr = treeDists['spr', , ]
  )
}
