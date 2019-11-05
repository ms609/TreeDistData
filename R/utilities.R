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

#' Distances between each pair of trees
#'
#' @param trees List of trees of class `phylo`.
#' @param Func distance function returning distance between two trees,
#' e.g. [phangorn::treedist][path.dist].
#' @return Matrix detailing distance between each pair of trees.
#' Identical trees are assumed to have zero distance.
#' @author Martin R. Smith
#' @family pairwise tree distances
#' @export
PairwiseDistances <- function (trees, Func) {
  ret <- matrix(0, length(trees), length(trees))
  for (i in seq_along(trees)) {
    trI <- trees[[i]]
    for (j in i + seq_len(length(trees) - i)) {
      ret[i, j] <- Func(trI, trees[[j]])
    }
  }
  ret[lower.tri(ret)] <- t(ret)[lower.tri(ret)]

  # Return:
  ret
}

#' All distances between each pair of trees
#' @param trees List of bifurcating trees of class `phylo`.
#' @author Martin R. Smith
#' @importFrom TreeTools as.Splits
#' @importFrom TreeDist VariationOfPhylogeneticInfo VariationOfMatchingSplitInfo
#' NyeTreeSimilarity MatchingSplitDistance
#' VariationOfClusteringInfo
#' @importFrom Quartet ManyToManyQuartetAgreement
#' @importFrom phangorn path.dist SPR.dist
#' @family pairwise tree distances
#' @export
CompareAllTrees <- function (trees) {
  elementStatus <- ManyToManyQuartetAgreement(trees)
  qd <- elementStatus[, , 'd'] / elementStatus[1, 1, 's']
  splits <- as.Splits(trees)

  list(
    vpi = VariationOfPhylogeneticInfo(splits, normalize=TRUE),
    vmsi = VariationOfMatchingSplitInfo(splits, normalize=TRUE),
    vci = VariationOfClusteringInfo(splits, normalize=TRUE),
    qd = qd,
    nts = 1 - NyeTreeSimilarity(splits, normalize=TRUE),
    msd = MatchingSplitDistance(splits),
    rf = RobinsonFoulds(splits),
    path = PairwiseDistances(trees, path.dist),
    spr = PairwiseDistances(trees, SPR.dist)
  )
}
