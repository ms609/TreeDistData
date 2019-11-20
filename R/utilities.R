#' All distances between a pair of trees
#' @param tr1,tr2 Phylogenetic trees of class `phylo`.
#' @author Martin R. Smith
#' @importFrom TreeDist VariationOfPhylogeneticInfo VariationOfMatchingSplitInfo
#' NyeTreeSimilarity MatchingSplitDistance
#' VariationOfClusteringInfo RobinsonFoulds
#' @importFrom Quartet QuartetDivergence QuartetStatus
#' @importFrom phangorn treedist mast SPR.dist
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
    length(mast(tr1, tr2, tree=FALSE)),
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
PairwiseDistances <- function (trees, Func, ...) {
  ret <- matrix(0, length(trees), length(trees))
  for (i in seq_along(trees)) {
    trI <- trees[[i]]
    for (j in i + seq_len(length(trees) - i)) {
      ret[i, j] <- Func(trI, trees[[j]], ...)
    }
  }
  ret[lower.tri(ret)] <- t(ret)[lower.tri(ret)]

  # Return:
  ret
}

#' All distances between each pair of trees
#'
#' @param trees List of bifurcating trees of class `phylo`.
#' @param exact Logical specifying whether to calculate exact rearrangement
#' distances.
#' @author Martin R. Smith
#' @importFrom TreeTools as.Splits Postorder
#' @importFrom TBRDist TBRDist SPRDist
#' @importFrom TreeDist VariationOfPhylogeneticInfo VariationOfMatchingSplitInfo
#' NyeTreeSimilarity MatchingSplitDistance
#' VariationOfClusteringInfo RobinsonFoulds
#' @importFrom phangorn path.dist SPR.dist mast
#' @importFrom Quartet ManyToManyQuartetAgreement
#'
#' @examples
#'   trees <- lapply(rep(8, 5), ape::rtree, br = NULL)
#'   CompareAllTrees(trees)
#'
#' @template MRS
#' @family pairwise tree distances
#' @export
CompareAllTrees <- function (trees, exact = FALSE) {
  elementStatus <- ManyToManyQuartetAgreement(trees)
  qd <- elementStatus[, , 'd'] / elementStatus[1, 1, 's']

  splits <- as.Splits(trees)
  if(!inherits(trees, 'multiPhylo')) {
    # Safest to re-order, as postordering avoids crash in path.dist
    trees <- structure(lapply(trees, Postorder), class='multiPhylo')
  }

  tbr <- TBRDist(trees, exact = exact)
  tbr <- if (exact) {
    list(tbr = as.matrix(tbr))
  } else {
    lapply(tbr, as.matrix)
  }

  c(list(
    vpi = VariationOfPhylogeneticInfo(splits, normalize=TRUE),
    vmsi = VariationOfMatchingSplitInfo(splits, normalize=TRUE),
    vci = VariationOfClusteringInfo(splits, normalize=TRUE),
    qd = qd,
    nts = 1 - NyeTreeSimilarity(splits, normalize=TRUE),
    msd = MatchingSplitDistance(splits),
    rf = RobinsonFoulds(splits),
    path = as.matrix(path.dist(trees)),
    mast = PairwiseDistances(trees, function (tree1, tree2)
      length(mast(tree1, tree2, tree = FALSE))),
    spr = as.matrix(SPR.dist(trees)),
    uspr = as.matrix(USPRDist(trees))
  ), tbr)
}
