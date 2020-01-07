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
#' @param valueLength Integer specifying expected length of the value returned
#' by `Func`.
#' @return Matrix detailing distance between each pair of trees.
#' Identical trees are assumed to have zero distance.
#' @template MRS
#' @family pairwise tree distances
#' @export
PairwiseDistances <- function (trees, Func, valueLength = 1L, ...) {
  ret <- array(0, c(length(trees), length(trees), valueLength))
  for (i in seq_along(trees)) {
    trI <- trees[[i]]
    for (j in i + seq_len(length(trees) - i)) {
      val <- Func(trI, trees[[j]], ...)
      ret[j, i, ] <- unlist(val)
    }
  }

  # Return:
  if (valueLength > 1L) {
    structure(lapply(seq_len(valueLength), function (i) {
      as.dist(ret[, , i], upper = TRUE)
    }), names = names(val))
  } else {
    as.dist(ret[, , 1], upper = TRUE)
  }
}

#' All distances between each pair of trees
#'
#' @param trees List of bifurcating trees of class `phylo`.
#' @param exact Logical specifying whether to calculate exact rearrangement
#' distances.
#' @author Martin R. Smith
#' @importFrom TreeTools as.Splits Postorder LnUnrooted
#' @importFrom TBRDist TBRDist USPRDist
#' @importFrom TreeDist VariationOfPhylogeneticInfo VariationOfMatchingSplitInfo
#' NyeTreeSimilarity MatchingSplitDistance MASTSize
#' VariationOfClusteringInfo RobinsonFoulds
#' @importFrom phangorn path.dist SPR.dist
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
    vpi = VariationOfPhylogeneticInfo(splits, normalize = TRUE),
    vmsi = VariationOfMatchingSplitInfo(splits, normalize = TRUE),
    vci = VariationOfClusteringInfo(splits, normalize = TRUE),
    qd = qd,
    nts = 1 - NyeTreeSimilarity(splits, normalize = TRUE),
    msd = MatchingSplitDistance(splits),
    rf = RobinsonFoulds(splits),
    path = as.matrix(path.dist(trees)),
    mast = PairwiseDistances(trees, MASTSize),
    nni = PairwiseDistances(trees, NNIDist, 3L),
    spr = as.matrix(SPR.dist(trees))
  ), tbr)
}
