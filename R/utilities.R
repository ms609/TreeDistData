#' All distances between a pair of trees
#' @param tr1,tr2 Phylogenetic trees of class `phylo`.
#' @param verbose Logical specifying whether to notify user of progress.
#'
#' @template MRS
#' @importFrom TreeDist MASTSize NNIDist
#' JaccardRobinsonFoulds
#' DifferentPhylogeneticInfo MatchingSplitInfoDistance
#' NyeTreeSimilarity MatchingSplitDistance
#' VariationOfClusteringInfo RobinsonFoulds
#' @importFrom Quartet QuartetDivergence QuartetStatus
#' @importFrom phangorn SPR.dist path.dist
#' @importFrom TBRDist TBRDist
#' @family pairwise tree distances
#' @export
AllDists <- function (tr1, tr2, verbose = FALSE) {
  if (verbose) cat('q')
  qd <- QuartetDivergence(QuartetStatus(tr1, tr2), similarity = FALSE)

  if (verbose) cat('m')
  mast <- MASTSize(tr1, tr2, rooted = FALSE)
  masti <-  LnUnrooted(mast) / log(2)
  attributes(masti) <- attributes(mast)

  if (verbose) cat('r')
  nni <- NNIDist(tr1, tr2)
  tbr <- TBRDist(tr1, tr2)

  if (verbose) cat('.')
  c(
    vpi = DifferentPhylogeneticInfo(tr1, tr2, normalize=TRUE),
    vmsi = MatchingSplitInfoDistance(tr1, tr2, normalize=TRUE),
    vci = VariationOfClusteringInfo(tr1, tr2, normalize=TRUE),
    qd = qd,
    nts = 1 - NyeTreeSimilarity(tr1, tr2, normalize=TRUE),

    ja2 = 1 - JaccardRobinsonFoulds(tr1, tr2, k = 2, arboreal = TRUE, normalize = TRUE),
    ja4 = 1 - JaccardRobinsonFoulds(tr1, tr2, k = 4, arboreal = TRUE, normalize = TRUE),
    jna2 =1 - JaccardRobinsonFoulds(tr1, tr2, k = 2, arboreal = FALSE, normalize = TRUE),
    jna4 =1 - JaccardRobinsonFoulds(tr1, tr2, k = 4, arboreal = FALSE, normalize = TRUE),

    msd = MatchingSplitDistance(tr1, tr2),
    mast = mast,
    masti = masti,
    nni_l = nni[['lower']],
    nni_t = nni[['tight_upper']],
    nni_u = nni[['loose_upper']],
    spr = SPR.dist(tr1, tr2)[['spr']],
    tbr_l = tbr$tbr_min,
    tbr_u = tbr$tbr_max,
    rf = RobinsonFoulds(tr1, tr2),
    path = path.dist(tr1, tr2)
  )
}

#' Distances between each pair of trees
#'
#' @param trees List of trees of class `phylo`.
#' @param Func distance function returning distance between two trees,
#' e.g. [phangorn::treedist][path.dist].
#' @param valueLength Integer specifying expected length of the value returned
#' by `Func`.
#' @param \dots Additional arguments to `Func`.
#'
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
#' @param slow Logical specifying whether to report distance measures that are
#' slow to calculate (quartet distance and maximum agreement subtree).
#' @param verbose Logical specifying whether to report which calculation
#' is presently being performed.
#'
#' @importFrom TreeTools as.Splits Postorder LnUnrooted
#' @importFrom TBRDist TBRDist USPRDist
#' @importFrom TreeDist DifferentPhylogeneticInfo MatchingSplitInfoDistance
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
CompareAllTrees <- function (trees, exact = FALSE, slow = TRUE,
                             verbose = FALSE) {
  MSG <- function (...) if (verbose) message(Sys.time(), ': ', ...)

  # Safest to re-order, as postordering avoids crash in path.dist and SPR.dist
  trees <- structure(lapply(trees, Postorder), class='multiPhylo')

  splits <- as.Splits(trees)

  if (slow) {
    MSG('QD')
    elementStatus <- ManyToManyQuartetAgreement(trees)
    qd <- elementStatus[, , 'd'] / elementStatus[1, 1, 's']

    MSG('MAST')
    mast <- as.matrix(PairwiseDistances(trees, MASTSize, rooted = FALSE))
    diag(mast) <- length(trees[[1]]$tip.label)
    masti <-  LnUnrooted(mast) / log(2)
    attributes(masti) <- attributes(mast)
  } else {
    qd <- NULL

    mast <- masti <- NULL
  }

  MSG('NNI')
  nni <- PairwiseDistances(trees, NNIDist, 3L)

  MSG('SPR')
  sprDist <- as.matrix(SPR.dist(trees))

  MSG('TBR')
  tbr <- TBRDist(trees, exact = exact)
  tbr <- if (exact) {
    list(tbr = as.matrix(tbr))
  } else {
    lapply(tbr, as.matrix)
  }

  MSG('path')
  pathDist <- as.matrix(path.dist(trees))

  MSG('VpI')
  vpi <- DifferentPhylogeneticInfo(splits, normalize = TRUE)

  MSG('VmsI')
  vmsi <- MatchingSplitInfoDistance(splits, normalize = TRUE)

  MSG('VcI')
  vci <- VariationOfClusteringInfo(splits, normalize = TRUE)

  MSG('Nye')
  nts <- 1 - NyeTreeSimilarity(splits, normalize = TRUE)

  MSG('MSD')
  msd <- MatchingSplitDistance(splits)

  MSG('RF')
  rf <- RobinsonFoulds(splits)

  MSG('Complete; listing.')
  list(
    vpi = vpi,
    vmsi = vmsi,
    vci = vci,
    qd = qd,
    nts = nts,

    ja2 =  JaccardRobinsonFoulds(splits, k = 2, arboreal = TRUE, normalize = TRUE),
    ja4 =  JaccardRobinsonFoulds(splits, k = 4, arboreal = TRUE, normalize = TRUE),
    jna2 = JaccardRobinsonFoulds(splits, k = 2, arboreal = FALSE, normalize = TRUE),
    jna4 = JaccardRobinsonFoulds(splits, k = 4, arboreal = FALSE, normalize = TRUE),

    msd = msd,
    mast = mast,
    masti = masti,
    nni_l = nni$lower,
    nni_t = nni$tight_upper,
    nni_u = nni$loose_upper,
    spr = sprDist,
    tbr_l = tbr$tbr_min,
    tbr_u = tbr$tbr_max,
    rf = rf,
    path = pathDist
  )
}
