#' All distances between a pair of trees
#'
#' `AllDists()` calculates the distances between two trees, using a suite
#' of distance measures.
#'
#' @templateVar nni_t TRUE
#' @template allDistMethods
#' @template methodRefs
#'
#' @param tr1,tr2 Phylogenetic trees of class `phylo`.
#' @param verbose Logical specifying whether to print messages allowing the
#' user to keep track of progress.
#'
#' @template MRS
#' @importFrom TreeDist MASTSize NNIDist SPRDist PathDist
#' JaccardRobinsonFoulds
#' DifferentPhylogeneticInfo MatchingSplitInfoDistance
#' NyeSimilarity MatchingSplitDistance
#' ClusteringInfoDistance RobinsonFoulds InfoRobinsonFoulds
#' @importFrom Quartet QuartetDivergence QuartetStatus
#' @importFrom TBRDist TBRDist
#' @family pairwise tree distances
#'
#' @examples
#' library('TreeTools', quietly = TRUE, warn.conflict = FALSE)
#' AllDists(BalancedTree(8), PectinateTree(8))
#' AllDists(list(BalancedTree(6), PectinateTree(6)), BalancedTree(6))
#'
#' @export
AllDists <- function (tr1, tr2, verbose = FALSE) {
  if (verbose) cat('q') # nocov
  qd <- QuartetDivergence(QuartetStatus(tr1, tr2), similarity = FALSE)

  if (verbose) cat('m') # nocov
  mast <- MASTSize(tr1, tr2, rooted = FALSE)
  masti <-  LnUnrooted(mast) / log(2)
  attributes(masti) <- attributes(mast)

  if (verbose) cat('n') # nocov
  nni <- NNIDist(tr1, tr2)
  if (verbose) cat('t') # nocov
  tbr <- TBRDist(tr1, tr2)

  NNIPart <- function (name) {
    if (is.null(names(nni))) nni[name, ] else nni[[name]]
  }
  if (verbose) cat('s') # nocov
  spr <- SPRDist(tr1, tr2)
  if (!is.null(names(spr))) spr <- spr[['spr']]

  if (verbose) cat('.') # nocov
  Bind <- if (is.null(names(nni))) rbind else c
  Bind(
    pid = DifferentPhylogeneticInfo(tr1, tr2, normalize = TRUE),
    msid = MatchingSplitInfoDistance(tr1, tr2, normalize = TRUE),
    cid = ClusteringInfoDistance(tr1, tr2, normalize = TRUE),
    qd = unname(qd),
    nye = NyeSimilarity(tr1, tr2, similarity = FALSE, normalize = TRUE),

    jnc2 = JaccardRobinsonFoulds(tr1, tr2, k = 2, allowConflict = FALSE,
                                 normalize = TRUE),
    jnc4 = JaccardRobinsonFoulds(tr1, tr2, k = 4, allowConflict = FALSE,
                                 normalize = TRUE),
    jco2 = JaccardRobinsonFoulds(tr1, tr2, k = 2, allowConflict = TRUE,
                                 normalize = TRUE),
    jco4 = JaccardRobinsonFoulds(tr1, tr2, k = 4, allowConflict = TRUE,
                                 normalize = TRUE),

    ms = MatchingSplitDistance(tr1, tr2),
    mast = mast,
    masti = masti,
    nni_l = NNIPart('lower'),
    nni_L = NNIPart('best_lower'),
    nni_t = NNIPart('tight_upper'),
    nni_U = NNIPart('best_upper'),
    nni_u = NNIPart('loose_upper'),
    spr = spr,
    tbr_l = tbr$tbr_min,
    tbr_u = tbr$tbr_max,
    rf = RobinsonFoulds(tr1, tr2),
    icrf = InfoRobinsonFoulds(tr1, tr2),
    path = PathDist(tr1, tr2)
  )
}

#' Distances between each pair of trees
#'
#' @param trees List of trees of class `phylo`.
#' @param Func Function returning distance between two trees.
#' @param valueLength Integer specifying expected length of the value returned
#' by `Func`.
#' @param \dots Additional arguments to `Func`.
#'
#' @return Matrix detailing distance between each pair of trees.
#' Identical trees are assumed to have zero distance.
#' @template MRS
#' @family pairwise tree distances
#' @importFrom stats as.dist
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
#' NyeSimilarity MatchingSplitDistance MASTSize SPRDist
#' ClusteringInfoDistance RobinsonFoulds InfoRobinsonFoulds
#' @importFrom Quartet ManyToManyQuartetAgreement
#'
#' @examples
#' trees <- lapply(rep(8, 5), TreeTools::RandomTree, root = TRUE)
#' CompareAllTrees(trees)
#' @template MRS
#' @family pairwise tree distances
#' @export
CompareAllTrees <- function (trees, exact = FALSE, slow = TRUE,
                             verbose = FALSE) {
  MSG <- function (...) if (verbose) message(Sys.time(), ': ', ...)

  # Re-order once; will happen when calling path.dist and SPR.dist
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
  nni <- PairwiseDistances(trees, NNIDist, 7L)

  MSG('SPR')
  sprDist <- as.matrix(SPRDist(trees))

  MSG('TBR')
  tbr <- TBRDist(trees, exact = exact)
  tbr <- if (exact) {
    list(tbr = as.matrix(tbr))
  } else {
    lapply(tbr, as.matrix)
  }

  MSG('path')
  pathDist <- as.matrix(PathDist(trees))

  MSG('PID')
  pid <- DifferentPhylogeneticInfo(splits, normalize = TRUE)

  MSG('msid')
  msid <- MatchingSplitInfoDistance(splits, normalize = TRUE)

  MSG('cid')
  cid <- ClusteringInfoDistance(splits, normalize = TRUE)

  MSG('Nye')
  nye <- 1 - NyeSimilarity(splits, normalize = TRUE)

  MSG('MSD')
  ms <- MatchingSplitDistance(splits)

  MSG('Complete; listing.')
  list(
    pid = pid,
    msid = msid,
    cid = cid,
    qd = qd,
    nye = nye,

    jnc2 =  JaccardRobinsonFoulds(splits, k = 2, allowConflict = FALSE,
                                  normalize = TRUE),
    jnc4 =  JaccardRobinsonFoulds(splits, k = 4, allowConflict = FALSE,
                                  normalize = TRUE),
    jco2 = JaccardRobinsonFoulds(splits, k = 2, allowConflict = TRUE,
                                 normalize = TRUE),
    jco4 = JaccardRobinsonFoulds(splits, k = 4, allowConflict = TRUE,
                                 normalize = TRUE),

    ms = ms,
    mast = mast,
    masti = masti,

    nni_l = nni$lower,
    nni_L = nni$best_lower,
    nni_t = nni$tight_upper,
    nni_U = nni$best_upper,
    nni_u = nni$loose_upper,

    spr = sprDist,
    tbr_l = tbr$tbr_min,
    tbr_u = tbr$tbr_max,

    rf = RobinsonFoulds(splits),
    icrf = InfoRobinsonFoulds(splits),
    path = pathDist
  )
}

#' Select colour from palette
#'
#' @param method Character specifying acronym for method: one of [`tdMethods`].
#' @param opacity Character specifying hex code for opacity; `"FF"` = opaque.
#'
#' @return `TreeDistCol()` returns a hex code for the colour matching the
#' specified method.
#'
#' @export
TreeDistCol <- function (method, opacity = '') {
  unset <- is.na(TreeDistData::tdCol[method])
  if (any(unset)) {
    warning("No colour set for ", method[unset])
  }
  paste0(TreeDistData::tdCol[method], opacity)
}

#' Tabulate method outputs
#'
#' Helper function to tabulate sortable data in package vignettes.
#'
#' @param Table `DT::DataTable()` or an equivalent tabulation function.
#' @param dat Matrix or vector of data to be tabulated, with (row) names
#' corresponding to the methods employed.
#'
#' @return `.TDDTable()` produces a table plotted using `Table()`.
#' @template MRS
#' @keywords internal
#' @export
.TDDTable <- function (Table, dat, ...) {
  method <- rownames(dat)
  if (is.null(method)) method <- names(dat)
  method[method == 'Nye _et al._'] <- 'Nye <i>et al.</i>'
  dat <- as.matrix(cbind(Method = method, dat))
  rownames(dat) <- NULL
  Table(dat, options = list(paging = FALSE, searching = FALSE, info = FALSE),
        escape = FALSE, ...)
}
