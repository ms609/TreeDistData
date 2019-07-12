#' Mean distances between random pairs of trees
#'
#' A three-dimensional array listing the normalized distances between 100 random
#' pairs of trees drawn from the uniform distribution using
#' `ape::rtree(nTip, br=NULL)`.
#'
#' Rows are named with an abbreviation of the tree comparison metric.
#'
#' Columns list the mean and standard deviation of calculated tree distances.
#'
#' The third dimension lists the number of tips in the trees compared.
#'
#' @keywords datasets
"randomTreeDistances"

#' Distances between random pairs of 25-tip trees
#'
#' A two-dimensional matrix listing the normalized distances between random
#' pairs of 25-tip trees drawn from the uniform distribution using
#' `ape::rtree(nTip, br=NULL)`.
#'
#' Rows are named with an abbreviation of the tree comparison metric.
#' * #TODO list these and state how normalized.
#'
#' Each column lists the calculated distances between each pair of trees.
#'
#' The pairs of trees are saved as data object [`randomTreePairs25`].
#'
#' @keywords datasets
"distanceDistribution25"

#' Pairs of random 25-tip trees
#'
#' A list of 10 000 25-tip trees drawn from the uniform distribution using
#' `ape::rtree(nTip, br=NULL)`.
#'
#' The distances between these pairs of trees are recorded in
#' the data object [`distanceDistribution25`].
#'
#' @keywords datasets
"randomTreePairs25"



#' Distances between unrooted seven-tip trees
#'
#' Distances between each possible pairing of the 945 unrooted seven-tip trees
#' (equivalent to rooted 6-tip trees).  Following Kendall and Colijn (2016).
#'
#' @references
#' \insertRef{Kendall2016}{TreeDistData}
#'
#' @keywords datasets
"sevenTipDistances"
