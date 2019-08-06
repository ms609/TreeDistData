#' Mean distances between random pairs of trees
#'
#' A three-dimensional array listing the normalized distances between 100 random
#' pairs of trees drawn from the uniform distribution using
#' `ape::rtree(nTip, br=NULL)`.
#'
#' Normalization is against the maximum possible value obtainable on a pair
#' of trees of the shapes given, with the exception of the SPR distance,
#' which is normalized against the upper bound (`spr`) and lower bound (`sprLB`)
#' of its diameter on a tree with _n_ tips.  The path and matching split
#' distances are not normalized.
#'
#' Rows are named with an abbreviation of the tree comparison metric:
#'
#' - `vpi`: Variation of Phylogenetic Information (Smith, forthcoming)
#' - `vmsi`: Variation of Matching Split Information (Smith, forthcoming)
#' - `vci`: Variation of Clustering Information (Smith, forthcoming)
#' - `qd`: Quartet divergence (Smith, 2019)
#' - `nts`: Nye _et al._ tree similarity (Nye _et al._ 2006)
#' - `msd`: Matching Split Distance (Bogdanowicz & Giaro, 2012)
#' - `rf`: Robinson-Foulds distance (Robinson & Foulds 1981)
#' - `path`: Path distance (Steel & Penny 1993)
#' - `spr`: SPR distance, normalized against upper bound of diameter (Allen & Steel, 2001)
#' - `sprLB`: SPR distance, normalized against lower bound of diameter (Allen & Steel, 2001)
#'
#' Columns list the mean and standard deviation of calculated tree distances.
#'
#' The third dimension lists the number of tips in the trees compared.
#'
#'
#' @references
#' \insertRef{Allen2001}{TreeSearch}
#' \insertRef{Bogdanowicz2012}{TreeDist}
#' \insertRef{Nye2006}{TreeDist}
#' \insertRef{Robinson1981}{TreeDist}
#' \insertRef{SmithTern}{TreeSearch}
#' \insertRef{SmithDist}{TreeDist}
#' \insertRef{Steel1993}{TreeDist}
#' @keywords datasets
"randomTreeDistances"

#' Distances between random pairs of 25-tip trees
#'
#' A two-dimensional matrix listing the normalized distances between random
#' pairs of 25-tip trees drawn from the uniform distribution using
#' `ape::rtree(nTip, br=NULL)`.
#'
#' Rows are named with an abbreviation of the tree comparison metric.
#' Variation of information measures are normalized against the maximum
#' possible variation of information for trees of the corresponding topologies.
#' The quartet distance and Nye _et al._ measures are normalized against their
#' maximum possible values.  The remaining measures are unnormalized.
#'
#' - `vpi`: Variation of Phylogenetic Information (Smith, forthcoming)
#' - `vmsi`: Variation of Matching Split Information (Smith, forthcoming)
#' - `vci`: Variation of Clustering Information (Smith, forthcoming)
#' - `qd`: Quartet divergence (Smith, 2019)
#' - `nts`: Nye _et al._ tree similarity (Nye _et al._ 2006)
#' - `msd`: Matching Split Distance (Bogdanowicz & Giaro, 2012), unnormalized
#' - `rf`: Robinson-Foulds distance (Robinson & Foulds 1985), unnormalized
#' - `path`: Path distance (Steel & Penny 1993), unnormalized
#' - `spr`: SPR distance, unnormalized
#'
#' Each column lists the calculated distances between each pair of trees.
#'
#' The pairs of trees are saved as data object [`randomTreePairs25`].
#'
#' @references
#' \insertRef{Bogdanowicz2012}{TreeDist}
#' \insertRef{Nye2006}{TreeDist}
#' \insertRef{Robinson1981}{TreeDist}
#' \insertRef{SmithTern}{TreeSearch}
#' \insertRef{SmithDist}{TreeDist}
#' \insertRef{Steel1993}{TreeDist}
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
## <!--Text copied from distanceDistribution25 above -->
#' The list entries are named with an abbreviation of the tree comparison metric.
#' Variation of information measures are normalized against the maximum
#' possible variation of information for trees of the corresponding topologies.
#' The quartet distance and Nye _et al._ measures are normalized against their
#' maximum possible values.  The remaining measures are unnormalized.
#'
#' - `vpi`: Variation of Phylogenetic Information (Smith, forthcoming)
#' - `vmsi`: Variation of Matching Split Information (Smith, forthcoming)
#' - `vci`: Variation of Clustering Information (Smith, forthcoming)
#' - `qd`: Quartet divergence (Smith, 2019)
#' - `nts`: Nye _et al._ tree similarity (Nye _et al._ 2006)
#' - `msd`: Matching Split Distance (Bogdanowicz & Giaro, 2012), unnormalized
#' - `rf`: Robinson-Foulds distance (Robinson & Foulds 1985), unnormalized
#' - `path`: Path distance (Steel & Penny 1993), unnormalized
#' - `spr`: SPR distance, unnormalized
#'
#' Each item in the list contains a 945&times;945 matrix reporting the distance
#' between each pair of seven-tip trees.
#'
#' The final entry of the list is `shapes`, whose value is an integer vector.
#' Each unique tree topology is represented by a distinct vector, whose digits
#' denote the number of nodes in the tree that possess 6, 5, 4, 3, 2 and 1
#' descendant nodes.
#'
#' @references
#' \insertRef{Bogdanowicz2012}{TreeDist}
#' \insertRef{Kendall2016}{TreeDistData}
#' \insertRef{Nye2006}{TreeDist}
#' \insertRef{Robinson1981}{TreeDist}
#' \insertRef{SmithTern}{TreeSearch}
#' \insertRef{SmithDist}{TreeDist}
#' \insertRef{Steel1993}{TreeDist}
#'
#' @keywords datasets
"sevenTipDistances"
