#' Mean distances between random pairs of trees
#'
#' A three-dimensional array listing the normalized distances between 1&nbsp;000
#' random pairs of trees drawn from the uniform distribution using
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
#' Columns list the summary statistics of calculated tree distances: the
#' minimum (`min`),
#' 1%, 5%, 10%, 25%, 50% (i.e. median), 75%, 90%, 95%, 99% percentiles,
#' maximum (`max`), mean (`mean`) and standard deviation (`sd`).
#'
#' The third dimension lists the number of tips in the trees compared.
#'
#'
#' @references
#' - \insertRef{Allen2001}{TreeSearch}
#'
#' - \insertRef{Bogdanowicz2012}{TreeDist}
#'
#' - \insertRef{Nye2006}{TreeDist}
#'
#' - \insertRef{Robinson1981}{TreeDist}
#'
#' - \insertRef{SmithTern}{TreeSearch}
#'
#' - \insertRef{SmithDist}{TreeDist}
#'
#' - \insertRef{Steel1993}{TreeDist}
#' @keywords datasets
"randomTreeDistances"

#' Distances between random pairs of trees
#'
#' Two-dimensional matrices listing the normalized distances between random
#' pairs of bifurcating 11- and 25-tip trees drawn from the uniform distribution
#' using `ape::rtree(nTip, br=NULL)`.
#'
#' `pectinateDistances11` reports distances between a pectinate 11-tip tree
#' and 100&nbsp;000 random bifurcating trees.
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
#' The pairs of 25-tip trees are saved as data object [`randomTreePairs25`].
#'
#' @references
#' \insertRef{Bogdanowicz2012}{TreeDist}
#' \insertRef{Nye2006}{TreeDist}
#' \insertRef{Robinson1981}{TreeDist}
#' \insertRef{SmithTern}{TreeSearch}
#' \insertRef{SmithDist}{TreeDist}
#' \insertRef{Steel1993}{TreeDist}
#'
#' @name distanceDistributions
#' @keywords datasets
NULL

#' @rdname distanceDistributions
"distanceDistribution25"
#' @rdname distanceDistributions
"distanceDistribution11"
#' @rdname distanceDistributions
"pectinateDistances11"

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
#'
#' \insertRef{Kendall2016}{TreeDistData}
#'
#' \insertRef{Nye2006}{TreeDist}
#'
#' \insertRef{Robinson1981}{TreeDist}
#'
#' \insertRef{SmithTern}{TreeSearch}
#'
#' \insertRef{SmithDist}{TreeDist}
#'
#' \insertRef{Steel1993}{TreeDist}
#'
#' @keywords datasets
"sevenTipDistances"

#' Bullseye test results
#'
#' Trees used to implement a 'Bullseye' test, after that proposed by Kuhner
#' and Yamato (2015).
#'
#' `bullseyeTrees` is a list with four elements, named `5 tips`, `10 tips`,
#' `20 tips` and `50 tips`.  Each element contains 1&nbsp;000 trees with _n_
#' tips, randomly sampled from the uniform distribution using `ape::rtree`.
#'
#' `bullXXInferred` is a list with four elements, named as in `bullseyeTrees`.
#' Each element contains 1&nbsp;000 subelements. Each subelement is a list of
#' ten trees, which have been inferred from progressively more degraded datasets,
#' originally simulated from the corresponding tree in `bullseyeTrees`.
#'
#' `bullXXScores` is a list with four elements, named as in `bullseyeTrees`.
#' Each element contains a three dimensional array, in which the first dimension
#' corresponds to the progressive degrees of degradation, labelled according to
#' the number of characters present or the percentage of tokens switched;
#' the second dimension is named with an abbreviation of the tree similarity /
#' distance metric used to score the trees, and the third dimension contains
#' 1&nbsp;000 entries corresponding to the trees in `bullseyeTrees`.
#' Each cell contains the distance between the inferred tree and the generative
#' tree under the stated tree distance metric.
#'
#' The `bullseyeMorph` prefix refers to the 'subsampling' experiment
#' described by Smith (forthcoming); the `bullMoDi` prefix refers to the
#' 'miscoding' experiment.
#'
#' `bullseyeDistances` contains four elements, each tabulating the distance
#' between each pair of _n_-tip trees in `bullseyeTrees`.  For details, see
#' \code{\link{sevenTipDistances}}.
#'
#'
#' @references
#' - \insertRef{Kuhner2015}{TreeDistData}
#'
#' - \insertRef{SmithDist}{TreeDist}
#'
#' @name bullseye

#' @rdname bullseye
'bullseyeTrees'
#' @rdname bullseye
'bullMoDiInferred'
#' @rdname bullseye
'bullMoDiScores'
#' @rdname bullseye
'bullseyeMorphInferred'
#' @rdname bullseye
'bullseyeMorphScores'
#' @rdname bullseye
'bullseyeDistances'
