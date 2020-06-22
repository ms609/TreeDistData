#' Mean distances between random pairs of trees
#'
#' A three-dimensional array listing the normalized distances between 1&nbsp;000
#' random pairs of trees drawn from the uniform distribution using
#' `ape::rtree(nTip, br = NULL)`.
#'
#' Normalization is against the maximum possible value obtainable on a pair
#' of trees of the shapes given, with the exception of the SPR distance,
#' which is normalized against the upper bound (`spr`) and lower bound (`sprLB`)
#' of its diameter on a tree with _n_ leaves.  The path and matching split
#' distances are not normalized.
#'
#' Rows are named with an abbreviation of the tree comparison metric:
#'
#' - `pid`: Phylogenetic Information Distance (Smith, forthcoming)
#' - `msid`: Matching Split Information Distance (Smith, forthcoming)
#' - `cid`: Clustering Information Distance (Smith, forthcoming)
#' - `qd`: Quartet divergence (Smith 2019)
#' - `nye`: Nye _et al._ tree similarity (Nye _et al._ 2006)
#' - `jnc2`, `jnc4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#'   conflicting pairings prohibited ('no-conflict')
#' - `jco2`, `jco4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#'   conflicting pairings permitted ('conflict-ok')
#' - `ms`: Matching Split Distance (Bogdanowicz & Giaro 2012)
#' - `mast`: Size of Maximum Agreement Subtree (Valiente 2009)
#' - `masti`: Information content of Maximum Agreement Subtree
#' - `nni_l`, `nni_u`: Lower and upper bounds for nearest-neighbour interchange
#'   distance (Li _et al._ 1996)
#' - `spr`: Approximate subtree prune and regraft distance
#' - `tbr_l`, `tbr_u`: Lower and upper bound for tree bisection and reconnection
#'   (TBR) distance, calculated using
#'   [TBRDist](https://ms609.github.io/TBRDist/)
#' - `rf`: Robinson-Foulds distance (Robinson & Foulds 1981)
#' - `icrf`: Robinson-Foulds distance, splits weighted by phylogenetic
#'   information content (Smith, forthcoming)
#' - `path`: Path distance (Steel & Penny 1993)
#'
#' Columns list the summary statistics of calculated tree distances: the
#' minimum (`min`),
#' 1%, 5%, 10%, 25%, 50% (i.e. median), 75%, 90%, 95%, 99% percentiles,
#' maximum (`max`), mean (`mean`) and standard deviation (`sd`).
#'
#' The third dimension lists the number of leaves in the trees compared.
#'
#' @templateVar vignette 09-expected-similarity
#' @template seeVignette
#' @template dataRaw
#'
#' @references
#' - \insertRef{Allen2001}{TreeDist}
#'
#' - \insertRef{Bocker2013}{TreeDist}
#'
#' - \insertRef{Bogdanowicz2012}{TreeDist}
#'
#'   \insertRef{Li1996}{TreeDist}
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
#'
#' @keywords datasets
"randomTreeDistances"

#' Distances between random pairs of trees
#'
#' Two-dimensional matrices listing the normalized distances between random
#' pairs of bifurcating trees with 11, 25 and 50 leaves drawn from the
#' uniform distribution using `ape::rtree(nTip, br=NULL)`.
#'
#' `pectinateDistances11` reports distances between a pectinate 11-leaf tree
#' and 100&nbsp;000 random bifurcating trees.
#'
#' Rows are named with an abbreviation of the tree comparison metric.
#' Information-based measures are normalized against the splitwise information
#' content/entropy for binary trees of the corresponding topologies.
#' <!--#TODO VERIFY-->
#' The quartet distance, Nye _et al._ and Jaccard-Robinson-Foulds measures
#' are normalized against their maximum possible values.
#' The remaining measures are unnormalized.
#'
#' - `pid`: Phylogenetic Information Distance (Smith, forthcoming)
#' - `msid`: Matching Split Information Distance (Smith, forthcoming)
#' - `cid`: Clustering Information Distance (Smith, forthcoming)
#' - `qd`: Quartet divergence (Smith 2019)
#' - `nye`: Nye _et al._ tree similarity (Nye _et al._ 2006)
#' - `jnc2`, `jnc4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#' conflicting pairings prohibited ('no-conflict')
#' - `jco2`, `jco4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#'  conflicting pairings permitted ('conflict-ok')
#' - `ms`: Matching Split Distance (Bogdanowicz & Giaro 2012), unnormalized
#' - `mast`: Size of Maximum Agreement Subtree (Valiente 2009)
#' - `masti`: Information content of Maximum Agreement Subtree
#' - `nni_l`, `nni_t`, `nni_u`: Lower bound, tight upper bound, and upper bound
#'         for  nearest-neighbour interchange distance (Li _et al._ 1996)
#' - `spr`: Approximate SPR distance, unnormalized
#' - `tbr_l`, `tbr_u`: Lower and upper bound for tree bisection and reconnection
#'          (TBR) distance
#' - `rf`: Robinson-Foulds distance (Robinson & Foulds 1981), unnormalized
#' - `icrf`: Robinson-Foulds distance, splits weighted by phylogenetic
#' information content (Smith, forthcoming)
#' - `path`: Path distance (Steel & Penny 1993), unnormalized
#'
#' Each column lists the calculated distances between each pair of trees.
#'
#' The pairs of 25-leaf trees are saved as data object [`randomTreePairs25`].
#'
#' @template dataRaw
#' @references
#' \insertRef{Bogdanowicz2012}{TreeDist}
#'
#' \insertRef{Li1996}{TreeDist}
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
#' \insertRef{Valiente2009}{TreeDist}
#'
#' @name distanceDistributions
#' @keywords datasets
NULL

#' @rdname distanceDistributions
"distanceDistribution50"
#' @rdname distanceDistributions
"distanceDistribution25"
#' @rdname distanceDistributions
"pectinateDistances11"

#' Pairs of random trees
#'
#' Lists of 10&nbsp;000 pairs of binary trees drawn from the uniform
#' distribution using `ape::rtree(nTip, br = NULL)`.
#'
#' The distances between these pairs of trees are recorded in
#' the data objects  [`distanceDistribution25`] and  [`distanceDistribution50`].
#'
#' @template dataRaw
#' @name randomTreePairs
#' @keywords datasets
NULL

#' @rdname randomTreePairs
"randomTreePairs25"
#' @rdname randomTreePairs
"randomTreePairs50"


#' Distances between unrooted seven-leaf trees
#'
#' Distances between each possible pairing of the 945 unrooted seven-leaf trees
#' (equivalent to rooted 6-leaf trees).  Following Kendall and Colijn (2016).
#'
## <!--Text copied from distanceDistribution25 above -->
#'
#' The list entries are named with an abbreviation of the tree comparison metric.
#' Information-based measures are normalized against the splitwise information
#' content/entropy for binary trees of the corresponding topologies.
#' <!--#TODO VERIFY-->
#' The quartet distance, Nye _et al._ and Jaccard-Robinson-Foulds measures
#' are normalized against their maximum possible values.
#' The remaining measures are unnormalized.
#'
#' - `pid`: Different Phylogenetic Information (Smith, forthcoming)
#' - `msid`: Matching Split Information Distance (Smith, forthcoming)
#' - `cid`: Clustering Information Distance (Smith, forthcoming)
#' - `qd`: Quartet divergence (Smith 2019)
#' - `nye`: Nye _et al._ tree similarity (Nye _et al._ 2006)
#' - `jnc2`, `jnc4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#' conflicting pairings prohibited ('no-conflict')
#' - `jco2`, `jco4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#'  conflicting pairings permitted ('conflict-ok')
#' - `ms`: Matching Split Distance (Bogdanowicz & Giaro 2012), unnormalized
#' - `mast`: Size of Maximum Agreement Subtree (Valiente 2009)
#' - `masti`: Information content of Maximum Agreement Subtree
#' - `nni_l`, `nni_t`, `nni_u`: Lower bound, tight upper bound, and upper bound
#'         for  nearest-neighbour interchange distance (Li _et al._ 1996)
#' - `spr`: Approximate SPR distance, unnormalized
#' - `tbr_l`, `tbr_u`: Lower and upper bound for tree bisection and reconnection
#'          (TBR) distance
#' - `rf`: Robinson-Foulds distance (Robinson & Foulds 1981), unnormalized
#' - `icrf`: Robinson-Foulds distance, splits weighted by phylogenetic
#' information content (Smith, forthcoming)
#' - `path`: Path distance (Steel & Penny 1993), unnormalized
#'
#' Each item in the list contains a 945&times;945 matrix reporting the distance
#' between each pair of seven-leaf trees.  The first 630 trees are pectinate
#' (tree shape 0), the final 315 are balanced (tree shape 1).
#'
#' @examples
#' library('TreeTools', quietly = TRUE, warn.conflicts = FALSE)
#'
#' # Pectinate unrooted tree shape:
#' plot(UnrootedTreeWithShape(0, 7))
#'
#' # Balanced unrooted tree shape:
#' plot(UnrootedTreeWithShape(1, 7))
#' @template dataRaw
#'
#' @references
#' \insertRef{Bogdanowicz2012}{TreeDist}
#'
#' \insertRef{Li1996}{TreeDist}
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
#' \insertRef{Valiente2009}{TreeDist}
#'
#' @keywords datasets
"sevenTipDistances"

#' Bullseye test results
#'
#' Implementation and results of a 'Bullseye' test, after that proposed by
#' Kuhner and Yamato (2015).
#'
#' `bullseyeTrees` is a list with four elements, named `5 leaves`, `10 leaves`,
#' `20 leaves` and `50 leaves`.  Each element contains 1&nbsp;000 trees with _n_
#' leaves, randomly sampled from the uniform distribution using `ape::rtree`.
#'
#' `bullXXInferred` is a list with four elements, named as in `bullseyeTrees`.
#' Each element contains 1&nbsp;000 sub-elements. Each sub-element is a list of
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
# <!--#TODO delete irrelevant statement -->
# `bullseyeDistances` contains four elements, each tabulating the distance
# between each pair of _n_-leaf trees in `bullseyeTrees`.  For details, see
# \code{\link{sevenTipDistances}}.
#
#' `bullseyeDistances` contains two elements, tabulating the distance
#' between each pair of 20- and 50-leaf trees in `bullseyeTrees`.
#' For details, see [`sevenTipDistances`].
#'
#'
#' @templateVar vignette 07-bullseye
#' @template seeVignette
#' @template dataRaw
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

#' Evaluating tree distance metrics by cluster recovery
#'
#' An effective measure of tree distance will recover clusters of similar
#' trees.  These datasets contain the results of tests modelled on those
#' in Lin _et al._ (2012).
#'
#'
#' I used three approaches to generate clusters of similar trees, and tested
#' each metric in its ability to recover these clusters (Lin _et al._, 2012).
#'
#'
#' For the first test, I generated 500 datasets of 100 binary trees with
#' _n_ = 40 leaves.
#' Each set of trees was created by randomly selecting two _k_-leaf
#' 'skeleton' trees, where _k_ ranges from 0.3 _n_ to 0.9 _n_.
#' From each skeleton, 50 trees were generated by adding each of the remaining
#' _n_ - _k_ leaves in turn at a uniformly selected point on the tree.
#'
#' For the second and third test, each dataset was constructed by selecting at
#' random two binary 40-leaf trees.
#' From each starting tree, I generated 50 binary trees by conducting _k_
#' leaf-label interchange (LLI) operations (test two) or _k_ subtree prune and
#' regraft (SPR) operations (test three) on the starting tree.
#' An LLI operation swaps the positions of two randomly selected leaves,
#' without affecting tree shape; an SPR operation moves a subtree to a new
#' location within the tree.
#'
#' For each dataset, I calculated the distance between each pair of trees.
#' Trees where then partitioned into clusters using five methods,
#' using the \pkg{stats} and \pkg{cluster}.
#' I define the success rate of each distance measure as the proportion of
#' datasets in which every tree generated from the same skeleton was placed
#' in the same cluster.
#'
#' @format A three-dimensional array.
#'
#' Rows correspond to the clustering methods:
#'
#' - `spc`: spectral clustering
#'
#' - `pam`: partitioning around medioids
#'
#' - `h...`: hierarchical clustering using:
#'  `h.cmp`, complete;
#'  `h.sng`, single; and
#'  `h.avg`, average linkage.
#'
#' Columns correspond to distance metrics; see 'Methods tested' below.
#'
#' Slices correspond to values of _k_:
#'
#' - `linTestOneResults`: _k_ = 30, 40, 50, 60, 70
#'
#' - `linTestTwoResults`: _k_ = 10, 20, 30, 40
#'
#' - `linTestSPRResults`: _k_ = 30, 40, 50, 60, 70
#'
#'
#' @templateVar vignette 06-lin-cluster-recovery
#' @template seeVignette
#'
#' @template methodAbbrevs
#'
#' @template dataRaw
#' @references \insertRef{Lin2012}{TreeDistData}
#' @name linTests
#'
#' @rdname linTests
'linTestOneResults'
#' @rdname linTests
'linTestTwoResults'
#' @rdname linTests
'linTestSPRResults'


#' Method parameters
#'
#' Metadata for methods examined in this package.
#'
#' `tdAbbrevs` lists abbreviations for each method, using expressions to allow
#' formatting of text when plotted.
#'
#' `tdPlotSequence` lists the 20 methods discussed in the main article,
#' in the sequence in which they are plotted in figures.
#'
#' `tdMdAbbrevs` uses markdown formatting.
#'
#' `tdBoxAbbrevs` uses line breaks to fit abbreviations in an approximately
#' square bounding box.
#'
#' `tdCol` provides each method with a suitable plotting colour.
#'
#' `TDFunctions` lists for each method a function that will calculate the
#' distance between two trees or lists of trees.
#'
#' `TDPair` lists for each method a function to calculate the distance
#' between one tree (`tr`) and another tree (`ref`).
#'
#' @template dataRaw
#' @name TreeDistMethods

#' @rdname TreeDistMethods
'tdAbbrevs'

#' @rdname TreeDistMethods
'tdPlotSequence'

#' @rdname TreeDistMethods
'tdMdAbbrevs'

#' @rdname TreeDistMethods
'tdBoxAbbrevs'

#' @rdname TreeDistMethods
'tdMethods'

#' @rdname TreeDistMethods
'tdCol'

#' @rdname TreeDistMethods
'TDFunctions'

#' @rdname TreeDistMethods
'TDPair'

#' Tree distance and SPR moves
#'
#' Datasets testing whether separating trees by increasingly many moves
#' results in a corresponding increase to their distance.
#'
#' I generated a chain of 100 50-leaf trees, starting from a pectinate tree
#' and deriving each tree in turn by performing an SPR operation on the previous
#' tree.
#' A consistent measure of tree similarity should correlate with the number of
#' SPR operations separating a pair of trees in this chain.
#' This said, because one SPR operation may counteract some of the difference
#' introduced by a previous one, perfect correlation is unlikely.
#'
#' @format A list of length 21.
#' Each entry is named according to the corresponding tree distance method; see
#' 'Methods tested' below.
#'
#' Each member of the list is a 100 x 100 matrix listing the distance
#' between each pair of trees in the SPR chain (see 'Details'),
#' numbered from 1 to 100.
#'
#'
#' @templateVar vignette 08-spr-walking
#' @template seeVignette
#' @template methodAbbrevs
#' @template dataRaw
'sprDistances'

#' Shape effect
#'
#' Results of tests exploring the influence of tree shape on reconstructed
#' tree distances.
#'
#' For each of the four binary unrooted tree shapes on eight leaves, I labelled
#' leaves at random until I had generated 100 distinct trees.
#'
#' I measured the distance from each tree to each of the other 399 trees.
#'
#' @templateVar vignette 05-tree-shape
#' @template seeVignette
#'
#' @format A list of length 21.
#' Each entry of the list is named according to the abbreviation of the
#' corresponding method (see 'Methods tested' below).
#'
#' Each entry is itself a list of ten elements.  Each element contains a numeric
#' vector listing the distances between each pair of trees with shape _x_ and
#' shape _y_, where:
#'
#' `x = 1, 1, 1, 1, 2, 2, 2, 3, 3, 4`
#' and
#' `y = 1, 2, 3, 4, 2, 3, 4, 3, 4, 4`.
#'
#' As trees are not compared with themselves (to avoid zero distances), elements
#' where _x_ = _y_ contain 4950 distances, whereas other elements contain 5050
#' distances.
#'
#'
#' @template methodAbbrevs
#' @template dataRaw
'shapeEffect'
