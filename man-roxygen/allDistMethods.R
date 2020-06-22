#' @details
#' # Methods tested
#'
#' - `pid`: Phylogenetic Information Distance (Smith, forthcoming), normalized
#' against the phylogenetic information content of the splits in the trees
#' being compared.
#' - `msid`: Matching Split Information Distance (Smith, forthcoming), normalized
#' against the phylogenetic information content of the splits in the trees
#' being compared.
#' - `cid`: Clustering Information Distance (Smith, forthcoming), normalized
#' against the entropy of the splits in the trees being compared.
#' - `qd`: Quartet divergence (Smith 2019), normalized against its maximum
#' possible value for _n_-leaf trees.
#' - `nye`: Nye _et al._ tree distance (Nye _et al._ 2006), normalized against
#' the total number of splits in the trees being compared.
#' - `jnc2`, `jnc4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#' conflicting pairings prohibited ('no-conflict'), normalized against
#' the total number of splits in the trees being compared.
#' - `jco2`, `jco4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#' conflicting pairings permitted ('conflict-ok'), normalized against
#' the total number of splits in the trees being compared.
#' - `ms`: Matching Split Distance (Bogdanowicz & Giaro 2012), unnormalized.
#' - `mast`: Size of Maximum Agreement Subtree (Valiente 2009), unnormalized.
#' - `masti`: Information content of Maximum Agreement Subtree, unnormalized.
#' - `nni_l`, `nni_t`, `nni_u`: Lower bound, tight upper bound, and upper bound
#' for nearest-neighbour interchange distance (Li _et al._ 1996), unnormalized.
#' - `spr`: Approximate subtree prune and regraft \acronym{SPR} distance,
#' unnormalized.
#' - `tbr_l`, `tbr_u`: Lower and upper bound for tree bisection and reconnection
#' (\acronym{TBR}) distance, calculated using
#'   [\pkg{TBRDist}](https://ms609.github.io/TBRDist/); unnormalized.
#' - `rf`: Robinson-Foulds distance (Robinson & Foulds 1981), unnormalized.
#' - `icrf`: Robinson-Foulds distance, splits weighted by phylogenetic
#' information content (Smith, forthcoming), unnormalized.
#' - `path`: Path distance (Steel & Penny 1993), unnormalized.
