#' @details
#' # Methods tested
#'
#' - `pid`: Phylogenetic Information Distance (Smith, forthcoming)
#' - `msid`: Matching Split Information Distance (Smith, forthcoming)
#' - `cid`: Clustering Information Distance (Smith, forthcoming)
#' - `qd`: Quartet divergence (Smith, 2019)
#' - `nye`: Nye _et al._ tree distance (Nye _et al._ 2006)
#' - `jnc2`, `jnc4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#' conflicting pairings prohibited ('no-conflict')
#' - `joc2`, `jco4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#'  conflicting pairings permitted ('conflict-ok')
#' - `ms`: Matching Split Distance (Bogdanowicz & Giaro 2012)
#' - `mast`: Size of Maximum Agreement Subtree (Valiente 2009)
#' - `masti`: Information content of Maximum Agreement Subtree
#' - `nni_l`,`r if (<%= nni_t %>) " ``nni_t``,"` `nni_u`: Lower
#'   `r ifelse(<%= nni_t %>, "bound, tight upper bound, and upper bound", "and upper bounds")`
#'   for nearest-neighbour interchange distance (Li _et al._ 1996)
#' - `spr`: Approximate SPR distance
#' - `tbr_l`, `tbr_u`: Lower and upper bound for tree bisection and reconnection
#' (TBR) distance, calculated using
#' [\pkg{TBRDist}](https://ms609.github.io/TBRDist/)
#' - `rf`: Robinson-Foulds distance (Robinson & Foulds 1981)
#' - `icrf`: Information-corrected Robinson-Foulds distance (Smith, forthcoming)
#' - `path`: Path distance (Steel & Penny 1993)
