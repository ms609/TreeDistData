library('TreeDist')

plotSequence <- c('dpi', 'cid', 'nts', 'qd',
                  'msd', 'path', 'spr', 'rf',
                  'rfi',
                  'ja2', 'ja4', 'jna2', 'jna4',
                  'nni_l', 'nni_t', 'nni_u', 'tbr_l',
                  'tbr_u',
                  'mast', 'masti')
legendSequence <- c('dpi', 'cid', 'nts', 'msd', 'qd', 'path', 'rf', 'rfi',
                    'ja2', 'ja4', 'jna2', 'jna4',
                    'nni_t', 'spr', 'tbr_u', 'mafi', 'mast', 'masti')

tdAbbrevs <- c(
  rf  = 'Robinson-Foulds',
  rfi = 'Robinson-Foulds Info',

  ja2 = 'JRF (k=2, arboreal)',
  ja4 = 'JRF (k=4, arboreal)',
  jna2 = 'JRF (k=2, non-arb.)',
  jna4 = 'JRF (k=4, non-arb.)',

  dpi = 'Diff. Phylog. Info',
  nni = 'NNI (approx.)',
  nni_u = 'NNI (upr bnd)',
  nni_t = 'NNI (ub tight)',
  nni_l = 'NNI (lwr bnd)',
  spr = 'SPR (approx.)',

  tbr = 'TBR (approx.)',
  tbr_l = 'TBR (lb approx.)',
  tbr_u = 'TBR (ub approx.)',
  mafi = 'MAF info',
  path = 'Path',
  mast = 'MAST size',
  masti = 'MAST info',
  mafi = 'MAF info',
  cid = 'Clust. Info. Dist.',

  nea = 'Nye et al.',
  nts = expression(paste(plain('Nye '), italic('et al.'))),
  qd  = 'Quartet',
  msid = 'Match. Split Info Dist',
  msd = 'Match. Split Dist.'
)

tdMethods <- names(tdAbbrevs)
tdMethods <- tdMethods[!tdMethods %in% c('nni', 'nea', 'tbr')]

JA2 <- function (...) TreeDist::JaccardRobinsonFoulds(..., k=2, arboreal=TRUE)
JA4 <- function (...) TreeDist::JaccardRobinsonFoulds(..., k=4, arboreal=TRUE)
JNA2 <- function (...) TreeDist::JaccardRobinsonFoulds(..., k=2, arboreal=FALSE)
JNA4 <- function (...) TreeDist::JaccardRobinsonFoulds(..., k=4, arboreal=FALSE)

TDFunctions <- list(
  rf = TreeDist::RobinsonFoulds,
  rfi = TreeDist::RobinsonFouldsInfo,
  ja2 = JA2,
  ja4 = JA4,
  jna2 = JNA2,
  jna4 = JNA4,

  dpi = DifferentPhylogeneticInfo,
  nni_u =  function(...) as.matrix(TreeDist::NNIDist(...)$loose_upper),
  nni_t =  function(...) as.matrix(TreeDist::NNIDist(...)$tight_upper),
  nni_l =  function(...) as.matrix(TreeDist::NNIDist(...)$lower),
  tbr_l =  function(...) as.matrix(TBRDist::TBRDist(...)$tbr_min),
  tbr_u =  function(...) as.matrix(TBRDist::TBRDist(...)$tbr_max),
  spr = phangorn::SPR.dist,
  tbr =  function(...) as.matrix(TBRDist::TBRDist(...)$tbr_max),
  path = phangorn::path.dist,
  mast =  function(...) TreeDist::MASTSize(..., rooted = FALSE),
  masti = function(...) TreeDist::MASTInfo(..., rooted = FALSE),
  cid = TreeDist::ClusteringInfoDistance,
  nts = TreeDist::NyeTreeSimilarity,
  msid = TreeDist::MatchingSplitInfoDistance,
  msd = TreeDist::MatchingSplitDistance,
  qd  = function (...) Quartet::QuartetDivergence(Quartet::ManyToManyQuartetAgreement(...)),
  mafi = TBRDist::MAFInfo
)

TDPair <- list(
  rf = function (tr, ref) TreeDist::RobinsonFoulds(tr, ref),
  rfi = function (tr, ref) TreeDist::RobinsonFouldsInfo(tr, ref),
  ja2 = function (tr, ref) JA2(tr, ref),
  jna2 = function (tr, ref) JNA2(tr, ref),
  ja4 = function (tr, ref) JA4(tr, ref),
  jna4 = function (tr, ref) JNA2(tr, ref),

  dpi = function (tr, ref) round(TreeDist::DifferentPhylogeneticInfo(tr, ref, normalize = TRUE), 4L),
  msid = function (tr, ref) round(TreeDist::MatchingSplitInfoDistance(tr, ref, normalize = TRUE), 4L),
  cid = function (tr, ref) round(TreeDist::ClusteringInfoDistance(tr, ref, normalize = TRUE), 4L),
  nts = function (tr, ref) round(1 - NyeTreeSimilarity(tr, ref, normalize = TRUE), 4L),
  tbr_u = function(tr, ref) TBRDist::TBRDist(tr, ref)$tbr_max,
  tbr_l = function(tr, ref) TBRDist::TBRDist(tr, ref)$tbr_min,
  nni_t = function(tr, ref) TreeDist::NNIDist(tr, ref)['tight_upper'],
  nni_l = function(tr, ref) TreeDist::NNIDist(tr, ref)['lower'],
  nni_u = function(tr, ref) TreeDist::NNIDist(tr, ref)['loose_upper'],
  mast = TreeDist::MASTSize, masti = TreeDist::MASTInfo, mafi = TBRDist::MAFInfo,
  msd = function (tr, ref) signif(TreeDist::MatchingSplitDistance(tr, ref), 4),
  qd = function (tr, ref) Quartet::QuartetStatus(list(tr, ref))[2, 'd'],
  path = function (tr, ref) signif(phangorn::path.dist(tr, ref), 4L),

  spr = function (tr, ref) phangorn::SPR.dist(tr, ref)
)




tableau10 <- c('#4e78a8',  '#f28e2c',  '#e15659',  '#75b7b2',  '#58a14e',
               '#edc949',  '#af7aa1',  '#fd9da7',  '#9d745f',  '#bab0ac')
# https://www.tableau.com/about/blog/2016/7/colors-upgrade-tableau-10-56782
tableau20 <- c('#4e78a8', '#A0CBE8',
               '#f28e2c', '#FFBE7D',
               '#59a14f', '#8CD17D',
               '#B6992D', '#F1CE63',
               '#499894', '#86BCB6',
               '#E15759', '#FF9D9A',
               '#79706E', '#BAB0AC',
               '#D37295', '#FABFD2',
               '#B07AA1', '#D4A6C8',
               '#9D7660', '#D7B5A6',
               '#000000','#000000','#000000','#000000','#000000','#000000','#000000',
               '#000000'
)
# https://jrnold.github.io/ggthemes/reference/tableau_color_pal.html

tdCol <- tableau20[c((1:10 * 2 - 1L), (seq_len(length(tdMethods) - 10L) * 2))]
names(tdCol) <- tdMethods
tdCol[c('nni', 'nea', 'tbr')] <- tdCol[c('nni_u', 'nts', 'tbr_u')]


usethis::use_data(tdAbbrevs, compress='xz', overwrite = TRUE)
usethis::use_data(tdMethods, compress='xz', overwrite = TRUE)
usethis::use_data(tdCol, compress='xz', overwrite = TRUE)
usethis::use_data(TDFunctions, compress='xz', overwrite = TRUE)
usethis::use_data(TDPair, compress='xz', overwrite = TRUE)
