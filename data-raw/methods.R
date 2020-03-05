library('TreeDist')

tdPlotSequence <- c("cid", "dpi",
                    "nts", "jna2", "jna4", "ja2", "ja4",
                    "rf", "rfi",
                    "path",
                    "msid", "msd",
                    "qd",
                    "mast", "masti",
                    "nni_u", "nni_l", "spr", "tbr_l", "tbr_u")
legendSequence <- c('dpi', 'cid', 'nts', 'msd', 'qd', 'path', 'rf', 'rfi',
                    'ja4', 'ja2', 'jna4', 'jna2',
                    'nni_t', 'spr', 'tbr_u', 'mafi', 'mast', 'masti')

tdAbbrevs <- c(
  rf  = 'Robinson-Foulds',
  rfi = 'Info. Corr. RF',

  ja4 = 'JRF (k=4, arboreal)',
  ja2 = 'JRF (k=2, arboreal)',
  jna4 = 'JRF (k=4, non-arb.)',
  jna2 = 'JRF (k=2, non-arb.)',

  nea = 'Nye et al.',
  nts = expression(paste(plain('Nye '), italic('et al.'))),

  dpi = 'Phylog. Info. Dist',
  cid = 'Clust. Info. Dist.',
  msid = 'Match. Split Info Dist',
  msd = 'Match. Split Dist.',

  nni = 'NNI (approx.)',
  nni_u = 'NNI (upr bnd)',
  nni_t = 'NNI (ub tight)',
  nni_l = 'NNI (lwr bnd)',
  spr = 'SPR (approx.)',
  tbr = 'TBR (approx.)',
  tbr_l = 'TBR (lwr bnd)',
  tbr_u = 'TBR (upr bnd)',

  mast = 'MAST size',
  masti = 'MAST info',
  mafi = 'MAF info',

  qd  = 'Quartet',
  path = 'Path'
)

tdMdAbbrevs <- tdAbbrevs
tdMdAbbrevs['nts'] <- 'Nye _et al._'

tdBoxAbbrevs <- c(
  rf  = 'Robinson\n-Foulds',
  rfi = 'Info.\nCorr.\nRF',

  ja4 = 'JRF\n(k=4,\narboreal)',
  ja2 = 'JRF\n(k=2,\narboreal)',
  jna4 = 'JRF\n(k=4,\nnon-arb.)',
  jna2 = 'JRF\n(k=2,\nnon-arb.)',

  nea = 'Nye\net al.',
  nts = 'Nye\net al.',#expression(paste(plain('Nye\n'), italic('et al.'))),

  dpi = 'Phylog.\nInfo.\nDist.',
  cid = 'Clust.\nInfo.\nDist.',
  msid = 'MS\nInfo\nDist',
  msd = 'Match.\nSplit\nDist.',

  nni = 'NNI\n(approx.)',
  nni_u = 'NNI\n(upr bnd)',
  nni_t = 'NNI\n(ub tight)',
  nni_l = 'NNI\n(lwr bnd)',
  spr = 'SPR\n(approx.)',
  tbr = 'TBR\n(approx.)',
  tbr_l = 'TBR\n(lwr bnd)',
  tbr_u = 'TBR\n(upr bnd)',

  mast = 'MAST\nsize',
  masti = 'MAST\ninfo',
  mafi = 'MAF\ninfo',

  qd  = 'Quartet',
  path = 'Path'
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
  ja2 = function(...) TreeDist::JaccardRobinsonFoulds(..., k=2, arboreal =TRUE),
  ja4 =  function(...) TreeDist::JaccardRobinsonFoulds(..., k=4, arboreal=TRUE),
  jna2 = function(...) TreeDist::JaccardRobinsonFoulds(..., k=2, arboreal=FALSE),
  jna4 = function(...) TreeDist::JaccardRobinsonFoulds(..., k=4, arboreal=FALSE),

  nts = function(...) TreeDist::NyeTreeSimilarity(..., similarity = FALSE),
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
  msid = TreeDist::MatchingSplitInfoDistance,
  msd = TreeDist::MatchingSplitDistance,
  qd  = function (...) Quartet::QuartetDivergence(Quartet::ManyToManyQuartetAgreement(...),
                                                  similarity = FALSE),
  mafi = TBRDist::MAFInfo
)

TDPair <- list(
  rf = function (tr, ref) TreeDist::RobinsonFoulds(tr, ref),
  rfi = function (tr, ref) TreeDist::RobinsonFouldsInfo(tr, ref),
  ja2 = function (tr, ref) round(TreeDist::JaccardRobinsonFoulds(
    tr, ref, k = 2, arboreal = TRUE, normalize = TRUE), 4),
  jna2 = function (tr, ref) round(TreeDist::JaccardRobinsonFoulds(
    tr, ref, k = 2, arboreal = FALSE, normalize = TRUE), 4),
  ja4 = function (tr, ref) round(TreeDist::JaccardRobinsonFoulds(
    tr, ref, k = 4, arboreal = TRUE, normalize = TRUE), 4),
  jna4 = function (tr, ref) round(TreeDist::JaccardRobinsonFoulds(
    tr, ref, k = 4, arboreal = FALSE, normalize = TRUE), 4),

  dpi = function (tr, ref) round(TreeDist::DifferentPhylogeneticInfo(
    tr, ref, normalize = TRUE), 4L),
  msid = function (tr, ref) round(TreeDist::MatchingSplitInfoDistance(
    tr, ref, normalize = TRUE), 4L),
  cid = function (tr, ref) round(TreeDist::ClusteringInfoDistance(
    tr, ref, normalize = TRUE), 4L),
  nts = function (tr, ref) round(NyeTreeSimilarity(
    tr, ref, similarity = FALSE, normalize = TRUE), 4L),
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
               '#9D7660', '#D7B5A6'
)

# https://jrnold.github.io/ggthemes/reference/tableau_color_pal.html
tab30 <- as.character(matrix(c(tableau20, tableau10), 3, byrow=TRUE))

# https://personal.sron.nl/~pault/data/colourschemes.pdf
# plot(inlmisc::GetColors(n = 22, scheme = 'discrete rainbow'))
dr22 <- c("#D9CCE3", "#CAACCB", "#BA8DB4", "#AA6F9E", "#994F88", "#882E72",
          "#1965B0", "#437DBF", "#6195CF", "#7BAFDE", "#4EB265", "#90C987",
          "#CAE0AB", "#F7F056", "#F7CB45", "#F4A736", "#EE8026", "#E65518",
          "#DC050C", "#A5170E", "#72190E", "#42150A")

#tdCol <- tab30[c((1:10 * 2 - 1L), (seq_len(length(tdMethods) - 10L) * 2))]
colOrder <- c(22, 21, 3, 4, 1, 2, nts = 10, 7, 11, 6:5,
              15:17, 14, 18:19, 8:9, 12, 20, 13)
if(any(duplicated(colOrder))) warning(ifelse(duplicated(colOrder), colOrder, 0))
if (any(which(!1:22 %in% colOrder))) warning(which(!1:22 %in% colOrder))
tdCol <- dr22[colOrder]
names(tdCol) <- tdMethods
tdCol[c('nni', 'nea', 'tbr')] <- tdCol[c('nni_u', 'nts', 'tbr_u')]


usethis::use_data(tdAbbrevs, compress='xz', overwrite = TRUE)
usethis::use_data(tdMdAbbrevs, compress='xz', overwrite = TRUE)
usethis::use_data(tdBoxAbbrevs, compress='xz', overwrite = TRUE)
usethis::use_data(tdMethods, compress='xz', overwrite = TRUE)
usethis::use_data(tdPlotSequence, compress='xz', overwrite = TRUE)
usethis::use_data(tdCol, compress='xz', overwrite = TRUE)
usethis::use_data(TDFunctions, compress='xz', overwrite = TRUE)
usethis::use_data(TDPair, compress='xz', overwrite = TRUE)
