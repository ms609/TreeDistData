library('TreeDist')

tdPlotSequence <- c("cid", "pid",
                    "nye", "jco2", "jco4", "jnc2", "jnc4",
                    "rf", "icrf",
                    "path",
                    "msid", "ms",
                    "qd",
                    "mast", "masti",
                    "nni_u", "nni_l", "spr", "tbr_l", "tbr_u")

tdAbbrevs <- c(
  pid = 'Phylog. Info. Dist',
  msid = 'Match. Split Info Dist',
  cid = 'Clust. Info. Dist.',
  qd  = 'Quartet',

  nea = 'Nye et al.',
  nye = expression(paste(plain('Nye '), italic('et al.'))),

  jnc2 = 'JRF (k = 2, no conflict)',
  jnc4 = 'JRF (k = 4, no conflict)',
  jco2 = 'JRF (k = 2, conflict ok)',
  jco4 = 'JRF (k = 4, conflict ok)',

  ms = 'Match. Split Dist.',
  mast = 'MAST size',
  masti = 'MAST info',

  nni = 'NNI (approx.)',
  nni_l = 'NNI (lwr bnd)',
  nni_L = 'NNI (tgt lwr bnd)',
  nni_t = 'NNI (ub tight)',
  nni_u = 'NNI (upr bnd)',
  nni_U = 'NNI (tgt upr bnd)',

  spr = 'SPR (approx.)',
  tbr = 'TBR (approx.)',
  tbr_l = 'TBR (lwr bnd)',
  tbr_u = 'TBR (upr bnd)',

  rf  = 'Robinson-Foulds',
  icrf = 'Info. Corr. RF',
  path = 'Path',

  mafi = 'MAF info'
)

tdMdAbbrevs <- tdAbbrevs
tdMdAbbrevs['nye'] <- 'Nye _et al._'

tdBoxAbbrevs <- c(
  pid = 'Phylog.\nInfo.\nDist.',
  msid = 'MS\nInfo\nDist',
  cid = 'Clust.\nInfo.\nDist.',
  qd  = 'Quartet',

  nea = 'Nye\net al.',
  nye = 'Nye\net al.',#expression(paste(plain('Nye\n'), italic('et al.'))),

  jnc2 = 'JRF\n(k = 2,\nno conf.)',
  jnc4 = 'JRF\n(k = 4,\nno conf.)',
  jco2 = 'JRF\n(k = 2,\nconf. ok)',
  jco4 = 'JRF\n(k = 4,\nconf. ok)',

  ms = 'Match.\nSplit\nDist.',
  mast = 'MAST\nsize',
  masti = 'MAST\ninfo',

  nni = 'NNI\n(approx.)',
  nni_l = 'NNI\n(lwr bnd)',
  nni_L = 'NNI\n(tgt lwr bnd)',
  nni_t = 'NNI\n(ub tight)',
  nni_u = 'NNI\n(upr bnd)',
  nni_U = 'NNI\n(tgt upr bnd)',
  spr = 'SPR\n(approx.)',
  tbr = 'TBR\n(approx.)',
  tbr_l = 'TBR\n(lwr bnd)',
  tbr_u = 'TBR\n(upr bnd)',

  rf  = 'Robins.\n-Foulds',
  icrf = 'Info.\nCorr.\nRF',
  path = 'Path',

  mafi = 'MAF\ninfo'
)

tdMethods <- names(tdAbbrevs)
tdMethods <- tdMethods[!tdMethods %in% c('nni', 'nea', 'tbr')]

JA2 <- function (...) TreeDist::JaccardRobinsonFoulds(..., k = 2,
                                                      allowConflict = FALSE)
JA4 <- function (...) TreeDist::JaccardRobinsonFoulds(..., k = 4,
                                                      allowConflict = FALSE)
JNA2 <- function (...) TreeDist::JaccardRobinsonFoulds(..., k = 2,
                                                       allowConflict = TRUE)
JNA4 <- function (...) TreeDist::JaccardRobinsonFoulds(..., k = 4,
                                                       allowConflict = TRUE)

TDFunctions <- list(
  pid = TreeDist::DifferentPhylogeneticInfo,
  msid = TreeDist::MatchingSplitInfoDistance,
  cid = TreeDist::ClusteringInfoDistance,
  qd  = function (...) Quartet::QuartetDivergence(
    Quartet::ManyToManyQuartetAgreement(...), similarity = FALSE),
  nye = function(...) TreeDist::NyeSimilarity(..., similarity = FALSE),

  jnc2 = function(...) TreeDist::JaccardRobinsonFoulds(..., k = 2,
                                                       allowConflict = FALSE),
  jnc4 =  function(...) TreeDist::JaccardRobinsonFoulds(..., k = 4,
                                                        allowConflict = FALSE),
  jco2 = function(...) TreeDist::JaccardRobinsonFoulds(..., k = 2,
                                                       allowConflict = TRUE),
  jco4 = function(...) TreeDist::JaccardRobinsonFoulds(..., k = 4,
                                                       allowConflict = TRUE),

  ms = TreeDist::MatchingSplitDistance,
  mast =  function(...) TreeDist::MASTSize(..., rooted = FALSE),
  masti = function(...) TreeDist::MASTInfo(..., rooted = FALSE),

  nni_u =  function(...) as.matrix(TreeDist::NNIDist(...)$loose_upper),
  nni_t =  function(...) as.matrix(TreeDist::NNIDist(...)$tight_upper),
  nni_T =  function(...) as.matrix(TreeDist::NNIDist(...)$best_upper),
  nni_l =  function(...) as.matrix(TreeDist::NNIDist(...)$lower),
  nni_L =  function(...) as.matrix(TreeDist::NNIDist(...)$best_lower),
  spr = phangorn::SPR.dist,
  tbr =  function(...) as.matrix(TBRDist::TBRDist(...)$tbr_max),
  tbr_l =  function(...) as.matrix(TBRDist::TBRDist(...)$tbr_min),
  tbr_u =  function(...) as.matrix(TBRDist::TBRDist(...)$tbr_max),
  rf = TreeDist::RobinsonFoulds,
  icrf = TreeDist::InfoRobinsonFoulds,
  path = phangorn::path.dist,
  mafi = TBRDist::MAFInfo
)

TDPair <- list(
  pid = function (tr, ref) round(TreeDist::DifferentPhylogeneticInfo(
    tr, ref, normalize = TRUE), 4L),
  msid = function (tr, ref) round(TreeDist::MatchingSplitInfoDistance(
    tr, ref, normalize = TRUE), 4L),
  cid = function (tr, ref) round(TreeDist::ClusteringInfoDistance(
    tr, ref, normalize = TRUE), 4L),
  nye = function (tr, ref) round(NyeSimilarity(
    tr, ref, similarity = FALSE, normalize = TRUE), 4L),
  jnc2 = function (tr, ref) round(TreeDist::JaccardRobinsonFoulds(
    tr, ref, k = 2, allowConflict = FALSE, normalize = TRUE), 4),
  jnc4 = function (tr, ref) round(TreeDist::JaccardRobinsonFoulds(
    tr, ref, k = 4, allowConflict = FALSE, normalize = TRUE), 4),
  jco2 = function (tr, ref) round(TreeDist::JaccardRobinsonFoulds(
    tr, ref, k = 2, allowConflict = TRUE, normalize = TRUE), 4),
  jco4 = function (tr, ref) round(TreeDist::JaccardRobinsonFoulds(
    tr, ref, k = 4, allowConflict = TRUE, normalize = TRUE), 4),

  ms = function (tr, ref) signif(TreeDist::MatchingSplitDistance(tr, ref), 4),
  qd = function (tr, ref) Quartet::QuartetStatus(list(tr, ref))[2, 'd'],
  nni_t = function(tr, ref) TreeDist::NNIDist(tr, ref)['tight_upper'],
  nni_l = function(tr, ref) TreeDist::NNIDist(tr, ref)['lower'],
  nni_L = function(tr, ref) TreeDist::NNIDist(tr, ref)['best_lower'],
  nni_u = function(tr, ref) TreeDist::NNIDist(tr, ref)['loose_upper'],
  nni_U = function(tr, ref) TreeDist::NNIDist(tr, ref)['best_upper'],
  spr = function (tr, ref) phangorn::SPR.dist(tr, ref),
  tbr_u = function(tr, ref) TBRDist::TBRDist(tr, ref)$tbr_max,
  tbr_l = function(tr, ref) TBRDist::TBRDist(tr, ref)$tbr_min,

  mast = function (...) TreeDist::MASTSize(...),
  masti = function (...) TreeDist::MASTInfo(...),
  rf = function (tr, ref) TreeDist::RobinsonFoulds(tr, ref),
  icrf = function (tr, ref) TreeDist::InfoRobinsonFoulds(tr, ref),
  path = function (tr, ref) signif(phangorn::path.dist(tr, ref), 4L),
  mafi = function (...) TBRDist::MAFInfo(...)
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
tab30 <- as.character(matrix(c(tableau20, tableau10), 3, byrow = TRUE))

# https://personal.sron.nl/~pault/data/colourschemes.pdf
# plot(inlmisc::GetColors(n = 22, scheme = 'discrete rainbow'))
dr22 <- c("#D9CCE3", "#CAACCB", "#BA8DB4", "#AA6F9E", "#994F88", "#882E72",
          "#1965B0", "#437DBF", "#6195CF", "#7BAFDE", "#4EB265", "#90C987",
          "#CAE0AB", "#F7F056", "#F7CB45", "#F4A736", "#EE8026", "#E65518",
          "#DC050C", "#A5170E", "#72190E", "#42150A")

#tdCol <- tab30[c((1:10 * 2 - 1L), (seq_len(length(tdMethods) - 10L) * 2))]
colOrder <- c(pid = 7, msid = 6, cid = 11, qd = 20, nye = 10,
              jnc2 = 3, jnc4 = 4, jco2 = 1, jco4 = 2,
              ms = 5, mast = 8, masti = 9,
              nni_l = 17, nni_t = 16, nni_u = 15,
              spr = 14, tbr_l = 18, tbr_u = 19,
              rf = 22, icrf = 21, path = 13, mafi = 12)
if(any(duplicated(colOrder))) warning(ifelse(duplicated(colOrder), colOrder, 0))
if (any(which(!1:22 %in% colOrder))) warning(which(!1:22 %in% colOrder))
tdCol <- dr22[colOrder[tdMethods]]
names(tdCol) <- tdMethods
tdCol[c('nni', 'nea', 'tbr', 'nni_L', 'nni_U')] <-
  tdCol[c('nni_u', 'nye', 'tbr_u', 'nni_l', 'nni_u')]


usethis::use_data(tdAbbrevs, compress = 'gzip', overwrite = TRUE)
usethis::use_data(tdMdAbbrevs, compress = 'gzip', overwrite = TRUE)
usethis::use_data(tdBoxAbbrevs, compress = 'gzip', overwrite = TRUE)
usethis::use_data(tdMethods, compress = 'gzip', overwrite = TRUE)
usethis::use_data(tdPlotSequence, compress = 'gzip', overwrite = TRUE)
usethis::use_data(tdCol, compress = 'gzip', overwrite = TRUE)
usethis::use_data(TDFunctions, compress = 'xz', overwrite = TRUE)
usethis::use_data(TDPair, compress = 'xz', overwrite = TRUE)
