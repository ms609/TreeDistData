library('phangorn')
library('TreeSearch')
library('TreeDist')
devtools::load_all() # necessary for correct path-to-inst/

data("bullseyeTrees", package = 'TreeDistData') # Generated in bullseyeTrees.R
tipsNames <- names(bullseyeTrees)
nTrees <- 1000L
subsamples <- 10:1 * 200

# Define functions:
WriteTNTData <- function (dataset, fileName) {
  index <- attr(dataset, 'index')
  write(paste0('nstates num ', attr(dataset, 'nc'), ';\n',
               'xread', '\n',
               length(index), ' ', length(dataset), '\n',
        paste(paste(names(dataset), lapply(dataset, function (x)
          paste0(x[index], collapse = '')), collapse = '\n')),
        '\n;'),
        fileName)
}


CacheFile <- function (name, ..., tmpDir = FALSE) {
  root <- if (tmpDir) {
    tempdir()
  } else {
    paste0(system.file(package = 'TreeDistData'), '/../data-raw/trees/')
  }
  paste0(root, name, ...)
}


data('bullseyeMorphInferred', package = "TreeDistData")
if (!exists('bullseyeMorphInferred')) {
  message("\n\n = == Infer trees = = =\n")
  bullseyeMorphInferred <- lapply(bullseyeTrees, function (x) NULL)

  for (tipName in names(bullseyeTrees)) {
    message('* ', tipName, ": Simulating sequences...")
    theseTrees <- bullseyeTrees[[tipName]][seq_len(nTrees)]
    seqs <- lapply(theseTrees, simSeq, l = 2000, type = 'USER', levels = 1:4)
    inferred <- vector(mode = 'list', nTrees)


    for (i in seq_along(seqs)) {
      if (i %% 100 == 1) message(i)
      seq00 <- formatC(i - 1L, width = 3, flag = '0')
      FilePattern <- function (n) {
        paste0(substr(tipName, 0, nchar(tipName) - 7),
               't-', seq00, '-k6-',
               formatC(n, width = 4, flag = '0'),
               '.tre')
      }

      if (!file.exists(CacheFile(FilePattern(200)))) {
        message("File not found at", CacheFile(FilePattern(200)), "; inferring:")
        seqFile <- CacheFile('seq-', seq00, '.tnt')
        WriteTNTData(seqs[[i]], file = seqFile)
        Line <- function (n) {
          paste0("piwe = 6 ;xmult;tsav *", FilePattern(n),
                 ";sav;tsav/;keep 0;hold 10000;\n",
          "ccode ] ", paste(seq_len(200) + n - 200L, collapse = ' '), ";\n"
          )
        }

        runRoot <- paste0(sample(letters, 8, replace = TRUE), collapse = '')
        runFile <- CacheFile(runRoot, '.run')
        file.create(runFile)
        write(paste("macro =;
         xmult:hits 1 level 4 chklevel 5 rat5 drift5;
         sect:slack 8;
         keep 0; hold 10000;\n",
         Line(2000),
         Line(1800),
         Line(1600),
         Line(1400),
         Line(1200),
         Line(1000),
         Line(0800),
         Line(0600),
         Line(0400),
         Line(0200),
         "quit;"), runFile)
        # Install TNT and add to the PATH environment variable before running:
        system(paste('tnt proc', seqFile, '; ', runRoot, ';'))
        file.remove(seqFile)
        file.remove(runFile)
      }

      inferred[[i]] <-
        lapply(formatC(subsamples, width = 4, flag = '0'),
               function (nChar) {
                 tr <- ReadTntTree(CacheFile(FilePattern(nChar)),
                                   relativePath = '.',
                    tipLabels = theseTrees[[i]]$tip.label)
                 # Return:
                 if (class(tr) == 'multiPhylo') tr[[1]] else tr
                 })
    }
    bullseyeMorphInferred[[tipName]] <- inferred
  }
  usethis::use_data(bullseyeMorphInferred, compress = 'bzip2', overwrite = TRUE)
}

message("\n\n = = = Calculate distances = = =\n")
bullseyeMorphScores <- setNames(vector('list', length(tipsNames)), tipsNames)

for (tipName in tipsNames) {
  cat('\u2714 Calculating tree distances:', tipName, ':\n')
  inferred <- bullseyeMorphInferred[[tipName]]
  trueTrees <- bullseyeTrees[[tipName]]
  theseScores <- vapply(seq_along(inferred), function (i) {
    cat('.')
    if (i %% 72 == 0) cat(' ', i, "\n")
    trueTree <- trueTrees[[i]]
    rootTip <- trueTree$tip.label[1]
    tr <- root(trueTree, rootTip, resolve.root = TRUE)
    tr$edge.length  <- NULL
    trs <- structure(lapply(inferred[[i]], root, rootTip, resolve.root = TRUE),
                     class = 'multiPhylo')

    mast <- vapply(trs, MASTSize, tr, rooted = FALSE, FUN.VALUE = 1L)
    masti <-  LnUnrooted(mast) / log(2)
    attributes(masti) <- attributes(mast)

    nni <- NNIDist(tr, trs)
    tbr <- TBRDist(tr, trs)

    cbind(
      # The order here MUST correspond to the dimnames template below!
      pid = DifferentPhylogeneticInfo(tr, trs, normalize = TRUE),
      msid = MatchingSplitInfoDistance(tr, trs, normalize = TRUE),
      cid = ClusteringInfoDistance(tr, trs, normalize = TRUE),
      nye = NyeSimilarity(tr, trs, similarity = FALSE, normalize = TRUE),
      qd = Quartet::QuartetDivergence(Quartet::QuartetStatus(trs, cf = tr),
                                      similarity = FALSE),

      jnc2 = JaccardRobinsonFoulds(tr, trs, k = 2, allowConflict = FALSE,
                                  normalize = TRUE),
      jnc4 = JaccardRobinsonFoulds(tr, trs, k = 4, allowConflict = FALSE,
                                  normalize = TRUE),
      jco2 = JaccardRobinsonFoulds(tr, trs, k = 2, allowConflict = TRUE,
                                   normalize = TRUE),
      jco4 = JaccardRobinsonFoulds(tr, trs, k = 4, allowConflict = TRUE,
                                   normalize = TRUE),

      ms = MatchingSplitDistance(tr, trs),
      mast = mast,
      masti = masti,

      nni_l = nni['lower', ],
      nni_L = nni['best_lower', ],
      nni_t = nni['tight_upper', ],
      nni_U = nni['best_upper', ],
      nni_u = nni['loose_upper', ],
      spr = SPR.dist(tr, trs),
      tbr_l = tbr$tbr_min,
      tbr_u = tbr$tbr_max,

      rf = RobinsonFoulds(tr, trs),
      icrf = InfoRobinsonFoulds(tr, trs),
      path = path.dist(tr, trs),

      kc = KendallColijn(tr, trs),
      es = KendallColijn(tr, trs, SplitVector)
    )
  }, matrix(0, nrow = 10L, ncol = 25L,
            dimnames = list(subsamples, c('pid', 'msid', 'cid', 'nye', 'qd',
                                          'jnc2', 'jnc4', 'jco2', 'jco4',
                                          'ms', 'mast', 'masti', 'nni_l',
                                          'nni_L', 'nni_t', 'nni_U', 'nni_u',
                                          'spr', 'tbr_l', 'tbr_u', 'rf',
                                          'icrf', 'path', 'kc', 'es')))
  )
  bullseyeMorphScores[[tipName]] <- theseScores
}
usethis::use_data(bullseyeMorphScores, compress = 'xz', overwrite = TRUE)
cat(" # # # COMPLETE # # # ")
