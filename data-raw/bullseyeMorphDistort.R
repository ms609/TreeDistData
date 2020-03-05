library('phangorn')
library('TreeSearch')
library('TreeDist')
devtools::load_all() # necessary for correct path-to-inst/

data("bullseyeTrees", package='TreeDistData') # Generated in bullseyeTrees.R
tipsNames <- names(bullseyeTrees)
subsamples <- 0:9 * 2 # Order in increasing dissimilarity, please

# Get results for a subset of trees:
nTrees <- 1000

# Define functions:
WriteTNTData <- function (dataset, fileName) {
  index <- attr(dataset, 'index')
  write(paste0('nstates num ', attr(dataset, 'nc'), ';\n',
               'xread', '\n',
               length(index), ' ', length(dataset), '\n',
               paste(paste(names(dataset), lapply(dataset, function (x)
                 paste0(x[index], collapse='')), collapse='\n')),
               '\n;'),
        fileName)
}

CacheFile <- function (..., tmpDir = FALSE) {
  root <- if (tmpDir) tempdir() else paste0(system.file(package='TreeDistData'),
                                            '/../data-raw/bullMoDi/')
  paste0(root, ...)
}

# Create subdirectories:
if (!dir.exists(CacheFile())) dir.create(CacheFile('../bullMoDi'))

bullMoDiInferred <- vector('list', length(tipsNames))
names(bullMoDiInferred) <- tipsNames

message("\n\n=== Infer trees ===\n")
for (tipName in names(bullseyeTrees)) {
  message('* ', tipName, ": Simulating sequences...")
  theseTrees <- bullseyeTrees[[tipName]][seq_len(nTrees)]
  simChar <- 2000L
  seqs <- lapply(theseTrees, simSeq, l = simChar, type='USER', levels=0:1)
  nData <- length(seqs[[1]]) * simChar
  inferred <- vector(mode='list', nTrees)


  for (i in seq_along(seqs)) {
    seq00 <- formatC(i - 1, width=3, flag='0')
    cat ("", seq00)
    FilePattern <- function (n) {
      CacheFile(substr(tipName, 0, nchar(tipName) - 7L),
             't-', seq00, '-k6-',
             formatC(n, width=2, flag='0'),
             '.tre')
    }

    if (!all(file.exists(FilePattern(subsamples)))) {
      cat("Generating missing files: ", FilePattern(subsamples)[!file.exists(FilePattern(subsamples))])
      seqFile <- paste0(tempdir(), '\\bullMoDi-', seq00, '.tnt')
      runRoot <- paste0(sample(letters, 8, replace = TRUE), collapse = '')
      runFile <- paste0(runRoot, '.run', collapse='')
      file.create(runFile)
      write(paste("macro =;
       xmult:hits 1 level 4 chklevel 5 rat5 drift5;
       sect:slack 8;
       keep 0; hold 10000;
       piwe=6 ;xmult;tsav *", FilePattern('%1'), ";sav;tsav/;keep 0;hold 10000;
       quit;"), runFile)

      sq <- PhyDatToMatrix(seqs[[i]])
      for (sub in subsamples) {
        switch <- sample(seq_len(nData), nData * sub / 100L)
        sq[switch] <- 1L - as.integer(sq[switch])
        WriteTNTData(MatrixToPhyDat(sq), file = seqFile)
        # Install TNT and add to the PATH environment variable before running:
        system(paste('tnt proc', seqFile, '; ', runRoot,
                     formatC(sub, width=2, flag='0'), ';'))
      }

      file.remove(seqFile)
      file.remove(runFile)
    }

    inferred[[i]] <-
      lapply(formatC(subsamples, width = 2, flag='0'),
             function (nChar) {
               tr <- ReadTntTree(FilePattern(nChar),
                                 relativePath = '.',
                                 tipLabels = theseTrees[[i]]$tip.label)
               # Return:
               if (class(tr) == 'multiPhylo') tr[[1]] else tr
             })
  }
  bullMoDiInferred[[tipName]] <- inferred
}
usethis::use_data(bullMoDiInferred, compress = 'bzip2', overwrite = TRUE)


message("\n\n === Calculate distances ===\n")
bullMoDiScores <- vector('list', length(tipsNames))
names(bullMoDiScores) <- tipsNames
sampledMethods <-
  c('dpi', 'msid', 'cid', 'nts', 'ja2', 'ja4', 'jna2', 'jna4',
    'msd', 'mast', 'masti', 'nni_l', 'nni_t', 'nni_u', 'spr', 'tbr_l', 'tbr_u',
    'rf', 'rfi', 'qd', 'path')
for (tipName in tipsNames) {
  inferred <- bullMoDiInferred[[tipName]]
  trueTrees <- bullseyeTrees[[tipName]]
  cat ("\n *** Scoring:", tipName, '***\n')
  theseScores <- vapply(seq_along(inferred), function (i) {
    cat('.')
    if (i %% 72 == 0) cat(' ', i, "\n")
    trueTree <- trueTrees[[i]]
    rootTip <- trueTree$tip.label[1]
    tr <- root(trueTree, rootTip, resolve.root=TRUE)
    tr$edge.length  <- NULL
    trs <- structure(lapply(inferred[[i]], root, rootTip, resolve.root=TRUE),
                     class = 'multiPhylo')

    mast <- vapply(trs, MASTSize, tr, rooted = FALSE, FUN.VALUE = 1L)
    masti <-  LnUnrooted(mast) / log(2)
    attributes(masti) <- attributes(mast)

    nni <- NNIDist(tr, trs)
    tbr <- TBRDist(tr, trs)

    normInfo <- SplitwiseInfo(tr)
    cbind(
      dpi = DifferentPhylogeneticInfo(tr, trs, normalize = TRUE),
      msid = MatchingSplitInfoDistance(tr, trs, normalize = TRUE),
      cid = ClusteringInfoDistance(tr, trs, normalize = TRUE),
      nts = NyeTreeSimilarity(tr, trs, similarity = FALSE, normalize = TRUE),

      ja2 = JaccardRobinsonFoulds(tr, trs, k = 2, arboreal = TRUE, normalize = TRUE),
      ja4 = JaccardRobinsonFoulds(tr, trs, k = 4, arboreal = TRUE, normalize = TRUE),
      jna2 =JaccardRobinsonFoulds(tr, trs, k = 2, arboreal = FALSE, normalize = TRUE),
      jna4 =JaccardRobinsonFoulds(tr, trs, k = 4, arboreal = FALSE, normalize = TRUE),

      msd = MatchingSplitDistance(tr, trs),
      mast = mast,
      masti = masti,

      nni_l = nni['lower', ],
      nni_t = nni['tight_upper', ],
      nni_u = nni['loose_upper', ],
      spr = SPR.dist(tr, trs),
      tbr_l = tbr$tbr_min,
      tbr_u = tbr$tbr_max,

      rf = RobinsonFoulds(tr, trs),
      rfi = InfoRobinsonFoulds(tr, trs),
      qd = Quartet::QuartetDivergence(Quartet::QuartetStatus(trs, cf=tr),
                                      similarity = FALSE),
      path = path.dist(tr, trs)
    )
  }, matrix(0, nrow = 10L, ncol = length(sampledMethods),
            dimnames=list(subsamples, sampledMethods))
  )
  bullMoDiScores[[tipName]] <- theseScores
}
usethis::use_data(bullMoDiScores, compress='xz', overwrite=TRUE)
cat('# # # BullseyeMorphDistort COMPLETE # # #')
