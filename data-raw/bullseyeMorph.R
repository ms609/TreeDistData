library('phangorn')
library('TreeSearch')
library('TreeDist')

data("bullseyeTrees", package='TreeDistData') # Generated in bullseye.R
tipsNames <- names(bullseyeTrees)
nTrees <- stop("nTrees Not set.")
subsamples <- 10:1 * 200

bullseyeMorphInferred <- vector('list', length(tipsNames))
names(bullseyeMorphInferred) <- tipsNames

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


CacheFile <- function (name, ..., tmpDir = FALSE) {
  root <- if (tmpDir) tempdir() else paste0(system.file(package='TreeDistData'),
                                            '/../data-raw/trees/')
  paste0(root, name, ...)
}

for (tipName in names(bullseyeTrees)) {
  theseTrees <- bullseyeTrees[[tipName]][seq_len(nTrees)]
  seqs <- lapply(theseTrees, simSeq, l = 2000, type='USER', levels=1:4)
  inferred <- vector(mode = 'list', nTrees)


  for (i in seq_along(seqs)) {
    seq00 <- formatC(i - 1L, width = 3, flag='0')
    FilePattern <- function (n) {
      paste0(substr(tipName, 0, nchar(tipName) - 5),
             't-', seq00, '-k6-',
             formatC(n, width=4, flag='0'),
             '.tre')
    }

    if (!file.exists(CacheFile(FilePattern(200)))) {
      seqFile <- CacheFile('seq-', seq00, '.tnt')
      WriteTNTData(seqs[[i]], file = seqFile)
      Line <- function (n) {
        paste0("piwe=6 ;xmult;tsav *", FilePattern(n),
               ";sav;tsav/;keep 0;hold 10000;\n",
        "ccode ] ", paste(seq_len(200) + n - 200L, collapse=' '), ";\n"
        )
      }

      runRoot <- paste0(sample(letters, 8, replace=TRUE), collapse='')
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
      lapply(formatC(subsamples, width=4, flag='0'),
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
usethis::use_data(bullseyeMorphInferred, compress='bzip2', overwrite=TRUE)


bullseyeMorphScores <- vector('list', length(tipsNames))
names(bullseyeMorphScores) <- tipsNames
for (tipName in tipsNames) {
  cat('\u2714 Calculating tree distances:', tipName, ':\n')
  inferred <- bullseyeMorphInferred[[tipName]]
  trueTrees <- bullseyeTrees[[tipName]]
  theseScores <- vapply(seq_along(inferred), function (i) {
    cat('.')
    trueTree <- trueTrees[[i]]
    rootTip <- trueTree$tip.label[1]
    tr <- root(trueTree, rootTip, resolve.root=TRUE)
    tr$edge.length  <- NULL
    trs <- lapply(inferred[[i]], root, rootTip, resolve.root=TRUE)

    mast <- lapply(trs, MASTSize, tr, rooted = FALSE)
    masti <-  LnUnrooted(mast) / log(2)
    attributes(masti) <- attributes(mast)

    nni <- NNIDist(tr, trs)
    tbr <- TBRDist(tr, trs)

    normInfo <- PartitionInfo(tr)
    cbind(
      spi = 1 - SharedPhylogeneticInfo(tr, trs, normalize=normInfo),
      dpi = DifferentPhylogeneticInfo(tr, trs, normalize=TRUE),
      msi = 1 - MatchingSplitInfo(tr, trs, normalize=normInfo),
      msid = MatchingSplitInfoDistance(tr, trs, normalize=TRUE),
      mci = 1 - MutualClusteringInfo(tr, trs, normalize=normInfo),
      cid = ClusteringInfoDistance(tr, trs, normalize=TRUE),
      qd = Quartet::QuartetDivergence(Quartet::QuartetStatus(trs, cf=tr), similarity = FALSE),

      ja2 = JaccardRobinsonFoulds(tr1, tr2, k = 2, arboreal = TRUE, normalize = TRUE),
      ja4 = JaccardRobinsonFoulds(tr1, tr2, k = 4, arboreal = TRUE, normalize = TRUE),
      jna2 =JaccardRobinsonFoulds(tr1, tr2, k = 2, arboreal = FALSE, normalize = TRUE),
      jna4 =JaccardRobinsonFoulds(tr1, tr2, k = 4, arboreal = FALSE, normalize = TRUE),

      mast = mast,
      masti = masti,
      nni_l = nni[['lower']],
      nni_t = nni[['tight_upper']],
      nni_u = nni[['loose_upper']],
      spr = SPR.dist(tr1, tr2)[['spr']],
      tbr_l = tbr$tbr_min,
      tbr_u = tbr$tbr_max,
      rf = RobinsonFoulds(tr1, tr2),
      path = path.dist(tr1, tr2)
      nts = NyeTreeSimilarity(tr1, tr2, similarity = FALSE, normalize = TRUE),
      msd = MatchingSplitDistance(tr, trs),
      t(vapply(trs, phangorn::treedist, tree2=tr, double(2))),
      spr = vapply(trs, phangorn::SPR.dist, tree2=tr, double(1))
    )
  }, matrix(0, nrow = 10L, ncol=12L,
            dimnames=list(subsamples,
                          c('spi', 'dpi', 'msi', 'msid', 'mci', 'cid',
                            'qd', 'nts', 'msd', 'rf', 'path', 'spr')
  )))
  bullseyeMorphScores[[tipName]] <- theseScores
}
usethis::use_data(bullseyeMorphScores, compress='xz', overwrite=TRUE)
cat(" # # # COMPLETE # # # ")
