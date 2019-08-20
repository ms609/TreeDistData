library('phangorn')
library('TreeSearch')
library('TreeDist')

data("bullseyeTrees", package='TreeDistData') # Generated in bullseye.R
tipsNames <- names(bullseyeTrees)
subsamples <- 0:9 * 2 # Order in increasing dissimilarity, please

wd <- getwd()
if (substr(wd, nchar(wd) - 7, nchar(wd)) != 'data-raw') setwd('data-raw')
###################

# Get results for a subset of trees:
tipsNames <- tipsNames[1:3]
bullseyeTrees <- bullseyeTrees[tipsNames]
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

# Create subdirectories:
if (!dir.exists('bullMoDi')) dir.create('bullMoDi')

bullMoDiInferred <- vector('list', length(tipsNames))
names(bullMoDiInferred) <- tipsNames

for (tipName in names(bullseyeTrees)) {
  theseTrees <- bullseyeTrees[[tipName]][seq_len(nTrees)]
  simChar <- 2000L
  seqs <- lapply(theseTrees, simSeq, l = simChar, type='USER', levels=0:1)
  nData <- length(seqs[[1]]) * simChar
  inferred <- vector(mode='list', nTrees)


  for (i in seq_along(seqs)) {
    seq00 <- formatC(i - 1, width=3, flag='0')
    FilePattern <- function (n) {
      paste0('bullMoDi/', substr(tipName, 0, nchar(tipName) - 5),
             't-', seq00, '-k6-',
             formatC(n, width=2, flag='0'),
             '.tre')
    }

    if (!file.exists(FilePattern(subsamples[(length(subsamples))]))) {
      seqFile <- paste0(tempdir(), '\\bullMoDi-', seq00, '.tnt')
      runRoot <- paste0(sample(letters, 8, replace=TRUE), collapse='')
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
      lapply(formatC(subsamples, width=2, flag='0'),
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
usethis::use_data(bullMoDiInferred, compress='xz', overwrite=TRUE)


bullMoDiScores <- vector('list', length(tipsNames))
names(bullMoDiScores) <- tipsNames
for (tipName in tipsNames) {
  inferred <- bullMoDiInferred[[tipName]]
  trueTrees <- bullseyeTrees[[tipName]]
  theseScores <- vapply(seq_along(inferred), function (i) {
    trueTree <- trueTrees[[i]]
    rootTip <- trueTree$tip.label[1]
    tr <- root(trueTree, rootTip, resolve.root=TRUE)
    tr$edge.length  <- NULL
    trs <- lapply(inferred[[i]], root, rootTip, resolve.root=TRUE)

    normInfo <- PartitionInfo(tr)
    cbind(
      mpi = 1 - MutualPhylogeneticInfo(tr, trs, normalize=normInfo),
      vpi = VariationOfPhylogeneticInfo(tr, trs, normalize=TRUE),
      mmsi = 1 - MutualMatchingSplitInfo(tr, trs, normalize=normInfo),
      vmsi = VariationOfMatchingSplitInfo(tr, trs, normalize=TRUE),
      mci = 1 - MutualClusteringInfo(tr, trs, normalize=normInfo),
      vci = VariationOfClusteringInfo(tr, trs, normalize=TRUE),
      qd = Quartet::QuartetDivergence(Quartet::QuartetStatus(trs, cf=tr), similarity = FALSE),
      nts = 1 - NyeTreeSimilarity(tr, trs, normalize=TRUE),
      msd = MatchingSplitDistance(tr, trs),
      t(vapply(trs, phangorn::treedist, tree2=tr, double(2))),
      spr = vapply(trs, phangorn::SPR.dist, tree2=tr, double(1))
    )
  }, matrix(0, nrow = 10L, ncol=12L,
            dimnames=list(subsamples,
                          c('mpi', 'vpi', 'mmsi', 'vmsi', 'mci', 'vci',
                            'qd', 'nts', 'msd', 'rf', 'path', 'spr')
  )))
  bullMoDiScores[[tipName]] <- theseScores
}
usethis::use_data(bullMoDiScores, compress='xz', overwrite=TRUE)
cat('# # # COMPLETE # # #')
