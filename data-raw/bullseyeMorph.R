library('phangorn')
library('TreeSearch')
library('TreeDist')


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

bullseyeMorphInferred <- vector('list', length(tipsNames))
names(bullseyeMorphInferred) <- tipsNames

for (tipName in names(bullseyeTrees)) {
  theseTrees <- bullseyeTrees[[tipName]][seq_len(nTrees)]
  seqs <- lapply(theseTrees, simSeq, l = 2000, type='USER', levels=1:4)
  inferred <- vector(mode='list', nTrees)


  for (i in seq_along(seqs)) {
    seq00 <- formatC(i - 1, width=3, flag='0')
    FilePattern <- function (n) {
      paste0(substr(tipName, 0, nchar(tipName) - 5),
             't-', seq00, '-k6-',
             formatC(n, width=4, flag='0'),
             '.tre')
    }

    if (!file.exists(FilePattern(200))) {
      seqFile <- paste0(tempdir(), '\\seq-', seq00, '.tnt')
      WriteTNTData(seqs[[i]], file = seqFile)
      Line <- function (n) {
        paste0("piwe=6 ;xmult;tsav *", FilePattern(n),
               ";sav;tsav/;keep 0;hold 10000;\n",
        "ccode ] ", paste(seq_len(200) + n - 200L, collapse=' '), ";\n"
        )
      }

      runRoot <- paste0(sample(letters, 8, replace=TRUE), collapse='')
      runFile <- paste0(runRoot, '.run', collapse='')
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
               tr <- ReadTntTree(FilePattern(nChar),
                                 relativePath = '.',
                  tipLabels = theseTrees[[i]]$tip.label)
               # Return:
               if (class(tr) == 'multiPhylo') tr[[1]] else tr
               })
  }
  bullseyeMorphInferred[[tipName]] <- inferred
}
usethis::use_data(bullseyeMorphInferred, compress='xz', overwrite=TRUE)


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
  bullseyeMorphScores[[tipName]] <- theseScores
}
usethis::use_data(bullseyeMorphScores, compress='xz', overwrite=TRUE)
cat(" # # # COMPLETE # # # ")
