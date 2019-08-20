library('phangorn')
library('TreeSearch')
library('TreeDist')

#data("bullseyeTrees", package='TreeDistData') # Generated in bullseye.R
nTips <- names(bullseyeTrees)

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

SubSample <- function (dat, n) {
  at <- attributes(dat)
  tokens <- t(vapply(seq_along(dat), function(taxon) {
    dat[[taxon]][at$index[seq_len(n)]]
  }, double(n)))
  rownames(tokens) <- names(dat)
  ret <- MatrixToPhyDat(tokens)
  attr(ret, 'levels') <- at$levels
  ret
}

bullseyeMorphInferred <- vector('list', length(nTips))
names(bullseyeMorphInferred) <- nTips
nTip = '5 tips'
for (nTip in names(bullseyeTrees)) {
  theseTrees <- bullseyeTrees[[nTip]]
  seqs <- lapply(theseTrees, simSeq, l = 2000, type='USER', levels=1:4)
  inferred <- vector(mode='list', length(theseTrees))


  for (i in seq_along(seqs)) {
    seq00 <- formatC(i - 1, width=3, flag='0')
    filePattern <- paste0(substr(nTip, 0, nchar(nTip) - 5),
                       't-', seq00, '-k6-%s.tre')

    if (!file.exists(sprintf(filePattern, 2000))) {
      seqFile <- paste0(tempdir(), '\\seq-', seq00, '.tnt')
      WriteTNTData(seqs[[i]], file = seqFile)
      Line <- function (n) {
        paste0("piwe=6 ;xmult;tsav *",
               sprintf(filePattern, formatC(n, width=4, flag='0')),
        ";sav;tsav/;keep 0;hold 10000;\n",
        "ccode ] ", paste(seq_len(200) + n - 200L, collapse=' '), ";\n"
        )
      }

      runFile <- 'doscript.run'
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
      system(paste('tnt proc', seqFile, '; doscript;'))
      file.remove(seqFile)
      file.remove(runFile)
    }

    inferred[[i]] <- iInferred <-
      lapply(formatC(1:10 * 200, width=4, flag='0'),
             function (nChar) {
               tr <- ReadTntTree(sprintf(filePattern, nChar),
                                 relativePath = '.',
                  tipLabels = theseTrees[[i]]$tip.label)
               # Return:
               if (class(tr) == 'multiPhylo') tr[[1]] else tr
               })
    bullseyeMorphInferred[[nTip]] <- inferred
  }
}
usethis::use_data(bullseyeMorphInferred, compress='xz', overwrite=TRUE)


bullseyeMorphScores <- vector('list', length(nTips))
names(bullseyeMorphScores) <- nTips
for (nTip in nTips) {
  inferred <- bullseyeMorphInferred[[nTip]]
  trueTrees <- bullseyeTrees[[nTip]]
  theseScores <- vapply(seq_along(inferred), function (i) {
    trueTree <- trueTrees[[i]]
    rootTip <- trueTree$tip.label[1]
    tr <- root(trueTree, rootTip, resolve.root=TRUE)
    tr$edge.length  <- NULL
    trs <- lapply(inferred[[i]], root, rootTip, resolve.root=TRUE)

    normInfo <- PartitionInfo(tr)
    cbind(
      vpi = VariationOfPhylogeneticInfo(tr, trs, normalize=normInfo),
      vmsi = VariationOfMatchingSplitInfo(tr, trs, normalize=normInfo),
      vci = VariationOfClusteringInfo(tr, trs, normalize=normInfo),
      qd = Quartet::QuartetDivergence(Quartet::QuartetStatus(trs, cf=tr), similarity = FALSE),
      nts = 1 - NyeTreeSimilarity(tr, trs, normalize=TRUE),
      msd = MatchingSplitDistance(tr, trs),
      t(vapply(trs, phangorn::treedist, tree2=tr, double(2))),
      spr = vapply(trs, phangorn::SPR.dist, tree2=tr, double(1))
    )
  }, matrix(0, nrow = 10L, ncol=9L, dimnames=list(
    10:1 * 200,
    c('vpi', 'vmsi', 'vci', 'qd', 'nts', 'msd', 'rf', 'path', 'spr')
  )))
  bullseyeMorphScores[[nTip]] <- theseScores
}
usethis::use_data(bullseyeMorphScores, compress='xz', overwrite=TRUE)
