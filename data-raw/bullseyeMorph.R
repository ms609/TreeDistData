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
    seq00 <- formatC(i, width=ceiling(log10(nTrees)), flag='0')
    seqFile <- paste0(tempdir(), '\\seq-', seq00, '.tnt')
    WriteTNTData(seqs[[i]], file = seqFile)

    runFile <- 'doscript.run'
    file.create(runFile)
    write(paste("macro =;
     xmult:hits 3 level 4 chklevel 5 rat10 drift10;
     sect:slack 8;
     keep 0; hold 10000;\n",
     "piwe=6 ;xmult;tsav *%1-k06-2000.tre;sav;tsav/;keep 0;hold 10000;\n",
     "ccode ] ", paste(1801:2000, collapse=' '), ";\n",
     "piwe=6 ;xmult;tsav *%1-k06-1800.tre;sav;tsav/;keep 0;hold 10000;\n",
     "ccode ] ", paste(1601:1800, collapse=' '), ";\n",
     "piwe=6 ;xmult;tsav *%1-k06-1600.tre;sav;tsav/;keep 0;hold 10000;\n",
     "ccode ] ", paste(1401:1600, collapse=' '), ";\n",
     "piwe=6 ;xmult;tsav *%1-k06-1400.tre;sav;tsav/;keep 0;hold 10000;\n",
     "ccode ] ", paste(1201:1400, collapse=' '), ";\n",
     "piwe=6 ;xmult;tsav *%1-k06-1200.tre;sav;tsav/;keep 0;hold 10000;\n",
     "ccode ] ", paste(1001:1200, collapse=' '), ";\n",
     "piwe=6 ;xmult;tsav *%1-k06-1000.tre;sav;tsav/;keep 0;hold 10000;\n",
     "ccode ] ", paste(801:1000, collapse=' '), ";\n",
     "piwe=6 ;xmult;tsav *%1-k06-0800.tre;sav;tsav/;keep 0;hold 10000;\n",
     "ccode ] ", paste(601:800, collapse=' '), ";\n",
     "piwe=6 ;xmult;tsav *%1-k06-0600.tre;sav;tsav/;keep 0;hold 10000;\n",
     "ccode ] ", paste(401:600, collapse=' '), ";\n",
     "piwe=6 ;xmult;tsav *%1-k06-0400.tre;sav;tsav/;keep 0;hold 10000;\n",
     "ccode ] ", paste(201:400, collapse=' '), ";\n",
     "piwe=6 ;xmult;tsav *%1-k06-0200.tre;sav;tsav/;keep 0;hold 10000;\n",
     "quit;"), runFile)

    # Install TNT and add to system path to make the next line work...
    system(paste('tnt proc', seqFile, '; doscript', seq00, ';'))
    file.remove(seqFile)
    file.remove(runFile)

    inferred[[i]] <- iInferred <- lapply(formatC(1:10 * 200, width=4, flag='0'), function (nChar) {
      treeFile <- paste0(seq00, '-k06-', nChar, '.tre')
      on.exit(file.remove(treeFile))
      ReadTntTree(treeFile, relativePath = '.',
                  tipLabels = theseTrees[[i]]$tip.label)
    })
    bullseyeMorphInferred[[nTip]] <- inferred
  }
}
usethis::use_data(bullseyeMorphInferred, compress='xz', overwrite=TRUE)


bullseyeMorphScores <- vector('list', length(nTips))
names(bullseyeMorphScores) <- nTips
for (trees in treesNames) {
  inferred <- bullseyeInferred[[trees]]
  theseTrees <- bullseyeTrees[[trees]]
  theseScores <- vapply(seq_along(inferred), function (i) {
    trueTree <- theseTrees[[i]]
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
  bullseyeScores[[trees]] <- theseScores
}
usethis::use_data(bullseyeScores, compress='xz', overwrite=TRUE)
