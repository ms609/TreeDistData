library('phangorn')
library('TreeSearch')
library('TreeDist')
set.seed(0)
nTrees <- 100L
nTip <- 15L

bullseyeTrees <- lapply(seq_len(nTrees), function (XX) ape::rtree(nTip))
usethis::use_data(bullseyeTrees, overwrite=TRUE)

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

seqs <- lapply(bullseyeTrees, simSeq, l = 2000, type='USER', levels=1:4)
bullseyeMorphInferred <- vector(mode='list', length = nTrees)

for (i in seq_along(seqs)) {
  seq00 <- formatC(i, width=ceiling(log10(nTrees)), flag='0')
  seqFile <- paste0(tempdir(), '\\seq-', seq00, '.tnt')
  WriteTNTData(seqs[[i]], file = seqFile)

  runFile <- paste0('doscript.run')
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

  bullseyeMorphInferred[[i]] <- iInferred <- lapply(formatC(1:10 * 200, width=4, flag='0'), function (nChar) {
    ReadTntTree(paste0(seq00, '-k06-', nChar, '.tre'), relativePath = '.',
                tipLabels = bullseyeTrees[[i]]$tip.label)
  })
}

usethis::use_data(bullseyeMorphInferred, compress='xz', overwrite=TRUE)

allScores <- vapply(seq_along(bullseyeMorphInferred), function (i) {
  trueTree <- bullseyeTrees[[i]]
  rootTip <- trueTree$tip.label[1]
  tr <- root(trueTree, rootTip, resolve.root=TRUE)
  tr$edge.length  <- NULL
  trs <- lapply(bullseyeMorphInferred[[i]], root, rootTip, resolve.root=TRUE)

  normInfo <- PartitionInfo(tr)
  scores <- cbind(
    vpi = VariationOfPhylogeneticInfo(tr, trs, normalize=normInfo),
    vmsi = VariationOfMatchingSplitInfo(tr, trs, normalize=normInfo),
    vci = VariationOfClusteringInfo(tr, trs, normalize=normInfo),
    qd = Quartet::QuartetDivergence(Quartet::QuartetStatus(trs, cf=tr), similarity = FALSE),
    nts = 1 - NyeTreeSimilarity(tr, trs, normalize=TRUE),
    msd = MatchingSplitDistance(tr, trs),
    t(vapply(trs, phangorn::treedist, tree2=tr, double(2))),
    spr = vapply(trs, phangorn::SPR.dist, tree2=tr, double(1))
  )
}, matrix(0, nrow = 10, ncol=9L))

rhos <- apply(allScores, 3L, function (scores) apply(scores, 2, function (x)
  cor.test(x,  y = seq_len(10L), method='spearman', exact=FALSE)$estimate))

rowMeans(rhos, na.rm=TRUE)

treeMeans <- apply(allScores, 1L, rowMeans)


par(mar=c(1,1, 1, 1))
plot(type='n', xlim=c(0, 10), ylim=0:1, axes=F, x=0, y=0)
lapply(1:9, function (i) {
  iScores <- treeMeans[i, ] / max(treeMeans[i, ])
  lines(iScores, col=i)
  }
)
legend('topright', legend=rownames(treeMeans), col=1:9, lty=1, bty='n')

