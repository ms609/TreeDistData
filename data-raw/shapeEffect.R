require('TreeTools')
require('TreeDist')
require('TreeDistData')
nSample <- 100L
nTip <- 8L

cat(paste0(seq_along(TDFunctions), ': ', names(TDFunctions), '\n'))


tipLabels <- paste0('t', seq_len(nTip))
nShapes <- NUnrootedShapes(nTip)
shapeKeys <- UnrootedKeys(nTip)
shapes <- lapply(shapeKeys, UnrootedTreeWithKey, nTip = nTip)
trees <- structure(unlist(lapply(shapes, function(skeleton) {
  unique(lapply(rep(0, nSample * 1.1), function(xx) {
    skeleton$tip.label <- sample(tipLabels)
    skeleton
  }))[seq_len(nSample)]
}), recursive = FALSE), class = 'multiPhylo')
shapeNumbers <- rep(seq_along(shapeKeys) - 1L, each = nSample)

dists <- CompareAllTrees(trees, verbose = TRUE)

CalcShapeEffect <- function (dist) {
  message('...')
  ShapeDists <- function (i, j) {
    d <- as.matrix(dist)[shapeNumbers == i, shapeNumbers == j]
    d[upper.tri(d, diag = (i != j))]
  }
  go <- matrix(seq_len(nShapes) - 1L, nShapes, nShapes)
  x <- t(go)[lower.tri(go, diag = TRUE)]
  y <- go[lower.tri(go, diag = TRUE)]
  ret <- mapply(ShapeDists, x, y)
}

shapeEffect <- lapply(dists, CalcShapeEffect)

usethis::use_data(shapeEffect, compress='xz', overwrite = TRUE)
