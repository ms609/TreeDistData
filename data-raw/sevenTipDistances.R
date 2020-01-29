library('TreeTools')
library('TreeDist')
nTip <- 7L
trees <- as.phylo(seq_len(945L) - 1, tipLabels = letters[seq_len(nTip)])
trees <- lapply(trees, ape::root, 'a', resolve.root = TRUE)
treeShapes <- vapply(trees, UnrootedTreeShape, 0L)
trees <- trees[order(treeShapes)]

cat(table(treeShapes))

sevenTipDistances <- TreeDistData::CompareAllTrees(trees, verbose = TRUE)

usethis::use_data(sevenTipDistances, compress='xz', overwrite = TRUE)
