library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
data("distanceDistribution25", package = "TreeDistData")
data("randomTreePairs25", package = "TreeDistData")
data("distanceDistribution50", package = "TreeDistData")
data("randomTreePairs50", package = "TreeDistData")

message(Sys.time(), ": Calculating TCI...")
tciDistribution25 <- vapply(randomTreePairs25, function (treePair) {
  abs(TotalCopheneticIndex(treePair[[1]]) - TotalCopheneticIndex(treePair[[2]]))
}, 0)

message(Sys.time(), ": Correlating...")
balance25 <- apply(distanceDistribution25, 1, cor, tciDistribution25) ^ 2
usethis::use_data(balance25, compress = 'xz', overwrite = TRUE)
message("25 Complete.")

message(Sys.time(), ": Calculating TCI...")
tciDistribution50 <- vapply(randomTreePairs50, function (treePair) {
  abs(TotalCopheneticIndex(treePair[[1]]) - TotalCopheneticIndex(treePair[[2]]))
}, 0)

message(Sys.time(), ": Correlating...")
balance50 <- apply(distanceDistribution25, 1, cor, tciDistribution50) ^ 2
usethis::use_data(balance50, compress = 'xz', overwrite = TRUE)
message("50 Complete.")
