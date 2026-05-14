# rebuild_acid_aphid.R -----------------------------------------------------
# Full re-sample of `cid` and `pid` cells for n in 8..28 using the *corrected*
# TreeTools::RandomTree(root = TRUE) sampler (acid-aphid branch of TreeTools).
#
# The original 1000-rep estimates and the earlier topup pool were drawn from
# the buggy non-uniform sampler. This script discards both and replaces the
# mean/sd cells (plus n_repls_* attributes and sidecar .rds) with fresh
# samples whose target is SE/mean <= 0.001.
#
# Run with the acid-aphid TreeTools worktree on the search path so the fix
# is in effect.

suppressPackageStartupMessages({
  devtools::load_all("C:/Users/pjjg18/GitHub/worktrees/acid-aphid-TreeTools",
                     quiet = TRUE)
  library("TreeDist")
  library("usethis")
})
RNGversion("3.6.0")

target <- 1e-3
safety <- 1.20
ns <- 8:28
absolute_cap <- 50000L
seed_base <- 2026L

paths <- fs::path("data", "randomTreeDistances", ext = "rda")
stopifnot(file.exists(proj_path(paths)))
load(proj_path(paths))

ns_chr <- dimnames(randomTreeDistances)[[3]]
n_repls_cid <- attr(randomTreeDistances, "n_repls_cid")
n_repls_pid <- attr(randomTreeDistances, "n_repls_pid")
if (is.null(n_repls_cid)) n_repls_cid <- setNames(rep(0L, length(ns_chr)), ns_chr)
if (is.null(n_repls_pid)) n_repls_pid <- setNames(rep(0L, length(ns_chr)), ns_chr)

DrawCidPid <- function(n, reps, seed) {
  set.seed(seed)
  cid <- numeric(reps); pid <- numeric(reps)
  for (i in seq_len(reps)) {
    t1 <- RandomTree(n, TRUE); t2 <- RandomTree(n, TRUE)
    cid[i] <- ClusteringInfoDistance(t1, t2, normalize = TRUE)
    pid[i] <- PhylogeneticInfoDistance(t1, t2, normalize = TRUE)
  }
  rbind(cid = cid, pid = pid)
}

# Pilot reps to estimate sd, then top up to hit target.
pilot <- 1000L
t0 <- Sys.time()
log_rows <- list()
for (n in ns) {
  tic <- Sys.time()
  raw <- DrawCidPid(n, pilot, seed_base + n)
  sd_cid <- sd(raw["cid", ]); m_cid <- mean(raw["cid", ])
  sd_pid <- sd(raw["pid", ]); m_pid <- mean(raw["pid", ])
  need <- max(ceiling((sd_cid / (target * m_cid))^2),
              ceiling((sd_pid / (target * m_pid))^2))
  need <- min(ceiling(need * safety), absolute_cap)
  more <- max(0L, need - pilot)
  if (more > 0L) {
    raw2 <- DrawCidPid(n, more, seed_base + n + 100000L)
    raw <- cbind(raw, raw2)
  }
  N <- ncol(raw)
  m_cid <- mean(raw["cid", ]); sd_cid <- sd(raw["cid", ])
  m_pid <- mean(raw["pid", ]); sd_pid <- sd(raw["pid", ])
  randomTreeDistances["cid", "mean", as.character(n)] <- m_cid
  randomTreeDistances["cid", "sd",   as.character(n)] <- sd_cid
  randomTreeDistances["pid", "mean", as.character(n)] <- m_pid
  randomTreeDistances["pid", "sd",   as.character(n)] <- sd_pid
  n_repls_cid[as.character(n)] <- N
  n_repls_pid[as.character(n)] <- N
  saveRDS(list(n = n, reps = N, raw = raw),
          proj_path("data-raw", sprintf("cid_pid_topup_n%02d.rds", n)))
  cid_se <- (sd_cid / sqrt(N)) / m_cid
  pid_se <- (sd_pid / sqrt(N)) / m_pid
  toc <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
  cat(sprintf("n=%2d  reps=%5d  cid mean=%.4f SE/m=%.5f  pid mean=%.4f SE/m=%.5f  (%.1fs)\n",
              n, N, m_cid, cid_se, m_pid, pid_se, toc))
  flush.console()
  log_rows[[as.character(n)]] <- data.frame(n = n, reps = N,
    cid_mean = m_cid, cid_se_mean = cid_se,
    pid_mean = m_pid, pid_se_mean = pid_se)
}

attr(randomTreeDistances, "n_repls_cid") <- n_repls_cid
attr(randomTreeDistances, "n_repls_pid") <- n_repls_pid
attr(randomTreeDistances, "topup_note") <- paste0(
  "cid and pid mean/sd cells for n in 8..28 re-sampled from scratch ",
  "against corrected TreeTools::RandomTree (acid-aphid branch) to ",
  "SE/mean <= ", target, " via data-raw/rebuild_acid_aphid.R. ",
  "Quantile cells NOT updated. Cells for n in 29..200 retain original ",
  "1000-rep estimates drawn under the buggy non-uniform sampler ",
  "(see PLAN.md Findings to follow up)."
)

use_data(randomTreeDistances, compress = "gzip", overwrite = TRUE)
cat(sprintf("\nDone. Wall %.1f s.\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(do.call(rbind, log_rows), digits = 4, row.names = FALSE)
