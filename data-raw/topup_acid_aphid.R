# topup_acid_aphid.R --------------------------------------------------------
# Top-up `cid` and `pid` cells of `randomTreeDistances` to SE/mean <= 0.001.
#
# Background
# ----------
# The original builder (`data-raw/randomTreeDistances.R`) ran 1000 replicates
# per `n` for all 24 methods. The audit
# (`acid-aphid/analysis/reference-data-audit.md`) showed this leaves SE/mean
# above the 0.001 target for small `n`. The ACID/APhID adjustment only needs
# `cid` and `pid`, so we top up those two columns only.
#
# Scope
# -----
# - n in 4..7 is handled by a parallel agent via explicit enumeration; this
#   script may touch those cells (it does not by default).
# - n where current 1000-rep SE/mean already meets target (n >= ~25) is
#   skipped automatically by the rep-count formula.
#
# What this script does
# ---------------------
# 1. Loads `data/randomTreeDistances.rda` from this worktree.
# 2. For each n in `topup_range`, computes how many ADDITIONAL replicates
#    are needed to bring SE/mean for cid and pid below `target`, using
#       reps_needed ~= (sd / (target * mean))^2
#    with a `safety` multiplier to absorb sd-estimation noise.
# 3. Generates that many additional random tree pairs via
#    `TreeTools::RandomTree(n, TRUE)` and computes only CID and PID via
#    `TreeDist::ClusteringInfoDistance` / `PhylogeneticInfoDistance`.
#    (This skips the expensive `AllDists()` machinery.)
# 4. Saves the raw new distances to `data-raw/cid_pid_topup_n<NN>.rds`.
# 5. Combines the new sample with the existing 1000-rep summary stats
#    via the pooled mean / pooled variance formulae (numerically stable
#    two-sample combination) and updates `randomTreeDistances['cid','mean',]`,
#    `['cid','sd',]`, `['pid','mean',]`, `['pid','sd',]` only.
# 6. Records the new effective replicate count per n in attributes
#    `n_repls_cid` and `n_repls_pid`.
# 7. NOTE: Does NOT update the quantile cells (`1%`, `5%`, ..., `99%`,
#    `min`, `max`) — these cannot be derived from summary stats alone.
#    The raw sidecar .rds files preserve the new data so future runs can
#    re-derive everything from a pooled raw distance vector.
# 8. Persists with `usethis::use_data(..., compress = 'gzip', ...)`
#    matching the existing builder.

suppressPackageStartupMessages({
  library("TreeTools")
  library("TreeDist")
  library("usethis")
})
RNGversion("3.6.0")

# -- configuration ---------------------------------------------------------

target <- 1e-3          # SE/mean target for cid and pid
existing_repls <- 1000  # The hard-coded `repls <- 1000` from the original
                        # builder; verified empirically in the audit.
safety <- 1.20          # 20 % over the analytic minimum to absorb sd noise.
topup_range <- 8:28     # n in 4..7 covered by enumeration elsewhere.
absolute_cap <- 100000  # Per-n cap on total reps (safety rail).

paths <- fs::path("data", "randomTreeDistances", ext = "rda")
stopifnot(file.exists(proj_path(paths)))
load(proj_path(paths))

# -- helpers ---------------------------------------------------------------

#' Reps required for SE/mean <= target, given sd and mean.
RepsNeeded <- function(sd, mean, target = 1e-3) {
  ceiling((sd / (target * mean))^2)
}

#' Two-sample pooled mean & sd. n1, m1, s1 from existing summary; n2 from
#' newly drawn raw vector x2.
PoolStats <- function(n1, m1, s1, x2) {
  n2 <- length(x2)
  if (n2 == 0L) return(list(n = n1, mean = m1, sd = s1))
  m2 <- mean(x2)
  v1 <- s1^2
  v2 <- if (n2 > 1L) var(x2) else 0
  N <- n1 + n2
  mN <- (n1 * m1 + n2 * m2) / N
  # Pooled sample variance (divisor N - 1):
  # SS = (n1 - 1) v1 + (n2 - 1) v2 + n1 (m1 - mN)^2 + n2 (m2 - mN)^2
  SS <- (n1 - 1) * v1 + (n2 - 1) * v2 +
    n1 * (m1 - mN)^2 + n2 * (m2 - mN)^2
  sdN <- sqrt(SS / (N - 1))
  list(n = N, mean = mN, sd = sdN)
}

#' Draw `reps` random pairs at n leaves and return a matrix with rows
#' "cid" and "pid".
DrawCidPid <- function(n, reps, seed) {
  set.seed(seed)
  cid <- numeric(reps)
  pid <- numeric(reps)
  for (i in seq_len(reps)) {
    t1 <- RandomTree(n, TRUE)
    t2 <- RandomTree(n, TRUE)
    cid[i] <- ClusteringInfoDistance(t1, t2, normalize = TRUE)
    pid[i] <- PhylogeneticInfoDistance(t1, t2, normalize = TRUE)
  }
  rbind(cid = cid, pid = pid)
}

# -- existing per-n rep counters (attribute) -------------------------------
# Carry forward if present; otherwise initialise from existing_repls.
ns_chr <- dimnames(randomTreeDistances)[[3]]
n_repls_cid <- attr(randomTreeDistances, "n_repls_cid")
n_repls_pid <- attr(randomTreeDistances, "n_repls_pid")
if (is.null(n_repls_cid)) {
  n_repls_cid <- setNames(rep(existing_repls, length(ns_chr)), ns_chr)
}
if (is.null(n_repls_pid)) {
  n_repls_pid <- setNames(rep(existing_repls, length(ns_chr)), ns_chr)
}

# -- main top-up loop ------------------------------------------------------

t0 <- Sys.time()
total_added <- 0L
log_rows <- list()

for (n in topup_range) {
  n_chr <- as.character(n)
  cur_n_cid <- as.integer(n_repls_cid[n_chr])
  cur_n_pid <- as.integer(n_repls_pid[n_chr])
  m_cid <- randomTreeDistances["cid", "mean", n_chr]
  s_cid <- randomTreeDistances["cid", "sd",   n_chr]
  m_pid <- randomTreeDistances["pid", "mean", n_chr]
  s_pid <- randomTreeDistances["pid", "sd",   n_chr]
  need_cid <- RepsNeeded(s_cid, m_cid, target)
  need_pid <- RepsNeeded(s_pid, m_pid, target)
  need_total <- max(need_cid, need_pid)
  # Apply safety multiplier and clip to absolute_cap; both cid and pid
  # are drawn from the same pair so add reps = max requirement.
  target_total <- min(ceiling(need_total * safety), absolute_cap)
  to_add <- max(0L, target_total - max(cur_n_cid, cur_n_pid))
  if (to_add == 0L) {
    log_rows[[n_chr]] <- data.frame(
      n = n, to_add = 0L, total_after = max(cur_n_cid, cur_n_pid),
      pid_se_mean_after = (s_pid / sqrt(cur_n_pid)) / m_pid,
      cid_se_mean_after = (s_cid / sqrt(cur_n_cid)) / m_cid
    )
    next
  }
  cat(sprintf("[%s] n=%d  need cid=%d pid=%d  adding=%d (safety x%.2f)\n",
              format(Sys.time(), "%H:%M:%S"), n, need_cid, need_pid,
              to_add, safety))
  flush.console()

  raw <- DrawCidPid(n, to_add, seed = 1000L + n)
  total_added <- total_added + to_add

  sidecar <- proj_path("data-raw", sprintf("cid_pid_topup_n%02d.rds", n))
  saveRDS(list(n = n, seed = 1000L + n, reps = to_add, raw = raw),
          sidecar)

  pooled_cid <- PoolStats(cur_n_cid, m_cid, s_cid, raw["cid", ])
  pooled_pid <- PoolStats(cur_n_pid, m_pid, s_pid, raw["pid", ])

  randomTreeDistances["cid", "mean", n_chr] <- pooled_cid$mean
  randomTreeDistances["cid", "sd",   n_chr] <- pooled_cid$sd
  randomTreeDistances["pid", "mean", n_chr] <- pooled_pid$mean
  randomTreeDistances["pid", "sd",   n_chr] <- pooled_pid$sd
  n_repls_cid[n_chr] <- pooled_cid$n
  n_repls_pid[n_chr] <- pooled_pid$n

  log_rows[[n_chr]] <- data.frame(
    n = n, to_add = to_add, total_after = pooled_cid$n,
    pid_se_mean_after = (pooled_pid$sd / sqrt(pooled_pid$n)) / pooled_pid$mean,
    cid_se_mean_after = (pooled_cid$sd / sqrt(pooled_cid$n)) / pooled_cid$mean
  )
}

# -- write attributes and persist -----------------------------------------

attr(randomTreeDistances, "n_repls_cid") <- n_repls_cid
attr(randomTreeDistances, "n_repls_pid") <- n_repls_pid
attr(randomTreeDistances, "topup_note") <- paste0(
  "cid and pid mean/sd cells updated to SE/mean <= ", target,
  " via data-raw/topup_acid_aphid.R. ",
  "Quantile cells (1% .. 99%, min, max) NOT updated — they remain ",
  "the original 1000-rep estimates. Raw new distances are in ",
  "data-raw/cid_pid_topup_n<NN>.rds for full re-derivation if needed."
)

use_data(randomTreeDistances, compress = "gzip", overwrite = TRUE)

cat(sprintf("\nTop-up complete. Wall time: %.1f s. Total reps added: %d.\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            total_added))

# -- verification: SE/mean at every n in the array -----------------------

cat("\nVerification: SE/mean for cid and pid across all n (showing n<=30 + worst)\n")
all_ns <- as.integer(dimnames(randomTreeDistances)[[3]])
chk <- data.frame(
  n = all_ns,
  pid_se_mean = NA_real_,
  cid_se_mean = NA_real_,
  pid_repls = NA_integer_,
  cid_repls = NA_integer_
)
for (i in seq_along(all_ns)) {
  nc <- as.character(all_ns[i])
  chk$pid_repls[i] <- as.integer(n_repls_pid[nc])
  chk$cid_repls[i] <- as.integer(n_repls_cid[nc])
  chk$pid_se_mean[i] <- (randomTreeDistances["pid", "sd", nc] /
                          sqrt(chk$pid_repls[i])) /
                         randomTreeDistances["pid", "mean", nc]
  chk$cid_se_mean[i] <- (randomTreeDistances["cid", "sd", nc] /
                          sqrt(chk$cid_repls[i])) /
                         randomTreeDistances["cid", "mean", nc]
}
print(head(chk, 30), digits = 4, row.names = FALSE)
cat("\nWorst SE/mean across all n (excluding n in 4:7 — enumeration):\n")
mask <- chk$n >= 8
worst_idx_pid <- which.max(chk$pid_se_mean[mask])
worst_idx_cid <- which.max(chk$cid_se_mean[mask])
chk_sub <- chk[mask, ]
cat(sprintf("  PID worst: n=%d  SE/mean = %.6f  (target <= %g)\n",
            chk_sub$n[worst_idx_pid], chk_sub$pid_se_mean[worst_idx_pid], target))
cat(sprintf("  CID worst: n=%d  SE/mean = %.6f  (target <= %g)\n",
            chk_sub$n[worst_idx_cid], chk_sub$cid_se_mean[worst_idx_cid], target))

cat("\nLog of top-ups performed:\n")
print(do.call(rbind, log_rows), digits = 4, row.names = FALSE)

if (any(chk_sub$pid_se_mean > target) || any(chk_sub$cid_se_mean > target)) {
  cat("\nNOTE: Some cells still exceed target — see table above.\n")
} else {
  cat("\nAll cells (n >= 8) meet SE/mean <= ", target, " target.\n", sep = "")
}
