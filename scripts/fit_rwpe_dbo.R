# rWPE (randomized weighted permutation entropy) applied to DBO time series.
# Sen et al. 2024, Methods Ecol Evol 15:1834-1846.
#
# Purpose here is NOT inference about whether a series is "truly" stochastic. It is to
# establish a PRACTICAL CEILING: if a series is indistinguishable from white noise at the
# length we actually have, no model can forecast it, and that bounds what the SDM analysis
# could ever have achieved. No multiple-comparison correction is applied -- no per-series
# claims are made, and correcting would only push more series into the "noise" verdict the
# conclusion already points toward.
#
# WPE: for each window of m CONSECUTIVE years, take its ordinal pattern and weight it by
# the window variance (so low-amplitude wiggles count less than real swings):
#   pw(pi) = sum(var(Xt) * 1[phi(Xt)=pi]) / sum(var(Xt));  WPE = -sum pw log2 pw / log2(m!)
# Normalized to [0,1]: 0 = perfectly predictable, 1 = maximally complex.
#
# rWPE test: shuffle the observations (breaking temporal order, approximating white noise),
# recompute WPE, and take p = P(WPE_null <= WPE_obs), one-sided. Small p => the ORDERING of
# the observations makes the series more predictable than white noise.
#
# GAPS: DBO series are broken -- a station-year is missing when unsampled, and a family-year
# when the family is absent (biomass = 0 is dropped upstream, matching the canonical
# analysis). Permutation entropy reads CONSECUTIVE ordinal patterns, so only windows whose
# m years are truly consecutive are used. The null reuses the same window positions, so the
# sampling structure is held fixed and only the temporal ordering is broken.
#
# Series tested: family log-biomass and environmental variables, each at station level and
# at regional level (DBO 1/2/3).
#
# Output: data/rwpe_dbo_results.rds

library(tidyverse)

M <- 3          # word length; m! = 6 patterns, appropriate for 16-year series
NPERM <- 999
MIN_N <- 10     # minimum series length (paper's lower bound)
MIN_W <- 6      # minimum number of valid consecutive windows

# --- WPE over a fixed set of windows (idx = precomputed window index matrix) ----------
# m=3 is vectorised: the ordinal pattern of (a,b,c) is uniquely encoded by the three
# pairwise comparisons, and the window variance has a closed form. This runs ~1.1M times
# across the full analysis, so the apply()-free path matters.
wpe_core <- function(v, idx, m) {
  W <- matrix(v[idx], ncol = m)
  if (m == 3L) {
    a <- W[, 1]; b <- W[, 2]; cc <- W[, 3]
    code <- (a < b) + 2L * (a < cc) + 4L * (b < cc)
    mu <- (a + b + cc) / 3
    wt <- ((a - mu)^2 + (b - mu)^2 + (cc - mu)^2) / 2
  } else {
    code <- apply(W, 1, function(w) paste(order(w), collapse = ""))
    wt <- apply(W, 1, stats::var)
  }
  tot <- sum(wt)
  if (!is.finite(tot) || tot <= 0) return(NA_real_)
  s <- rowsum(wt, code, reorder = FALSE)
  p <- s[s > 0] / tot
  -sum(p * log2(p)) / log2(factorial(m))
}

# --- randomized test -----------------------------------------------------------------
rwpe <- function(vals, years, m = M, nperm = NPERM) {
  o <- order(years); vals <- vals[o]; years <- years[o]
  keep <- is.finite(vals); vals <- vals[keep]; years <- years[keep]
  n <- length(vals)
  if (n < MIN_N) return(NULL)
  starts <- which(vapply(seq_len(n - m + 1),
    function(i) all(diff(years[i:(i + m - 1)]) == 1), logical(1)))
  if (length(starts) < MIN_W) return(NULL)
  idx <- outer(starts, 0:(m - 1), "+")          # precompute once, reuse every permutation
  obs <- wpe_core(vals, idx, m)
  if (!is.finite(obs)) return(NULL)
  null <- replicate(nperm, wpe_core(sample(vals), idx, m))
  null <- null[is.finite(null)]
  if (!length(null)) return(NULL)
  # +1/+1 permutation p-value (never exactly 0)
  tibble(n = n, n_win = length(starts), wpe = obs,
    wpe_null = median(null), p = (sum(null <= obs) + 1) / (length(null) + 1))
}

set.seed(1)
d <- readRDS("data/data_dbo.rds") %>% filter(biomass > 0) %>% mutate(logN = log(biomass))
base <- c("Temp", "Salinity", "integchla", "sedchla", "Ammonia", "Phosphate",
          "NiTriTra", "Silicate", "phigte5", "TOC", "cn")

# ---- 1. family biomass, STATION level ------------------------------------------------
fam_stn <- d %>% select(family, DBOreg, StationNme, DataYear, logN) %>%
  group_by(family, StationNme) %>% group_split() %>%
  map_dfr(function(g) {
    r <- rwpe(g$logN, g$DataYear); if (is.null(r)) return(NULL)
    bind_cols(tibble(group = "Family biomass", level = "Station",
      series = paste(g$family[1], g$StationNme[1]), var = g$family[1],
      DBOreg = g$DBOreg[1]), r)
  })

# ---- 2. family biomass, REGIONAL level (>=2 stations per region-year) -----------------
fam_reg <- d %>% group_by(family, DBOreg, DataYear) %>%
  filter(n() >= 2) %>% summarise(v = mean(logN), .groups = "drop") %>%
  group_by(family, DBOreg) %>% group_split() %>%
  map_dfr(function(g) {
    r <- rwpe(g$v, g$DataYear); if (is.null(r)) return(NULL)
    bind_cols(tibble(group = "Family biomass", level = "Regional",
      series = paste(g$family[1], g$DBOreg[1]), var = g$family[1],
      DBOreg = g$DBOreg[1]), r)
  })

# ---- 3. environment, STATION level ----------------------------------------------------
env <- d %>% distinct(StationNme, DataYear, DBOreg, across(all_of(base))) %>%
  pivot_longer(all_of(base), names_to = "var", values_to = "v")
env_stn <- env %>% group_by(var, StationNme) %>% group_split() %>%
  map_dfr(function(g) {
    r <- rwpe(g$v, g$DataYear); if (is.null(r)) return(NULL)
    bind_cols(tibble(group = "Environment", level = "Station",
      series = paste(g$var[1], g$StationNme[1]), var = g$var[1],
      DBOreg = g$DBOreg[1]), r)
  })

# ---- 4. environment, REGIONAL level ---------------------------------------------------
env_reg <- env %>% group_by(var, DBOreg, DataYear) %>%
  summarise(v = mean(v), .groups = "drop") %>%
  group_by(var, DBOreg) %>% group_split() %>%
  map_dfr(function(g) {
    r <- rwpe(g$v, g$DataYear); if (is.null(r)) return(NULL)
    bind_cols(tibble(group = "Environment", level = "Regional",
      series = paste(g$var[1], g$DBOreg[1]), var = g$var[1],
      DBOreg = g$DBOreg[1]), r)
  })

res <- bind_rows(fam_stn, fam_reg, env_stn, env_reg) %>%
  mutate(sig = p < 0.05,
    group = factor(group, c("Family biomass", "Environment")),
    level = factor(level, c("Station", "Regional")))
saveRDS(res, "data/rwpe_dbo_results.rds")

cat(sprintf("\n=== rWPE, DBO (m=%d, %d permutations, min n=%d, min windows=%d) ===\n",
  M, NPERM, MIN_N, MIN_W))
cat("WPE normalized to [0,1]; 0 = perfectly predictable, 1 = maximally complex.\n")
cat("'signif' = more predictable than white noise at p<0.05 (5% expected by chance).\n\n")
cat(sprintf("%-16s %-10s %7s %8s %10s %10s %9s\n",
  "group", "level", "series", "med n", "med WPE", "med null", "signif"))
res %>% group_by(group, level) %>%
  summarise(n_series = n(), mn = median(n), mw = median(wpe), mnull = median(wpe_null),
    sg = mean(sig), .groups = "drop") %>%
  pwalk(function(group, level, n_series, mn, mw, mnull, sg)
    cat(sprintf("%-16s %-10s %7d %8.0f %10.3f %10.3f %8.0f%%\n",
      group, level, n_series, mn, mw, mnull, 100 * sg)))

cat("\n--- environment, regional, by variable ---\n")
env_reg %>% mutate(sig = p < 0.05) %>% group_by(var) %>%
  summarise(n = n(), wpe = round(median(wpe), 3), minp = round(min(p), 3),
    sig = sum(sig), .groups = "drop") %>% arrange(wpe) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nsaved -> data/rwpe_dbo_results.rds\n")
