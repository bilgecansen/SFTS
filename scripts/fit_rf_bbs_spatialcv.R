# Spatio-temporal block cross-validation (RF). On TOP of the temporal split
# (train < 2011 / test 2011-2020) we add a SPATIAL holdout: each species'
# populations are k-means clustered by (lat, long) into K=5 contiguous blocks.
# For each block: train on the OTHER blocks' pre-2011 data (the held-out block's
# pre-2011 data is discarded -- those populations never inform the model), and
# predict the held-out block's 2011-2020. Every population is thus predicted once,
# out-of-block and in the future. Predictions are POOLED across folds and a single
# species-wide and population-level correlation is computed per species.
#
# Three models: static / dynamic / decomposed (no SVC). Only species with
# n_strata >= 25 (enough populations to block, train, and give a reliable pooled
# species-wide). Same 8 bioclim vars, elevs + party_hours, log(abundance).

library(tidyverse)
library(ranger)
library(foreach)

LIMIT <- as.numeric(Sys.getenv("LIMIT", "Inf"))
K <- as.numeric(Sys.getenv("K", "8"))
min_strata_cv <- 25     # populations required to enter this analysis
min_pool_strata <- 15   # pooled held-out populations required for a valid species-wide r
min_test_years <- 5
n_trees <- 1000

vars <- c("bio2", "bio3", "bio5", "bio8", "bio9", "bio15", "bio16", "bio18")
comp_terms <- as.vector(t(outer(vars, c("spatial", "temporal", "residual"), paste, sep = "_")))

bbs <- readRDS("data/data_bbs_nozero.rds")
bbs_train <- filter(bbs, year < 2011)
bbs_test  <- filter(bbs, year >= 2011 & year < 2021)

year_n <- bbs_train %>%
  group_by(species_id, site_id) %>%
  summarise(n = length(abundance[abundance > 0]), .groups = "drop") %>% filter(n > 4)
bbs_train$species_site <- paste(bbs_train$species_id, bbs_train$site_id)
bbs_train2 <- filter(bbs_train, species_site %in% paste(year_n$species_id, year_n$site_id))
site_n <- bbs_train2 %>%
  group_by(species_id) %>% summarize(n = length(unique(site_id)), .groups = "drop") %>% filter(n > 4)
bbs_train3 <- filter(bbs_train2, species_id %in% site_n$species_id)
species <- unique(bbs_train3$species_id)

# per-stratum centroid for spatial blocking
coord <- bbs_train3 %>% group_by(species_id, strata) %>%
  summarise(lat = mean(lat), long = mean(long), .groups = "drop")
rm(bbs, bbs_train, bbs_train2)

keep <- c("abundance", "elevs", "party_hours", vars, comp_terms)
agg <- function(df, sp, strata_keep) {
  filter(df, species_id == sp, strata %in% strata_keep) %>%
    group_by(strata, year) %>% summarise(across(all_of(keep), mean), .groups = "drop")
}

f_static <- as.formula(paste("log(abundance) ~",
  paste(c("elevs", "party_hours", paste0(vars, "_spatial")), collapse = " + ")))
f_dynamic <- as.formula(paste("log(abundance) ~",
  paste(c("elevs", "party_hours", vars), collapse = " + ")))
f_decomp <- as.formula(paste("log(abundance) ~",
  paste(c("elevs", "party_hours", comp_terms), collapse = " + ")))

plcor <- function(strata, y, pred) tibble(strata, y, pred) %>% group_by(strata) %>%
  summarise(r = suppressWarnings(cor(y, pred)), .groups = "drop") %>% pull(r) %>% median(na.rm = TRUE)
swcor <- function(strata, y, pred) { s <- tibble(strata, y, pred) %>% group_by(strata) %>%
  summarise(o = mean(y), p = mean(pred), .groups = "drop"); suppressWarnings(cor(s$o, s$p)) }

rf <- function(f, dat) ranger(f, data = dat, num.trees = n_trees, seed = 1, num.threads = 0)
pr <- function(m, newd) predict(m, data = as.data.frame(newd))$predictions

n_run <- min(LIMIT, length(species))
cat(sprintf("RF BBS spatial block CV: K=%d, n_strata>=%d, %d candidate species\n",
  K, min_strata_cv, length(species))); flush.console()
t0 <- Sys.time()

res <- foreach(i = 1:n_run) %do% {
  sp <- species[i]
  co <- filter(coord, species_id == sp)
  if (i %% 25 == 0 || i == n_run) {
    cat(sprintf("  %d/%d (%.0fs)\n", i, n_run,
      as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
  }
  if (nrow(co) < min_strata_cv) return(NULL)

  # K contiguous spatial blocks (longitude scaled to ~equal distance)
  xy <- cbind(x = co$long * cos(mean(co$lat) * pi / 180), y = co$lat)
  km <- tryCatch({ set.seed(1); kmeans(xy, centers = K, nstart = 10) }, error = function(e) NULL)
  if (is.null(km)) return(NULL)
  co$block <- km$cluster

  # leave-one-block-out: pool out-of-block future predictions
  pooled <- foreach(k = 1:K, .combine = bind_rows) %do% {
    tr <- agg(bbs_train3, sp, co$strata[co$block != k])
    te <- agg(bbs_test, sp, co$strata[co$block == k]) %>%
      group_by(strata) %>% filter(n() >= min_test_years) %>% ungroup()
    if (n_distinct(tr$strata) < 5 || nrow(te) < min_test_years) return(NULL)
    te_static <- te; te_static[paste0(vars, "_spatial")] <- te[vars]
    tryCatch(tibble(strata = te$strata, year = te$year, y = log(te$abundance),
      p_static = pr(rf(f_static, tr), te_static),
      p_dynamic = pr(rf(f_dynamic, tr), te),
      p_decomp = pr(rf(f_decomp, tr), te)), error = function(e) NULL)
  }
  if (is.null(pooled) || n_distinct(pooled$strata) < min_pool_strata) return(NULL)

  tibble(species_id = sp, n_strata = nrow(co), n_pool = n_distinct(pooled$strata),
    static_sw = swcor(pooled$strata, pooled$y, pooled$p_static),
    static_pl = plcor(pooled$strata, pooled$y, pooled$p_static),
    dynamic_sw = swcor(pooled$strata, pooled$y, pooled$p_dynamic),
    dynamic_pl = plcor(pooled$strata, pooled$y, pooled$p_dynamic),
    decomp_sw = swcor(pooled$strata, pooled$y, pooled$p_decomp),
    decomp_pl = plcor(pooled$strata, pooled$y, pooled$p_decomp))
}
res <- bind_rows(res)
out_path <- sprintf("data/rf_bbs_spatialcv_k%d_results.rds", K)
saveRDS(res, out_path)

# Summary vs the standard (non-spatial) RF on the SAME n_strata>=25 species -----
q <- function(x) sprintf("%.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
  quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))
std <- tryCatch(readRDS("data/rf_bbs_results.rds") %>% filter(n_strata >= min_strata_cv),
  error = function(e) NULL)

cat(sprintf("\n\n=== RF spatial block CV (%d species, n_strata>=%d, K=%d): pooled out-of-block ===\n",
  nrow(res), min_strata_cv, K))
cat(sprintf("median pooled populations/species: %.0f\n\n", median(res$n_pool)))
cat(sprintf("%-14s %-22s %-22s | %-16s\n", "Model", "Species-wide (CV)", "Pop-level (CV)", "SW std (same sp.)"))
for (m in c("static", "dynamic", "decomp")) {
  sv <- if (is.null(std)) "-" else sprintf("%.2f", median(std[[paste0(m, "_sw")]], na.rm = TRUE))
  cat(sprintf("%-14s %-22s %-22s | %-16s\n", m,
    q(res[[paste0(m, "_sw")]]), q(res[[paste0(m, "_pl")]]), sv))
}
cat(sprintf("\nsaved -> %s\n", out_path))
