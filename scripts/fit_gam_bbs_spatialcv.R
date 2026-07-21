# GAM spatio-temporal block CV (main 4 model types). Same design as the RF
# version: K spatial blocks per species (k-means on lat/long), leave-one-block-out
# (train on other blocks pre-2011, predict held-out block's 2011-2020), pool the
# out-of-block future predictions, compute species-wide and population-level r.
#
# Models:
#   static     : 2 vars only (bio1 AMT + bio12 precip) -- did better than 8 vars
#   dynamic    : 8 raw vars
#   decomposed : 8 vars x {spatial, temporal, residual}
#   svc        : 8 vars spatial+temporal + per-population random intercept &
#                random residual slopes. For SPATIALLY HELD-OUT (unseen)
#                populations the random effects are undefined, so we predict with
#                them set to zero (exclude=), i.e. SVC falls back to its fixed
#                spatial+temporal part -- its per-population machinery cannot
#                transfer to new locations.
#
# K=8, n_strata >= 25. train < 2011 / test 2011-2020. log(abundance).

library(tidyverse)
library(mgcv)
library(foreach)

LIMIT <- as.numeric(Sys.getenv("LIMIT", "Inf"))
K <- as.numeric(Sys.getenv("K", "8"))
min_strata_cv <- 25
min_pool_strata <- 15
min_test_years <- 5

static_vars <- c("bio1", "bio12")
static_comp <- paste0(static_vars, "_spatial")
vars <- c("bio2", "bio3", "bio5", "bio8", "bio9", "bio15", "bio16", "bio18")
comp_terms <- as.vector(t(outer(vars, c("spatial", "temporal", "residual"), paste, sep = "_")))
svc_terms <- paste0(vars, "_residual")

bbs <- readRDS("data/data_bbs_nozero.rds")
bbs_train <- filter(bbs, year < 2011)
bbs_test  <- filter(bbs, year >= 2011 & year < 2021)

year_n <- bbs_train %>% group_by(species_id, site_id) %>%
  summarise(n = length(abundance[abundance > 0]), .groups = "drop") %>% filter(n > 4)
bbs_train$species_site <- paste(bbs_train$species_id, bbs_train$site_id)
bbs_train2 <- filter(bbs_train, species_site %in% paste(year_n$species_id, year_n$site_id))
site_n <- bbs_train2 %>% group_by(species_id) %>%
  summarize(n = length(unique(site_id)), .groups = "drop") %>% filter(n > 4)
bbs_train3 <- filter(bbs_train2, species_id %in% site_n$species_id)
species <- unique(bbs_train3$species_id)
coord <- bbs_train3 %>% group_by(species_id, strata) %>%
  summarise(lat = mean(lat), long = mean(long), .groups = "drop")
rm(bbs, bbs_train, bbs_train2)

keep <- c("abundance", "elevs", "party_hours", vars, static_vars, comp_terms, static_comp)
agg <- function(df, sp, strata_keep) filter(df, species_id == sp, strata %in% strata_keep) %>%
  group_by(strata, year) %>% summarise(across(all_of(keep), mean), .groups = "drop")

sm <- function(x) paste(sprintf("s(%s)", x), collapse = " + ")
f_static  <- as.formula(paste("log(abundance) ~", sm(static_comp), "+ elevs + party_hours"))
f_dynamic <- as.formula(paste("log(abundance) ~", sm(vars), "+ elevs + party_hours"))
f_decomp  <- as.formula(paste("log(abundance) ~", sm(comp_terms), "+ elevs + party_hours"))
f_svc <- as.formula(paste("log(abundance) ~",
  sm(c(paste0(vars, "_spatial"), paste0(vars, "_temporal"))), "+ elevs + party_hours",
  "+ s(strata, bs='re') +", paste(sprintf("s(strata, %s, bs='re')", svc_terms), collapse = " + ")))

plcor <- function(strata, y, pred) tibble(strata, y, pred) %>% group_by(strata) %>%
  summarise(r = suppressWarnings(cor(y, pred)), .groups = "drop") %>% pull(r) %>% median(na.rm = TRUE)
swcor <- function(strata, y, pred) { s <- tibble(strata, y, pred) %>% group_by(strata) %>%
  summarise(o = mean(y), p = mean(pred), .groups = "drop"); suppressWarnings(cor(s$o, s$p)) }

gfit <- function(f, dat) bam(f, data = dat, discrete = TRUE, method = "fREML")
gpred <- function(m, newd) as.numeric(predict(m, newd))
# SVC on unseen populations: zero the random effects (dummy training level + exclude)
gpred_svc <- function(m, te, tr) {
  te2 <- te; te2$strata <- factor(levels(tr$strata)[1], levels = levels(tr$strata))
  re <- vapply(m$smooth, function(s) if (inherits(s, "random.effect")) s$label else NA_character_, character(1))
  as.numeric(predict(m, te2, exclude = re[!is.na(re)]))
}

n_run <- min(LIMIT, length(species))
cat(sprintf("GAM BBS spatial block CV: K=%d, n_strata>=%d, %d candidate species\n",
  K, min_strata_cv, length(species))); flush.console()
t0 <- Sys.time()

res <- foreach(i = 1:n_run) %do% {
  sp <- species[i]
  co <- filter(coord, species_id == sp)
  if (i %% 10 == 0 || i == n_run) {
    cat(sprintf("  %d/%d (%.0fs)\n", i, n_run,
      as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
  }
  if (nrow(co) < min_strata_cv) return(NULL)
  xy <- cbind(x = co$long * cos(mean(co$lat) * pi / 180), y = co$lat)
  km <- tryCatch({ set.seed(1); kmeans(xy, centers = K, nstart = 10) }, error = function(e) NULL)
  if (is.null(km)) return(NULL)
  co$block <- km$cluster

  pooled <- foreach(k = 1:K, .combine = bind_rows) %do% {
    tr <- agg(bbs_train3, sp, co$strata[co$block != k]) %>% mutate(strata = factor(strata))
    te <- agg(bbs_test, sp, co$strata[co$block == k]) %>%
      group_by(strata) %>% filter(n() >= min_test_years) %>% ungroup()
    if (n_distinct(tr$strata) < 5 || nrow(te) < min_test_years) return(NULL)
    te_static <- te; te_static[static_comp] <- te[static_vars]
    tryCatch(tibble(strata = te$strata, year = te$year, y = log(te$abundance),
      p_static  = gpred(gfit(f_static, tr), te_static),
      p_dynamic = gpred(gfit(f_dynamic, tr), te),
      p_decomp  = gpred(gfit(f_decomp, tr), te),
      p_svc     = gpred_svc(gfit(f_svc, tr), te, tr)), error = function(e) NULL)
  }
  if (is.null(pooled) || n_distinct(pooled$strata) < min_pool_strata) return(NULL)

  tibble(species_id = sp, n_strata = nrow(co), n_pool = n_distinct(pooled$strata),
    static_sw = swcor(pooled$strata, pooled$y, pooled$p_static), static_pl = plcor(pooled$strata, pooled$y, pooled$p_static),
    dynamic_sw = swcor(pooled$strata, pooled$y, pooled$p_dynamic), dynamic_pl = plcor(pooled$strata, pooled$y, pooled$p_dynamic),
    decomp_sw = swcor(pooled$strata, pooled$y, pooled$p_decomp), decomp_pl = plcor(pooled$strata, pooled$y, pooled$p_decomp),
    svc_sw = swcor(pooled$strata, pooled$y, pooled$p_svc), svc_pl = plcor(pooled$strata, pooled$y, pooled$p_svc))
}
res <- bind_rows(res)
saveRDS(res, sprintf("data/gam_bbs_spatialcv_k%d_results.rds", K))

q <- function(x) sprintf("%.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
  quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))
cat(sprintf("\n\n=== GAM spatial block CV (%d species, n_strata>=%d, K=%d) ===\n", nrow(res), min_strata_cv, K))
cat(sprintf("median pooled populations/species: %.0f\n\n", median(res$n_pool)))
cat(sprintf("%-24s %-22s %-22s\n", "Model", "Species-wide (CV)", "Pop-level (CV)"))
cat(sprintf("%-24s %-22s %-22s\n", "static (2 vars)", q(res$static_sw), q(res$static_pl)))
cat(sprintf("%-24s %-22s %-22s\n", "dynamic (8)", q(res$dynamic_sw), q(res$dynamic_pl)))
cat(sprintf("%-24s %-22s %-22s\n", "decomposed (8)", q(res$decomp_sw), q(res$decomp_pl)))
cat(sprintf("%-24s %-22s %-22s\n", "svc (8, REs->0 for new)", q(res$svc_sw), q(res$svc_pl)))
cat(sprintf("\nsaved -> data/gam_bbs_spatialcv_k%d_results.rds\n", K))
