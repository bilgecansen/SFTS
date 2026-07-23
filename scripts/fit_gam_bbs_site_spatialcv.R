# SITE-LEVEL GAM, SPATIO-TEMPORAL BLOCK CV. Same block-CV design as
# fit_gam_bbs_spatialcv.R (K spatial blocks per species, leave-one-block-out, pool
# out-of-block future predictions) but at the site level with the smooth-SVC form
# (see fit_gam_bbs_site.R). The key difference from the strata SVC: the smooth
# s(long,lat) fields are defined over continuous space, so for a HELD-OUT block the
# SVC is simply evaluated at the unseen sites' coordinates (extrapolation toward
# the block, not undefined) -- no random effects to zero. K=8, n_sites >= 25.

library(tidyverse)
library(mgcv)
library(foreach)

LIMIT <- as.numeric(Sys.getenv("LIMIT", "Inf"))
K <- as.numeric(Sys.getenv("K", "8"))
min_sites_cv <- 25
min_pool_sites <- 15
min_test_years <- 5

static_vars <- c("bio1", "bio12")
static_comp <- paste0(static_vars, "_spatial")
vars <- c("bio2", "bio3", "bio5", "bio8", "bio9", "bio15", "bio16", "bio18")
comp_terms <- as.vector(t(outer(vars, c("spatial", "temporal", "residual"), paste, sep = "_")))
resid_terms <- paste0(vars, "_residual")

bbs <- readRDS("data/data_bbs_nozero.rds")
bbs_train <- filter(bbs, year < 2011)
bbs_test  <- filter(bbs, year >= 2011 & year < 2021)

year_n <- bbs_train %>% group_by(species_id, site_id) %>%
  summarise(n = sum(abundance > 0), .groups = "drop") %>% filter(n > 4)
bbs_train$species_site <- paste(bbs_train$species_id, bbs_train$site_id)
bbs_train2 <- filter(bbs_train, species_site %in% paste(year_n$species_id, year_n$site_id))
site_n <- bbs_train2 %>% group_by(species_id) %>%
  summarize(n = length(unique(site_id)), .groups = "drop") %>% filter(n > 4)
bbs_train3 <- filter(bbs_train2, species_id %in% site_n$species_id)
species <- unique(bbs_train3$species_id)
coord <- bbs_train3 %>% group_by(species_id, site_id) %>%
  summarise(lat = mean(lat), long = mean(long), .groups = "drop")
rm(bbs, bbs_train, bbs_train2)

keep <- c("abundance", "elevs", "party_hours", "lat", "long",
  vars, static_vars, comp_terms, static_comp)
agg <- function(df, sp, sites_keep) filter(df, species_id == sp, site_id %in% sites_keep) %>%
  group_by(site_id, year) %>% summarise(across(all_of(keep), mean), .groups = "drop")

sm <- function(x) paste(sprintf("s(%s)", x), collapse = " + ")
ctrl <- "elevs + party_hours"
rhs_static  <- paste(sm(static_comp), "+", ctrl)
rhs_dynamic <- paste(sm(vars), "+", ctrl)
rhs_decomp  <- paste(sm(comp_terms), "+", ctrl)
rhs_svc_fix <- paste(sm(c(paste0(vars, "_spatial"), paste0(vars, "_temporal"))), "+", ctrl)

plcor <- function(site, y, pred) tibble(site, y, pred) %>% group_by(site) %>%
  summarise(r = suppressWarnings(cor(y, pred)), .groups = "drop") %>% pull(r) %>% median(na.rm = TRUE)
swcor <- function(site, y, pred) { s <- tibble(site, y, pred) %>% group_by(site) %>%
  summarise(o = mean(y), p = mean(pred), .groups = "drop"); suppressWarnings(cor(s$o, s$p)) }

gfit <- function(f, dat) bam(f, data = dat, discrete = TRUE, method = "fREML")
gp <- function(m, newd) as.numeric(predict(m, newd))

n_run <- min(LIMIT, length(species))
cat(sprintf("GAM SITE spatial block CV: K=%d, n_sites>=%d, %d candidate species\n",
  K, min_sites_cv, length(species))); flush.console()
t0 <- Sys.time()

res <- foreach(i = 1:n_run) %do% {
  sp <- species[i]
  co <- filter(coord, species_id == sp)
  if (i %% 10 == 0 || i == n_run) {
    cat(sprintf("  %d/%d (%.0fs)\n", i, n_run,
      as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
  }
  if (nrow(co) < min_sites_cv) return(NULL)
  xy <- cbind(x = co$long * cos(mean(co$lat) * pi / 180), y = co$lat)
  km <- tryCatch({ set.seed(1); kmeans(xy, centers = K, nstart = 10) }, error = function(e) NULL)
  if (is.null(km)) return(NULL)
  co$block <- km$cluster

  pooled <- foreach(k = 1:K, .combine = bind_rows) %do% {
    tr <- agg(bbs_train3, sp, co$site_id[co$block != k])
    te <- agg(bbs_test, sp, co$site_id[co$block == k]) %>%
      group_by(site_id) %>% filter(n() >= min_test_years) %>% ungroup()
    nsite <- n_distinct(tr$site_id)
    if (nsite < 5 || nrow(te) < min_test_years) return(NULL)
    k_sp <- min(60, nsite - 1); k_by <- min(25, nsite - 1)
    svc_int <- sprintf("s(long, lat, k=%d)", k_sp)
    svc_by  <- paste(sprintf("s(long, lat, by=%s, k=%d)", resid_terms, k_by), collapse = " + ")
    f_static <- as.formula(paste("log(abundance) ~", rhs_static))
    f_dynamic <- as.formula(paste("log(abundance) ~", rhs_dynamic))
    f_decomp  <- as.formula(paste("log(abundance) ~", rhs_decomp))
    f_svc <- as.formula(paste("log(abundance) ~", rhs_svc_fix, "+", svc_int, "+", svc_by))
    te_static <- te; te_static[static_comp] <- te[static_vars]
    tryCatch(tibble(site_id = te$site_id, year = te$year, y = log(te$abundance),
      p_static  = gp(gfit(f_static, tr), te_static),
      p_dynamic = gp(gfit(f_dynamic, tr), te),
      p_decomp  = gp(gfit(f_decomp, tr), te),
      p_svc     = gp(gfit(f_svc, tr), te)), error = function(e) NULL)
  }
  if (is.null(pooled) || n_distinct(pooled$site_id) < min_pool_sites) return(NULL)

  tibble(species_id = sp, n_sites = nrow(co), n_pool = n_distinct(pooled$site_id),
    static_sw = swcor(pooled$site_id, pooled$y, pooled$p_static), static_pl = plcor(pooled$site_id, pooled$y, pooled$p_static),
    dynamic_sw = swcor(pooled$site_id, pooled$y, pooled$p_dynamic), dynamic_pl = plcor(pooled$site_id, pooled$y, pooled$p_dynamic),
    decomp_sw = swcor(pooled$site_id, pooled$y, pooled$p_decomp), decomp_pl = plcor(pooled$site_id, pooled$y, pooled$p_decomp),
    svc_sw = swcor(pooled$site_id, pooled$y, pooled$p_svc), svc_pl = plcor(pooled$site_id, pooled$y, pooled$p_svc))
}
res <- bind_rows(res)
saveRDS(res, sprintf("data/gam_bbs_site_spatialcv_k%d_results.rds", K))

q <- function(x) sprintf("%.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
  quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))
cat(sprintf("\n\n=== GAM SITE spatial block CV (%d species, n_sites>=%d, K=%d) ===\n", nrow(res), min_sites_cv, K))
cat(sprintf("median pooled sites/species: %.0f\n\n", median(res$n_pool)))
cat(sprintf("%-16s %-22s %-22s\n", "Model", "Species-wide (CV)", "Pop-level (CV)"))
cat(sprintf("%-16s %-22s %-22s\n", "static (2)", q(res$static_sw), q(res$static_pl)))
cat(sprintf("%-16s %-22s %-22s\n", "dynamic (8)", q(res$dynamic_sw), q(res$dynamic_pl)))
cat(sprintf("%-16s %-22s %-22s\n", "decomposed (8)", q(res$decomp_sw), q(res$decomp_pl)))
cat(sprintf("%-16s %-22s %-22s\n", "svc (smooth)", q(res$svc_sw), q(res$svc_pl)))
cat(sprintf("\nsaved -> data/gam_bbs_site_spatialcv_k%d_results.rds\n", K))
