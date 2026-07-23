# SUPPLEMENTARY: LINEAR GAM, SPATIO-TEMPORAL BLOCK CV. Same block-CV design as
# fit_gam_bbs_spatialcv.R (K spatial blocks, leave-one-block-out, pool out-of-block
# future predictions) but all climate effects are linear and the static model uses
# 8 variables (see fit_lin_bbs.R). static/dynamic/decomposed are lm(); svc is
# mgcv::bam and, for spatially held-out (unseen) strata, its random effects are
# zeroed (dummy training level + exclude) so it falls back to its fixed
# spatial+temporal linear part. K=8, n_strata >= 25, log(abundance).

library(tidyverse)
library(mgcv)
library(foreach)

LIMIT <- as.numeric(Sys.getenv("LIMIT", "Inf"))
K <- as.numeric(Sys.getenv("K", "8"))
min_strata_cv <- 25
min_pool_strata <- 15
min_test_years <- 5

vars <- c("bio2", "bio3", "bio5", "bio8", "bio9", "bio15", "bio16", "bio18")
static_comp <- paste0(vars, "_spatial")
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

keep <- c("abundance", "elevs", "party_hours", vars, comp_terms)
agg <- function(df, sp, strata_keep) filter(df, species_id == sp, strata %in% strata_keep) %>%
  group_by(strata, year) %>% summarise(across(all_of(keep), mean), .groups = "drop")

lin <- function(x) paste(x, collapse = " + ")
ctrl <- "elevs + party_hours"
f_static  <- as.formula(paste("log(abundance) ~", lin(static_comp), "+", ctrl))
f_dynamic <- as.formula(paste("log(abundance) ~", lin(vars), "+", ctrl))
f_decomp  <- as.formula(paste("log(abundance) ~", lin(comp_terms), "+", ctrl))
f_svc <- as.formula(paste("log(abundance) ~",
  lin(c(paste0(vars, "_spatial"), paste0(vars, "_temporal"))), "+", ctrl,
  "+ s(strata, bs='re') +", paste(sprintf("s(strata, %s, bs='re')", svc_terms), collapse = " + ")))

plcor <- function(strata, y, pred) tibble(strata, y, pred) %>% group_by(strata) %>%
  summarise(r = suppressWarnings(cor(y, pred)), .groups = "drop") %>% pull(r) %>% median(na.rm = TRUE)
swcor <- function(strata, y, pred) { s <- tibble(strata, y, pred) %>% group_by(strata) %>%
  summarise(o = mean(y), p = mean(pred), .groups = "drop"); suppressWarnings(cor(s$o, s$p)) }

lp <- function(f, dat, newd) as.numeric(predict(lm(f, data = dat), newd))
# SVC on unseen populations: zero the random effects (dummy training level + exclude)
svc_pred <- function(tr, te) {
  m <- bam(f_svc, data = tr, discrete = TRUE, method = "fREML")
  te2 <- te; te2$strata <- factor(levels(tr$strata)[1], levels = levels(tr$strata))
  re <- vapply(m$smooth, function(s) if (inherits(s, "random.effect")) s$label else NA_character_, character(1))
  as.numeric(predict(m, te2, exclude = re[!is.na(re)]))
}

n_run <- min(LIMIT, length(species))
cat(sprintf("LINEAR BBS spatial block CV: K=%d, n_strata>=%d, %d candidate species\n",
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
    te_static <- te; te_static[static_comp] <- te[vars]
    tryCatch(tibble(strata = te$strata, year = te$year, y = log(te$abundance),
      p_static  = lp(f_static, tr, te_static),
      p_dynamic = lp(f_dynamic, tr, te),
      p_decomp  = lp(f_decomp, tr, te),
      p_svc     = svc_pred(tr, te)), error = function(e) NULL)
  }
  if (is.null(pooled) || n_distinct(pooled$strata) < min_pool_strata) return(NULL)

  tibble(species_id = sp, n_strata = nrow(co), n_pool = n_distinct(pooled$strata),
    static_sw = swcor(pooled$strata, pooled$y, pooled$p_static), static_pl = plcor(pooled$strata, pooled$y, pooled$p_static),
    dynamic_sw = swcor(pooled$strata, pooled$y, pooled$p_dynamic), dynamic_pl = plcor(pooled$strata, pooled$y, pooled$p_dynamic),
    decomp_sw = swcor(pooled$strata, pooled$y, pooled$p_decomp), decomp_pl = plcor(pooled$strata, pooled$y, pooled$p_decomp),
    svc_sw = swcor(pooled$strata, pooled$y, pooled$p_svc), svc_pl = plcor(pooled$strata, pooled$y, pooled$p_svc))
}
res <- bind_rows(res)
saveRDS(res, sprintf("data/lin_bbs_spatialcv_k%d_results.rds", K))

q <- function(x) sprintf("%.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
  quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))
cat(sprintf("\n\n=== LINEAR spatial block CV (%d species, n_strata>=%d, K=%d) ===\n", nrow(res), min_strata_cv, K))
cat(sprintf("median pooled populations/species: %.0f\n\n", median(res$n_pool)))
cat(sprintf("%-24s %-22s %-22s\n", "Model", "Species-wide (CV)", "Pop-level (CV)"))
cat(sprintf("%-24s %-22s %-22s\n", "static (8 vars)", q(res$static_sw), q(res$static_pl)))
cat(sprintf("%-24s %-22s %-22s\n", "dynamic (8)", q(res$dynamic_sw), q(res$dynamic_pl)))
cat(sprintf("%-24s %-22s %-22s\n", "decomposed (8)", q(res$decomp_sw), q(res$decomp_pl)))
cat(sprintf("%-24s %-22s %-22s\n", "svc (8, REs->0 for new)", q(res$svc_sw), q(res$svc_pl)))
cat(sprintf("\nsaved -> data/lin_bbs_spatialcv_k%d_results.rds\n", K))
