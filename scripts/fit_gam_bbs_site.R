# SITE-LEVEL GAM (population unit = individual BBS route, not stratum). Four model
# types, same as fit_gam_bbs.R, EXCEPT the SVC uses smooth spatially-varying
# coefficients instead of per-site random slopes -- at the site level the RE form
# (one intercept + 8 slopes per site, thousands of coefficients) is intractable,
# so we replace
#     s(strata, bs="re")                 -> s(long, lat, k=k_sp)
#     s(strata, X_residual, bs="re")     -> s(long, lat, by=X_residual, k=k_by)
# a fixed-basis spatial field whose cost is independent of the number of sites and
# which (unlike the RE version) can be evaluated at unseen sites. k adapts per
# species: k_sp = min(60, n_sites-1), k_by = min(25, n_sites-1).
#
# Models (response log(abundance), controls elevs + party_hours):
#   static     : 2 vars (bio1 + bio12), spatial components, projected onto the
#                test year's actual climate
#   dynamic    : 8 raw vars
#   decomposed : 8 vars x {spatial, temporal, residual}
#   svc        : 8 spatial + 8 temporal smooths + smooth spatial intercept +
#                smooth spatially-varying residual coefficients
#
# NOTE: s(long,lat) and the s(bio_spatial) terms are both functions of location,
# so they are concurved -- fine for prediction, but the individual SVC smooths are
# not interpretable as pure climate effects.
#
# Test 2011-2020. Temporal (year < 2011) and Buffered (year < 2001, gap dropped).
# The Spatiotemporal block-CV regime is in fit_gam_bbs_site_spatialcv.R.

library(tidyverse)
library(mgcv)
library(foreach)

LIMIT <- as.numeric(Sys.getenv("LIMIT", "Inf"))
min_sites <- 15
min_test_years <- 5

static_vars <- c("bio1", "bio12")
static_comp <- paste0(static_vars, "_spatial")
vars <- c("bio2", "bio3", "bio5", "bio8", "bio9", "bio15", "bio16", "bio18")
comp_terms <- as.vector(t(outer(vars, c("spatial", "temporal", "residual"), paste, sep = "_")))
resid_terms <- paste0(vars, "_residual")

bbs <- readRDS("data/data_bbs_nozero.rds")

keep <- c("abundance", "elevs", "party_hours", "lat", "long",
  vars, static_vars, comp_terms, static_comp)
agg <- function(df, sp) {
  filter(df, species_id == sp) %>%
    group_by(site_id, year) %>%
    summarise(across(all_of(keep), mean), .groups = "drop")
}

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


# Fit the four models for a training regime -------------------------------

run_gam_bbs_site <- function(train_cutoff, out_path, tag) {
  bbs_train <- filter(bbs, year < train_cutoff)
  bbs_test <- filter(bbs, year >= 2011 & year < 2021)

  year_n <- bbs_train %>% group_by(species_id, site_id) %>%
    summarise(n = sum(abundance > 0), .groups = "drop") %>% filter(n > 4)
  bbs_train$species_site <- paste(bbs_train$species_id, bbs_train$site_id)
  bbs_train2 <- filter(bbs_train, species_site %in% paste(year_n$species_id, year_n$site_id))
  site_n <- bbs_train2 %>% group_by(species_id) %>%
    summarize(n = length(unique(site_id)), .groups = "drop") %>% filter(n > 4)
  bbs_train3 <- filter(bbs_train2, species_id %in% site_n$species_id)
  species <- unique(bbs_train3$species_id)

  n_run <- min(LIMIT, length(species))
  cat(sprintf("\nGAM SITE [%s]: %d species (of %d), 4 models (smooth SVC), train<%d test 2011-2020\n",
    tag, n_run, length(species), train_cutoff)); flush.console()
  t0 <- Sys.time()

  res <- foreach(i = 1:n_run) %do% {
    sp <- species[i]
    tr <- agg(bbs_train3, sp)
    te <- agg(bbs_test, sp) %>%
      filter(site_id %in% unique(tr$site_id)) %>%
      group_by(site_id) %>% filter(n() >= min_test_years) %>% ungroup()

    if (i %% 25 == 0 || i == n_run) {
      cat(sprintf("  %d/%d (%.0fs)\n", i, n_run,
        as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
    }
    nsite <- n_distinct(tr$site_id)
    if (nsite < min_sites || nrow(te) < min_test_years) return(NULL)

    k_sp <- min(60, nsite - 1); k_by <- min(25, nsite - 1)
    svc_int <- sprintf("s(long, lat, k=%d)", k_sp)
    svc_by  <- paste(sprintf("s(long, lat, by=%s, k=%d)", resid_terms, k_by), collapse = " + ")
    f_static <- as.formula(paste("log(abundance) ~", rhs_static))
    f_dynamic <- as.formula(paste("log(abundance) ~", rhs_dynamic))
    f_decomp  <- as.formula(paste("log(abundance) ~", rhs_decomp))
    f_svc <- as.formula(paste("log(abundance) ~", rhs_svc_fix, "+", svc_int, "+", svc_by))

    te_static <- te
    te_static[static_comp] <- te[static_vars]
    y <- log(te$abundance)
    fit <- function(f, newd) as.numeric(predict(bam(f, data = tr, discrete = TRUE, method = "fREML"), newd))

    tryCatch({
      ps <- fit(f_static, te_static); pd <- fit(f_dynamic, te)
      pc <- fit(f_decomp, te); pv <- fit(f_svc, te)
      tibble(species_id = sp, n_sites = nsite, k_sp = k_sp,
        static_sw = swcor(te$site_id, y, ps), static_pl = plcor(te$site_id, y, ps),
        dynamic_sw = swcor(te$site_id, y, pd), dynamic_pl = plcor(te$site_id, y, pd),
        decomp_sw = swcor(te$site_id, y, pc), decomp_pl = plcor(te$site_id, y, pc),
        svc_sw = swcor(te$site_id, y, pv), svc_pl = plcor(te$site_id, y, pv))
    }, error = function(e) NULL)
  }
  res <- bind_rows(res)
  saveRDS(res, out_path)

  q <- function(x) sprintf("%.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
    quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))
  cat(sprintf("\n=== BBS GAM SITE [%s] (%d species) ===\n", tag, nrow(res)))
  cat(sprintf("%-16s %-22s %-22s\n", "Model", "Species-wide", "Pop-level"))
  cat(sprintf("%-16s %-22s %-22s\n", "static (2)", q(res$static_sw), q(res$static_pl)))
  cat(sprintf("%-16s %-22s %-22s\n", "dynamic (8)", q(res$dynamic_sw), q(res$dynamic_pl)))
  cat(sprintf("%-16s %-22s %-22s\n", "decomposed (8)", q(res$decomp_sw), q(res$decomp_pl)))
  cat(sprintf("%-16s %-22s %-22s\n", "svc (smooth)", q(res$svc_sw), q(res$svc_pl)))
  cat(sprintf("saved -> %s\n", out_path))
  invisible(res)
}

run_gam_bbs_site(2011, "data/gam_bbs_site_results.rds", "Temporal")
run_gam_bbs_site(2001, "data/gam_bbs_site_buffer_results.rds", "Buffered")
