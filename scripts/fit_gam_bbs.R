# GAM counterpart to fit_rf_bbs.R -- same data, filters, split and metrics, but
# mgcv GAMs (smooth terms). Four model types:
#   static     : 2 vars only (bio1 AMT + bio12 precip), spatial (climatological)
#                components, projected onto the test year's actual climate. The
#                2-variable static outperformed the 8-variable version, so the
#                nonlinear static GAM uses 2 vars throughout.
#   dynamic    : raw, undecomposed climate (8 vars)
#   decomposed : spatial + temporal + residual (8 vars)
#   svc        : decomposed with a per-population random intercept and random
#                residual slopes (spatially-varying coefficients)
#
# Controls: elevs + party_hours. No lat/long. Site identity (strata) is only a
# predictor in svc (its random effects); elsewhere it just filters species to
# those occupying >= min_strata populations. Response is log(abundance).
#
# Test is always 2011-2020. Two training regimes ("treatments") are run:
#   Temporal : train year < 2011  (future years, same populations)
#   Buffered : train year < 2001, 2001-2010 dropped -- a 10-year gap forcing
#              genuine temporal extrapolation.
# (The third treatment, Spatiotemporal block CV, is in fit_gam_bbs_spatialcv.R.)

library(tidyverse)
library(mgcv)
library(foreach)

LIMIT <- as.numeric(Sys.getenv("LIMIT", "Inf"))   # set small for a timing probe
min_strata <- 15
min_test_years <- 5

static_vars <- c("bio1", "bio12")
static_comp <- paste0(static_vars, "_spatial")
vars <- c("bio2", "bio3", "bio5", "bio8", "bio9", "bio15", "bio16", "bio18")
comp_terms <- as.vector(t(outer(vars, c("spatial", "temporal", "residual"), paste, sep = "_")))
svc_terms <- paste0(vars, "_residual")

bbs <- readRDS("data/data_bbs_nozero.rds")

keep <- c("abundance", "elevs", "party_hours", vars, static_vars, comp_terms, static_comp)
agg <- function(df, sp) {
  filter(df, species_id == sp) %>%
    group_by(strata, year) %>%
    summarise(across(all_of(keep), mean), .groups = "drop")
}

sm <- function(x) paste(sprintf("s(%s)", x), collapse = " + ")
ctrl <- "elevs + party_hours"
f_static  <- as.formula(paste("log(abundance) ~", sm(static_comp), "+", ctrl))
f_dynamic <- as.formula(paste("log(abundance) ~", sm(vars), "+", ctrl))
f_decomp  <- as.formula(paste("log(abundance) ~", sm(comp_terms), "+", ctrl))
f_svc <- as.formula(paste("log(abundance) ~",
  sm(c(paste0(vars, "_spatial"), paste0(vars, "_temporal"))), "+", ctrl,
  "+ s(strata, bs='re') +", paste(sprintf("s(strata, %s, bs='re')", svc_terms), collapse = " + ")))

plcor <- function(strata, y, pred) {
  tibble(strata, y, pred) %>% group_by(strata) %>%
    summarise(r = suppressWarnings(cor(y, pred)), .groups = "drop") %>%
    pull(r) %>% median(na.rm = TRUE)
}
swcor <- function(strata, y, pred) {
  s <- tibble(strata, y, pred) %>% group_by(strata) %>%
    summarise(o = mean(y), p = mean(pred), .groups = "drop")
  suppressWarnings(cor(s$o, s$p))
}


# Fit the four models for a given training regime -------------------------

run_gam_bbs <- function(train_cutoff, out_path, tag) {
  bbs_train <- filter(bbs, year < train_cutoff)
  bbs_test <- filter(bbs, year >= 2011 & year < 2021)

  year_n <- bbs_train %>%
    group_by(species_id, site_id) %>%
    summarise(n = length(abundance[abundance > 0]), .groups = "drop") %>% filter(n > 4)
  bbs_train$species_site <- paste(bbs_train$species_id, bbs_train$site_id)
  bbs_train2 <- filter(bbs_train, species_site %in% paste(year_n$species_id, year_n$site_id))
  site_n <- bbs_train2 %>%
    group_by(species_id) %>% summarize(n = length(unique(site_id)), .groups = "drop") %>% filter(n > 4)
  bbs_train3 <- filter(bbs_train2, species_id %in% site_n$species_id)
  species <- unique(bbs_train3$species_id)

  n_run <- min(LIMIT, length(species))
  cat(sprintf("\nGAM BBS [%s]: %d species (of %d), 4 models, train<%d test 2011-2020\n",
    tag, n_run, length(species), train_cutoff))
  flush.console()
  t0 <- Sys.time()

  res <- foreach(i = 1:n_run) %do% {
    sp <- species[i]
    tr <- agg(bbs_train3, sp) %>% mutate(strata = factor(strata))
    te <- agg(bbs_test, sp) %>%
      filter(strata %in% levels(tr$strata)) %>%
      group_by(strata) %>% filter(n() >= min_test_years) %>% ungroup() %>%
      mutate(strata = factor(strata, levels = levels(tr$strata)))

    if (i %% 25 == 0 || i == n_run) {
      cat(sprintf("  %d/%d (%.0fs)\n", i, n_run,
        as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
    }
    if (n_distinct(tr$strata) < min_strata || nrow(te) < min_test_years) return(NULL)

    te_static <- te
    te_static[static_comp] <- te[static_vars]   # project raw climate -> _spatial
    y <- log(te$abundance)
    fit <- function(f, newd) as.numeric(predict(
      bam(f, data = tr, discrete = TRUE, method = "fREML"), newd))

    tryCatch({
      ps <- fit(f_static, te_static); pd <- fit(f_dynamic, te)
      pc <- fit(f_decomp, te); pv <- fit(f_svc, te)
      tibble(species_id = sp, n_strata = n_distinct(tr$strata),
        static_sw = swcor(te$strata, y, ps), static_pl = plcor(te$strata, y, ps),
        dynamic_sw = swcor(te$strata, y, pd), dynamic_pl = plcor(te$strata, y, pd),
        decomp_sw = swcor(te$strata, y, pc), decomp_pl = plcor(te$strata, y, pc),
        svc_sw = swcor(te$strata, y, pv), svc_pl = plcor(te$strata, y, pv))
    }, error = function(e) NULL)
  }
  res <- bind_rows(res)
  saveRDS(res, out_path)

  q <- function(x) sprintf("%.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
    quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))
  cat(sprintf("\n=== BBS GAMs [%s] (%d species): test 2011-2020 ===\n", tag, nrow(res)))
  cat(sprintf("%-22s %-22s %-22s\n", "Model", "Species-wide", "Pop-level"))
  cat(sprintf("%-22s %-22s %-22s\n", "static (2 vars)", q(res$static_sw), q(res$static_pl)))
  cat(sprintf("%-22s %-22s %-22s\n", "dynamic (8)", q(res$dynamic_sw), q(res$dynamic_pl)))
  cat(sprintf("%-22s %-22s %-22s\n", "decomposed (8)", q(res$decomp_sw), q(res$decomp_pl)))
  cat(sprintf("%-22s %-22s %-22s\n", "svc (8)", q(res$svc_sw), q(res$svc_pl)))
  cat(sprintf("saved -> %s\n", out_path))
  invisible(res)
}


# Run both treatments -----------------------------------------------------

run_gam_bbs(2011, "data/gam_bbs_results.rds", "Temporal")
run_gam_bbs(2001, "data/gam_bbs_buffer_results.rds", "Buffered")
