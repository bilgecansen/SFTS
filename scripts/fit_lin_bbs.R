# SUPPLEMENTARY: LINEAR GAM (all climate effects enter linearly, no smooths).
# Parallels fit_gam_bbs.R exactly -- same data, filters, split, metrics, four
# model types -- but each s(x) smooth is replaced by a linear term x. Because
# linear terms cannot overfit, the static model uses the full 8 variables here
# (the nonlinear static was trimmed to 2 to avoid wiggly overfitting; that
# concern does not apply to a linear fit).
#
# Models (response log(abundance), controls elevs + party_hours):
#   static     : 8 spatial (climatological) terms, projected onto the test year's
#                actual climate (raw renamed to _spatial)
#   dynamic    : 8 raw climate terms
#   decomposed : 8 vars x {spatial, temporal, residual} = 24 linear terms
#   svc        : linear spatial + temporal terms + per-population random intercept
#                and random residual slopes (identical RE structure to the
#                nonlinear SVC -- those were already linear random effects)
# static / dynamic / decomposed are pure linear models (lm); svc needs the random
# effects so it is fit with mgcv::bam.
#
# Test 2011-2020. Temporal (year < 2011) and Buffered (year < 2001, gap dropped).
# The Spatiotemporal block-CV regime is in fit_lin_bbs_spatialcv.R.

library(tidyverse)
library(mgcv)
library(foreach)

LIMIT <- as.numeric(Sys.getenv("LIMIT", "Inf"))
min_strata <- 15
min_test_years <- 5

vars <- c("bio2", "bio3", "bio5", "bio8", "bio9", "bio15", "bio16", "bio18")
static_comp <- paste0(vars, "_spatial")   # 8-var static (linear)
comp_terms <- as.vector(t(outer(vars, c("spatial", "temporal", "residual"), paste, sep = "_")))
svc_terms <- paste0(vars, "_residual")

bbs <- readRDS("data/data_bbs_nozero.rds")

keep <- c("abundance", "elevs", "party_hours", vars, comp_terms)
agg <- function(df, sp) {
  filter(df, species_id == sp) %>%
    group_by(strata, year) %>%
    summarise(across(all_of(keep), mean), .groups = "drop")
}

lin <- function(x) paste(x, collapse = " + ")
ctrl <- "elevs + party_hours"
f_static  <- as.formula(paste("log(abundance) ~", lin(static_comp), "+", ctrl))
f_dynamic <- as.formula(paste("log(abundance) ~", lin(vars), "+", ctrl))
f_decomp  <- as.formula(paste("log(abundance) ~", lin(comp_terms), "+", ctrl))
f_svc <- as.formula(paste("log(abundance) ~",
  lin(c(paste0(vars, "_spatial"), paste0(vars, "_temporal"))), "+", ctrl,
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

# Fit the four models for a training regime -------------------------------

run_lin_bbs <- function(train_cutoff, out_path, tag) {
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
  cat(sprintf("\nLINEAR BBS [%s]: %d species (of %d), 4 models (8-var static), train<%d test 2011-2020\n",
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
    te_static[static_comp] <- te[vars]   # project raw climate -> _spatial
    y <- log(te$abundance)
    lp <- function(f, newd) as.numeric(predict(lm(f, data = tr), newd))

    tryCatch({
      ps <- lp(f_static, te_static); pd <- lp(f_dynamic, te); pc <- lp(f_decomp, te)
      pv <- as.numeric(predict(bam(f_svc, data = tr, discrete = TRUE, method = "fREML"), te))
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
  cat(sprintf("\n=== BBS LINEAR [%s] (%d species): test 2011-2020 ===\n", tag, nrow(res)))
  cat(sprintf("%-22s %-22s %-22s\n", "Model", "Species-wide", "Pop-level"))
  cat(sprintf("%-22s %-22s %-22s\n", "static (8 vars)", q(res$static_sw), q(res$static_pl)))
  cat(sprintf("%-22s %-22s %-22s\n", "dynamic (8)", q(res$dynamic_sw), q(res$dynamic_pl)))
  cat(sprintf("%-22s %-22s %-22s\n", "decomposed (8)", q(res$decomp_sw), q(res$decomp_pl)))
  cat(sprintf("%-22s %-22s %-22s\n", "svc (8)", q(res$svc_sw), q(res$svc_pl)))
  cat(sprintf("saved -> %s\n", out_path))
  invisible(res)
}

run_lin_bbs(2011, "data/lin_bbs_results.rds", "Temporal")
run_lin_bbs(2001, "data/lin_bbs_buffer_results.rds", "Buffered")
