# SUPPLEMENTARY: GAM with an added SPATIOTEMPORAL SMOOTHER. Identical to
# fit_gam_bbs.R (same data, filters, 2-var static, split, metrics, four model
# types) except every model gains a space-time tensor term
#     te(long, lat, year, d = c(2, 1), k = c(k_sp, k_yr))
# a 2-D spatial marginal (strata centroid) crossed with a 1-D year marginal, with
# per-species adaptive bases k_sp <- min(10, max(3, n_strata-1)) and
# k_yr <- min(5, max(3, n_years-1)). This soaks up smooth spatial + temporal +
# space-time structure beyond the climate terms, testing whether accounting for
# spatiotemporal autocorrelation changes the species-wide / population-level
# result. Note the year marginal EXTRAPOLATES into the 2011-2020 test window (and,
# under spatial block CV, the spatial marginal extrapolates to unseen strata) --
# that is intrinsic to a forecasting test and is the point: a flexible space-time
# smoother still cannot forecast unseen dynamics.
#
# Test is 2011-2020. Two training regimes: Temporal (year < 2011) and Buffered
# (year < 2001, 2001-2010 dropped). The Spatiotemporal block-CV regime is in
# fit_gam_bbs_st_spatialcv.R.

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

keep <- c("abundance", "elevs", "party_hours", "lat", "long",
  vars, static_vars, comp_terms, static_comp)
agg <- function(df, sp) {
  filter(df, species_id == sp) %>%
    group_by(strata, year) %>%
    summarise(across(all_of(keep), mean), .groups = "drop")
}

sm <- function(x) paste(sprintf("s(%s)", x), collapse = " + ")
ctrl <- "elevs + party_hours"
rhs_static  <- paste(sm(static_comp), "+", ctrl)
rhs_dynamic <- paste(sm(vars), "+", ctrl)
rhs_decomp  <- paste(sm(comp_terms), "+", ctrl)
rhs_svc <- paste(sm(c(paste0(vars, "_spatial"), paste0(vars, "_temporal"))), "+", ctrl,
  "+ s(strata, bs='re') +", paste(sprintf("s(strata, %s, bs='re')", svc_terms), collapse = " + "))

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


# Fit the four models (+ spatiotemporal smoother) for a training regime --------

run_gam_bbs_st <- function(train_cutoff, out_path, tag) {
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
  cat(sprintf("\nGAM-ST BBS [%s]: %d species (of %d), 4 models + te(long,lat,year), train<%d test 2011-2020\n",
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

    k_sp <- min(10, max(3, n_distinct(tr$strata) - 1))
    k_yr <- min(5, max(3, n_distinct(tr$year) - 1))
    st <- sprintf("te(long, lat, year, d = c(2,1), k = c(%d, %d))", k_sp, k_yr)
    f_static  <- as.formula(paste("log(abundance) ~", rhs_static,  "+", st))
    f_dynamic <- as.formula(paste("log(abundance) ~", rhs_dynamic, "+", st))
    f_decomp  <- as.formula(paste("log(abundance) ~", rhs_decomp,  "+", st))
    f_svc     <- as.formula(paste("log(abundance) ~", rhs_svc,     "+", st))

    te_static <- te
    te_static[static_comp] <- te[static_vars]   # project raw climate -> _spatial
    y <- log(te$abundance)
    fit <- function(f, newd) as.numeric(predict(
      bam(f, data = tr, discrete = TRUE, method = "fREML"), newd))

    tryCatch({
      ps <- fit(f_static, te_static); pd <- fit(f_dynamic, te)
      pc <- fit(f_decomp, te); pv <- fit(f_svc, te)
      tibble(species_id = sp, n_strata = n_distinct(tr$strata), k_sp = k_sp, k_yr = k_yr,
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
  cat(sprintf("\n=== BBS GAM-ST [%s] (%d species): test 2011-2020 ===\n", tag, nrow(res)))
  cat(sprintf("%-22s %-22s %-22s\n", "Model", "Species-wide", "Pop-level"))
  cat(sprintf("%-22s %-22s %-22s\n", "static (2 vars)", q(res$static_sw), q(res$static_pl)))
  cat(sprintf("%-22s %-22s %-22s\n", "dynamic (8)", q(res$dynamic_sw), q(res$dynamic_pl)))
  cat(sprintf("%-22s %-22s %-22s\n", "decomposed (8)", q(res$decomp_sw), q(res$decomp_pl)))
  cat(sprintf("%-22s %-22s %-22s\n", "svc (8)", q(res$svc_sw), q(res$svc_pl)))
  cat(sprintf("saved -> %s\n", out_path))
  invisible(res)
}

run_gam_bbs_st(2011, "data/gam_bbs_st_results.rds", "Temporal")
run_gam_bbs_st(2001, "data/gam_bbs_st_buffer_results.rds", "Buffered")
