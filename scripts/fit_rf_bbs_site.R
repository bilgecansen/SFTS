# SITE-LEVEL counterpart to fit_rf_bbs.R. Identical data, models, split and
# metrics, but the population unit is the individual BBS route (site_id) rather
# than the state x BCR stratum -- i.e. NO averaging of sites within a stratum.
# This is the finest-grain ("robustness") version of the main RF analysis:
# strata-level stays the main result, and this shows whether the same pattern
# holds at the route level. GAMs are not run at site level (RF only).
#
# The underlying file is already at species x site x year grain and the
# spatial/temporal/residual climate decomposition is already computed per site,
# so the only change from fit_rf_bbs.R is the aggregation key (strata -> site_id)
# and the population-count filter (min_strata -> min_sites). Response is
# log(abundance); zeros are excluded (nozero file), so a site "population" is the
# run of years that route detected the species.
#
# Three model types (no SVC):
#   static     : spatial (climatological) components only, projected onto the
#                test year's actual climate (raw bio renamed to _spatial)
#   dynamic    : raw, undecomposed climate
#   decomposed : spatial + temporal + residual
#
# Controls are elevs + party_hours. Test is always 2011-2020. Two training
# regimes are run:
#   standard : train year < 2011
#   buffer   : train year < 2001, 2001-2010 dropped entirely.

library(tidyverse)
library(ranger)
library(foreach)

LIMIT <- as.numeric(Sys.getenv("LIMIT", "Inf"))   # set small for a timing probe
n_trees <- 1000
min_test_years <- 5
min_sites <- 15   # species must occupy >= this many populations (sites/routes)

vars <- c("bio2", "bio3", "bio5", "bio8", "bio9", "bio15", "bio16", "bio18")
comp_terms <- as.vector(t(outer(vars,
  c("spatial", "temporal", "residual"), paste, sep = "_")))

bbs <- readRDS("data/data_bbs_nozero.rds")

keep <- c("abundance", "elevs", "party_hours", vars, comp_terms)
# each (site_id, year) is already a single row, so the group/summarise is an
# identity -- kept in this form so the only difference from fit_rf_bbs.R's agg()
# is the grouping key (strata -> site_id).
agg <- function(df, sp) {
  filter(df, species_id == sp) %>%
    group_by(site_id, year) %>%
    summarise(across(all_of(keep), mean), .groups = "drop")
}


# Model predictor sets (response = log(abundance)) ------------------------

f_static <- as.formula(paste("log(abundance) ~",
  paste(c("elevs", "party_hours", paste0(vars, "_spatial")), collapse = " + ")))
f_dynamic <- as.formula(paste("log(abundance) ~",
  paste(c("elevs", "party_hours", vars), collapse = " + ")))
f_decomp <- as.formula(paste("log(abundance) ~",
  paste(c("elevs", "party_hours", comp_terms), collapse = " + ")))


# Metrics -- population unit is now the site --------------------------------

plcor <- function(site, y, pred) {
  tibble(site, y, pred) %>%
    group_by(site) %>%
    summarise(r = suppressWarnings(cor(y, pred)), .groups = "drop") %>%
    pull(r) %>%
    median(na.rm = TRUE)
}
swcor <- function(site, y, pred) {
  s <- tibble(site, y, pred) %>%
    group_by(site) %>%
    summarise(o = mean(y), p = mean(pred), .groups = "drop")
  suppressWarnings(cor(s$o, s$p))
}


# Fit the three models for a given training regime ------------------------
# Filters (>=5 non-zero years/site, >=5 sites/species) are applied to the
# training set of that regime, exactly as in fit_rf_bbs.R.

run_rf_bbs_site <- function(train_cutoff, out_path, tag) {
  bbs_train <- filter(bbs, year < train_cutoff)
  bbs_test <- filter(bbs, year >= 2011 & year < 2021)

  year_n <- bbs_train %>%
    group_by(species_id, site_id) %>%
    summarise(n = length(abundance[abundance > 0]), .groups = "drop") %>%
    filter(n > 4)
  bbs_train$species_site <- paste(bbs_train$species_id, bbs_train$site_id)
  bbs_train2 <- filter(bbs_train, species_site %in% paste(year_n$species_id, year_n$site_id))
  site_n <- bbs_train2 %>%
    group_by(species_id) %>%
    summarize(n = length(unique(site_id)), .groups = "drop") %>%
    filter(n > 4)
  bbs_train3 <- filter(bbs_train2, species_id %in% site_n$species_id)
  species <- unique(bbs_train3$species_id)

  n_run <- min(LIMIT, length(species))
  cat(sprintf("\nRF BBS SITE [%s]: %d species (of %d), 3 models, %d trees, train<%d test 2011-2020\n",
    tag, n_run, length(species), n_trees, train_cutoff))
  flush.console()
  t0 <- Sys.time()

  res <- foreach(i = 1:n_run) %do% {
    sp <- species[i]
    tr <- agg(bbs_train3, sp)
    te <- agg(bbs_test, sp) %>%
      filter(site_id %in% unique(tr$site_id)) %>%
      group_by(site_id) %>%
      filter(n() >= min_test_years) %>%
      ungroup()

    if (i %% 25 == 0 || i == n_run) {
      cat(sprintf("  %d/%d (%.0fs)\n", i, n_run,
        as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      flush.console()
    }
    if (n_distinct(tr$site_id) < min_sites || nrow(te) < min_test_years) return(NULL)

    # static projected onto the test year's actual climate (raw bio fed in under
    # the _spatial names), as in fit_rf_bbs.R
    te_static <- te
    te_static[paste0(vars, "_spatial")] <- te[vars]

    y <- log(te$abundance)
    rf <- function(f, dat) ranger(f, data = dat, num.trees = n_trees, seed = 1, num.threads = 0)
    pr <- function(m, newd) predict(m, data = as.data.frame(newd))$predictions

    tryCatch({
      p_static <- pr(rf(f_static, tr), te_static)
      p_dynamic <- pr(rf(f_dynamic, tr), te)
      p_decomp <- pr(rf(f_decomp, tr), te)
      tibble(
        species_id = sp, n_sites = n_distinct(tr$site_id),
        static_sw = swcor(te$site_id, y, p_static), static_pl = plcor(te$site_id, y, p_static),
        dynamic_sw = swcor(te$site_id, y, p_dynamic), dynamic_pl = plcor(te$site_id, y, p_dynamic),
        decomp_sw = swcor(te$site_id, y, p_decomp), decomp_pl = plcor(te$site_id, y, p_decomp)
      )
    }, error = function(e) NULL)
  }
  res <- bind_rows(res)
  saveRDS(res, out_path)

  q <- function(x) sprintf("%.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
    quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))
  cat(sprintf("\n=== BBS RFs SITE [%s] (%d species): test 2011-2020 ===\n", tag, nrow(res)))
  cat(sprintf("median sites/species: %.0f\n", median(res$n_sites)))
  cat(sprintf("%-16s %-22s %-22s\n", "Model", "Species-wide (test)", "Pop-level (test)"))
  cat(sprintf("%-16s %-22s %-22s\n", "RF static", q(res$static_sw), q(res$static_pl)))
  cat(sprintf("%-16s %-22s %-22s\n", "RF dynamic", q(res$dynamic_sw), q(res$dynamic_pl)))
  cat(sprintf("%-16s %-22s %-22s\n", "RF decomposed", q(res$decomp_sw), q(res$decomp_pl)))
  cat(sprintf("saved -> %s\n", out_path))
  invisible(res)
}


# Run both regimes --------------------------------------------------------

run_rf_bbs_site(2011, "data/rf_bbs_site_results.rds", "standard")
run_rf_bbs_site(2001, "data/rf_bbs_site_buffer_results.rds", "buffer: 2001-2010 dropped")
