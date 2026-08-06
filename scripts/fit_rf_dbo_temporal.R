# DBO RF, leave-one-year-out, TEMPORAL prediction levels (regional + population-level).
# The spatial (species-wide) prediction level is not computed here.
#
# Models -- the standard four. Response = log biomass, always spatiotemporal; what
# differs is how the ENVIRONMENTAL predictors are represented:
#   Static     : TRAIN on per-station means (*_spatial). PREDICT with the held-out year's
#                RAW env fed into those same slots -> a spatial relationship used to
#                predict across time (the space-for-time test).
#   Dynamic    : TRAIN and TEST both raw spatiotemporal env.
#   Decomposed : TRAIN and TEST both decomposed, SPATIAL BLOCK DROPPED -> temporal + residual.
#   SVC        : Decomposed (as above) + Latitude + Longitude.
# Depth is NOT included here (LOYO only; the LOSO pipeline still uses it as the control).
# Depth is constant within a station, so under LOYO it is one more station-identifying
# feature; dropping it removes that route to a per-station intercept entirely.
#
# Why Decomposed/SVC drop the spatial block here. Under LOYO the eleven *_spatial features
# are constant within a station, so they act as a station identifier: the forest fits a
# per-station intercept, and under LOYO that intercept is the leave-one-out MEAN of the
# station's response, which is mechanically anti-correlated with the held-out year. The
# pure station-LOO-mean predictor scores population-level r = -1.00, and the full
# decomposed model collapsed onto it (within-station SD(pred)/SD(obs) = 0.22), giving
# r = -0.33. Dropping the spatial block removes the intercept. Depth on its own was not
# the culprit (keeping vs dropping it moved nothing: -0.037 vs -0.042), but it is dropped
# here anyway so that no constant-per-station feature remains except SVC's coordinates.
# NOTE: SVC's Latitude/Longitude are also constant per station and may reintroduce the
# same intercept -- that is a result to read off this run, not an assumption.
#
# Decomposition is the GLOBAL full-period one from data_dbo.rds, never recomputed per
# fold, so the spatial component is a fixed constant in every fold.
#
# Prediction levels (both temporal):
#   regional         : average obs & pred across stations within a DBO region per year,
#                      correlate the regional trajectory over years, median across the
#                      family's regions
#   population-level : within-station correlation over years, median across stations
#
# Input : data/data_dbo.rds
# Output: data/rf_dbo_temporal_results.rds (<model>_rg / <model>_pl per family)
#         data/rf_dbo_temporal_preds.rds   (station-level LOYO predictions)

library(tidyverse)
library(ranger)

set.seed(1)
NTREE <- 2000
min_years <- 5 # station (or region) needs >= this many years for a correlation
min_stn <- 10  # family needs >= this many stations, matching the other DBO scripts

d <- readRDS("data/data_dbo.rds") %>% filter(biomass > 0)

base <- c(
  "Temp", "Salinity", "integchla", "sedchla", "Ammonia", "Phosphate",
  "NiTriTra", "Silicate", "phigte5", "TOC", "cn"
)
td_terms <- c(paste0(base, "_temporal"), paste0(base, "_residual"))  # spatial block dropped
coords <- c("Latitude", "Longitude")

fit_pred <- function(train, newd) {
  m <- ranger(log(biomass) ~ ., data = train, num.trees = NTREE)
  as.numeric(predict(m, data = newd)$predictions)
}

families <- unique(d$family)
cat(sprintf("DBO RF temporal: %d families, LOYO, %d trees, 4 models\n",
  length(families), NTREE))
flush.console()
t0 <- Sys.time()

preds <- map_dfr(seq_along(families), function(i) {
  fm <- families[i]
  df <- filter(d, family == fm)
  if (n_distinct(df$StationNme) < min_stn) return(NULL)
  if (i %% 10 == 0) {
    cat(sprintf("  %d/%d (%.0fs)\n", i, length(families),
      as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
  }
  map_dfr(sort(unique(df$DataYear)), function(y) {
    tr <- filter(df, DataYear != y); te <- filter(df, DataYear == y)
    if (nrow(tr) < 10 || nrow(te) < 1) return(NULL)

    # static: fit on climatological means, predict with the held-out year's raw env
    tr_st <- select(tr, biomass, all_of(paste0(base, "_spatial")))
    te_st <- te %>% select(all_of(base)) %>%
      rename_with(~ paste0(., "_spatial"), all_of(base))

    tibble(
      family = fm, DBOreg = te$DBOreg, stn = te$StationNme, year = y,
      y = log(te$biomass),
      p_static  = fit_pred(tr_st, te_st),
      p_dynamic = fit_pred(select(tr, biomass, all_of(base)),
                           select(te, all_of(base))),
      p_decomp  = fit_pred(select(tr, biomass, all_of(td_terms)),
                           select(te, all_of(td_terms))),
      p_svc     = fit_pred(select(tr, biomass, all_of(td_terms), all_of(coords)),
                           select(te, all_of(td_terms), all_of(coords)))
    )
  })
})
saveRDS(preds, "data/rf_dbo_temporal_preds.rds")

plcor <- function(dd, col) dd %>%
  group_by(stn) %>% filter(n() >= min_years) %>%
  summarise(r = suppressWarnings(cor(y, .data[[col]])), .groups = "drop") %>%
  pull(r) %>% median(na.rm = TRUE)

rgcor <- function(dd, col) dd %>%
  group_by(DBOreg, year) %>% filter(n() >= 2) %>%          # >= 2 stations per regional mean
  summarise(o = mean(y), p = mean(.data[[col]]), .groups = "drop") %>%
  group_by(DBOreg) %>% filter(n() >= min_years) %>%        # >= 5 region-years to correlate
  summarise(r = suppressWarnings(cor(o, p)), .groups = "drop") %>%
  pull(r) %>% median(na.rm = TRUE)

mods <- c(static = "p_static", dynamic = "p_dynamic", decomp = "p_decomp", svc = "p_svc")
res <- preds %>% group_split(family) %>% map_dfr(function(dd) {
  out <- tibble(family = dd$family[1], n_stn = n_distinct(dd$stn))
  for (m in names(mods)) {
    out[[paste0(m, "_rg")]] <- rgcor(dd, mods[[m]])
    out[[paste0(m, "_pl")]] <- plcor(dd, mods[[m]])
  }
  out
})
saveRDS(res, "data/rf_dbo_temporal_results.rds")

# collapse diagnostic: how much of the observed within-station variation does each
# model's prediction actually reproduce? (low = collapsed onto a per-station intercept)
wsd <- function(col) preds %>% group_by(family, stn) %>% filter(n() >= min_years) %>%
  summarise(x = sd(.data[[col]]) / sd(y), .groups = "drop") %>% pull(x) %>% median(na.rm = TRUE)

q <- function(x) sprintf("%+.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
  quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))
cat(sprintf("\n\n=== DBO RF, LOYO, temporal levels (%d families) ===\n", nrow(res)))
cat(sprintf("%-14s %-22s %-22s %s\n", "Model", "Regional", "Population-level", "SD(pred)/SD(obs)"))
for (m in names(mods))
  cat(sprintf("%-14s %-22s %-22s %.2f\n", m, q(res[[paste0(m, "_rg")]]),
    q(res[[paste0(m, "_pl")]]), wsd(mods[[m]])))
cat("\nsaved -> data/rf_dbo_temporal_results.rds, data/rf_dbo_temporal_preds.rds\n")
