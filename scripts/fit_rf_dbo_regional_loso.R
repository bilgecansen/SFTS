# DBO RF, LEAVE-ONE-STATION-OUT (LOSO), all three prediction levels in one place.
# Canonical source for figures/dbo_rf_loso_3row.png. Previously the regional row was
# computed in a scratchpad script that no longer exists; this rebuilds it and, as a check,
# also recomputes the spatial and station-level rows so they can be verified against the
# existing data/rf_dbo_loso_results.rds and data/rf_dbo_svc_loso_results.rds.
#
# Holdout: an entire station (all its years) is held out; the model trains on the OTHER
# stations (all years) and predicts the held-out station. Spatial prediction therefore
# becomes honest extrapolation to an unseen location. Decomposition is the GLOBAL
# full-period one baked into data_dbo.rds -- correct here because the environment is
# exogenous and the held-out station's own climatology is a legitimate input.
#
# Models (response = log biomass, always spatiotemporal; only the ENV representation differs):
#   Static     : TRAIN on per-station means (*_spatial). PREDICT with spatiotemporal env fed
#                into those slots (regional/station level); species-wide uses the held-out
#                station's MEAN env -> one prediction per station.
#   Dynamic    : TRAIN and TEST both raw spatiotemporal env.
#   Decomposed : TRAIN and TEST both decomposed (spatial + temporal + residual).
#   SVC        : Decomposed + Latitude + Longitude.
# Depth is the control in all four.
#
# Prediction levels:
#   spatial  : station-mean observed vs predicted, correlated across the family's stations
#   regional : average observed & predicted across stations within a DBO region per year,
#              correlate the regional trajectory over years -> one r per family x region
#   station  : within-station correlation over years, median across stations
#
# REGIONAL AGGREGATION IS COMPUTED TWO WAYS, for comparison:
#   _log : mean(log biomass)      -- geometric mean; each station weighted equally in
#          proportional terms. This is what every other script in the project does.
#   _raw : log(mean(biomass))     -- arithmetic mean; the highest-biomass station dominates
#          the regional trajectory. Predictions are back-transformed with exp() first, which
#          is a median (not mean) retransformation, so _raw carries a known downward bias in
#          level; the bias is roughly common to obs and pred and largely drops out of a
#          correlation, but it is not exactly neutral.
#
# Output: data/rf_dbo_regional_loso_results.rds  (family x region, both aggregations)
#         data/rf_dbo_loso_preds.rds             (out-of-fold station-year predictions)
#         data/rf_dbo_loso_levels.rds            (family-level spatial + station-level check)

library(tidyverse)
library(ranger)

set.seed(1)
NTREE <- 2000
min_years <- 5     # station needs >= this many years for a within-station correlation
min_stn <- 10      # family needs >= this many stations (matches fit_rf_dbo_loso.R)
min_reg_stn <- 3   # family needs >= this many stations in a region for a regional cell
min_reg_years <- 5 # and >= this many region-years

d <- readRDS("data/data_dbo.rds") %>% filter(biomass > 0)
base <- c("Temp", "Salinity", "integchla", "sedchla", "Ammonia", "Phosphate",
          "NiTriTra", "Silicate", "phigte5", "TOC", "cn")
comp_terms <- as.vector(t(outer(base, c("spatial", "temporal", "residual"), paste, sep = "_")))
ctrl <- "Depth"

fit_pred <- function(train, newd) {
  m <- ranger(log(biomass) ~ ., data = train, num.trees = NTREE)
  as.numeric(predict(m, data = newd)$predictions)
}

families <- unique(d$family)
cat(sprintf("DBO RF LOSO: %d families, %d trees, 4 models\n", length(families), NTREE))
flush.console(); t0 <- Sys.time()

out <- map(seq_along(families), function(i) {
  fm <- families[i]
  df <- filter(d, family == fm)
  if (n_distinct(df$StationNme) < min_stn) return(NULL)
  if (i %% 10 == 0) {
    cat(sprintf("  %d/%d (%.0fs)\n", i, length(families),
      as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
  }
  per <- map(sort(unique(df$StationNme)), function(s) {
    tr <- filter(df, StationNme != s); te <- filter(df, StationNme == s)
    if (nrow(tr) < 10 || nrow(te) < 1) return(NULL)

    # static: trained on climatological means of the OTHER stations
    m_st <- ranger(log(biomass) ~ .,
      data = select(tr, biomass, all_of(ctrl), all_of(paste0(base, "_spatial"))),
      num.trees = NTREE)
    te_st <- te %>% select(all_of(ctrl), all_of(base)) %>%
      rename_with(~ paste0(., "_spatial"), all_of(base))
    p_static_yr <- as.numeric(predict(m_st, te_st)$predictions)              # yearly -> regional/station
    p_static_sw <- as.numeric(predict(m_st,
      summarise(te_st, across(everything(), mean)))$predictions)             # mean env -> spatial

    list(
      yr = tibble(family = fm, DBOreg = te$DBOreg, stn = te$StationNme, year = te$DataYear,
        y = log(te$biomass), biomass = te$biomass,
        p_static  = p_static_yr,
        p_dynamic = fit_pred(select(tr, biomass, all_of(ctrl), all_of(base)),
                             select(te, all_of(ctrl), all_of(base))),
        p_decomp  = fit_pred(select(tr, biomass, all_of(ctrl), all_of(comp_terms)),
                             select(te, all_of(ctrl), all_of(comp_terms))),
        p_svc     = fit_pred(select(tr, biomass, all_of(ctrl), all_of(comp_terms), Latitude, Longitude),
                             select(te, all_of(ctrl), all_of(comp_terms), Latitude, Longitude))),
      sw = tibble(family = fm, stn = s, o = mean(log(te$biomass)), p_static_sw = p_static_sw))
  })
  per <- per[!map_lgl(per, is.null)]
  if (length(per) < 2) return(NULL)
  list(yr = bind_rows(map(per, "yr")), sw = bind_rows(map(per, "sw")))
})
out <- out[!map_lgl(out, is.null)]
preds <- map_dfr(out, "yr"); sw_static <- map_dfr(out, "sw")
saveRDS(preds, "data/rf_dbo_loso_preds.rds")

mods <- c(static = "p_static", dynamic = "p_dynamic", decomp = "p_decomp", svc = "p_svc")

# ---- regional: one r per family x region, computed BOTH ways -------------------------
reg <- map_dfr(names(mods), function(m) {
  pc <- mods[[m]]
  preds %>%
    group_by(family, DBOreg, year) %>%
    summarise(nst = n_distinct(stn),
      o_log = mean(y),                     # mean(log biomass)  -> geometric
      p_log = mean(.data[[pc]]),
      o_raw = log(mean(biomass)),          # log(mean biomass)  -> arithmetic
      p_raw = log(mean(exp(.data[[pc]]))),
      .groups = "drop") %>%
    group_by(family, DBOreg) %>%
    filter(max(nst) >= min_reg_stn, n() >= min_reg_years) %>%
    summarise(n_stn = max(nst), n_yr = n(),
      r_log = suppressWarnings(cor(o_log, p_log)),
      r_raw = suppressWarnings(cor(o_raw, p_raw)), .groups = "drop") %>%
    mutate(model = m)
})
reg_w <- reg %>%
  pivot_wider(names_from = model, values_from = c(r_log, r_raw), names_glue = "{model}_{.value}")
saveRDS(reg_w, "data/rf_dbo_regional_loso_results.rds")

# ---- spatial + station level (check against the existing LOSO files) ------------------
plcor <- function(dd, col) dd %>% group_by(stn) %>% filter(n() >= min_years) %>%
  summarise(r = suppressWarnings(cor(y, .data[[col]])), .groups = "drop") %>%
  pull(r) %>% median(na.rm = TRUE)
swcor <- function(dd, col) { q <- dd %>% group_by(stn) %>%
    summarise(o = mean(y), p = mean(.data[[col]]), .groups = "drop")
  suppressWarnings(cor(q$o, q$p)) }

lev <- preds %>% group_split(family) %>% map_dfr(function(dd) {
  fm <- dd$family[1]
  s <- filter(sw_static, family == fm)
  o <- tibble(family = fm,
    static_sw = suppressWarnings(cor(s$o, s$p_static_sw)),   # static spatial: mean-env prediction
    static_pl = plcor(dd, "p_static"))
  for (m in c("dynamic", "decomp", "svc")) {
    o[[paste0(m, "_sw")]] <- swcor(dd, mods[[m]])
    o[[paste0(m, "_pl")]] <- plcor(dd, mods[[m]])
  }
  o
})
saveRDS(lev, "data/rf_dbo_loso_levels.rds")

# ---- report ---------------------------------------------------------------------------
md <- function(x) sprintf("%+.2f", median(x, na.rm = TRUE))
cat(sprintf("\n\n=== DBO RF LOSO (%d families) ===\n", nrow(lev)))
cat("\n-- spatial and station level (compare to existing rf_dbo_loso_results.rds) --\n")
cat(sprintf("%-10s %-10s %-10s\n", "model", "spatial", "station"))
for (m in names(mods))
  cat(sprintf("%-10s %-10s %-10s\n", m, md(lev[[paste0(m, "_sw")]]), md(lev[[paste0(m, "_pl")]])))

cat("\n-- REGIONAL, mean(log biomass) vs log(mean biomass) --\n")
cat(sprintf("%-10s %-8s %-6s %-12s %-12s %-8s\n",
  "model", "region", "n_fam", "mean(log)", "log(mean)", "diff"))
reg %>% group_by(model, DBOreg) %>%
  summarise(n = n(), a = median(r_log, na.rm = TRUE), b = median(r_raw, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(model = factor(model, names(mods))) %>% arrange(model, DBOreg) %>%
  pwalk(function(model, DBOreg, n, a, b)
    cat(sprintf("%-10s DBO %-4s %-6d %-12s %-12s %+.3f\n", model, DBOreg, n,
      sprintf("%+.3f", a), sprintf("%+.3f", b), b - a)))
cat("\nsaved -> data/rf_dbo_regional_loso_results.rds, data/rf_dbo_loso_preds.rds, data/rf_dbo_loso_levels.rds\n")
