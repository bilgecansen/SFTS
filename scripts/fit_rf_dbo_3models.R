# DBO random forests, following the same logic as scripts/fit_gam_bbs.R.
#
# Three model types (no SVC -- DBO's per-family sample cannot support the random
# effect structure: at the median 14 stations, a random intercept plus one random
# slope already costs ~27 parameters against ~113 training rows):
#   static     : spatial (climatological) components only
#   dynamic    : raw, undecomposed environment
#   decomposed : spatial + temporal + residual
#
# Control variable is Depth only (DBO's elevation analog; there is no party_hours
# equivalent). Latitude, Longitude, DBOreg and StationNme are NOT used -- no
# site-identifying predictors, as in the BBS GAMs. Lags are not used.
#
# "Taxon" = macroinvertebrate family, population = station, response = log biomass.
# Evaluation is leave-one-year-out CV, as in fit_rf_dbo.R. Decomposition is global
# (full-period, from data_dbo.rds) and is never recomputed per fold, so the spatial
# component is a fixed constant in every fold; folds differ only in which years'
# temporal/residual rows train the model. Reports three prediction levels on the
# out-of-sample predictions: species-wide (across-station), regional (average stations
# within a DBO region per year, correlate the regional trajectory over years), and
# population-level (within-station temporal). Regional is the temporal analog aggregated
# to the DBO latitudinal band; population-level is the pure within-station temporal test.
#
# Input: data/data_dbo.rds (from wrangle_dbo_data.R).

library(tidyverse)
library(ranger)

set.seed(1)
NTREE <- 2000
min_years <- 5 # station needs >= this many years for a within-station correlation
min_stn <- 10  # family needs >= this many stations for a reliable species-wide r
               # (below ~8-10 stations the across-station correlation is not
               # estimable and swings spuriously negative under leave-one-out)

d <- readRDS("data/data_dbo.rds") %>% filter(biomass > 0)

base <- c(
  "Temp", "Salinity", "integchla", "sedchla", "Ammonia", "Phosphate",
  "NiTriTra", "Silicate", "phigte5", "TOC", "cn"
)
comp_terms <- as.vector(t(outer(base,
  c("spatial", "temporal", "residual"), paste, sep = "_")))
ctrl <- "Depth"

plcor <- function(stn, y, pred) {
  tibble(stn, y, pred) %>%
    group_by(stn) %>%
    filter(n() >= min_years) %>%
    summarise(r = suppressWarnings(cor(y, pred)), .groups = "drop") %>%
    pull(r) %>%
    median(na.rm = TRUE)
}
swcor <- function(stn, y, pred) {
  s <- tibble(stn, y, pred) %>%
    group_by(stn) %>%
    summarise(o = mean(y), p = mean(pred), .groups = "drop")
  suppressWarnings(cor(s$o, s$p))
}
# regional (within-region temporal): average obs & pred across stations within a DBO
# region per year -> one regional trajectory, correlate over years, median across the
# family's regions. No family-pooling needed (17-year series per family x region).
rgcor <- function(reg, year, y, pred) {
  tibble(reg, year, y, pred) %>%
    group_by(reg, year) %>%
    filter(n() >= 2) %>%                     # >= 2 stations to form a regional mean
    summarise(o = mean(y), p = mean(pred), .groups = "drop") %>%
    group_by(reg) %>%
    filter(n() >= min_years) %>%             # >= 5 region-years to correlate
    summarise(r = suppressWarnings(cor(o, p)), .groups = "drop") %>%
    pull(r) %>%
    median(na.rm = TRUE)
}

fit_pred <- function(train, newd) {
  m <- ranger(log(biomass) ~ ., data = train, num.trees = NTREE)
  as.numeric(predict(m, data = newd)$predictions)
}

families <- unique(d$family)
cat(sprintf("DBO RF: %d families, LOYO CV, %d trees\n", length(families), NTREE))
flush.console()
t0 <- Sys.time()

out <- map(seq_along(families), function(i) {
  fm <- families[i]
  df <- filter(d, family == fm)
  if (n_distinct(df$StationNme) < min_stn) return(NULL)
  years <- sort(unique(df$DataYear))

  oos <- map_dfr(years, function(y) {
    tr <- filter(df, DataYear != y)
    te <- filter(df, DataYear == y)
    if (nrow(tr) < 10 || nrow(te) < 1) return(NULL)

    # static: trained on climatological means, projected onto the held-out
    # year's actual environment (raw values renamed to _spatial), as fit_rf_dbo.R does
    tr_st <- select(tr, biomass, all_of(ctrl), all_of(paste0(base, "_spatial")))
    te_st <- te %>%
      select(all_of(ctrl), all_of(base)) %>%
      rename_with(~ paste0(., "_spatial"), all_of(base))

    tr_dy <- select(tr, biomass, all_of(ctrl), all_of(base))
    te_dy <- select(te, all_of(ctrl), all_of(base))

    tr_dc <- select(tr, biomass, all_of(ctrl), all_of(comp_terms))
    te_dc <- select(te, all_of(ctrl), all_of(comp_terms))

    tibble(
      family = fm, DBOreg = te$DBOreg, stn = te$StationNme, year = y,
      y = log(te$biomass),
      p_static = fit_pred(tr_st, te_st),
      p_dynamic = fit_pred(tr_dy, te_dy),
      p_decomp = fit_pred(tr_dc, te_dc)
    )
  })
  if (nrow(oos) < min_years) return(NULL)

  if (i %% 10 == 0 || i == length(families)) {
    cat(sprintf("  %d/%d (%.0fs)\n", i, length(families),
      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    flush.console()
  }

  summ <- tibble(
    family = fm, n_stn = n_distinct(df$StationNme), rows = nrow(df),
    static_sw = swcor(oos$stn, oos$y, oos$p_static),
    static_rg = rgcor(oos$DBOreg, oos$year, oos$y, oos$p_static),
    static_pl = plcor(oos$stn, oos$y, oos$p_static),
    dynamic_sw = swcor(oos$stn, oos$y, oos$p_dynamic),
    dynamic_rg = rgcor(oos$DBOreg, oos$year, oos$y, oos$p_dynamic),
    dynamic_pl = plcor(oos$stn, oos$y, oos$p_dynamic),
    decomp_sw = swcor(oos$stn, oos$y, oos$p_decomp),
    decomp_rg = rgcor(oos$DBOreg, oos$year, oos$y, oos$p_decomp),
    decomp_pl = plcor(oos$stn, oos$y, oos$p_decomp)
  )
  list(summary = summ, preds = oos)
})
out <- out[!map_lgl(out, is.null)]
res <- map_dfr(out, "summary")
preds <- map_dfr(out, "preds")   # full station-level LOYO predictions (for region breakdowns / re-metrics)
saveRDS(res, "data/rf_dbo_results.rds")
saveRDS(preds, "data/rf_dbo_loyo_preds.rds")

q <- function(x) sprintf("%.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
  quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))

cat(sprintf("\n\n=== DBO RF, LOYO CV (%d of %d families) ===\n", nrow(res), length(families)))
cat(sprintf("%-16s %-20s %-20s %-20s\n", "Model", "Species-wide", "Regional", "Pop-level"))
cat(sprintf("%-16s %-20s %-20s %-20s\n", "RF static", q(res$static_sw), q(res$static_rg), q(res$static_pl)))
cat(sprintf("%-16s %-20s %-20s %-20s\n", "RF dynamic", q(res$dynamic_sw), q(res$dynamic_rg), q(res$dynamic_pl)))
cat(sprintf("%-16s %-20s %-20s %-20s\n", "RF decomposed", q(res$decomp_sw), q(res$decomp_rg), q(res$decomp_pl)))

# regional temporal r broken down by DBO region (per family x region, then median)
cat("\n--- Regional temporal r by DBO region (median across families) ---\n")
rg_break <- preds %>%
  pivot_longer(c(p_static, p_dynamic, p_decomp), names_to = "model", values_to = "pred") %>%
  group_by(model, DBOreg, family, year) %>%
  filter(n() >= 2) %>%
  summarise(o = mean(y), p = mean(pred), .groups = "drop") %>%
  group_by(model, DBOreg, family) %>%
  filter(n() >= min_years) %>%
  summarise(r = suppressWarnings(cor(o, p)), .groups = "drop") %>%
  group_by(model, DBOreg) %>%
  summarise(med_r = round(median(r, na.rm = TRUE), 3), n_fam = n(), .groups = "drop") %>%
  mutate(model = recode(model, p_static = "static", p_dynamic = "dynamic", p_decomp = "decomp"),
         model = factor(model, c("static", "dynamic", "decomp"))) %>%
  arrange(DBOreg, model)
print(as.data.frame(rg_break), row.names = FALSE)
cat("\nsaved -> data/rf_dbo_results.rds, data/rf_dbo_loyo_preds.rds\n")
