# DBO random forests, LEAVE-ONE-STATION-OUT (LOSO) cross-validation.
#
# Companion to fit_rf_dbo_3models.R (which uses leave-one-YEAR-out, LOYO). LOYO
# trains on the same 23 fixed stations in every other year, so the spatial
# gradient is memorised before the held-out year is scored -- the species-wide
# correlation is spatial INTERPOLATION among known locations, not prediction,
# which is why LOYO species-wide sits near ceiling. LOSO instead holds out an
# entire station (all its years), trains on the OTHER stations (all years), and
# predicts the held-out station. Species-wide then becomes honest spatial
# EXTRAPOLATION (place an unseen location on the gradient from environment
# alone); population-level also becomes stricter because the model has never seen
# the held-out station's own history -- it must transfer year-to-year signal from
# other places. The DBO analog of BBS's "unseen routes" holdout, spatial-only (a
# temporal train-early/test-late split is not affordable with 17 years).
#
# Three model types (no SVC), Depth as the only control; Latitude, Longitude,
# DBOreg, StationNme are NOT predictors. "Taxon" = family, population = station,
# response = log biomass. Input: data/data_dbo.rds (from wrangle_dbo_data.R).

library(tidyverse)
library(ranger)

set.seed(1)
NTREE <- 2000
min_years <- 5 # station needs >= this many years for a within-station correlation
min_stn <- 10  # family needs >= this many stations: LOSO holds one out, and an
               # across-station correlation on < ~8-10 stations is not estimable
               # (it swings spuriously negative, e.g. the 5-station families)

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

fit_pred <- function(train, newd) {
  m <- ranger(log(biomass) ~ ., data = train, num.trees = NTREE)
  as.numeric(predict(m, data = newd)$predictions)
}

families <- unique(d$family)
cat(sprintf("DBO RF: %d families, LOSO CV (leave-one-station-out), %d trees\n",
  length(families), NTREE))
flush.console()
t0 <- Sys.time()

res <- map_dfr(seq_along(families), function(i) {
  fm <- families[i]
  df <- filter(d, family == fm)
  if (n_distinct(df$StationNme) < min_stn) return(NULL)
  stations <- sort(unique(df$StationNme))

  # leave-one-station-out: hold out an entire station (all its years)
  oos <- map_dfr(stations, function(s) {
    tr <- filter(df, StationNme != s)
    te <- filter(df, StationNme == s)
    if (nrow(tr) < 10 || nrow(te) < 1) return(NULL)

    # static: trained on climatological (_spatial) means of the OTHER stations,
    # projected onto the held-out station's actual per-year environment (raw
    # values renamed to _spatial), as in fit_rf_dbo_3models.R
    tr_st <- select(tr, biomass, all_of(ctrl), all_of(paste0(base, "_spatial")))
    te_st <- te %>%
      select(all_of(ctrl), all_of(base)) %>%
      rename_with(~ paste0(., "_spatial"), all_of(base))

    tr_dy <- select(tr, biomass, all_of(ctrl), all_of(base))
    te_dy <- select(te, all_of(ctrl), all_of(base))

    tr_dc <- select(tr, biomass, all_of(ctrl), all_of(comp_terms))
    te_dc <- select(te, all_of(ctrl), all_of(comp_terms))

    tibble(
      stn = te$StationNme, year = te$DataYear, y = log(te$biomass),
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

  tibble(
    family = fm, n_stn = n_distinct(df$StationNme), rows = nrow(df),
    static_sw = swcor(oos$stn, oos$y, oos$p_static),
    static_pl = plcor(oos$stn, oos$y, oos$p_static),
    dynamic_sw = swcor(oos$stn, oos$y, oos$p_dynamic),
    dynamic_pl = plcor(oos$stn, oos$y, oos$p_dynamic),
    decomp_sw = swcor(oos$stn, oos$y, oos$p_decomp),
    decomp_pl = plcor(oos$stn, oos$y, oos$p_decomp)
  )
})
saveRDS(res, "data/rf_dbo_loso_results.rds")

q <- function(x) sprintf("%.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
  quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))

# LOYO reference (same families) so the species-wide collapse is visible --------
loyo <- tryCatch(readRDS("data/rf_dbo_results.rds"), error = function(e) NULL)

cat(sprintf("\n\n=== DBO RF, LOSO CV (%d of %d families) ===\n", nrow(res), length(families)))
cat(sprintf("%-14s %-22s %-22s | %-16s\n", "Model", "Species-wide (LOSO)", "Pop-level (LOSO)", "SW LOYO (ref)"))
for (m in c("static", "dynamic", "decomp")) {
  lv <- if (is.null(loyo)) "-" else sprintf("%.2f", median(loyo[[paste0(m, "_sw")]], na.rm = TRUE))
  cat(sprintf("%-14s %-22s %-22s | %-16s\n", m,
    q(res[[paste0(m, "_sw")]]), q(res[[paste0(m, "_pl")]]), lv))
}
cat("\nsaved -> data/rf_dbo_loso_results.rds\n")
