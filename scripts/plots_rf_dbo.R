# DBO RF main-text figure (median-band style), three models
# (static / dynamic / decomposed). One point per family. Each metric is drawn
# from the holdout that honestly tests it:
#   species-wide     <- leave-one-STATION-out (spatial extrapolation to an unseen
#                        station); avoids the leave-one-year-out inflation where
#                        the fixed 23-station gradient is memorised.
#   population-level  <- leave-one-YEAR-out (temporal extrapolation to an unseen
#                        year); leave-one-station-out would leak the year effect.
# Both files are filtered to families with >= 10 stations, so they cover the same
# 45 taxa. The population-level panel auto-ranges, showing its negative values.

source("scripts/plot_helpers.R")

make_medband_split(
  sw_path = "data/rf_dbo_loso_results.rds",   # species-wide: leave-one-station-out
  pl_path = "data/rf_dbo_results.rds",        # population-level: leave-one-year-out
  # SVC candidate (decomposed + Latitude/Longitude), same two holdouts; merged in as
  # the 4th model so DBO matches the BBS figures (from fit_rf_dbo_svc.R).
  sw_svc_path = "data/rf_dbo_svc_loso_results.rds",
  pl_svc_path = "data/rf_dbo_svc_results.rds",
  out_stem = "dbo_rf",
  title = "DBO — RF (species-wide: leave-one-station-out; pop-level: leave-one-year-out)",
  model_lv = c("Static", "Dynamic", "Decomposed", "SVC"),
  fig_width = 175)
