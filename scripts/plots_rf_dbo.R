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
  sw_path = "data/rf_dbo_loso_results.rds",           # species-wide: leave-one-station-out
  pl_path = "data/rf_dbo_traintest_results.rds",      # population-level: train 2001-2010 / test 2011-2019
  # Species-wide SVC comes from the separate LOSO svc file; population-level SVC
  # (svc_pl) already lives in the train/test file, so no separate pl_svc_path.
  # Pop-level moved off leave-one-year-out: LOYO injected a regression-to-mean artifact
  # (see fit_rf_dbo_traintest.R); the train/test split removes it. Species-wide stays
  # LOSO, which has no such artifact.
  sw_svc_path = "data/rf_dbo_svc_loso_results.rds",
  pl_svc_path = NULL,
  out_stem = "dbo_rf",
  title = "DBO — RF (species-wide: LOSO; pop-level: train/test 2001-10 → 2011-19)",
  model_lv = c("Static", "Dynamic", "Decomposed", "SVC"),
  fig_width = 175)
