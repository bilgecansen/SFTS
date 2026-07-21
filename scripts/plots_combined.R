# Combined BBS main-text figures -- the DEFAULT BBS figure structure. Each figure
# is a 2x3 grid: species-wide over population-level (rows), Temporal / Buffered /
# Spatiotemporal (columns), with the y-scale shared across treatments within each
# row and the population-level row 1/3 the height of species-wide. One figure per
# model family (see make_medband_grid in plot_helpers.R). Journal full width.
#   gam_bbs_combined      - GAM, strata (4 models incl. SVC)      [main]
#   rf_bbs_combined       - RF,  strata (3 models)                [main]
#   rf_bbs_site_combined  - RF,  site-level / route (3 models)    [SI robustness]

source("scripts/plot_helpers.R")

treatments <- c("Temporal", "Buffered", "Spatiotemporal")

# GAM (strata) -- four models including SVC
make_medband_grid(
  paths = c("data/gam_bbs_results.rds", "data/gam_bbs_buffer_results.rds",
            "data/gam_bbs_spatialcv_k8_results.rds"),
  treatments = treatments, out_stem = "gam_bbs_combined",
  model_lv = c("Static", "Dynamic", "Decomposed", "SVC"))

# RF (strata) -- three models (no SVC)
make_medband_grid(
  paths = c("data/rf_bbs_results.rds", "data/rf_bbs_buffer_results.rds",
            "data/rf_bbs_spatialcv_k8_results.rds"),
  treatments = treatments, out_stem = "rf_bbs_combined",
  model_lv = c("Static", "Dynamic", "Decomposed"))

# RF (site-level / route) -- three models, SI robustness companion
make_medband_grid(
  paths = c("data/rf_bbs_site_results.rds", "data/rf_bbs_site_buffer_results.rds",
            "data/rf_bbs_site_spatialcv_k8_results.rds"),
  treatments = treatments, out_stem = "rf_bbs_site_combined",
  model_lv = c("Static", "Dynamic", "Decomposed"))

# GAM + spatiotemporal smoother (strata, 4 models) -- SI. Only built once its
# three result files exist (see fit_gam_bbs_st.R / fit_gam_bbs_st_spatialcv.R).
st_paths <- c("data/gam_bbs_st_results.rds", "data/gam_bbs_st_buffer_results.rds",
              "data/gam_bbs_st_spatialcv_k8_results.rds")
if (all(file.exists(st_paths)))
  make_medband_grid(paths = st_paths, treatments = treatments,
    out_stem = "gam_bbs_st_combined", model_lv = c("Static", "Dynamic", "Decomposed", "SVC"))

# Linear GAM (strata, 4 models, 8-var static) -- SI. Built once its three result
# files exist (see fit_lin_bbs.R / fit_lin_bbs_spatialcv.R).
lin_paths <- c("data/lin_bbs_results.rds", "data/lin_bbs_buffer_results.rds",
               "data/lin_bbs_spatialcv_k8_results.rds")
if (all(file.exists(lin_paths)))
  make_medband_grid(paths = lin_paths, treatments = treatments,
    out_stem = "lin_bbs_combined", model_lv = c("Static", "Dynamic", "Decomposed", "SVC"))
