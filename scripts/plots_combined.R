# SUPPLEMENTARY BBS grid figures (2x3 median-band grid: species-wide over
# population-level rows x Temporal/Buffered/Spatiotemporal columns, one figure per
# model family). See make_medband_grid in plot_helpers.R. Writes to figures/.
#
# The MAIN paper figures are now the "two-axis" figures in plots_2axis.R (GAM+RF
# combined). This script keeps the per-family grids for the SI comparisons only:
# un-clamped GAM (vs the clamped main) and the 2-variable (bio1+bio12) robustness.

source("scripts/plot_helpers.R")
treatments <- c("Temporal", "Buffered", "Spatiotemporal")
rf_lv  <- c("Static", "Dynamic", "Decomposed")
gam_lv <- c("Static", "Dynamic", "Decomposed", "SVC")

# stem -> (result triple prefix, model levels, title)
specs <- list(
  # --- SI: un-clamped GAM (vs the clamped main analysis) ---
  list("gam_bbs_noclamp_combined",     "gam_bbs_noclamp",     gam_lv, "GAM (strata) — un-clamped"),
  list("gam_bbs_site_noclamp_combined","gam_bbs_site_noclamp",gam_lv, "GAM (route) — un-clamped"),
  # --- SI: 2 variables (bio1 + bio12) ---
  list("rf_bbs_2var_combined",         "rf_bbs_2var",         rf_lv,  "RF (strata) — 2 vars"),
  list("rf_bbs_site_2var_combined",    "rf_bbs_site_2var",    rf_lv,  "RF (route) — 2 vars"),
  list("gam_bbs_2var_combined",        "gam_bbs_2var",        gam_lv, "GAM (strata) — 2 vars"),
  list("gam_bbs_site_2var_combined",   "gam_bbs_site_2var",   gam_lv, "GAM (route) — 2 vars")
)

for (s in specs) {
  stem <- s[[1]]; pre <- s[[2]]; lv <- s[[3]]; title <- s[[4]]
  paths <- c(sprintf("data/%s_results.rds", pre),
             sprintf("data/%s_buffer_results.rds", pre),
             sprintf("data/%s_spatialcv_k8_results.rds", pre))
  if (!all(file.exists(paths))) { cat(sprintf("skip %s (missing results)\n", stem)); next }
  make_medband_grid(paths, treatments, stem, lv, title = title)
}
