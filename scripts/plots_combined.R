# Combined BBS figures (CANONICAL). Each figure is the 2x3 median-band grid:
# species-wide over population-level (rows) x Temporal / Buffered / Spatiotemporal
# (columns), from the per-block-decomposition results. See make_medband_grid in
# plot_helpers.R. Writes to figures/.
#
# MAIN:  RF strata, RF site, GAM strata (clamped), GAM site (clamped).
# SI:    GAM strata/site un-clamped ; RF & GAM at 2 vars (bio1+bio12), strata/site.
# (GAM main analysis clamps predictors to the training range at predict time; the
#  un-clamped versions are the SI extrapolation comparison.)

source("scripts/plot_helpers.R")
treatments <- c("Temporal", "Buffered", "Spatiotemporal")
rf_lv  <- c("Static", "Dynamic", "Decomposed")
gam_lv <- c("Static", "Dynamic", "Decomposed", "SVC")

# stem -> (result triple prefix, model levels, title)
specs <- list(
  # --- main ---
  list("rf_bbs_combined",              "rf_bbs",              rf_lv,  "RF (strata)"),
  list("rf_bbs_site_combined",         "rf_bbs_site",         rf_lv,  "RF (route)"),
  list("gam_bbs_combined",             "gam_bbs",             gam_lv, "GAM (strata)"),
  list("gam_bbs_site_combined",        "gam_bbs_site",        gam_lv, "GAM (route)"),
  # --- SI: un-clamped GAM ---
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
