# Combined BBS main-text figures -- the DEFAULT BBS figure structure. Each figure
# is a 2x3 grid: species-wide over population-level (rows), Temporal / Buffered /
# Spatiotemporal (columns), with the y-scale shared across treatments within each
# row and the population-level row 1/3 the height of species-wide. One figure per
# model family, titled by model, y-axes labelled as correlation (r).
#
# For each family two versions are written:
#   <stem>        - pooled across all species
#   <stem>_mig2   - split by migratory status (Resident vs Migrant, coloured;
#                   AVONET Migration collapsed to sedentary=Resident vs 2/3=Migrant,
#                   see data/species_migration.rds). Grouped bands are dodged side
#                   by side with no median-connecting line.
#
#   gam_bbs_combined       - GAM,  strata (4 models incl. SVC)
#   rf_bbs_combined        - RF,   strata (3 models)
#   rf_bbs_site_combined   - RF,   route/site (3 models)
#   gam_bbs_site_combined  - GAM,  route/site (4 models incl. SVC)

source("scripts/plot_helpers.R")

treatments <- c("Temporal", "Buffered", "Spatiotemporal")

specs <- list(
  list(stem = "gam_bbs_combined", title = "GAM (strata)",
       paths = c("data/gam_bbs_results.rds", "data/gam_bbs_buffer_results.rds",
                 "data/gam_bbs_spatialcv_k8_results.rds"),
       models = c("Static", "Dynamic", "Decomposed", "SVC")),
  list(stem = "rf_bbs_combined", title = "Random Forest (strata)",
       paths = c("data/rf_bbs_results.rds", "data/rf_bbs_buffer_results.rds",
                 "data/rf_bbs_spatialcv_k8_results.rds"),
       models = c("Static", "Dynamic", "Decomposed")),
  list(stem = "rf_bbs_site_combined", title = "Random Forest (route)",
       paths = c("data/rf_bbs_site_results.rds", "data/rf_bbs_site_buffer_results.rds",
                 "data/rf_bbs_site_spatialcv_k8_results.rds"),
       models = c("Static", "Dynamic", "Decomposed")),
  list(stem = "gam_bbs_site_combined", title = "GAM (route)",
       paths = c("data/gam_bbs_site_results.rds", "data/gam_bbs_site_buffer_results.rds",
                 "data/gam_bbs_site_spatialcv_k8_results.rds"),
       models = c("Static", "Dynamic", "Decomposed", "SVC"))
)

for (s in specs) {
  if (!all(file.exists(s$paths))) {
    cat(sprintf("skip %s (missing result files)\n", s$stem)); next
  }
  make_medband_grid(s$paths, treatments, s$stem, s$models, title = s$title)
  make_medband_grid(s$paths, treatments, paste0(s$stem, "_mig2"), s$models,
                    title = s$title, split = "mig2")
}
