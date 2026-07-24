# Effect-size tables that accompany the combined figures: for each non-static
# model, the paired per-species improvement over the STATIC model -- median delta
# (Fisher-z is unnecessary for a median of differences) and the % of species that
# improve. Comparisons are within-species paired (same held-out data), so delta is
# taken row-wise (per species) and summarised across species. Written per model
# family to tables/ as CSV, and printed.

library(tidyverse)

model_lab <- c(dynamic = "Dynamic", decomp = "Decomposed", svc = "SVC")

# per-metric paired improvement over static, for the given results tibble
one_metric <- function(r, suf, models) {
  base <- r[[paste0("static_", suf)]]
  map_dfr(models, function(m) {
    col <- paste0(m, "_", suf)
    if (!col %in% names(r)) return(NULL)
    d <- r[[col]] - base
    d <- d[is.finite(d)]
    tibble(model = model_lab[[m]], n = length(d),
      median_delta = median(d), q05_delta = quantile(d, 0.05), q95_delta = quantile(d, 0.95),
      pct_improving = 100 * mean(d > 0))
  })
}

# BBS families: three treatment files sharing the same column layout
delta_bbs <- function(paths, treatments, models) {
  map2_dfr(paths, treatments, function(p, tt) {
    r <- readRDS(p)
    bind_rows(
      mutate(one_metric(r, "sw", models), metric = "Species-wide"),
      mutate(one_metric(r, "pl", models), metric = "Population-level")
    ) %>% mutate(treatment = tt)
  }) %>% mutate(treatment = factor(treatment, levels = treatments))
}

treatments <- c("Temporal", "Buffered", "Spatiotemporal")
fam <- list(
  `GAM (strata)` = list(
    paths = c("data/gam_bbs_results.rds", "data/gam_bbs_buffer_results.rds", "data/gam_bbs_spatialcv_k8_results.rds"),
    models = c("dynamic", "decomp", "svc")),
  `RF (strata)` = list(
    paths = c("data/rf_bbs_results.rds", "data/rf_bbs_buffer_results.rds", "data/rf_bbs_spatialcv_k8_results.rds"),
    models = c("dynamic", "decomp")),
  `RF (site-level)` = list(
    paths = c("data/rf_bbs_site_results.rds", "data/rf_bbs_site_buffer_results.rds", "data/rf_bbs_site_spatialcv_k8_results.rds"),
    models = c("dynamic", "decomp")),
  `GAM + spatiotemporal smoother` = list(
    paths = c("data/gam_bbs_st_results.rds", "data/gam_bbs_st_buffer_results.rds", "data/gam_bbs_st_spatialcv_k8_results.rds"),
    models = c("dynamic", "decomp", "svc")),
  `Linear GAM` = list(
    paths = c("data/lin_bbs_results.rds", "data/lin_bbs_buffer_results.rds", "data/lin_bbs_spatialcv_k8_results.rds"),
    models = c("dynamic", "decomp", "svc")))

dir.create("tables", showWarnings = FALSE)

pretty_print <- function(tab, title) {
  cat(sprintf("\n=== %s -- improvement over Static ===\n", title))
  tab %>%
    mutate(metric = factor(metric, levels = c("Species-wide", "Population-level"))) %>%
    arrange(metric, treatment, factor(model, levels = model_lab)) %>%
    mutate(`median Δ` = sprintf("%+.2f", median_delta),
      `5-95%` = sprintf("%+.2f to %+.2f", q05_delta, q95_delta),
      `% improving` = sprintf("%.0f%%", pct_improving)) %>%
    select(metric, treatment, model, n, `median Δ`, `5-95%`, `% improving`) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

all_tab <- imap_dfr(fam, function(f, nm) {
  if (!all(file.exists(f$paths))) { cat(sprintf("skip (missing files): %s\n", nm)); return(NULL) }
  tab <- delta_bbs(f$paths, treatments, f$models) %>% mutate(family = nm)
  write_csv(select(tab, family, metric, treatment, model, n, median_delta, q05_delta, q95_delta, pct_improving),
    sprintf("tables/delta_%s.csv", str_replace_all(tolower(nm), "[^a-z0-9]+", "_")))
  pretty_print(tab, nm)
  tab
})

# DBO: species-wide from LOSO, population-level from LOYO (its two honest holdouts)
if (all(file.exists(c("data/rf_dbo_loso_results.rds", "data/rf_dbo_results.rds")))) {
  loso <- readRDS("data/rf_dbo_loso_results.rds"); loyo <- readRDS("data/rf_dbo_results.rds")
  dbo <- bind_rows(
    mutate(one_metric(loso, "sw", c("dynamic", "decomp")), metric = "Species-wide (LOSO)"),
    mutate(one_metric(loyo, "pl", c("dynamic", "decomp")), metric = "Population-level (LOYO)")) %>%
    mutate(family = "DBO (RF)")
  write_csv(select(dbo, family, metric, model, n, median_delta, q05_delta, q95_delta, pct_improving), "tables/delta_dbo.csv")
  cat("\n=== DBO (RF) -- improvement over Static ===\n")
  dbo %>% mutate(`median Δ` = sprintf("%+.2f", median_delta),
      `5-95%` = sprintf("%+.2f to %+.2f", q05_delta, q95_delta),
      `% improving` = sprintf("%.0f%%", pct_improving)) %>%
    select(metric, model, n, `median Δ`, `5-95%`, `% improving`) %>% as.data.frame() %>% print(row.names = FALSE)
}

write_csv(all_tab, "tables/delta_all_bbs.csv")
cat("\nsaved -> tables/delta_*.csv\n")
