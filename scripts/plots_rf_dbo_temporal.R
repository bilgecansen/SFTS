# DBO temporal figure (house median-band style). The spatial (species-wide) prediction
# level is not shown -- this figure is about temporal predictive ability only, at the
# two temporal levels:
#   regional          (top)    : stations averaged within a DBO region per year, then
#                                the regional trajectory correlated over years
#   population-level  (bottom) : within-station correlation over years
# One box per model = spread ACROSS FAMILIES (box = IQR, whiskers = Tukey 1.5*IQR,
# median on a white chip), not a confidence interval. Both panels anchored at 0.
#
# Models are the standard four; Decomposed and SVC drop the *_spatial block (see
# scripts/fit_rf_dbo_temporal.R for why -- station-intercept artifact under LOYO).
#
# Input: data/rf_dbo_temporal_results.rds (from fit_rf_dbo_temporal.R)

source("scripts/plot_helpers.R")

model_lv <- c("Static", "Dynamic", "Decomposed", "SVC")
cols <- c(rg = "#1B7837", pl = "#254F9A")   # regional green, population-level blue

res <- readRDS("data/rf_dbo_temporal_results.rds")
summ <- res %>%
  select(ends_with("_rg"), ends_with("_pl")) %>%
  pivot_longer(everything(), names_to = "key", values_to = "r") %>%
  filter(!is.na(r)) %>%
  mutate(metric = if_else(str_ends(key, "_rg"), "rg", "pl"),
    model = recode(str_remove(key, "_(rg|pl)$"),
      static = "Static", dynamic = "Dynamic", decomp = "Decomposed", svc = "SVC"),
    model = factor(model, levels = model_lv)) %>%
  group_by(metric, model) %>%
  summarise(med = median(r), q25 = quantile(r, .25), q75 = quantile(r, .75),
    ymin = min(r[r >= quantile(r, .25) - 1.5 * IQR(r)]),
    ymax = max(r[r <= quantile(r, .75) + 1.5 * IQR(r)]),
    .groups = "drop")

g_rg <- .sfts_panel(filter(summ, metric == "rg"), cols["rg"], "Regional", FALSE, zero_line = TRUE) +
  ggtitle(sprintf("DBO - RF temporal prediction, leave-one-year-out (%d families)", nrow(res))) +
  theme(plot.title = element_text(size = 10.5))
g_pl <- .sfts_panel(filter(summ, metric == "pl"), cols["pl"], "Population-level", TRUE, zero_line = TRUE)

p_main <- g_rg / g_pl + plot_layout(heights = c(1, 1))
ylab_col <- patchwork::wrap_elements(grid::textGrob("Correlation (r)", rot = 90,
  gp = grid::gpar(fontsize = 13)))
p <- patchwork::wrap_plots(ylab_col, p_main, widths = c(1, 26))

ggsave("figures/dbo_rf_temporal.pdf", p, width = 175, height = 165, units = "mm", dpi = 600)
ggsave("figures/dbo_rf_temporal.png", p, width = 175, height = 165, units = "mm", dpi = 600)
cat("saved -> figures/dbo_rf_temporal.png\n")
print(as.data.frame(summ %>% transmute(metric, model, med = round(med, 3))), row.names = FALSE)
