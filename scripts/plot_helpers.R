# Shared median-band plot for the SFTS BBS/DBO figures.
# Two stacked panels (species-wide on top, population-level below), each with a
# bold IQR (25-75%, no caps) over a light capped 5-95% band, medians connected by
# a line, and white-backed value labels. Each panel auto-ranges to its own data
# (species-wide has no zero line; population-level is anchored at 0). The
# population-level panel is shorter, matching its smaller spread.
# NOTE: bands are the SPREAD ACROSS SPECIES/FAMILIES, not a confidence interval.

library(tidyverse)
library(patchwork)

.sfts_cols <- c(sw = "#792C1F", pl = "#485B7C")

.sfts_summ <- function(res_path, model_lv) {
  readRDS(res_path) %>%
    select(ends_with("_sw"), ends_with("_pl")) %>%
    pivot_longer(everything(), names_to = "key", values_to = "r") %>%
    mutate(metric = if_else(str_ends(key, "_sw"), "sw", "pl"),
      model = recode(str_remove(key, "_(sw|pl)$"),
        static = "Static", dynamic = "Dynamic", decomp = "Decomposed", svc = "SVC")) %>%
    filter(model %in% model_lv) %>%
    group_by(metric, model) %>%
    summarise(med = median(r, na.rm = TRUE), q05 = quantile(r, .05, na.rm = TRUE),
      q25 = quantile(r, .25, na.rm = TRUE), q75 = quantile(r, .75, na.rm = TRUE),
      q95 = quantile(r, .95, na.rm = TRUE), .groups = "drop") %>%
    mutate(model = factor(model, levels = model_lv))
}

.sfts_panel <- function(dd, colr, ylab, show_x, zero_line) {
  g <- ggplot(dd, aes(model, med, group = 1))
  if (zero_line) g <- g + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey65")
  g <- g +
    geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.12, linewidth = 0.8, alpha = 0.45, colour = colr) +
    geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 1.9, colour = colr) +
    geom_line(linewidth = 0.9, colour = colr) +
    geom_point(size = 3, colour = colr) +
    geom_label(aes(label = sprintf("%.2f", med)), colour = colr, fill = "white",
      label.size = 0, label.padding = unit(0.1, "lines"), nudge_x = 0.17, size = 3.3) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = expansion(mult = c(0.06, 0.10))) +
    coord_cartesian(clip = "off") +
    labs(y = ylab) +
    theme_bw(base_size = 12) +
    theme(axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 11.5, margin = margin(r = 8)),
      panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
      plot.margin = margin(t = 6, r = 14, b = 6, l = 4))
  if (show_x) g + theme(axis.text.x = element_text(size = 13))
  else g + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

# res_path : rds with <model>_sw / <model>_pl columns for model in
#            {static, dynamic, decomp, svc}. model_lv selects/orders the models.
make_medband <- function(res_path, out_stem, title, model_lv, fig_width = 175) {
  d <- .sfts_summ(res_path, model_lv)
  g_sw <- .sfts_panel(filter(d, metric == "sw"), .sfts_cols["sw"], "Species-wide", FALSE, zero_line = FALSE) +
    ggtitle(title) + theme(plot.title = element_text(size = 12))
  g_pl <- .sfts_panel(filter(d, metric == "pl"), .sfts_cols["pl"], "Population-level", TRUE, zero_line = TRUE)
  p <- g_sw / g_pl + plot_layout(heights = c(1.6, 1))
  ggsave(sprintf("figures/%s.pdf", out_stem), p, width = fig_width, height = 190, units = "mm", dpi = 600)
  ggsave(sprintf("figures/%s.png", out_stem), p, width = fig_width, height = 190, units = "mm", dpi = 600)
  cat(sprintf("saved -> figures/%s.png\n", out_stem))
}

# Combined grid: all treatments in one figure, metric as rows (species-wide over
# population-level) and treatment as columns. Built as two facet_wrap(~treatment)
# rows stacked with patchwork; plot_layout(heights) makes the population-level row
# shorter (default 1/3 of species-wide), matching its smaller spread. The y-scale
# is shared across treatments WITHIN each row, so the species-wide decline from
# Temporal -> Spatiotemporal reads as the points stepping down. Treatment strips
# sit on the top row only; model labels on the bottom row only. Journal full width.
#   paths      : character vector of result rds files, one per treatment
#   treatments : column labels, same length/order as paths
#   heights    : c(top, bottom) relative row heights
make_medband_grid <- function(paths, treatments, out_stem, model_lv,
                              fig_width = 180, fig_height = 120, title = NULL,
                              heights = c(3, 1), labels = FALSE) {
  d <- purrr::imap_dfr(paths, function(p, i)
    mutate(.sfts_summ(p, model_lv), treatment = treatments[i])) %>%
    mutate(treatment = factor(treatment, levels = treatments))

  row_panel <- function(dd, colr, ylab, is_top, zero_line, nbreaks) {
    g <- ggplot(dd, aes(model, med, group = 1))
    if (zero_line)
      g <- g + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey65")
    g <- g +
      geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.12, linewidth = 0.7, alpha = 0.45, colour = colr) +
      geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 1.7, colour = colr) +
      geom_line(linewidth = 0.8, colour = colr) +
      geom_point(size = 2.3, colour = colr) +
      facet_wrap(~ treatment, nrow = 1) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = nbreaks),
        expand = expansion(mult = c(0.10, 0.12))) +
      coord_cartesian(clip = "off") +
      labs(x = NULL, y = ylab) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
        panel.spacing.x = unit(7, "pt"),
        axis.title.y = element_text(size = 11, margin = margin(r = 6)),
        axis.text.y = element_text(size = 9.5))
    if (labels)
      g <- g + geom_label(aes(label = sprintf("%.2f", med)), colour = colr, fill = "white",
        label.size = 0, label.padding = unit(0.08, "lines"), nudge_x = 0.16, size = 2.6)
    if (is_top)
      g + theme(strip.background = element_rect(fill = "grey95", colour = "grey80", linewidth = 0.5),
        strip.text = element_text(size = 10.5),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())
    else
      g + theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x = element_text(size = 9.5, angle = 30, hjust = 1))
  }

  top <- row_panel(filter(d, metric == "sw"), unname(.sfts_cols["sw"]),
    "Species-wide", TRUE, FALSE, nbreaks = 5)
  bot <- row_panel(filter(d, metric == "pl"), unname(.sfts_cols["pl"]),
    "Population-level", FALSE, TRUE, nbreaks = 3)
  if (!is.null(title))
    top <- top + ggtitle(title) + theme(plot.title = element_text(size = 12))

  p <- top / bot + plot_layout(heights = heights)
  ggsave(sprintf("figures/%s.pdf", out_stem), p, width = fig_width, height = fig_height, units = "mm", dpi = 600)
  ggsave(sprintf("figures/%s.png", out_stem), p, width = fig_width, height = fig_height, units = "mm", dpi = 600)
  cat(sprintf("saved -> figures/%s.png\n", out_stem))
}

# Two-source variant: the species-wide panel is drawn from sw_path and the
# population-level panel from pl_path. Used for DBO, where each metric is taken
# from the holdout that honestly tests it -- species-wide from leave-one-STATION-
# out (spatial extrapolation) and population-level from leave-one-YEAR-out
# (temporal extrapolation). Both files must cover the same taxa for the panels to
# be comparable.
make_medband_split <- function(sw_path, pl_path, out_stem, title, model_lv, fig_width = 175) {
  d_sw <- filter(.sfts_summ(sw_path, model_lv), metric == "sw")
  d_pl <- filter(.sfts_summ(pl_path, model_lv), metric == "pl")
  g_sw <- .sfts_panel(d_sw, .sfts_cols["sw"], "Species-wide", FALSE, zero_line = FALSE) +
    ggtitle(title) + theme(plot.title = element_text(size = 12))
  g_pl <- .sfts_panel(d_pl, .sfts_cols["pl"], "Population-level", TRUE, zero_line = TRUE)
  p <- g_sw / g_pl + plot_layout(heights = c(1.6, 1))
  ggsave(sprintf("figures/%s.pdf", out_stem), p, width = fig_width, height = 190, units = "mm", dpi = 600)
  ggsave(sprintf("figures/%s.png", out_stem), p, width = fig_width, height = 190, units = "mm", dpi = 600)
  cat(sprintf("saved -> figures/%s.png\n", out_stem))
}
