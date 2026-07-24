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

# Colour palettes for the migratory/resident split (see data/species_migration.rds).
# Colour encodes GROUP; the metric is already encoded by row.
.sfts_grp_cols <- list(
  mig2 = c(Resident = "#2C7FB8", Migrant = "#D95F02"),
  mig3 = c(Sedentary = "#2C7FB8", Partial = "#8856A7", Migratory = "#D95F02"))

# split = NULL -> pooled (one band per model, current behaviour).
# split = "mig2" or "mig3" -> join data/species_migration.rds and summarise per
# group, dropping species with no trait. Adds a `grp` column (factor).
.sfts_summ <- function(res_path, model_lv, split = NULL) {
  r <- readRDS(res_path)
  if (is.null(split)) {
    long <- r %>% select(ends_with("_sw"), ends_with("_pl")) %>%
      pivot_longer(everything(), names_to = "key", values_to = "r") %>%
      mutate(grp = factor("all"))
  } else {
    tr <- readRDS("data/species_migration.rds") %>% select(species_id, grp = all_of(split))
    long <- r %>% left_join(tr, by = "species_id") %>% filter(!is.na(grp)) %>%
      select(grp, ends_with("_sw"), ends_with("_pl")) %>%
      pivot_longer(c(ends_with("_sw"), ends_with("_pl")), names_to = "key", values_to = "r") %>%
      mutate(grp = factor(grp, levels = names(.sfts_grp_cols[[split]])))
  }
  long %>%
    mutate(metric = if_else(str_ends(key, "_sw"), "sw", "pl"),
      model = recode(str_remove(key, "_(sw|pl)$"),
        static = "Static", dynamic = "Dynamic", decomp = "Decomposed", svc = "SVC")) %>%
    filter(model %in% model_lv) %>%
    group_by(metric, model, grp) %>%
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
  g_sw <- .sfts_panel(filter(d, metric == "sw"), .sfts_cols["sw"], "Species-wide\ncorrelation (r)", FALSE, zero_line = FALSE) +
    ggtitle(title) + theme(plot.title = element_text(size = 12))
  g_pl <- .sfts_panel(filter(d, metric == "pl"), .sfts_cols["pl"], "Population-level\ncorrelation (r)", TRUE, zero_line = TRUE)
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
#   split      : NULL (pooled) or "mig2"/"mig3". When set, each model shows one
#                vertical band PER GROUP, dodged side by side and coloured by
#                group; there is no median-connecting line (the bands stand alone).
make_medband_grid <- function(paths, treatments, out_stem, model_lv,
                              fig_width = 180, fig_height = 120, title = NULL,
                              heights = c(3, 1), labels = TRUE, split = NULL) {
  d <- purrr::imap_dfr(paths, function(p, i)
    mutate(.sfts_summ(p, model_lv, split = split), treatment = treatments[i])) %>%
    mutate(treatment = factor(treatment, levels = treatments))
  grouped <- !is.null(split)
  pal <- if (grouped) .sfts_grp_cols[[split]] else NULL
  pd <- position_dodge(width = 0.6)

  row_panel <- function(dd, colr, ylab, is_top, zero_line, nbreaks) {
    if (grouped) g <- ggplot(dd, aes(model, med, colour = grp, group = grp))
    else         g <- ggplot(dd, aes(model, med, group = 1))
    if (zero_line)
      g <- g + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey65")
    if (grouped) {
      g <- g +
        geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.14, linewidth = 0.7, alpha = 0.5, position = pd) +
        geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 1.7, position = pd) +
        geom_point(size = 2.1, position = pd) +
        scale_colour_manual(values = pal, name = NULL)
    } else {
      g <- g +
        geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.12, linewidth = 0.7, alpha = 0.45, colour = colr) +
        geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 1.7, colour = colr) +
        geom_line(linewidth = 0.8, colour = colr) +
        geom_point(size = 2.3, colour = colr)
    }
    g <- g +
      facet_wrap(~ treatment, nrow = 1) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = nbreaks),
        expand = expansion(mult = if (is_top) c(0.08, 0.20) else c(0.14, 0.14))) +
      coord_cartesian(clip = "off") +
      labs(x = NULL, y = ylab) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
        panel.spacing.x = unit(7, "pt"),
        axis.title.y = element_text(size = 9, margin = margin(r = 2)),
        axis.text.y = element_text(size = 9.5))
    if (labels && !grouped) {
      if (is_top)  # species-wide: float the number above the 95% cap -> no overlap
        g <- g + geom_label(aes(y = q95, label = sprintf("%.2f", med)), colour = colr,
          fill = "white", label.size = 0, label.padding = unit(0.05, "lines"),
          vjust = -0.35, size = 2.3)
      else         # population-level: nudge beside the point (space is tight, overlap ok)
        g <- g + geom_label(aes(label = sprintf("%.2f", med)), colour = colr, fill = "white",
          label.size = 0, label.padding = unit(0.05, "lines"), nudge_x = 0.24, size = 2.2)
    }
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

  p_main <- top / bot + plot_layout(heights = heights, guides = "collect")
  if (grouped) p_main <- p_main & theme(legend.position = "bottom", legend.margin = margin(t = 0, b = 0))
  # outer, larger y-axis title spanning both rows ("Correlation (r)"), set a bit
  # further out than the small per-panel subtitles ("Species-wide"/"Population-level").
  ylab_col <- patchwork::wrap_elements(grid::textGrob("Correlation (r)", rot = 90,
    gp = grid::gpar(fontsize = 13)))
  p <- patchwork::wrap_plots(ylab_col, p_main, widths = c(1, 26))
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
  g_sw <- .sfts_panel(d_sw, .sfts_cols["sw"], "Species-wide\ncorrelation (r)", FALSE, zero_line = FALSE) +
    ggtitle(title) + theme(plot.title = element_text(size = 12))
  g_pl <- .sfts_panel(d_pl, .sfts_cols["pl"], "Population-level\ncorrelation (r)", TRUE, zero_line = TRUE)
  p <- g_sw / g_pl + plot_layout(heights = c(1.6, 1))
  ggsave(sprintf("figures/%s.pdf", out_stem), p, width = fig_width, height = 190, units = "mm", dpi = 600)
  ggsave(sprintf("figures/%s.png", out_stem), p, width = fig_width, height = 190, units = "mm", dpi = 600)
  cat(sprintf("saved -> figures/%s.png\n", out_stem))
}


# Two-axis MAIN figure (CANONICAL). One figure per metric x scale. The two axes of
# changing performance are made explicit:
#   x  = model type (Static -> SVC)                = relaxing the STE assumption ->
#   rows (facet) = treatment (Temporal/Buffered/Spatiotemporal) = harder to transfer (down)
# GAM vs RF are dodged and coloured within each model type (GAM centres at SVC,
# which RF lacks). Median point + IQR (thick) + 5-95% (thin), NO connecting line;
# the median value sits beside each point (GAM to the left, RF to the right). y is
# free-scaled per row, gridlines every 0.2 (species-wide) / 0.1 (population-level).
# Two thick labelled arrows frame the axes. Reads data/{rf,gam}_bbs[_site]_*.rds.
#   metric : "sw" (species-wide) or "pl" (population-level)
#   scale  : "strata" or "site"
#   out    : figures/<out>.{png,pdf}
make_2axis <- function(metric, scale, out, width = 165, height = 178) {
  model_lv <- c("Static","Dynamic","Decomposed","SVC")
  trt_lv   <- c("Temporal","Buffered","Spatiotemporal")
  fam_cols <- c(GAM = "#DD9100", RF = "#254F9A")
  pre  <- if (scale == "strata") c(rf="rf_bbs", gam="gam_bbs") else c(rf="rf_bbs_site", gam="gam_bbs_site")
  trt  <- c(Temporal="results", Buffered="buffer_results", Spatiotemporal="spatialcv_k8_results")
  suff <- paste0("_", metric)
  mlab <- if (metric == "sw") "Species-wide" else "Population-level"
  slab <- if (scale == "strata") "Strata" else "Site"
  ystep <- if (metric == "sw") 0.2 else 0.1

  summ <- function(path) readRDS(path) %>% select(ends_with(suff)) %>%
    pivot_longer(everything(), names_to = "key", values_to = "r") %>%
    mutate(model = recode(str_remove(key, paste0(suff, "$")),
                          static = "Static", dynamic = "Dynamic", decomp = "Decomposed", svc = "SVC")) %>%
    filter(model %in% model_lv) %>% group_by(model) %>%
    summarise(med = median(r, na.rm = TRUE), q05 = quantile(r, .05, na.rm = TRUE),
              q25 = quantile(r, .25, na.rm = TRUE), q75 = quantile(r, .75, na.rm = TRUE),
              q95 = quantile(r, .95, na.rm = TRUE), .groups = "drop")
  fam_paths <- list(
    RF  = setNames(sprintf("data/%s_%s.rds", pre["rf"],  trt), names(trt)),
    GAM = setNames(sprintf("data/%s_%s.rds", pre["gam"], trt), names(trt)))
  d <- purrr::imap_dfr(fam_paths, function(paths, fam)
         purrr::imap_dfr(paths, function(p, tr) mutate(summ(p), family = fam, treatment = tr))) %>%
    mutate(family = factor(family, levels = c("GAM","RF")), model = factor(model, levels = model_lv),
           treatment = factor(treatment, levels = trt_lv))   # no complete() -> GAM centres at SVC

  pd <- position_dodge(width = 0.55)
  main <- ggplot(d, aes(model, med, colour = family, group = family)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey70") +
    geom_linerange(aes(ymin = q05, ymax = q95), linewidth = 0.6, alpha = 0.5, position = pd, na.rm = TRUE) +
    geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 1.9, position = pd, na.rm = TRUE) +
    geom_point(size = 2.1, position = pd, na.rm = TRUE) +
    geom_text(aes(y = med, label = sprintf("%.2f", med), hjust = ifelse(family == "GAM", 1.45, -0.45)),
              position = pd, vjust = 0.4, size = 2.0, show.legend = FALSE, na.rm = TRUE) +
    facet_grid(treatment ~ ., scales = "free_y") +
    scale_colour_manual(values = fam_cols, name = NULL) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(breaks = seq(-1, 1, ystep), expand = expansion(mult = c(0.05, 0.07))) +
    labs(x = NULL, y = paste0(mlab, " correlation (r)")) + theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
          strip.background = element_rect(fill = "grey95", colour = "grey80", linewidth = 0.5),
          strip.text.y = element_text(size = 10), legend.position = "bottom",
          axis.text.x = element_text(size = 10.5), axis.title.y = element_text(size = 10.5, margin = margin(r = 3)))
  top_arrow <- ggplot() +
    geom_segment(aes(x = 0.08, xend = 0.97, y = 0.32, yend = 0.32),
      arrow = arrow(length = unit(3.2, "mm"), type = "closed"), linewidth = 1.8, colour = "grey20") +
    annotate("text", x = 0.52, y = 0.82, label = "STE assumption relaxed", size = 3.9, fontface = 2) +
    xlim(0, 1) + ylim(0, 1) + theme_void()
  left_arrow <- ggplot() +
    geom_segment(aes(x = 0.72, xend = 0.72, y = 0.96, yend = 0.04),
      arrow = arrow(length = unit(3.2, "mm"), type = "closed"), linewidth = 1.8, colour = "grey20") +
    annotate("text", x = 0.32, y = 0.5, label = "Harder to transfer", angle = 90, size = 3.9, fontface = 2) +
    xlim(0, 1) + ylim(0, 1) + theme_void()
  p <- wrap_plots(A = top_arrow, B = left_arrow, C = main, design = "#A\nBC",
                  widths = c(0.07, 1), heights = c(0.08, 1)) +
    plot_annotation(title = paste0(mlab, " — ", slab),
                    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5)))
  ggsave(sprintf("figures/%s.png", out), p, width = width, height = height, units = "mm", dpi = 300)
  ggsave(sprintf("figures/%s.pdf", out), p, width = width, height = height, units = "mm", dpi = 300)
  cat(sprintf("saved -> figures/%s.png\n", out))
}
