# Main random-forest figure for the manuscript.
#
# A condensed view of the four per-component supplementary figures
# (plot_dbo_ladder_{spatial,global,regional,station}.R). Those show every model type;
# this one drops that axis and reports, for each cell, the skill of the BEST-PERFORMING
# model among static, dynamic and decomposed + Lat/Long. Which model that is varies by
# cell and is not shown; it is printed to the console and belongs in the supplement.
#
# Reads data/rf_ladder_bu_preds.rds (fit_rf_ladder_bu.R) -- nothing is refitted, those
# predictions repartitioned.
#
# ICE. The decomposed + Lat/Long model additionally carries breakup_global,
# breakup_regional and breakup_residual. Static and dynamic are unchanged and ice-free.
# breakup only, and its temporal components only: between stations the ice variables
# correlate 0.77-0.97 with each other and 0.99 with latitude, so a spatial ice term
# re-labels position that Lat/Long already supply. The four supplementary per-component
# figures (plot_dbo_ladder_*.R) still read the ice-free data/rf_dbo_ladder_preds.rds. Every part is a piece of the same four-way decomposition (station ->
# year -> region:year -> residual) and each row's parts recombine to the canonical
# component, asserted below.
#
# TWO FIGURES, one per holdout, panels side by side and colour the response definition:
#
#   dbo_main_loso        Abundance distribution | Station trajectory
#   dbo_main_traintest   Abundance distribution | Global and regional trajectory |
#                        Station trajectory
#
# The trajectory panel is omitted under LOSO by design: the network-wide year effect and
# the bands' departures from it are reconstructed from stations sampled in the same years,
# and the user's judgement is that this does not address the question the paper asks. The
# numbers remain in the supplementary per-component figures.
#
#   Abundance distribution        across regions, then within each band. The station
#                                 effect split into its between- and within-band parts.
#   Global and regional           the network-wide year effect (leftmost), then each
#     trajectory                  band's departure from it.
#   Station trajectory            what is left at one station in one year, by band.
#
# Both figures share one y range so the two holdouts can be read against each other.
#
# SELECTION CAVEAT. Taking the maximum over three models is selection on the test data
# and is optimistic relative to naming a model in advance. The console output reports how
# much the maximum exceeds the decomposed + Lat/Long model, which is the natural
# pre-specified choice, so the size of that optimism is on record.
#
# Output: figures/dbo_main_loso.{png,pdf}, figures/dbo_main_traintest.{png,pdf}

suppressMessages(library(tidyverse))

MIN_RESID_DF <- 15
HOLD <- c("LOSO", "Train/test")
RUNG_LAB <- c(family = "Family (57)", class = "Class (7)", total = "Total (1)")
COLS <- c("Family (57)" = "#9ECAE1", "Class (7)" = "#4292C6", "Total (1)" = "#08519C")
PT <- 2.8; FLOOR <- -0.5; GUT <- -0.72; HIDE <- 99; TOP <- 1.04

p <- readRDS("data/rf_ladder_bu_preds.rds") %>%
  filter(holdout %in% HOLD) %>% mutate(ry = paste(DBOreg, year))
regs <- sort(unique(p$DBOreg))

# ---- every part of the decomposition, per taxon --------------------------------------
parts_of <- function(t) {
  g <- function(v) {
    f_reg <- fitted(lm(v ~ factor(t$DBOreg)))
    f_stn <- fitted(lm(v ~ factor(t$DBOreg) + factor(t$stn)))
    f1 <- fitted(lm(v ~ factor(t$stn)))
    f2 <- fitted(lm(v ~ factor(t$stn) + factor(t$year)))
    f3 <- fitted(lm(v ~ factor(t$stn) + factor(t$year) + factor(t$ry)))
    by_reg <- function(x) as.numeric(tapply(x, t$DBOreg, sum)[regs])
    c(sum((f_reg - mean(v))^2), by_reg((f_stn - f_reg)^2),    # abundance: across, within
      sum((f2 - f1)^2), by_reg((f3 - f2)^2),                  # global, regional
      by_reg((v - f3)^2))                                     # station
  }
  key <- c("ab_across", paste0("ab_", regs), "tr_global", paste0("tr_", regs),
           paste0("st_", regs))
  tibble(part = key, e = g(t$y - t$p), o = g(t$y)) %>% replace_na(list(e = 0, o = 0))
}

acc <- p %>% group_split(rung, holdout, model) %>% map_dfr(function(x) {
  r <- x %>% group_split(taxon) %>% map_dfr(function(t) {
    if (nrow(t) - (n_distinct(t$stn) + n_distinct(t$ry)) < MIN_RESID_DF) return(NULL)
    parts_of(t) })
  if (!nrow(r)) return(NULL)
  r %>% group_by(part) %>% summarise(SS_obs = sum(o), SSE = sum(e), .groups = "drop") %>%
    mutate(rung = x$rung[1], holdout = x$holdout[1], model = x$model[1],
           skill = 1 - SSE / SS_obs)
})

# the parts must recombine to the canonical components
old <- readRDS("data/rf_ladder_bu_skill.rds")
chk <- function(pre, comp) {
  a <- acc %>% filter(rung == "family", holdout == "LOSO", model == "svc",
                      str_starts(part, pre)) %>%
    summarise(s = 1 - sum(SSE) / sum(SS_obs)) %>% pull(s)
  b <- old %>% filter(rung == "family", holdout == "LOSO", model == "svc",
                      component == comp) %>% pull(skill)
  cat(sprintf("  %-10s recombined %+.3f vs canonical %+.3f\n", comp, a, b))
  stopifnot(abs(a - b) < 0.005)
}
cat("family / LOSO / decomposed + Lat/Long:\n")
chk("ab_", "abundance"); chk("tr_g", "global"); chk("tr_[123]", "regional")
chk("st_", "station")

# ---- best model per cell ---------------------------------------------------------------
best <- acc %>% group_by(rung, holdout, part) %>%
  slice_max(skill, n = 1, with_ties = FALSE) %>% ungroup()

svc <- acc %>% filter(model == "svc") %>% select(rung, holdout, part, svc = skill)
cat(sprintf("\nbest model beats decomposed + Lat/Long by a median of %+.3f (max %+.3f)\n",
            median(left_join(best, svc, by = c("rung","holdout","part"))$skill -
                   left_join(best, svc, by = c("rung","holdout","part"))$svc),
            max(left_join(best, svc, by = c("rung","holdout","part"))$skill -
                left_join(best, svc, by = c("rung","holdout","part"))$svc)))
cat("\nwhich model wins, by cell\n")
print(as.data.frame(best %>% count(holdout, model) %>%
        pivot_wider(names_from = model, values_from = n, values_fill = 0)),
      row.names = FALSE)

# ---- assemble ---------------------------------------------------------------------------
LAB <- c(ab_across = "Across\nregions", ab_1 = "DBO 1", ab_2 = "DBO 2", ab_3 = "DBO 3",
         tr_global = "Global", tr_1 = "DBO 1", tr_2 = "DBO 2", tr_3 = "DBO 3",
         st_1 = "DBO 1", st_2 = "DBO 2", st_3 = "DBO 3")
# One shared level order across panels, with the two aggregate categories ahead of the
# bands. free_x drops the unused level in each panel, so "Across regions" and "Global"
# each sit leftmost in the panel where they apply. Ordering by first appearance instead
# would put "Global" last, since "DBO 1" is already used by the abundance panel.
XLEV <- c("Across\nregions", "Global", "DBO 1", "DBO 2", "DBO 3")
PANEL <- list(
  ab = list(lev = c("ab_across","ab_1","ab_2","ab_3"), title = "Abundance distribution"),
  tr = list(lev = c("tr_global","tr_1","tr_2","tr_3"),
            title = "Global and regional trajectory"),
  st = list(lev = c("st_1","st_2","st_3"), title = "Station trajectory"))

dodge <- position_dodge(width = 0.68, preserve = "single")

# facet_grid with a single row: scales = "free_x" gives each panel its own categories and
# space = "free_x" makes panel width proportional to how many it has
make_fig <- function(hold, keys, stem, width) {
  lev <- unlist(map(PANEL[keys], "lev"))
  ttl <- map_chr(PANEL[keys], "title")
  d <- best %>% filter(holdout == hold, part %in% lev) %>%
    mutate(panel = factor(ttl[str_sub(part, 1, 2)], ttl),
           x = factor(LAB[part], XLEV),
           rung = factor(RUNG_LAB[rung], RUNG_LAB),
           off = skill < FLOOR,
           y_pt = ifelse(off, HIDE, skill), y_tri = ifelse(off, FLOOR, HIDE),
           y_lab = ifelse(off, FLOOR - 0.055, HIDE),
           lab = ifelse(off, sprintf("%.2f", skill), ""))
  stopifnot(!any(is.na(d$panel)), !any(is.na(d$x)))

  g <- ggplot(d, aes(x, skill, colour = rung, group = rung)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = GUT, ymax = FLOOR,
             fill = "grey94", colour = NA) +
    geom_hline(yintercept = FLOOR, colour = "grey70", linewidth = 0.3, linetype = "13") +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_linerange(aes(ymin = 0, ymax = pmax(skill, FLOOR)), position = dodge,
                   linewidth = 0.55) +
    geom_point(aes(y = y_pt), size = PT, position = dodge) +
    geom_point(aes(y = y_tri), shape = 17, size = 1.5, position = dodge,
               show.legend = FALSE) +
    geom_text(aes(y = y_lab, label = lab), position = dodge, angle = 90, hjust = 1,
              size = 1.95, show.legend = FALSE) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
    scale_colour_manual(values = COLS, name = "Response") +
    scale_y_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       labels = function(z) sprintf("%.1f", z)) +
    coord_cartesian(ylim = c(GUT, TOP)) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(x = NULL, y = "Skill  (1 - SSE / SS)", title = hold) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = "grey80"),
          plot.title = element_text(size = 11),
          axis.text.x = element_text(size = 7.5, lineheight = 0.9),
          legend.position = "bottom")

  ggsave(sprintf("figures/%s.png", stem), g, width = width, height = 105,
         units = "mm", dpi = 300)
  ggsave(sprintf("figures/%s.pdf", stem), g, width = width, height = 105, units = "mm")
  cat(sprintf("saved -> figures/%s.png\n", stem))
}

make_fig("LOSO",       c("ab", "st"),       "dbo_main_loso",      140)
make_fig("Train/test", c("ab", "tr", "st"), "dbo_main_traintest", 190)
cat("\n")

print(as.data.frame(best %>%
        transmute(holdout, rung, part, skill = round(skill, 3)) %>%
        pivot_wider(names_from = part, values_from = skill) %>%
        arrange(holdout, rung)), row.names = FALSE)
