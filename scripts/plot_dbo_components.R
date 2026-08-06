# DBO component figure -- the real-data analogue of the simulation panels in
# figures/figure2_plots.R, but built on the orthogonal error decomposition rather than on
# four different aggregations of the same points.
#
# Each panel scatters the OBSERVED against the PREDICTED projection onto one subspace of
#     station  ->  year  ->  region:year  ->  residual
# mapped to the four prediction levels:
#     Abundance distribution : station effect      (spatial, incl. the latitudinal gradient)
#     Global trajectory      : year effect         (the signal common to all three bands)
#     Regional trajectory    : region:year effect  (how the bands' trajectories diverge)
#     Station trajectory     : residual            (station-specific deviation)
# Because those projections are mutually orthogonal, the four panels partition the total
# error exactly: no point's error is counted twice, and the SSE shares sum to 100%.
#
# THE REFERENCE LINE IS 1:1, NOT A REGRESSION FIT. The simulation panels drew lm fits
# because a correlation is invariant to slope. Skill is not -- it charges for every
# departure from the identity -- so the diagonal is the correct reference. Three line types,
# drawn bottom to top so each stays legible against the next:
#   faint thin solid  = 1:1 identity (the skill reference)
#   heavy solid       = pooled fit across all points, the panel summary
#   dashed            = one fit per family, so family spread is not averaged away
# The gap between the heavy solid and the faint 1:1 is the amplitude error that skill
# charges for and correlation ignores.
#
# Predictions carry only ~40% of the observed amplitude, which shows up as clouds that are
# NARROW IN X (little predicted spread) with family lines STEEPER than the diagonal -- if
# pred ~ a*obs with a < 1, regressing obs on pred returns slope 1/a > 1. The extreme case is
# the station-trajectory panels, where clouds are essentially vertical: observations vary
# widely, predictions barely at all.
#
# Rows are the two holdouts. The global-trajectory panel is the sharpest contrast: tight on
# the diagonal under LOSO (skill +0.87), collapsed to a blob under train/test (-0.00).
#
# Decomposed model only; raw component values (so the annotated POOLED skill and the figure
# agree). Points are deduplicated per panel -- a station effect is constant within a
# station, a year effect within a year, so plotting one point per station-year would
# replicate identical values.
#
# Input : data/rf_dbo_4rt_{loso,traintest}_preds.rds
# Output: figures/dbo_components.{png,pdf}

library(tidyverse)

PC <- "p_decomp"
lev <- c("Abundance distribution", "Global trajectory", "Regional trajectory",
         "Station trajectory")
cols <- setNames(c("#DD9100", "#5AAE61", "#1B7837", "#254F9A"), lev)

# sequential orthogonal projections, returned as VECTORS (not just sums of squares)
components <- function(v, stn, yr, ry) {
  f0 <- rep(mean(v), length(v))
  f1 <- fitted(lm(v ~ factor(stn)))
  f2 <- fitted(lm(v ~ factor(stn) + factor(yr)))
  f3 <- fitted(lm(v ~ factor(stn) + factor(yr) + factor(ry)))
  tibble(ad = f1 - f0, gt = f2 - f1, rt = f3 - f2, st = v - f3)
}

build <- function(path, lab) {
  p <- readRDS(path) %>% mutate(ry = paste(DBOreg, year))
  p %>% group_split(family) %>% map_dfr(function(d) {
    if (nrow(d) - (n_distinct(d$stn) + n_distinct(d$ry)) < 15) return(NULL)
    O <- components(d$y, d$stn, d$year, d$ry)
    P <- components(d[[PC]], d$stn, d$year, d$ry)
    tibble(family = d$family[1], stn = d$stn, year = d$year, ry = d$ry,
           o_ad = O$ad, p_ad = P$ad, o_gt = O$gt, p_gt = P$gt,
           o_rt = O$rt, p_rt = P$rt, o_st = O$st, p_st = P$st)
  }) %>% mutate(holdout = lab)
}

cmp <- bind_rows(build("data/rf_dbo_4rt_loso_preds.rds",      "LOSO"),
                 build("data/rf_dbo_4rt_traintest_preds.rds", "Train/test"))

# deduplicate: each component is constant within its own grouping
long <- bind_rows(
  cmp %>% distinct(holdout, family, stn,  .keep_all = TRUE) %>%
    transmute(holdout, family, obs = o_ad, pred = p_ad, component = lev[1]),
  cmp %>% distinct(holdout, family, year, .keep_all = TRUE) %>%
    transmute(holdout, family, obs = o_gt, pred = p_gt, component = lev[2]),
  cmp %>% distinct(holdout, family, ry,   .keep_all = TRUE) %>%
    transmute(holdout, family, obs = o_rt, pred = p_rt, component = lev[3]),
  cmp %>% transmute(holdout, family, obs = o_st, pred = p_st, component = lev[4])) %>%
  mutate(component = factor(component, lev),
         holdout = factor(holdout, c("LOSO", "Train/test")))

# pooled skill + share of total SSE, on the FULL (undeduplicated) cells
sse_tab <- cmp %>% group_by(holdout) %>% summarise(
  !!lev[1] := sum((o_ad - p_ad)^2), !!lev[2] := sum((o_gt - p_gt)^2),
  !!lev[3] := sum((o_rt - p_rt)^2), !!lev[4] := sum((o_st - p_st)^2), .groups = "drop") %>%
  pivot_longer(-holdout, names_to = "component", values_to = "sse") %>%
  group_by(holdout) %>% mutate(share = sse / sum(sse)) %>% ungroup()
skl_tab <- cmp %>% group_by(holdout) %>% summarise(
  !!lev[1] := 1 - sum((o_ad - p_ad)^2)/sum(o_ad^2),
  !!lev[2] := 1 - sum((o_gt - p_gt)^2)/sum(o_gt^2),
  !!lev[3] := 1 - sum((o_rt - p_rt)^2)/sum(o_rt^2),
  !!lev[4] := 1 - sum((o_st - p_st)^2)/sum(o_st^2), .groups = "drop") %>%
  pivot_longer(-holdout, names_to = "component", values_to = "skill")
ann <- left_join(skl_tab, sse_tab, by = c("holdout", "component")) %>%
  mutate(component = factor(component, lev),
         holdout = factor(holdout, c("LOSO", "Train/test")),
         label = sprintf("skill %+.2f\n%.0f%% of SSE", skill, 100 * share))

# square, symmetric limits per panel so the 1:1 line sits at 45 degrees
lim <- long %>% group_by(holdout, component) %>%
  summarise(m = max(abs(c(obs, pred))), .groups = "drop")
blank <- bind_rows(lim %>% transmute(holdout, component, obs = m, pred = m),
                   lim %>% transmute(holdout, component, obs = -m, pred = -m))

p <- ggplot(long, aes(pred, obs)) +
  geom_blank(data = blank) +
  geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.3) +
  geom_point(aes(colour = component), size = 0.55, alpha = 0.18, show.legend = FALSE) +
  # 1:1 reference drawn UNDER the fits and faintest of all, so it reads as a background
  # reference rather than competing with the family lines. Alpha is baked into the colour
  # because geom_smooth(se = FALSE) does not reliably honour an alpha parameter on the line.
  geom_abline(slope = 1, intercept = 0,
              colour = scales::alpha("black", 0.18), linewidth = 0.4) +
  # pooled fit across all points in the panel, drawn UNDER the family fits: heavy enough to
  # read as the summary, but translucent so it does not mask the fan it summarises. Its
  # slope against the 1:1 line is the amplitude error that skill charges for.
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x,
              colour = scales::alpha("black", 0.5), linewidth = 0.9) +
  # one dashed fit per family, alpha-blended so all are legible together
  geom_smooth(aes(group = family), method = "lm", se = FALSE, formula = y ~ x,
              colour = scales::alpha("grey15", 0.30), linetype = 2, linewidth = 0.5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.12, vjust = 1.15, size = 2.9, fontface = 2, lineheight = 0.95) +
  scale_colour_manual(values = cols) +
  facet_grid(holdout ~ component, scales = "free") +
  labs(x = "Predicted component", y = "Observed component") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        strip.text = element_text(size = 9),
        aspect.ratio = 1)

ggsave("figures/dbo_components.png", p, width = 250, height = 140, units = "mm", dpi = 300)
ggsave("figures/dbo_components.pdf", p, width = 250, height = 140, units = "mm", dpi = 300)
cat("saved -> figures/dbo_components.png\n")
print(as.data.frame(ann %>% transmute(holdout, component, skill = round(skill, 3),
                                      sse_share = sprintf("%.0f%%", 100 * share))),
      row.names = FALSE)
