# Prediction-skill figure for the DBO random-forest analysis.
#
# Reads data/rf_dbo_ladder_skill.rds. Rows are the three model types (static / dynamic /
# decomposed + coordinates), columns the three holdouts, and every panel shares one y
# axis.
#
# WITHIN a panel the x axis is the four prediction components, and the three response
# definitions -- 57 families, 7 classes, one total -- are drawn side by side at each
# component. They are NOT connected: taxonomic aggregation has no consistent effect on
# skill, so the three are replicate estimates of the same quantity rather than a
# gradient, and a connecting line would assert a trend that is not there. Their
# agreement is the point.
#
# Each estimate is drawn as a stem from zero, because skill = 1 - SSE/SS_obs is defined
# against zero: at zero the model does no better than predicting no variation along that
# axis, and below zero it does worse.
#
# THE AXIS IS CUT AT -0.5. Eleven estimates fall below it, nine of them the LORO
# abundance distribution, and at full extent they would take two thirds of the panel
# height to display values that all mean the same thing -- worse than predicting no
# variation -- while compressing the range where the results actually differ. Skill
# values are NOT transformed or rescaled. The stem is drawn to the cut, a filled triangle
# marks the estimate as off-scale, and its exact value is printed in the shaded gutter
# below. Everything at or above -0.5 is drawn on a plain linear scale.
#
# POINT SIZE is that component's share of the observed variance within its rung. Without
# it the figure over-reads small components: the global trajectory carries 13% of the
# observed variance at the family rung but only 3-4% at class and total, so a high skill
# there is worth less than the same number on the abundance distribution, which carries
# 52-61%.
#
# GDM rungs are excluded: their skill is scored on pairwise Bray-Curtis rather than on
# biomass, and mean dissimilarity itself changes with the number of taxa (0.461 at 57
# families, 0.300 at 8 classes), so the two are not on a common footing. GDM results are
# in data/gdm_dbo_ladder_results.rds and in the figures from scripts/plot_gdm_dbo.R.
#
# Output: figures/dbo_ladder.{png,pdf}

suppressMessages(library(tidyverse))

RUNG_LAB <- c(family = "Family (57)", class = "Class (7)", total = "Total (1)")
HOLD_LEV <- c("LOSO", "LORO", "Train/test")
MOD_LAB  <- c(static = "Static", dynamic = "Dynamic", svc = "Decomposed + Lat/Long")
COMP_LAB <- c(abundance = "Abundance\ndistribution", global = "Global\ntrajectory",
              regional  = "Regional\ntrajectory",    station = "Station\ntrajectory")
# one hue, light to dark: the rungs are the same measurement at three grains
COLS <- c("Family (57)" = "#9ECAE1", "Class (7)" = "#4292C6", "Total (1)" = "#08519C")

s <- readRDS("data/rf_dbo_ladder_skill.rds")

# share of observed variance, computed within rung x holdout x model over the four
# components (the "pooled" row is the undecomposed total, not a fifth component)
d <- s %>% filter(component != "pooled") %>%
  group_by(rung, holdout, model) %>% mutate(share = SS_obs / sum(SS_obs)) %>% ungroup() %>%
  transmute(component = factor(COMP_LAB[component], COMP_LAB),
            rung      = factor(RUNG_LAB[rung], RUNG_LAB),
            holdout   = factor(holdout, HOLD_LEV),
            model     = factor(MOD_LAB[model], MOD_LAB),
            skill, share)
stopifnot(!any(is.na(d$component)), !any(is.na(d$rung)), !any(is.na(d$model)))

FLOOR <- -0.5                       # axis cut; estimates below are shown in the gutter
GUT   <- -0.72                      # bottom of the panel, leaving room for the labels
HIDE  <- 99                         # parked above the view, so coord_cartesian clips it

# EVERY layer must carry ALL THREE rungs at EVERY component, or the dodge breaks.
# position_dodge allocates slots from the groups present in that layer's data at that x
# position: filtering the off-scale estimates out of the point layer left two groups
# where there had been three, and the survivors were re-spaced across the full width,
# landing off their own stems. So no layer is subset -- rows that should not show are
# parked at y = HIDE and clipped by coord_cartesian, and preserve = "single" fixes the
# slot width independently of how many elements are actually visible.
dodge <- position_dodge(width = 0.68, preserve = "single")

d <- d %>% mutate(
  off      = skill < FLOOR,
  y_pt     = ifelse(off, HIDE, skill),          # on-scale points only
  y_tri    = ifelse(off, FLOOR, HIDE),
  y_lab    = ifelse(off, FLOOR - 0.055, HIDE),
  lab      = ifelse(off, sprintf("%.2f", skill), ""))

p <- ggplot(d, aes(component, skill, colour = rung, group = rung)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = GUT, ymax = FLOOR,
           fill = "grey94", colour = NA) +
  geom_hline(yintercept = FLOOR, colour = "grey70", linewidth = 0.3, linetype = "13") +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_linerange(aes(ymin = 0, ymax = pmax(skill, FLOOR)),
                 position = dodge, linewidth = 0.55) +
  geom_point(aes(y = y_pt, size = share), position = dodge) +
  # off-scale: triangle at the cut, exact value printed sideways in the gutter
  geom_point(aes(y = y_tri), shape = 17, size = 1.5, position = dodge,
             show.legend = FALSE) +
  geom_text(aes(y = y_lab, label = lab), position = dodge, angle = 90, hjust = 1,
            size = 1.95, show.legend = FALSE) +
  facet_grid(model ~ holdout) +
  scale_colour_manual(values = COLS, name = "Response") +
  scale_size_area(max_size = 3.4, name = "Share of observed variance",
                  breaks = c(0.05, 0.25, 0.5), labels = scales::percent_format(1)) +
  # literal breaks, not seq(): seq() lands on 2.2e-16 instead of 0 and ggplot then
  # labels the whole axis in scientific notation
  scale_y_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8),
                     labels = function(x) sprintf("%.1f", x)) +
  coord_cartesian(ylim = c(GUT, 0.92)) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2.6)),
         size = guide_legend(order = 2)) +
  labs(x = NULL, y = "Skill  (1 - SSE / SS)",
       title = "Prediction skill by component, model type and holdout",
       subtitle = paste("Response always spatiotemporal; only the environmental",
                        "representation differs. Skill 0 = no better than predicting no",
                        "variation along that axis.\nThe three response definitions are",
                        "shown side by side, not joined: aggregation has no consistent",
                        "effect. Triangles mark estimates\nbelow the axis cut at -0.5;",
                        "their values are printed in the grey band and are not",
                        "rescaled.")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x = element_text(size = 7.5, lineheight = 0.9),
        legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(size = 8.5), legend.text = element_text(size = 8))

ggsave("figures/dbo_ladder.png", p, width = 200, height = 175, units = "mm", dpi = 300)
ggsave("figures/dbo_ladder.pdf", p, width = 200, height = 175, units = "mm")
cat("saved -> figures/dbo_ladder.png\n\n")

print(as.data.frame(d %>% select(component, rung, holdout, model, skill) %>%
        pivot_wider(names_from = rung, values_from = skill) %>%
        arrange(holdout, model, component) %>%
        mutate(across(all_of(unname(RUNG_LAB)), ~ round(.x, 3)))), row.names = FALSE)

cat("\nshare of observed variance (%):\n")
print(as.data.frame(d %>% filter(model == MOD_LAB[["svc"]]) %>%
        transmute(holdout, rung, component, share = round(100 * share, 1)) %>%
        pivot_wider(names_from = component, values_from = share)), row.names = FALSE)
