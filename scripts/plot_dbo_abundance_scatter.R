# Observed against predicted station means -- the abundance distribution component, in the
# scatter form that dbo_main_loso.png summarises as skill bars.
#
# TOTAL BIOMASS, DECOMPOSED + Lat/Long, LOSO. One point per station: the station's mean
# fourth-root biomass over its sampled years against the mean of its held-out predictions.
# That IS the abundance distribution component -- the decomposition takes the station effect
# as fitted(lm(v ~ station)), which is the station mean, centred on the grand mean. Plotting
# the uncentred means shifts both axes by the same constant and leaves the picture unchanged.
#
# The two relationships the main figure separates are both visible here:
#   across regions   the spread BETWEEN the three colours, scored as ab_across
#   within a region  the spread WITHIN each colour, scored as ab_1, ab_2, ab_3
# A heavy grey line fits all 16 stations; the coloured lines fit each region on its own.
# Skill is 1 - SSE/SS on the corresponding projection, matching plot_dbo_main.R exactly, so
# the numbers annotated here are the ones plotted as bars there.
#
# Reads data/rf_ladder_bu_preds.rds. Nothing refitted.
#
# Output: figures/dbo_abundance_scatter.{png,pdf}

suppressMessages(library(tidyverse))

COLR <- c("DBO 1" = "#4E79A7", "DBO 2" = "#F28E2B", "DBO 3" = "#59A14F")

p <- readRDS("data/rf_ladder_bu_preds.rds") %>%
  filter(rung == "total", holdout == "LOSO", model == "svc")
stopifnot(nrow(p) == 226)

# the component skills, computed the way plot_dbo_main.R computes them
regs <- sort(unique(p$DBOreg))
f_reg <- fitted(lm(p$y ~ factor(p$DBOreg)))
f_stn <- fitted(lm(p$y ~ factor(p$DBOreg) + factor(p$stn)))
e <- p$y - p$p
e_reg <- fitted(lm(e ~ factor(p$DBOreg)))
e_stn <- fitted(lm(e ~ factor(p$DBOreg) + factor(p$stn)))
sk_across <- 1 - sum((e_reg - mean(e))^2) / sum((f_reg - mean(p$y))^2)
sk_within <- map_dbl(regs, function(r) {
  k <- p$DBOreg == r
  1 - sum((e_stn - e_reg)[k]^2) / sum((f_stn - f_reg)[k]^2)
})
names(sk_within) <- paste("DBO", regs)

st <- p %>% group_by(stn, DBOreg) %>%
  summarise(obs = mean(y), pred = mean(p), n = n(), .groups = "drop") %>%
  mutate(region = paste("DBO", DBOreg))

# Six stations sit in a tight cluster around (1.92-2.06, 2.00-2.15), right where the fitted
# line passes, and ggrepel re-solves the layout whenever the panel geometry changes -- so
# left to itself a label lands on the line, and fixing one moves another onto it. Each of
# the six gets an explicit direction instead, which also makes the figure reproducible.
# The other ten are unconstrained and repulsion places them.
NUDGE <- tribble(
  ~stn,     ~nx,    ~ny,
  "SLIP3", -0.100,  0.075,
  "UTBS2",  0.025,  0.085,
  "SLIP2", -0.135, -0.005,
  "SLIP1", -0.055, -0.080,
  "SLIP4",  0.075, -0.055,
  "SLIP5",  0.105,  0.015)
st <- st %>% left_join(NUDGE, by = "stn") %>%
  mutate(nx = replace_na(nx, 0), ny = replace_na(ny, 0))

rng <- range(c(st$obs, st$pred))
g <- ggplot(st, aes(obs, pred, colour = region)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey30",
              linewidth = 1.3, linetype = "22", alpha = 0.5, aes(group = 1)) +
  geom_point(size = 4.2, alpha = 0.75) +
  ggrepel::geom_text_repel(aes(label = stn), size = 2.8, colour = "grey35",
                           nudge_x = st$nx, nudge_y = st$ny,
                           max.overlaps = Inf, seed = 1, show.legend = FALSE,
                           point.padding = 0.5, box.padding = 0.6,
                           min.segment.length = 0.2, segment.colour = "grey65",
                           segment.size = 0.3) +
  scale_colour_manual(values = COLR, name = NULL) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = "Observed station mean, fourth-root total biomass",
       y = "Predicted station mean") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        axis.text = element_text(size = 11), axis.title = element_text(size = 13),
        legend.text = element_text(size = 10))

ggsave("figures/dbo_abundance_scatter.png", g, width = 165, height = 175,
       units = "mm", dpi = 300)
ggsave("figures/dbo_abundance_scatter.pdf", g, width = 165, height = 175, units = "mm")
cat("saved -> figures/dbo_abundance_scatter.png\n\n")

cat(sprintf("across regions   skill %+.3f\n", sk_across))
for (nm in names(sk_within)) cat(sprintf("within %-8s  skill %+.3f\n", nm, sk_within[[nm]]))
cat("\n")
print(as.data.frame(st %>% arrange(region, obs) %>%
  transmute(station = stn, region, n_years = n,
            observed = round(obs, 3), predicted = round(pred, 3),
            error = round(pred - obs, 3))), row.names = FALSE)
