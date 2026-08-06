# All GDM figures for the DBO community-turnover analysis.
# Reads data/gdm_dbo_preds.rds (written by scripts/fit_gdm_dbo.R) and produces the three
# canonical figures:
#
#   1. figures/gdm_dbo_region_loso.{png,pdf}       spatial turnover by region pair, LOSO
#   2. figures/gdm_dbo_region_traintest.{png,pdf}  spatial turnover by region pair, train/test
#   3. figures/gdm_dbo_temporal.{png,pdf}          temporal turnover, one panel per holdout
#
# Why the spatial figures are split by region pair. The pooled spatial scatter is bimodal:
# within-region pairs sit near Bray-Curtis 0.32, cross-region pairs near 0.49, and that
# split alone carries much of the pooled correlation -- centring within region-pair group
# drops r from 0.754 to 0.500 (LOSO) and 0.755 to 0.482 (train/test). Conditioning on the
# group shows what fine-grained ability remains: essentially all of it is within DBO 3, the
# best-sampled band, while every cross-region group has near-zero or negative skill.
#
# The temporal figure is NOT split by region: temporal pairs are same-station by
# construction, so they always fall inside one band and no cross-region pairs exist.
#
# All three figures share ONE axis range, computed over both models' predictions, so they
# read as a set. That is deliberate -- it shows temporal dissimilarities are both SMALLER
# than spatial ones and less predictable, rather than rescaling and hiding the first of those.
#
# In every panel: faint line is 1:1 (the skill reference), heavy translucent line is the
# pooled fit within that panel. Temporal points sit systematically ABOVE the 1:1 -- observed
# turnover exceeds predicted -- which is the sampling floor: two grabs from one station in
# different years differ for reasons no environmental model can anticipate. Skill is
# therefore scored against the mean observed dissimilarity, asking whether the model
# explains VARIATION in turnover rather than its absolute level.

suppressMessages(library(tidyverse))

glev <- c("Within DBO 1", "Within DBO 2", "Within DBO 3",
          "DBO 1 - DBO 2", "DBO 1 - DBO 3", "DBO 2 - DBO 3")
hlev <- c("LOSO", "Train/test")
kcols <- c("Within region" = "#DD9100", "Between regions" = "#254F9A")
TEMP_COL <- "#254F9A"

allp <- readRDS("data/gdm_dbo_preds.rds") %>% mutate(holdout = factor(holdout, hlev))
rng <- range(c(allp$obs, allp$pred))          # shared by all three figures

skill_of <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)

base_layers <- function() list(
  geom_abline(slope = 1, intercept = 0,
              colour = scales::alpha("black", 0.25), linewidth = 0.4),
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x,
              colour = scales::alpha("black", 0.55), linewidth = 0.9),
  coord_cartesian(xlim = rng, ylim = rng),
  labs(x = "Predicted dissimilarity (Bray-Curtis)",
       y = "Observed dissimilarity (Bray-Curtis)"),
  theme_bw(base_size = 10),
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        strip.text = element_text(size = 9.5),
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        aspect.ratio = 1))

lab_layer <- function(ann, size = 2.8) list(
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.12, vjust = 1.12, size = size, fontface = 2, lineheight = 0.95))

# ---- figures 1 & 2: spatial turnover by region pair ----------------------------------
spat <- allp %>% filter(model == "Spatial turnover") %>%
  mutate(regpair = factor(regpair, glev))

region_fig <- function(hold, stem, title) {
  d <- filter(spat, holdout == hold)
  ann <- d %>% group_by(regpair, kind) %>%
    summarise(n = n(), r = cor(obs, pred), skill = skill_of(obs, pred), .groups = "drop") %>%
    mutate(label = sprintf("r %+.2f\nskill %+.2f\nn = %d", r, skill, n))
  cw <- d %>% group_by(regpair) %>%
    mutate(o = obs - mean(obs), p = pred - mean(pred)) %>% ungroup()
  p <- ggplot(d, aes(pred, obs)) + base_layers()[1] +
    geom_point(aes(colour = kind), size = 0.8, alpha = 0.3, show.legend = FALSE) +
    base_layers()[-1] + lab_layer(ann) +
    scale_colour_manual(values = kcols) +
    facet_wrap(~ regpair, nrow = 2) +
    labs(title = title,
         subtitle = sprintf("pooled across all pairs: r %+.2f, skill %+.2f  |  within-group centred: r %+.2f",
                            cor(d$obs, d$pred), skill_of(d$obs, d$pred), cor(cw$o, cw$p)))
  ggsave(sprintf("figures/%s.png", stem), p, width = 190, height = 150, units = "mm", dpi = 300)
  ggsave(sprintf("figures/%s.pdf", stem), p, width = 190, height = 150, units = "mm", dpi = 300)
  cat(sprintf("saved -> figures/%s.png\n", stem))
  ann %>% transmute(figure = stem, regpair, n, r = round(r, 3), skill = round(skill, 3))
}

out <- bind_rows(
  region_fig("LOSO",       "gdm_dbo_region_loso",
             "GDM spatial turnover by region pair - LOSO"),
  region_fig("Train/test", "gdm_dbo_region_traintest",
             "GDM spatial turnover by region pair - train/test"))

# ---- figure 3: temporal turnover, one panel per holdout -------------------------------
temp <- allp %>% filter(model == "Temporal turnover")
ann_t <- temp %>% group_by(holdout) %>%
  summarise(n = n(), r = cor(obs, pred), skill = skill_of(obs, pred), .groups = "drop") %>%
  mutate(label = sprintf("r %+.2f\nskill %+.2f\nn = %s", r, skill, format(n, big.mark = ",")))

p3 <- ggplot(temp, aes(pred, obs)) + base_layers()[1] +
  geom_point(colour = TEMP_COL, size = 0.8, alpha = 0.28) +
  base_layers()[-1] + lab_layer(ann_t, size = 3.0) +
  facet_wrap(~ holdout, nrow = 1) +
  labs(title = "GDM temporal turnover (same station, different year)")
ggsave("figures/gdm_dbo_temporal.png", p3, width = 150, height = 92, units = "mm", dpi = 300)
ggsave("figures/gdm_dbo_temporal.pdf", p3, width = 150, height = 92, units = "mm", dpi = 300)
cat("saved -> figures/gdm_dbo_temporal.png\n\n")

print(as.data.frame(out), row.names = FALSE)
cat("\n")
print(as.data.frame(ann_t %>% transmute(figure = "gdm_dbo_temporal", holdout, n,
      r = round(r, 3), skill = round(skill, 3))), row.names = FALSE)
