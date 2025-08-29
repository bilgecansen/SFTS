library(tidyverse)
library(patchwork)

plots_sim <- readRDS("plots_sim.rds")
plots_bbs <- readRDS("plots_bbs.rds")
plots_dbo <- readRDS("plots_dbo.rds")

theme_set(theme_bw())

# Include only species_wide r > 0.5
filter_cor <- function(p, l) {
  z <- filter(p$data, type == "Species-wide")
  idx <- which(z$r >= 0.5)
  idx <- c(idx, idx + l)
  p$data <- p$data[idx, ]

  p
}

plots_bbs$g_pred1 <- filter_cor(plots_bbs$g_pred, l = 293)
plots_bbs$g_pred2 <- filter_cor(plots_bbs$g_pred2, l = 293)
plots_bbs$g_pred3 <- filter_cor(plots_bbs$g_pred3, l = 293)

plots_dbo$g_pred1 <- filter_cor(plots_dbo$g_pred1, l = 61)
plots_dbo$g_pred2 <- filter_cor(plots_dbo$g_pred2, l = 61)
plots_dbo$g_pred3 <- filter_cor(plots_dbo$g_pred3, l = 61)

# Figure 2
g0a <- plots_sim$gm_uneven$spatial +
  theme(
    axis.title = element_text(size = 14),
    element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(1050, 1250),
    breaks = c(1100, 1150, 1200, 1250)
  )

k <- g0a$layers[[1]]
g0a$layers[[1]] <- g0a$layers[[2]]
g0a$layers[[2]] <- k

g0a$layers[[1]]$aes_params$colour <- "#792C1F"
g0a$layers[[1]]$aes_params$alpha <- 0.8

g0a$layers[[2]]$aes_params$colour <- "#485B7C"
g0a$layers[[2]]$aes_params$alpha <- 0.9
g0a$layers[[2]]$aes_params$linewidth <- 1.5

g0b <- plots_sim$g_pred_uneven$spatial +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(1050, 1250),
    breaks = c(1100, 1150, 1200, 1250)
  )
g0b$layers[[2]]$aes_params$colour <- "#792C1F"
g0b$layers[[2]]$aes_params$alpha <- 1
g0b$layers[[2]]$aes_params$size <- 3
g0b$layers[[1]]$aes_params$colour <- "#485B7C"
g0b$layers[[1]]$aes_params$alpha <- 0.8

g0c <- plots_sim$g_pred_uneven_multi$spatial +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(-0.3, 1),
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.9)
  ) +
  scale_color_manual(
    values = c(
      "Species-wide" = "#792C1F",
      "Population-level" = "#485B7C"
    )
  )
g0c$layers[[2]]$aes_params$size <- 2
g0c$layers[[2]]$aes_params$alpha <- 0.5

layout <- "
AA##
BBCC
"

g0a +
  g0b +
  g0c +
  #plot_annotation(
  #tag_levels = "a",
  #title = "Simulation: Spatial (+), Temporal (0)"
  #) +
  plot_layout(design = layout)

ggsave("figures/fig2.pdf", width = 220, height = 180, units = "mm", dpi = 600)
ggsave("figures/fig2.jpeg", width = 220, height = 180, units = "mm", dpi = 600)

# Figure 3
th <- theme(
  axis.title = element_blank(),
  plot.title = element_blank(),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 16),
  legend.position = "none"
)

j1 <- plots_bbs$g1 + th
j1$layers[[2]]$aes_params$colour <- "#792C1F"
j1$layers[[2]]$aes_params$alpha <- 0.7
j1$layers[[1]]$aes_params$colour <- "#485B7C"
j1$layers[[1]]$aes_params$alpha <- 0.8

j2 <- plots_dbo$g1 + th
j2$layers[[2]]$aes_params$colour <- "#792C1F"
j2$layers[[2]]$aes_params$alpha <- 0.7
j2$layers[[2]]$aes_params$size <- 2
j2$layers[[1]]$aes_params$colour <- "#485B7C"
j2$layers[[1]]$aes_params$alpha <- 0.8

j3 <- plots_bbs$g2 + th
j3$layers[[2]]$aes_params$colour <- "#792C1F"
j3$layers[[2]]$aes_params$alpha <- 0.7
j3$layers[[1]]$aes_params$colour <- "#485B7C"
j3$layers[[1]]$aes_params$alpha <- 0.8

j4 <- plots_dbo$g2 + th
j4$layers[[2]]$aes_params$colour <- "#792C1F"
j4$layers[[2]]$aes_params$alpha <- 0.7
j4$layers[[1]]$aes_params$colour <- "#485B7C"
j4$layers[[1]]$aes_params$alpha <- 0.8

j5 <- plots_bbs$g3 + th
j5$layers[[2]]$aes_params$colour <- "#792C1F"
j5$layers[[2]]$aes_params$alpha <- 0.7
j5$layers[[1]]$aes_params$colour <- "#485B7C"
j5$layers[[1]]$aes_params$alpha <- 0.8

j6 <- plots_dbo$g3 + th
j6$layers[[2]]$aes_params$colour <- "#792C1F"
j6$layers[[2]]$aes_params$alpha <- 0.7
j6$layers[[1]]$aes_params$colour <- "#485B7C"
j6$layers[[1]]$aes_params$alpha <- 0.8

(j1 + j2 + j3 + j4 + j5 + j6) +
  plot_layout(nrow = 3, ncol = 2)

ggsave("figures/fig3.pdf", width = 220, height = 220, units = "mm", dpi = 600)
ggsave("figures/fig3.jpeg", width = 220, height = 220, units = "mm", dpi = 600)

# Figure 4
mn <- scale_color_manual(
  values = c(
    "Species-wide" = "#792C1F",
    "Population-level" = "#485B7C"
  )
)

## Static SDM
d1 <- plots_bbs$g_pred1 +
  #labs(title = "Breeding Bird Survey") +
  scale_y_continuous(
    limits = c(-0.65, 1),
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.9)
  ) +
  th +
  mn
d1$layers[[2]]$aes_params$alpha <- 0.4
d1$layers[[2]]$aes_params$size <- 2

d2 <- plots_dbo$g_pred1 +
  #labs(title = "Distributed Biological Observatory") +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
  ) +
  th +
  mn
d2$layers[[2]]$aes_params$alpha <- 0.4
d2$layers[[2]]$aes_params$size <- 2

## Dynamic SDM
d3 <- plots_bbs$g_pred2 +
  #labs(title = "Breeding Bird Survey") +
  scale_y_continuous(
    limits = c(-0.65, 1),
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.9)
  ) +
  th +
  mn
d3$layers[[2]]$aes_params$alpha <- 0.4
d3$layers[[2]]$aes_params$size <- 2

d4 <- plots_dbo$g_pred2 +
  #labs(title = "Distributued Biological Observatory") +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
  ) +
  th +
  mn
d4$layers[[2]]$aes_params$alpha <- 0.4
d4$layers[[2]]$aes_params$size <- 2

## Decomposed SDM
d5 <- plots_bbs$g_pred3 +
  #labs(title = "Breeding Bird Survey") +
  scale_y_continuous(
    limits = c(-0.65, 1),
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.9)
  ) +
  th +
  mn
d5$layers[[2]]$aes_params$alpha <- 0.4
d5$layers[[2]]$aes_params$size <- 2

d6 <- plots_dbo$g_pred3 +
  #labs(title = "Distributed Biological Observatory") +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
  ) +
  th +
  mn
d6$layers[[2]]$aes_params$alpha <- 0.4
d6$layers[[2]]$aes_params$size <- 2

(d1 + d2 + d3 + d4 + d5 + d6) +
  plot_layout(nrow = 3, ncol = 2)

ggsave("figures/fig4.pdf", width = 220, height = 220, units = "mm", dpi = 600)
ggsave("figures/fig4.png", width = 220, height = 220, units = "mm", dpi = 600)


# Supplemental figures
# Figure S1
g1 <- plots_sim$g_imp_uneven$`spatio-temporal` +
  theme(plot.title = element_text(size = 10))

g2 <- plots_sim$g_imp_uneven_multi$`spatio-temporal` +
  scale_y_continuous(breaks = 1:12) +
  theme(plot.title = element_text(size = 10), legend.position = "none") +
  scale_color_manual(
    values = c(
      "Residual" = "#907047",
      "Spatial" = "#792C1F",
      "Temporal" = "#485B7C"
    )
  )
g2$layers[[2]]$aes_params$alpha <- 0.8

g3 <- plots_sim$g_pred_uneven$`spatio-temporal` +
  theme(plot.title = element_text(size = 10))
g3$layers[[2]]$aes_params$colour <- "#792C1F"
g3$layers[[2]]$aes_params$alpha <- 1
g3$layers[[2]]$aes_params$size <- 2
g3$layers[[1]]$aes_params$colour <- "#485B7C"
g3$layers[[1]]$aes_params$alpha <- 0.8

g4 <- plots_sim$g_pred_uneven_multi$`spatio-temporal` +
  theme(plot.title = element_text(size = 10), legend.position = "none") +
  scale_color_manual(
    values = c(
      "Species-wide" = "#792C1F",
      "Population-level" = "#485B7C"
    )
  )
g4$layers[[2]]$aes_params$alpha <- 0.8

(g1 + g2) / free(g3 + g4) + plot_annotation(tag_levels = "a")

ggsave("figures/fig2.pdf", width = 220, height = 180, units = "mm", dpi = 600)
ggsave("figures/fig2.jpeg", width = 220, height = 180, units = "mm", dpi = 600)

# Fig S1
s1 <- plots_bbs$g3 +
  labs(title = "Breeding Bird Survey")
s2 <- plots_dbo$g4 +
  labs(title = "Disributed Biological Observatory")

s1 +
  s2 +
  plot_annotation(tag_levels = "a")

ggsave(
  "figures/figS1.pdf",
  width = 220,
  height = 100,
  units = "mm",
  dpi = 600
)

ggsave(
  "figures/figS1.jpeg",
  width = 220,
  height = 100,
  units = "mm",
  dpi = 600
)

# Fig S2

d3 +
  d5 +
  plot_annotation(tag_levels = "a") +
  plot_annotation(title = "Dynamic SDM")

ggsave("figures/figS2.pdf", width = 220, height = 100, units = "mm", dpi = 600)
ggsave("figures/figS2.jpeg", width = 220, height = 100, units = "mm", dpi = 600)

d4 +
  d6 +
  plot_annotation(tag_levels = "a") +
  plot_annotation(title = "Decomposed SDM")

ggsave("figures/figS3.pdf", width = 220, height = 100, units = "mm", dpi = 600)
ggsave("figures/figS3.jpeg", width = 220, height = 100, units = "mm", dpi = 600)

# Fig S?
# Figure 3
h1 <- plots_sim$g_imp_uneven_multi$spatial +
  theme(
    plot.title = element_text(size = 10),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(
    values = c(
      "Residual" = "#907047",
      "Spatial" = "#792C1F",
      "Temporal" = "#485B7C"
    )
  )
h1$layers[[2]]$aes_params$alpha <- 0.7

h2 <- plots_bbs$g_imp +
  theme(
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = c(1, 10, 20, 30, 40)) +
  labs(title = "Breeding Bird Survey") +
  scale_color_manual(
    values = c(
      "Residual" = "#907047",
      "Spatial" = "#792C1F",
      "Temporal" = "#485B7C"
    )
  )
h2$layers[[2]]$aes_params$alpha <- 0.7

h3 <- plots_dbo$g_imp +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 25)) +
  labs(title = "Distributed Biological Observatory") +
  scale_y_continuous(breaks = c(1, 10, 20, 30, 40)) +
  scale_color_manual(
    values = c(
      "Residual" = "#907047",
      "Spatial" = "#792C1F",
      "Temporal" = "#485B7C"
    )
  )
h3$layers[[2]]$aes_params$alpha <- 0.8

h1 +
  h2 +
  h3 +
  plot_annotation(tag_levels = "a") +
  plot_layout(design = layout)

ggsave("figures/fig3.pdf", width = 220, height = 180, units = "mm", dpi = 600)
ggsave("figures/fig3.jpeg", width = 220, height = 180, units = "mm", dpi = 600)
