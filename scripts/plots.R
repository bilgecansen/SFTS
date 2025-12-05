library(tidyverse)
library(patchwork)

plots_sim <- readRDS("figures/plots_sim.rds")
plots_bbs <- readRDS("figures/plots_bbs.rds")
plots_dbo <- readRDS("figures/plots_dbo.rds")

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

# Figure 3
th <- theme(
  axis.title = element_blank(),
  plot.title = element_blank(),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 16),
  legend.position = "none"
)

g1 <- plots_sim$gm_uneven$pos + th
g1$layers[[2]]$aes_params$colour <- "#792C1F"
g1$layers[[2]]$aes_params$alpha <- 0.7
g1$layers[[1]]$aes_params$colour <- "#485B7C"
g1$layers[[1]]$aes_params$alpha <- 0.8

g2 <- plots_sim$g_pred_uneven$pos + th
g2$layers[[2]]$aes_params$colour <- "#792C1F"
g2$layers[[2]]$aes_params$alpha <- 0.7
g2$layers[[1]]$aes_params$colour <- "#485B7C"
g2$layers[[1]]$aes_params$alpha <- 0.8

g3 <- plots_sim$gm_uneven$neut + th
g3$layers[[2]]$aes_params$colour <- "#792C1F"
g3$layers[[2]]$aes_params$alpha <- 0.7
g3$layers[[1]]$aes_params$colour <- "#485B7C"
g3$layers[[1]]$aes_params$alpha <- 0.8

g4 <- plots_sim$g_pred_uneven$neut + th
g4$layers[[2]]$aes_params$colour <- "#792C1F"
g4$layers[[2]]$aes_params$alpha <- 0.7
g4$layers[[1]]$aes_params$colour <- "#485B7C"
g4$layers[[1]]$aes_params$alpha <- 0.8

g5 <- plots_sim$gm_uneven$neg + th
g5$layers[[2]]$aes_params$colour <- "#792C1F"
g5$layers[[2]]$aes_params$alpha <- 0.7
g5$layers[[1]]$aes_params$colour <- "#485B7C"
g5$layers[[1]]$aes_params$alpha <- 0.8

g6 <- plots_sim$g_pred_uneven$neg + th
g6$layers[[2]]$aes_params$colour <- "#792C1F"
g6$layers[[2]]$aes_params$alpha <- 0.7
g6$layers[[1]]$aes_params$colour <- "#485B7C"
g6$layers[[1]]$aes_params$alpha <- 0.8

(g1 + g2 + g3 + g4 + g5 + g6) +
  plot_layout(nrow = 3, ncol = 2)

ggsave("figures/fig3.pdf", width = 220, height = 220, units = "mm", dpi = 600)
ggsave("figures/fig3.jpeg", width = 220, height = 220, units = "mm", dpi = 600)

# Figure 4
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

ggsave("figures/fig4.pdf", width = 220, height = 220, units = "mm", dpi = 600)
ggsave("figures/fig4.jpeg", width = 220, height = 220, units = "mm", dpi = 600)

# Figure 5
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

ggsave("figures/fig5.pdf", width = 220, height = 220, units = "mm", dpi = 600)
ggsave("figures/fig5.png", width = 220, height = 220, units = "mm", dpi = 600)
