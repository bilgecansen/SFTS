library(tidyverse)
library(patchwork)

plots_sim <- readRDS("plots_sim.rds")
plots_bbs <- readRDS("plots_bbs.rds")
plots_dbo <- readRDS("plots_dbo.rds")

theme_set(theme_bw())

# Figure 1
g1 <- plots_sim$g_imp_uneven$`spatio-temporal` +
  theme(plot.title = element_text(size = 10))

g2 <- plots_sim$g_imp_uneven_multi$`spatio-temporal` +
  scale_y_continuous(breaks = 1:12) +
  theme(plot.title = element_text(size = 10), legend.position = "none")

g3 <- plots_sim$g_pred_uneven$`spatio-temporal` +
  theme(plot.title = element_text(size = 10))

g4 <- plots_sim$g_pred_uneven_multi$`spatio-temporal` +
  theme(plot.title = element_text(size = 10), legend.position = "none")

(g1 + g2) / free(g3 + g4) + plot_annotation(tag_levels = "a")
ggsave("figures/fig1.pdf", width = 180, height = 180, units = "mm", dpi = 600)
ggsave("figures/fig1.jpeg", width = 180, height = 180, units = "mm", dpi = 600)

# Figure 2
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

h2 <- plots_sim$g_imp_uneven_multi$`spatial + temporal-opposite` +
  theme(
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(
    values = c(
      "Residual" = "#907047",
      "Spatial" = "#792C1F",
      "Temporal" = "#485B7C"
    )
  )

plots_bbs$g_imp$layers[[2]]$aes_params$alpha <- 0.3
h3 <- plots_bbs$g_imp +
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

h4 <- plots_dbo$g1 +
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

h1 +
  h2 +
  h3 +
  h4 +
  plot_annotation(tag_levels = "a") +
  plot_layout(ncol = 2, nrow = 2)

ggsave("figures/fig2.pdf", width = 220, height = 180, units = "mm", dpi = 600)
ggsave("figures/fig2.jpeg", width = 220, height = 180, units = "mm", dpi = 600)

# Figure 3
j1 <- plots_sim$g_pred_uneven$spatial +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    #axis.title.x = element_blank(),
    legend.position = "none"
  )
j1$layers[[1]]$aes_params$colour <- "#792C1F"
j1$layers[[1]]$aes_params$alpha <- 0.4
j1$layers[[2]]$aes_params$colour <- "#485B7C"
j1$layers[[2]]$aes_params$alpha <- 0.6

j2 <- plots_bbs$g1 +
  theme(
    plot.title = element_text(size = 10, hjust = 0),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    #axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    legend.position = "none"
  )
j2$layers[[1]]$aes_params$colour <- "#792C1F"
j2$layers[[1]]$aes_params$alpha <- 0.4
j2$layers[[2]]$aes_params$colour <- "#485B7C"
j2$layers[[2]]$aes_params$alpha <- 0.4

plots_sim$g_pred_uneven_multi$spatial$layers[[1]]$aes_params$size <- 2
plots_sim$g_pred_uneven_multi$spatial$layers[[1]]$aes_params$alpha <- 0.2
d1 <- plots_sim$g_pred_uneven_multi$spatial +
  theme(
    plot.title = element_text(size = 10),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    #axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(-0.3, 0.9),
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.9)
  ) +
  scale_color_manual(
    values = c(
      "Spatio-temporal" = "#792C1F",
      "Temporal" = "#485B7C"
    )
  )

plots_bbs$g_pred$layers[[1]]$aes_params$alpha <- 0.4
d2 <- plots_bbs$g_pred +
  theme(
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    #axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  labs(title = "Breeding Bird Survey") +
  scale_y_continuous(
    limits = c(-0.65, 0.85),
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.9)
  ) +
  scale_color_manual(
    values = c(
      "Spatio-temporal" = "#792C1F",
      "Temporal" = "#485B7C"
    )
  )

j1 +
  d1 +
  j2 +
  d2 +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 2, ncol = 2)

ggsave("figures/fig3.pdf", width = 220, height = 180, units = "mm", dpi = 600)

# Figure 4
j4 <- plots_dbo$g2 +
  labs(title = "Naididae") +
  theme(
    plot.title = element_text(size = 10, hjust = 0),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    #axis.title.y = element_blank(),
    legend.position = "none"
  )
j4$layers[[1]]$aes_params$colour <- "#792C1F"
j4$layers[[2]]$aes_params$colour <- "#485B7C"

j5 <- plots_dbo$g3 +
  labs(title = "Pleustidae") +
  theme(
    plot.title = element_text(size = 10, hjust = 0),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    #axis.title.y = element_blank(),
    legend.position = "none"
  )
j5$layers[[1]]$aes_params$colour <- "#792C1F"
j5$layers[[2]]$aes_params$colour <- "#485B7C"

j6 <- plots_dbo$g4 +
  labs(title = "Tellinidae") +
  theme(
    plot.title = element_text(size = 10, hjust = 0),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    #axis.title.y = element_blank(),
    legend.position = "none"
  )
j6$layers[[1]]$aes_params$colour <- "#792C1F"
j6$layers[[2]]$aes_params$colour <- "#485B7C"

plots_dbo$g5$layers[[1]]$aes_params$alpha <- 0.5
d4 <- plots_dbo$g5 +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    #axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  labs(title = "Distributed Biological Observatory") +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
  ) +
  scale_color_manual(
    values = c(
      "Spatio-temporal" = "#792C1F",
      "Temporal" = "#485B7C"
    )
  )

d4 +
  j6 +
  j4 +
  j5 +
  plot_annotation(tag_levels = "a") +
  plot_layout(ncol = 2, nrow = 2)

ggsave("figures/fig4.pdf", width = 220, height = 180, units = "mm", dpi = 600)
