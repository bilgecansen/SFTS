library(tidyverse)
library(ranger)
library(foreach)
library(patchwork)

theme_set(theme_bw())

bbs <- readRDS("data/data_bbs_nozero.rds")

bbs_train <- filter(bbs, year < 2011)
bbs_test <- filter(bbs, year >= 2011 & year < 2021)

# only sites with at least 5 years of observations
year_n <- bbs_train %>%
  group_by(species_id, site_id) %>%
  summarise(n = length(abundance[abundance > 0])) %>%
  filter(n > 4)

idx_ss <- paste(year_n$species_id, year_n$site_id)
bbs_train$species_site <- paste(bbs_train$species_id, bbs_train$site_id)

bbs_train2 <- filter(bbs_train, species_site %in% idx_ss)

# only species with at least 5 sites
site_n <- bbs_train2 %>%
  group_by(species_id) %>%
  summarize(n = length(unique(site_id))) %>%
  filter(n > 4)

bbs_train3 <- filter(bbs_train2, species_id %in% site_n$species_id)

species <- unique(bbs_train3$species_id)

# remove unused large objects
rm(bbs, bbs_train, bbs_train2)


# Abundance RF models on strata scale -------------------------------------

# aggregate to strata and decompose data for RF
data_rf_str <- foreach(i = 1:length(species)) %do%
  {
    filter(bbs_train3, species_id == species[i]) %>%
      dplyr::select(-site_id, -species_site) %>%
      group_by(species_id, strata, year) %>%
      summarize(across(everything(), mean)) %>%
      ungroup() %>%
      select(
        abundance,
        elevs,
        lat,
        long,
        party_hours,
        amo,
        amo_lag,
        enso,
        enso_lag,
        nao,
        nao_lag,
        pdo,
        pdo_lag,
        contains(c("spatial", "temporal", "residual")),
        -contains("lag_spatial"),
        strata
      )
  }
names(data_rf_str) <- species

## Remove species with too few strata
idx_str <- which(map_dbl(data_rf_str, function(x) length(unique(x$strata))) > 9)
data_rf_str <- data_rf_str[idx_str]

pb <- txtProgressBar(max = length(data_rf_str), style = 3)
var_imp_str <- list()
R2_str <- c()

# All variables
for (i in 1:length(data_rf_str)) {
  res_rf <- ranger(
    log(abundance) ~ .,
    data = data_rf_str[[i]],
    num.trees = 500,
    importance = "permutation"
  )

  R2_str[i] <- res_rf$r.squared
  var_imp_str[[i]] <- res_rf$variable.importance

  setTxtProgressBar(pb, i)
}

hist(R2_str)
var_imp_str_sel <- var_imp_str[which(R2_str >= 0.25)]
dat_imp <- var_imp_str_sel[[1]][order(var_imp_str_sel[[1]], decreasing = T)][
  1:40
]
names_imp <- factor(names(dat_imp)[40:1], levels = names(dat_imp)[40:1])

g_imp1 <- ggplot() +
  geom_segment(
    aes(y = names_imp, x = 0, xend = dat_imp[40:1]),
    color = "#333D79FF",
    linewidth = 1.1
  ) +
  labs(
    title = "Northern Bobwhite",
    y = "Variables",
    x = "Permutation Importance"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    #axis.title.x = element_text(margin = margin(t = 10)),
    #axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

med_imp_str1 <- foreach(i = 1:length(var_imp_str_sel), .combine = "rbind") %do%
  {
    z <- var_imp_str_sel[[i]][order(var_imp_str_sel[[i]], decreasing = T)]

    data.frame(
      best_imp = c(
        min(str_which(names(z), "spatial")),
        min(str_which(names(z), "temporal")),
        min(str_which(names(z), "residual")),
        min(str_which(names(z), "amo")),
        min(str_which(names(z), "enso")),
        min(str_which(names(z), "nao")),
        min(str_which(names(z), "pdo"))
      ),
      #median(str_which(names(z), "party_hours"))),
      type = c("spatial", "temporal", "residual", "amo", "enso", "nao", "pdo")
    ) #, "elevs",
    #"lat", "lon", "party_hours"))
  }

med_imp_str2 <- foreach(i = 1:length(var_imp_str_sel), .combine = "rbind") %do%
  {
    z <- var_imp_str_sel[[i]][order(var_imp_str_sel[[i]], decreasing = T)]

    data.frame(
      best_imp = c(
        min(str_which(names(z), "spatial")),
        min(str_which(names(z), "temporal")),
        min(str_which(names(z), "residual"))
      ),
      #median(str_which(names(z), "party_hours"))),
      type = c("Spatial", "Temporal", "Residual")
    ) #, "elevs",
    #"lat", "lon", "party_hours"))
  }

ggplot() +
  geom_boxplot(data = med_imp_str1, aes(x = type, y = best_imp))

g_imp2 <- ggplot() +
  geom_jitter(
    data = med_imp_str2,
    aes(x = type, y = best_imp),
    color = "#333D79FF",
    alpha = 0.5
  ) +
  labs(y = "Rank of Top Variable", x = "Variable Composition Type") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 10)
  )

g_imp1 + g_imp2 + plot_annotation(tag_levels = 'a')

ggsave("fig2.pdf", width = 180, height = 150, units = "mm", dpi = 600)


# Predict test set with space only models ---------------------------------

bbs_test$species_site <- paste(bbs_test$species_id, bbs_test$site_id)

#data_rf_base <- data_rf_str[which(R2_str >= 0.25)] %>%
#map(function(x) select(x, abundance, elevs, lat, long, party_hours))

#data_rf_sp <- data_rf_str[which(R2_str >= 0.25)]
data_rf_sp <- map(
  data_rf_str,
  function(x)
    select(
      x,
      abundance,
      elevs,
      #lat,
      #long,
      #strata,
      party_hours,
      contains("spatial")
    )
)

# spatial model
pb <- txtProgressBar(max = length(data_rf_sp), style = 3)
res_rf_sp <- list()
R2_sp <- c()

for (i in 1:length(data_rf_sp)) {
  res_rf_sp[[i]] <- ranger(
    log(abundance) ~ .,
    data = data_rf_sp[[i]],
    num.trees = 500
  )

  R2_sp[i] <- res_rf_sp[[i]]$r.squared

  setTxtProgressBar(pb, i)
}

hist(R2_sp)

# base model
#pb <- txtProgressBar(max = length(data_rf_base), style = 3)
#res_rf_base <- list()
#R2_base <- c()

#for (i in 1:length(data_rf_base)) {
#res_rf_base[[i]] <- ranger(
#log(abundance) ~ .,
#data = data_rf_base[[i]],
#num.trees = 500
#)

#R2_base[i] <- res_rf_base[[i]]$r.squared

#setTxtProgressBar(pb, i)
#}

#hist(R2_base)

# Predictions with SFTS
species_sp <- species[idx_str]

data_pred_r <- foreach(i = 1:length(species_sp)) %do%
  {
    dat <- filter(bbs_test, species_id == as.numeric(species_sp[i])) %>%
      dplyr::select(-site_id, -species_site) %>%
      group_by(species_id, strata, year) %>%
      summarize(across(everything(), mean)) %>%
      ungroup()

    sites <- unique(dat$strata)

    dat_pred <- foreach(h = 1:length(sites), .combine = "rbind") %do%
      {
        dat_site <- filter(dat, strata == sites[h]) %>%
          select(elevs, party_hours, contains("bio")) %>%
          select(-contains("lag")) %>%
          rename_with(~ paste(., "spatial", sep = "_"), contains("bio")) %>%
          as.data.frame()

        #dat_site_base <- select(dat_site, -contains("bio"))

        if (nrow(dat_site) < 5) return(NA)

        y_pred_sp <- predict(res_rf_sp[[i]], data = dat_site)
        y_pred_sp <- (y_pred_sp$predictions -
          min(log(data_rf_sp[[i]]$abundance))) /
          (max(log(data_rf_sp[[i]]$abundance)) -
            min(log(data_rf_sp[[i]]$abundance)))
        #y_pred_base <- y_pred_base$predictions - min(data_rf_base[[i]]$abundance)/
        #max(data_rf_base[[i]]$abundance) - min(data_rf_base[[i]]$abundance)
        y <- (log(filter(dat, strata == sites[h])$abundance) -
          min(log(data_rf_sp[[i]]$abundance))) /
          (max(log(data_rf_sp[[i]]$abundance)) -
            min(log(data_rf_sp[[i]]$abundance)))

        data.frame(
          y_pred_sp = y_pred_sp,
          #y_pred_base = y_pred_base,
          y = y,
          species_id = species_sp[i],
          strata = sites[h]
        )
      }

    dat_pred[!is.na(dat_pred$y), ]
  }

pred_r <- map_dbl(data_pred_r, function(x) cor(x$y_pred_sp, x$y))
idx_pred <- which(pred_r > 0.5)
data_pred_r <- data_pred_r[idx_pred]
pred_r <- pred_r[idx_pred]
species_sp_sel <- species_sp[idx_pred]

order(pred_r, decreasing = T)[1:10]

g1 <- ggplot() +
  geom_point(
    aes(
      x = data_pred_r[[32]]$y_pred_sp,
      y = data_pred_r[[32]]$y,
      group = data_pred_r[[32]]$strata
    ),
    color = "pink3",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[32]]$y_pred_sp,
      y = data_pred_r[[32]]$y,
      group = data_pred_r[[32]]$strata
    ),
    linewidth = 1.25,
    color = "#333D79FF",
    alpha = 0.6,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "American Robin",
    y = "Observed Abudance (log)",
    x = "Predicted Abundance (log)"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

g2 <- ggplot() +
  geom_point(
    aes(
      x = data_pred_r[[52]]$y_pred_sp,
      y = data_pred_r[[52]]$y,
      group = data_pred_r[[52]]$strata
    ),
    color = "pink3",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[52]]$y_pred_sp,
      y = data_pred_r[[52]]$y,
      group = data_pred_r[[52]]$strata
    ),
    linewidth = 1.25,
    color = "#333D79FF",
    alpha = 0.6,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Western Meadowlark",
    y = "Observed Abudance (log)",
    x = "Predicted Abundance (log)"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

g3 <- ggplot() +
  geom_point(
    aes(
      x = data_pred_r[[26]]$y_pred_sp,
      y = data_pred_r[[26]]$y,
      group = data_pred_r[[26]]$strata
    ),
    color = "pink3",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[26]]$y_pred_sp,
      y = data_pred_r[[26]]$y,
      group = data_pred_r[[26]]$strata
    ),
    linewidth = 1.25,
    color = "#333D79FF",
    alpha = 0.6,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Northern Mockingbird",
    y = "Observed Abudance (log)",
    x = "Predicted Abundance (log)"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

g4 <- ggplot() +
  geom_point(
    aes(
      x = data_pred_r[[17]]$y_pred_sp,
      y = data_pred_r[[17]]$y,
      group = data_pred_r[[17]]$strata
    ),
    color = "pink3",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[17]]$y_pred_sp,
      y = data_pred_r[[17]]$y,
      group = data_pred_r[[17]]$strata
    ),
    linewidth = 1.25,
    color = "#333D79FF",
    alpha = 0.6,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Northern Cardinal",
    y = "Observed Abudance (log)",
    x = "Predicted Abundance (log)"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

g1 +
  g2 +
  g3 +
  g4 +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag.position = c(0.07, 1))

ggsave("fig3.pdf", width = 180, height = 140, units = "mm", dpi = 600)

# Correlation comparison
data_pred_r2 <- do.call(rbind, data_pred_r)
data_pred_r2$gr <- paste(
  data_pred_r2$species_id,
  data_pred_r2$strata,
  sep = "-"
)

cor_temp <- data_pred_r2 %>%
  group_by(species_id, strata) %>%
  summarise(r = cor(y, y_pred_sp)) %>%
  ungroup() %>%
  group_by(species_id) %>%
  summarise(r = median(r[!is.na(r)]))

cor_all <- data_pred_r2 %>%
  group_by(species_id) %>%
  summarise(r = cor(y, y_pred_sp))

ggplot() +
  geom_point(
    aes(x = cor_all$r, y = cor_temp$r),
    color = "#333D79FF",
    alpha = 0.75,
    size = 2
  ) +
  labs(y = "Median Temporal Correlation", x = "Spatiotemporal Correlation") +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    #axis.title.x = element_text(margin = margin(t = 10)),
    #axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 10)
  )

ggsave("fig4.pdf", width = 140, height = 100, units = "mm", dpi = 600)

# Poster plots
ggplot() +
  geom_jitter(
    data = med_imp_str2,
    aes(x = type, y = best_imp, color = type, ),
    alpha = 0.60,
    size = 4
  ) +
  labs(y = "Rank of Top Variable", x = "Variable Composition Type") +
  scale_color_manual(
    values = c(
      "Residual" = "#2A3132",
      "Spatial" = "#763626",
      "Temporal" = "#90AFC5"
    )
  ) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "#EADBCB"),
    panel.background = element_rect(fill = "#EADBCB"),
    text = element_text(color = "#2A3132"),
    panel.grid.major.y = element_line(colour = "#2A3132", linetype = 3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

ggsave("poster_fig1.jpeg", width = 180, height = 140, units = "mm", dpi = 600)
ggsave("poster_fig1.svg", width = 180, height = 140, units = "mm", dpi = 600)

gp1 <- ggplot() +
  geom_point(
    aes(
      x = data_pred_r[[32]]$y_pred_sp,
      y = data_pred_r[[32]]$y,
      group = data_pred_r[[32]]$strata
    ),
    color = "#763626",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[32]]$y_pred_sp,
      y = data_pred_r[[32]]$y,
      group = data_pred_r[[32]]$strata
    ),
    linewidth = 1.25,
    color = "#90AFC5",
    alpha = 0.5,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  #scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "American Robin",
    y = "Observed Abudance",
    x = "Predicted Abundance"
  ) +
  theme(
    plot.background = element_rect(fill = "#EADBCB"),
    panel.background = element_rect(fill = "#EADBCB"),
    text = element_text(color = "#2A3132"),
    panel.grid.major = element_line(colour = "#2A3132", linetype = 3),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16)
  )

gp2 <- ggplot() +
  geom_point(
    aes(
      x = data_pred_r[[52]]$y_pred_sp,
      y = data_pred_r[[52]]$y,
      group = data_pred_r[[52]]$strata
    ),
    color = "#763626",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[52]]$y_pred_sp,
      y = data_pred_r[[52]]$y,
      group = data_pred_r[[52]]$strata
    ),
    linewidth = 1.25,
    color = "#90AFC5",
    alpha = 0.5,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  #scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Western Meadowlark",
    y = "Observed Abudance",
    x = "Predicted Abundance"
  ) +
  theme(
    plot.background = element_rect(fill = "#EADBCB"),
    panel.background = element_rect(fill = "#EADBCB"),
    text = element_text(color = "#2A3132"),
    panel.grid.major = element_line(colour = "#2A3132", linetype = 3),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16)
  )

gp3 <- ggplot() +
  geom_point(
    aes(
      x = data_pred_r[[26]]$y_pred_sp,
      y = data_pred_r[[26]]$y,
      group = data_pred_r[[26]]$strata
    ),
    color = "#763626",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[26]]$y_pred_sp,
      y = data_pred_r[[26]]$y,
      group = data_pred_r[[26]]$strata
    ),
    linewidth = 1.25,
    color = "#90AFC5",
    alpha = 0.5,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  #scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Northern Mockingbird",
    y = "Observed Abudance",
    x = "Predicted Abundance"
  ) +
  theme(
    plot.background = element_rect(fill = "#EADBCB"),
    panel.background = element_rect(fill = "#EADBCB"),
    text = element_text(color = "#2A3132"),
    panel.grid.major = element_line(colour = "#2A3132", linetype = 3),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16)
  )

gp4 <- ggplot() +
  geom_point(
    aes(
      x = data_pred_r[[17]]$y_pred_sp,
      y = data_pred_r[[17]]$y,
      group = data_pred_r[[17]]$strata
    ),
    color = "#763626",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[17]]$y_pred_sp,
      y = data_pred_r[[17]]$y,
      group = data_pred_r[[17]]$strata
    ),
    linewidth = 1.25,
    color = "#90AFC5",
    alpha = 0.5,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  #scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Northern Cardinal",
    y = "Observed Abudance",
    x = "Predicted Abundance"
  ) +
  theme(
    plot.background = element_rect(fill = "#EADBCB"),
    panel.background = element_rect(fill = "#EADBCB"),
    text = element_text(color = "#2A3132"),
    panel.grid.major = element_line(colour = "#2A3132", linetype = 3),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16)
  )

gp1 +
  gp2 +
  gp3 +
  gp4 +
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag.position = "topright")

ggsave("poster_fig2al.svg", width = 180, height = 140, units = "mm", dpi = 600)

ggplot() +
  geom_point(
    aes(x = cor_all$r, y = cor_temp$r),
    color = "#763626",
    #alpha = 0.75,
    size = 4
  ) +
  labs(y = "Median Temporal Correlation", x = "Spatiotemporal Correlation") +
  theme(
    plot.background = element_rect(fill = "#EADBCB"),
    panel.background = element_rect(fill = "#EADBCB"),
    text = element_text(color = "#2A3132"),
    panel.grid.major = element_line(colour = "#2A3132", linetype = 3),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16)
  )

ggsave("poster_fig3.svg", width = 180, height = 140, units = "mm", dpi = 600)

map1 <- load_map(stratify_by = "bbs_usgs") %>%
  filter(str_detect(strata_name, "US")) %>%
  filter(!str_detect(strata_name, "AK"))

ggplot(data = map1) +
  geom_sf(fill = "#763626", col = "#2A3132", linewidth = 1.5, alpha = 0.9) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    panel.grid = element_blank()
  )

ggsave("poster_map1.svg", width = 180, height = 140, units = "mm", dpi = 600)
