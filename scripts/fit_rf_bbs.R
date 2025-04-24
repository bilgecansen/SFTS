library(tidyverse)
library(ranger)
library(foreach)
library(patchwork)

theme_set(theme_bw())

bbs <- readRDS("data/data_bbs_nozero.rds")

bbs_train <- filter(bbs, year < 2011)
bbs_test <- filter(bbs, year >= 2011 & year < 2021)

# only sites with at least 10 years of observations
year_n <- bbs_train %>%
  group_by(species_id, site_id) %>%
  summarise(n = length(abundance[abundance > 0])) %>%
  filter(n > 9)

idx_ss <- paste(year_n$species_id, year_n$site_id)
bbs_train$species_site <- paste(bbs_train$species_id, bbs_train$site_id)

bbs_train2 <- filter(bbs_train, species_site %in% idx_ss)

# only species with at least 10 sites
site_n <- bbs_train2 %>%
  group_by(species_id) %>%
  summarize(n = length(unique(site_id))) %>%
  filter(n > 9)

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
    alpha = 0.75,
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
      x = data_pred_r[[50]]$y_pred_sp,
      y = data_pred_r[[50]]$y,
      group = data_pred_r[[50]]$strata
    ),
    color = "pink3",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[50]]$y_pred_sp,
      y = data_pred_r[[50]]$y,
      group = data_pred_r[[50]]$strata
    ),
    linewidth = 1.25,
    color = "#333D79FF",
    alpha = 0.75,
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
    alpha = 0.75,
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
      x = data_pred_r[[89]]$y_pred_sp,
      y = data_pred_r[[89]]$y,
      group = data_pred_r[[89]]$strata
    ),
    color = "pink3",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[89]]$y_pred_sp,
      y = data_pred_r[[89]]$y,
      group = data_pred_r[[89]]$strata
    ),
    linewidth = 1.25,
    color = "#333D79FF",
    alpha = 0.75,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Yellow-rumped Warbler",
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


# Simulations -------------------------------------------------------------

a <- 10
b <- 10
mu_x <- 1:10
mu_ab <- a + b * mu_x

m_ab <- lapply(mu_ab, function(x) rnorm(20, x, 10))
m_ab <- do.call(rbind, m_ab)
m_ab <- data.frame(m_ab)
colnames(m_ab) <- 1:20
m_ab$site <- 1:10

m_x <- lapply(mu_x, function(x) rnorm(20, x, 2))
m_x <- do.call(rbind, m_x)
m_x <- data.frame(m_x)
colnames(m_x) <- 1:20
m_x$site <- 1:10

m_ab <- pivot_longer(m_ab, !site, names_to = "time", values_to = "abundance")

m_x <- pivot_longer(m_x, !site, names_to = "time", values_to = "env")

m <- left_join(m_ab, m_x, by = c("site", "time"))

gm1 <- ggplot(data = m) +
  geom_point(
    aes(x = env, y = abundance),
    color = "pink3",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = env,
      y = abundance,
      group = site
    ),
    linewidth = 1.25,
    color = "#333D79FF",
    alpha = 0.75,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  labs(y = "Simulated Abundance", x = "Environmental Variable") +
  scale_y_continuous(breaks = seq(0, 125, 25), limits = c(0, 135)) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    #axis.title.x = element_text(margin = margin(t = 10)),
    #axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 10)
  )

# Predict a test set
m_ab_test <- lapply(mu_ab, function(x) rnorm(20, x, 10))
m_ab_test <- do.call(rbind, m_ab_test)
m_ab_test <- data.frame(m_ab_test)
colnames(m_ab_test) <- 1:20
m_ab_test$site <- 1:10

m_x_test <- lapply(mu_x, function(x) rnorm(20, x, 2))
m_x_test <- do.call(rbind, m_x_test)
m_x_test <- data.frame(m_x_test)
colnames(m_x_test) <- 1:20
m_x_test$site <- 1:10

m_ab_test <- pivot_longer(
  m_ab_test,
  !site,
  names_to = "time",
  values_to = "abundance"
)

m_x_test <- pivot_longer(
  m_x_test,
  !site,
  names_to = "time",
  values_to = "env_spatial"
)

m_test <- left_join(m_ab_test, m_x_test, by = c("site", "time"))

# decompose simulated env data
m_train <- m %>%
  # center without scaling
  dplyr::mutate(across(env, ~ scale(., scale = FALSE)[, 1])) %>%

  # compute the spatial component
  # (mean by pixel across all years of centered variables)
  dplyr::group_by(site) %>%
  dplyr::mutate(across(
    env,
    ~ mean(., na.rm = TRUE),
    .names = "{.col}_spatial"
  )) %>%
  dplyr::ungroup() %>%
  # compute the temporal component
  # (mean by year across all pixels of centered variables)
  dplyr::group_by(time) %>%
  dplyr::mutate(across(
    env,
    ~ mean(., na.rm = TRUE),
    .names = "{.col}_temporal"
  )) %>%
  dplyr::ungroup() %>%

  # compute residual for each site i and year j as centered variable value
  # i,j - spatial mean i - temporal mean j
  dplyr::mutate(across(
    env,
    ~ . -
      get(paste0(cur_column(), "_spatial")) -
      get(paste0(cur_column(), "_temporal")),
    .names = "{.col}_residual"
  )) %>%
  dplyr::arrange(site, time)

res_sim <- lm(abundance ~ env_spatial, data = m_train)
pred_sim <- predict.lm(res_sim, as.data.frame(m_test[, 4]))

m_test$pred <- pred_sim

gm2 <- ggplot(data = m_test) +
  geom_point(
    aes(x = pred, y = abundance),
    color = "pink3",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = pred,
      y = abundance,
      group = site
    ),
    linewidth = 1.25,
    color = "#333D79FF",
    alpha = 0.75,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  labs(y = "Simulated Abundance", x = "Predicted Abundance") +
  scale_y_continuous(breaks = seq(0, 125, 25), limits = c(0, 135)) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    #axis.title.x = element_text(margin = margin(t = 10)),
    #axis.title.y = element_text(margin = margin(r = 10)),
    axis.title = element_text(size = 10)
  )

gm1 + gm2 + plot_layout(axes = "collect") + plot_annotation(tag_levels = "a")

ggsave("fig1.pdf", width = 180, height = 90, units = "mm", dpi = 600)
