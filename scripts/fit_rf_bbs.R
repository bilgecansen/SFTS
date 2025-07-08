library(tidyverse)
library(ranger)
library(foreach)
library(patchwork)
library(ggbeeswarm)

theme_set(theme_bw())

bbs <- readRDS("data/data_bbs_nozero.rds")

bbs_train <- filter(bbs, year < 2011)
bbs_test <- filter(bbs, year >= 2011 & year < 2021)

# only sites with at least 5 years of non-zero observations
year_n <- bbs_train %>%
  group_by(species_id, site_id) %>%
  summarise(n = length(abundance[abundance > 0])) %>%
  filter(n > 4)

idx_ss <- paste(year_n$species_id, year_n$site_id)
bbs_train$species_site <- paste(bbs_train$species_id, bbs_train$site_id)

bbs_train2 <- filter(bbs_train, species_site %in% idx_ss)

# only species with at least 5 sites of occurrence
site_n <- bbs_train2 %>%
  group_by(species_id) %>%
  summarize(n = length(unique(site_id))) %>%
  filter(n > 4)

bbs_train3 <- filter(bbs_train2, species_id %in% site_n$species_id)

species <- unique(bbs_train3$species_id)

# remove unused large objects
rm(bbs, bbs_train, bbs_train2)


# Abundance RF models on strata scale -------------------------------------

# aggregate to strata and select decomposed variables
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
        contains(c(
          "bio2",
          "bio3",
          "bio5",
          "bio8",
          "bio9",
          "bio15",
          "bio16",
          "bio18"
        )),
        -ends_with(c(
          "bio2",
          "bio3",
          "bio5",
          "bio8",
          "bio9",
          "bio15",
          "bio16",
          "bio18"
        )),
        -contains("lag_spatial"),
        strata
      )
  }
names(data_rf_str) <- species

## Remove species with too few strata
idx_str <- which(map_dbl(data_rf_str, function(x) length(unique(x$strata))) > 4)
data_rf_str <- data_rf_str[idx_str]

# Random forests and variable importance
pb <- txtProgressBar(max = length(data_rf_str), style = 3)
var_imp_str <- list()
R2_str <- c()

for (i in 1:length(data_rf_str)) {
  res_rf <- ranger(
    log(abundance) ~ .,
    data = data_rf_str[[i]],
    num.trees = 2000,
    importance = "permutation"
  )

  R2_str[i] <- res_rf$r.squared
  var_imp_str[[i]] <- res_rf$variable.importance

  setTxtProgressBar(pb, i)
}

hist(R2_str)
var_imp_str_sel <- var_imp_str[which(R2_str >= 0.25)]

med_imp_str <- foreach(i = 1:length(var_imp_str_sel), .combine = "rbind") %do%
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

g_imp <- ggplot(med_imp_str) +
  geom_violin(aes(x = type, y = best_imp), linewidth = 1.5, scale = "width") +
  geom_quasirandom(
    aes(x = type, y = best_imp + rnorm(length(best_imp), 0, 0.5), color = type),
    alpha = 0.5,
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "Residual" = "#2A3132",
      "Spatial" = "#763626",
      "Temporal" = "#90AFC5"
    )
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

g_imp


# Predict test set with space only models ---------------------------------

data_rf_sp <- map(
  data_rf_str,
  function(x)
    select(
      x,
      abundance,
      elevs,
      lat,
      long,
      strata,
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
    num.trees = 2000
  )

  R2_sp[i] <- res_rf_sp[[i]]$r.squared

  setTxtProgressBar(pb, i)
}

hist(R2_sp)

# Predictions with SFTS
species_sp <- species[idx_str]

data_pred_r <- foreach(i = 1:length(species_sp)) %do%
  {
    dat <- filter(bbs_test, species_id == as.numeric(species_sp[i])) %>%
      dplyr::select(-site_id) %>%
      group_by(species_id, strata, year) %>%
      summarize(across(everything(), mean)) %>%
      ungroup()

    sites <- unique(dat$strata)

    dat_pred <- foreach(h = 1:length(sites), .combine = "rbind") %do%
      {
        dat_site <- filter(dat, strata == sites[h]) %>%
          select(elevs, lat, long, strata, party_hours, contains("bio")) %>%
          select(-contains("lag")) %>%
          rename_with(~ paste(., "spatial", sep = "_"), contains("bio")) %>%
          as.data.frame()

        if (nrow(dat_site) < 5) return(NA)

        y_pred_sp <- predict(res_rf_sp[[i]], data = dat_site)

        y_pred_sp <- y_pred_sp$predictions #-
        #min(log(data_rf_sp[[i]]$abundance))) /
        #(max(log(data_rf_sp[[i]]$abundance)) -
        #min(log(data_rf_sp[[i]]$abundance)))

        y <- log(filter(dat, strata == sites[h])$abundance) #-
        #min(log(data_rf_sp[[i]]$abundance))) /
        #(max(log(data_rf_sp[[i]]$abundance)) -
        #min(log(data_rf_sp[[i]]$abundance)))

        data.frame(
          y_pred_sp = y_pred_sp,
          y = y,
          species_id = species_sp[i],
          strata = sites[h]
        )
      }

    dat_pred[!is.na(dat_pred$y), ]
  }

pred_r <- map_dbl(data_pred_r, function(x) cor(x$y_pred_sp, x$y))
order(pred_r, decreasing = T)[1:10]

g1 <- ggplot(data_pred_r[[47]]) +
  geom_point(
    aes(
      x = y_pred_sp,
      y = y,
      group = strata
    ),
    color = "#763626",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = y_pred_sp,
      y = y,
      group = strata
    ),
    linewidth = 1.25,
    color = "#90AFC5",
    alpha = 0.6,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  labs(
    title = "American Robin (Turdus migratorius)",
    y = "Observed Abudance",
    x = "Predicted Abundance"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

# Correlation comparison
data_pred_r2 <- do.call(rbind, data_pred_r)

cor_temp <- data_pred_r2 %>%
  group_by(species_id, strata) %>%
  summarise(r = cor(y, y_pred_sp)) %>%
  ungroup() %>%
  group_by(species_id) %>%
  summarise(
    r_med = median(r[!is.na(r)]),
    r_min = quantile(r[!is.na(r)], 0.05),
    r_max = quantile(r[!is.na(r)], 0.95)
  )

cor_sptemp <- data_pred_r2 %>%
  group_by(species_id) %>%
  summarise(r = cor(y, y_pred_sp))

dat_cor <- data.frame(
  r = c(cor_sptemp$r, cor_temp$r_med),
  type = rep(
    c("Spatio-temporal", "Temporal"),
    each = length(cor_sptemp$r)
  )
)

g_pred <- ggplot(dat_cor) +
  geom_violin(aes(x = type, y = r), linewidth = 1.5) +
  geom_quasirandom(aes(x = type, y = r, color = type), alpha = 0.5, size = 2) +
  scale_color_manual(
    values = c(
      "Spatio-temporal" = "#763626",
      "Temporal" = "#90AFC5"
    )
  ) +
  labs(y = "Prediction Correlation", x = "Prediction Type") +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 10)
  ) +
  scale_y_continuous(limits = c(-0.55, 0.9), breaks = c(-0.4, 0, 0.4, 0.8))

## Sptemp vs temp across sites
z <- left_join(cor_sptemp, cor_temp, by = "species_id")

g2 <- ggplot(data = z) +
  geom_errorbar(
    aes(x = r, ymin = r_min, ymax = r_max),
    color = "#763626",
    alpha = 0.5
  ) +
  geom_point(aes(y = r_med, x = r), color = "#763626", alpha = 0.9, size = 2) +
  labs(x = "Spatio-temporal correlation", y = "Temporal correlation") +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 10)
  )

# Proportion of sites with temporal cor > 0.5
d <- data_pred_r2 %>%
  group_by(species_id, strata) %>%
  summarise(r = cor(y, y_pred_sp)) %>%
  ungroup() %>%
  group_by(species_id) %>%
  summarise(n = length(which(r >= 0.5)) / length(r))

g3 <- ggplot(data = d) +
  geom_histogram(aes(n), color = "#90AFC5", fill = "#763626", alpha = 0.8) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 10)
  ) +
  labs(y = "Frequency", x = "Proportion of sites with r > 0.5")

plots_bbs <- list(
  g1 = g1,
  g2 = g2,
  g3 = g3,
  g_imp = g_imp,
  g_pred = g_pred
)

saveRDS(plots_bbs, "plots_bbs.rds")
