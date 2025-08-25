library(tidyverse)
library(readxl)
library(foreach)
library(patchwork)
library(ranger)
library(ggbeeswarm)

theme_set(theme_bw())


# Prep data ---------------------------------------------------------------

# Env data for sites except UTN
data_env <-
  read_excel(
    "data/EnvtlVariables_UniqueID_0-5nm_DBO_063023_v3_abwk122023.xlsx",
    sheet = 4
  ) %>%
  select(-orignme, -UnofficialDBO)
names(data_env)[3] <- "StationNme"

# Env data for UTN
data_env_utn <-
  read_excel("data/EnvironmentalVariables_UTN_2000-2019.xlsx", sheet = 2) %>%
  filter(StationName_Standard != "UTN5=SEC1") %>%
  select(-grabnum)
names(data_env_utn)[3] <- "StationNme"
data_env_utn$uniqueID <- paste(
  data_env_utn$CruiseID,
  data_env_utn$StationNum,
  data_env_utn$StationNme,
  sep = ""
)
data_env_utn <- relocate(data_env_utn, uniqueID, .before = DBOnme)

#--
data_env <- rbind(data_env, data_env_utn)

data_env$integchla <- as.numeric(data_env$integchla)
data_env$sedchla <- as.numeric(data_env$sedchla)
data_env$NiTriTra <- as.numeric(data_env$NiTriTra)
data_env$Ammonia <- as.numeric(data_env$Ammonia)
data_env$Phosphate <- as.numeric(data_env$Phosphate)
data_env$Silicate <- as.numeric(data_env$Silicate)

# Family biomass data
data_fam <- readRDS("data/data_family_bio.rds")

data_env <- data_env[data_env$uniqueID %in% data_fam$uniqueID, ]
data_fam <- data_fam[data_fam$uniqueID %in% data_env$uniqueID, ]

data_env <- data_env[order(data_env$uniqueID), ]
data_fam <- data_fam[order(data_fam$uniqueID), ]

# Select relevant variables and stations
data_fam_gc <- data_fam %>%
  filter(str_detect(CruiseID, pattern = "SWL")) %>%
  filter(!StationNme %in% c("DBO2.7", "BCL6C", "UTBS2A"))

data_env_sel <- select(
  data_env,
  -DBOnme,
  -Abundance,
  -BiomWW,
  -BiomGC,
  -Taxanum,
  -del13C,
  -del15N,
  -delO18,
  -SCOC,
  -DepthBOT,
  -pressure,
  -phi1,
  -phi2,
  -phi3,
  -phi4,
  -ModalSz,
  -TON,
  -philt0,
  -starts_with("Resp"),
  -phi1to4,
  -bwChla,
  -swe,
  -swi
) %>%
  filter(str_detect(CruiseID, pattern = "SWL")) %>%
  filter(!StationNme %in% c("DBO2.7", "BCL6C", "UTBS2A"))

## Subset to DBO 1,2,3
idx_dbo <- which(data_env_sel$DBOreg %in% c(1, 2, 3))
data_fam_gc <- data_fam_gc[idx_dbo, ]
data_env_sel <- data_env_sel[idx_dbo, ]

## remove NAs
idx_na <- foreach(i = 1:ncol(data_env_sel), .combine = "c") %do%
  {
    which(is.na(data_env_sel[, i]))
  }
idx_na <- unique(idx_na)

data_fam_gc <- data_fam_gc[-idx_na, ]
data_env_sel <- data_env_sel[-idx_na, ]

## Remove families not observed in any stations
idx_fam <- map_lgl(data_fam_gc, function(x) any(x > 0))
idx_fam[c(7, 8)] <- TRUE
data_fam_gc <- data_fam_gc[, which(idx_fam)]

## Remove sites with less than 5 years of observations,
## separately for each family
data_fam_gc <- data_fam_gc %>%
  pivot_longer(
    cols = ends_with("AE"),
    names_to = "family",
    values_to = "biomass"
  )

year_n <- data_fam_gc %>%
  group_by(family, StationNme) %>%
  summarise(n = length(biomass[biomass > 0])) %>%
  filter(n > 4)

idx_ss <- paste(year_n$family, year_n$StationNme)
data_fam_gc$family_site <- paste(data_fam_gc$family, data_fam_gc$StationNme)

data_fam_gc2 <- filter(data_fam_gc, family_site %in% idx_ss)

## Remove families with less than 5 stations
site_n <- data_fam_gc2 %>%
  group_by(family) %>%
  summarize(n = length(unique(StationNme))) %>%
  filter(n > 4)

data_fam_gc3 <- filter(data_fam_gc2, family %in% site_n$family) %>%
  select(
    -CruiseID,
    -StationNum,
    -uniqueID,
    -DataDate,
    -family_site
  )

## remove 0s
idx0 <- which(data_fam_gc3$biomass == 0)
data_fam_gc4 <- data_fam_gc3[-idx0, ]

# Add lags to env data
z <- select(
  data_env_sel,
  -CruiseID,
  -StationNum,
  -uniqueID,
  -DataDate,
  -Latitude,
  -Longitude,
  -Depth,
  -DBOreg
)

stations <- unique(z$StationNme)

z2 <- foreach(i = 1:length(stations), .combine = "rbind") %do%
  {
    filter(z, StationNme == stations[i]) %>%
      filter(DataYear != max(DataYear))
  }
z2 <- z2[order(z2$DataYear), ]

z3 <- foreach(i = 1:length(stations), .combine = "rbind") %do%
  {
    filter(data_env_sel, StationNme == stations[i]) %>%
      filter(DataYear != min(DataYear))
  }
z3 <- z3[order(z3$DataYear), ]

z2$DataYear <- z3$DataYear
colnames(z2)[-(1:2)] <- paste(colnames(z2)[-(1:2)], "lag", sep = "_")

data_env_sel2 <- cbind(z3, z2[, -(1:2)])

# decompose data for RF
## these covariates vary across space and time and we should decompose them
decomp_vars <- colnames(data_env_sel2)[-(1:10)]

data_env_dec <- data_env_sel2 %>%
  # center without scaling
  dplyr::mutate(across(all_of(decomp_vars), ~ scale(., scale = FALSE)[, 1])) %>%

  # compute the spatial component
  # (mean by pixel across all years of centered variables)
  dplyr::group_by(StationNme) %>%
  dplyr::mutate(across(
    all_of(decomp_vars),
    ~ mean(., na.rm = TRUE),
    .names = "{.col}_spatial"
  )) %>%
  dplyr::ungroup() %>%
  # compute the temporal component
  # (mean by year across all pixels of centered variables)
  dplyr::group_by(DataYear) %>%
  dplyr::mutate(across(
    all_of(decomp_vars),
    ~ mean(., na.rm = TRUE),
    .names = "{.col}_temporal"
  )) %>%
  dplyr::ungroup() %>%

  # compute residual for each site i and year j as centered variable value
  # i,j - spatial mean i - temporal mean j
  dplyr::mutate(across(
    tidyselect::all_of(decomp_vars),
    ~ . -
      get(paste0(cur_column(), "_spatial")) -
      get(paste0(cur_column(), "_temporal")),
    .names = "{.col}_residual"
  )) %>%
  dplyr::arrange(StationNme, DataYear) %>%
  select(
    -CruiseID,
    -StationNum,
    -uniqueID,
    -DataDate,
    -Latitude,
    -Longitude
  )


# Biomass RF models -------------------------------------------------------

data_rf <- inner_join(
  data_fam_gc4,
  data_env_dec,
  by = c("StationNme", "DataYear")
)

families <- unique(data_rf$family)

data_rf <- foreach(i = 1:length(families)) %do%
  {
    filter(data_rf, family == families[i])
  }

names(data_rf) <- families

pb <- txtProgressBar(max = length(data_rf), style = 3)
var_imp_fam <- list()
R2_fam <- c()

# All variables
for (i in 1:length(data_rf)) {
  res_rf <- ranger(
    log(biomass) ~ .,
    data = select(
      data_rf[[i]],
      biomass,
      DBOreg,
      Latitude,
      Longitude,
      Depth,
      contains(c("spatial", "temporal", "residual")),
    ),
    num.trees = 2000,
    importance = "permutation"
  )

  R2_fam[i] <- res_rf$r.squared
  var_imp_fam[[i]] <- res_rf$variable.importance

  setTxtProgressBar(pb, i)
}

hist(R2_fam)
var_imp_fam_sel <- var_imp_fam[which(R2_fam >= 0.25)]

med_imp_fam <- foreach(i = 1:length(var_imp_fam_sel), .combine = "rbind") %do%
  {
    z <- var_imp_fam_sel[[i]][order(var_imp_fam_sel[[i]], decreasing = T)]

    data.frame(
      best_imp = c(
        min(str_which(names(z), "spatial")),
        min(str_which(names(z), "temporal")),
        min(str_which(names(z), "residual"))
      ),
      type = c("Spatial", "Temporal", "Residual")
    )
  }

g_imp <- ggplot(med_imp_fam) +
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


# Cross-validation with space only models ---------------------------------

data_rf_sp <- map(
  data_rf,
  function(x)
    select(
      x,
      DataYear,
      StationNme,
      biomass,
      DBOreg,
      Latitude,
      Longitude,
      Depth,
      contains("spatial"),
      -contains("lag")
    )
)

data_rf_decomp <- map(
  data_rf,
  function(x)
    select(
      x,
      DataYear,
      StationNme,
      biomass,
      DBOreg,
      Latitude,
      Longitude,
      Depth,
      contains(c("spatial", "temporal", "residual")),
      -contains("lag_spatial")
    )
)

data_rf_sptemp <- map(
  data_rf,
  function(x)
    select(
      x,
      -contains(c("spatial", "temporal", "residual")),
      -family
    )
)

data_rf_test <- data_rf

# Predictions with SFTS
data_pred <- foreach(i = 1:length(families)) %do%
  {
    years <- unique(data_rf_sp[[i]]$DataYear)

    foreach(t = 1:length(years), .combine = "rbind") %do%
      {
        dat_train_sp <- filter(
          data_rf_sp[[i]],
          DataYear != years[t]
        ) %>%
          select(-StationNme, -DataYear)

        dat_train_sptemp <- filter(
          data_rf_sptemp[[i]],
          DataYear != years[t]
        ) %>%
          select(-StationNme, -DataYear)

        dat_train_decomp <- filter(
          data_rf_decomp[[i]],
          DataYear != years[t]
        ) %>%
          select(-StationNme, -DataYear)

        dat_test_sp <- filter(
          data_rf_test[[i]],
          DataYear == years[t]
        ) %>%
          select(
            -family,
            -contains(c("spatial", "temporal", "residual")),
            -contains("lag")
          )

        colnames(dat_test_sp)[-(1:7)] <- paste(
          colnames(dat_test_sp)[-(1:7)],
          "spatial",
          sep = "_"
        )

        dat_test_sptemp <- filter(
          data_rf_test[[i]],
          DataYear == years[t]
        ) %>%
          as.data.frame()

        res_rf_sp <- ranger(
          log(biomass) ~ .,
          data = dat_train_sp,
          num.trees = 2000
        )

        res_rf_decomp <- ranger(
          log(biomass) ~ .,
          data = dat_train_decomp,
          num.trees = 2000
        )

        res_rf_sptemp <- ranger(
          log(biomass) ~ .,
          data = dat_train_sptemp,
          num.trees = 2000
        )

        y_pred_sp <- predict(
          res_rf_sp,
          data = select(dat_test_sp, -StationNme, -DataYear)
        )

        y_pred_decomp <- predict(
          res_rf_decomp,
          data = select(dat_test_sptemp, -StationNme, -DataYear)
        )

        y_pred_sptemp <- predict(
          res_rf_sptemp,
          data = select(dat_test_sptemp, -StationNme, -DataYear)
        )

        y_pred_sp <- (y_pred_sp$predictions)
        y_pred_sptemp <- (y_pred_sptemp$predictions)
        y_pred_decomp <- (y_pred_decomp$predictions)

        y <- (log(dat_test_sp$biomass))

        data.frame(
          y_pred_sp = y_pred_sp,
          y_pred_decomp = y_pred_decomp,
          y_pred_sptemp = y_pred_sptemp,
          y = y,
          species_id = families[i],
          site = dat_test_sp$StationNme,
          year = dat_test_sp$DataYear
        )
      }
  }

plot_pred <- function(i, y_pred) {
  z <- data_pred[[i]]

  ggplot(z) +
    geom_line(
      aes(
        x = {{ y_pred }},
        y = y,
        group = site
      ),
      linewidth = 1.25,
      color = "#90AFC5",
      alpha = 0.75,
      se = F,
      method = "lm",
      stat = "smooth"
    ) +
    geom_point(
      data = group_by(z, site) %>%
        summarise(pred = mean({{ y_pred }}), y = mean(y)),
      aes(
        x = pred,
        y = y,
        group = site
      ),
      color = "#763626",
      ,
      alpha = 0.5,
      size = 2
    ) +
    labs(
      title = families[i],
      y = "Observed Biomass",
      x = "Predicted Biomass"
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
}

g1 <- plot_pred(11, y_pred_sp)
g2 <- plot_pred(11, y_pred_sptemp)
g3 <- plot_pred(11, y_pred_decomp)

# Correlation comparison for SFTS
plot_sp_multi <- function(y_pred) {
  data_pred_sp <- do.call(rbind, data_pred)

  cor_temp <- data_pred_sp %>%
    group_by(species_id, site) %>%
    summarise(r = cor(y, {{ y_pred }})) %>%
    ungroup() %>%
    group_by(species_id) %>%
    summarise(
      r_med = median(r[!is.na(r)]),
      r_min = quantile(r[!is.na(r)], 0.05),
      r_max = quantile(r[!is.na(r)], 0.95)
    )

  cor_sp <- data_pred_sp %>%
    group_by(species_id, site) %>%
    summarise(y = mean(y), y_pred_sp = mean({{ y_pred }})) %>%
    ungroup() %>%
    group_by(species_id) %>%
    summarise(r = cor(y, y_pred_sp))

  dat_cor <- data.frame(
    r = c(cor_sp$r, cor_temp$r_med),
    type = rep(
      c("Species-wide", "Population-level"),
      each = length(cor_sp$r)
    )
  )

  ggplot(dat_cor) +
    geom_violin(aes(x = type, y = r), linewidth = 1.5) +
    geom_quasirandom(
      aes(x = type, y = r, color = type),
      alpha = 0.5,
      size = 2
    ) +
    scale_color_manual(
      values = c(
        "Species-wide" = "#763626",
        "Population-level" = "#90AFC5"
      )
    ) +
    labs(y = "Prediction Correlation", x = "Prediction Type") +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(size = 10),
      legend.position = "none"
    ) +
    scale_y_continuous(
      limits = c(-1.01, 1),
      breaks = c(-0.6, -0.3, 0, 0.3, 0.6, 0.9)
    )
}

g_pred1 <- plot_sp_multi(y_pred_sp)
g_pred2 <- plot_sp_multi(y_pred_sptemp)
g_pred3 <- plot_sp_multi(y_pred_decomp)

# Component variability plot
v <- distinct(
  data_env_dec,
  StationNme,
  DataYear,
  Temp_spatial,
  Temp_temporal,
  Temp_residual
) %>%
  pivot_longer(cols = c(Temp_spatial, Temp_temporal, Temp_residual))

g4 <- ggplot(data = v) +
  geom_violin(
    aes(x = name, y = value),
    fill = "#485B7C",
    color = "#792C1F",
    alpha = 0.5,
    linewidth = 1.05
  ) +
  labs(x = "Variable Components", y = "Standardized Temperature") +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 10)
  ) +
  scale_x_discrete(
    labels = c(
      "Temp_temporal" = "Temporal",
      "Temp_spatial" = "Spatial",
      "Temp_residual" = "Residual"
    )
  )

plots_dbo <- list(
  g1 = g1,
  g2 = g2,
  g3 = g3,
  g4 = g4,
  g_imp = g_imp,
  g_pred1 = g_pred1,
  g_pred2 = g_pred2,
  g_pred3 = g_pred3
)

saveRDS(plots_dbo, "plots_dbo.rds")
