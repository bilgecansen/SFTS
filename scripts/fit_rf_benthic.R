library(tidyverse)
library(readxl)
library(foreach)
library(patchwork)
library(ranger)

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

# Sea ice persistence
data_per <- readRDS("data/data_sic_per.rds")

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

# add sic and per data
data_env_sel <- left_join(
  data_env_sel,
  data_per,
  by = c("StationNme", "DataYear")
)

# decompose data for RF
## these covariates vary across space and time and we should decompose them
decomp_vars <- colnames(data_env_sel)[-(1:10)]

data_env_dec <- data_env_sel %>%
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

data_rf <- left_join(
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

med_imp_fam1 <- foreach(i = 1:length(var_imp_fam_sel), .combine = "rbind") %do%
  {
    z <- var_imp_fam_sel[[i]][order(var_imp_fam_sel[[i]], decreasing = T)]

    data.frame(
      best_imp = c(
        min(str_which(names(z), "spatial")),
        min(str_which(names(z), "temporal")),
        min(str_which(names(z), "residual"))
      ),
      type = c("spatial", "temporal", "residual")
    )
  }

ggplot() +
  geom_jitter(
    data = med_imp_fam1,
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


# Cross-validation with space only models ---------------------------------

families_sp <- families[which(R2_fam >= 0.25)]

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
      contains("spatial")
    )
)

data_rf_sp <- data_rf_sp[which(R2_fam >= 0.25)]
data_rf_test <- data_rf[which(R2_fam >= 0.25)]

# Predictions with SFTS
data_pred_r <- foreach(i = 1:length(families_sp)) %do%
  {
    years <- unique(data_rf_sp[[i]]$DataYear)

    foreach(t = 1:length(years), .combine = "rbind") %do%
      {
        dat_train <- filter(data_rf_sp[[i]], DataYear != years[t]) %>%
          select(-StationNme, -DataYear)

        dat_test <- filter(data_rf_test[[i]], DataYear == years[t]) %>%
          select(
            -family,
            -contains(c("spatial", "temporal", "residual")),
          )

        colnames(dat_test)[-(1:7)] <- paste(
          colnames(dat_test)[-(1:7)],
          "spatial",
          sep = "_"
        )

        res_rf_sp <- ranger(
          biomass ~ .,
          data = dat_train,
          num.trees = 2000
        )

        y_pred <- predict(
          res_rf_sp,
          data = select(dat_test, -StationNme, -DataYear)
        )

        y_pred <- (y_pred$predictions -
          min(dat_train$biomass)) /
          (max(dat_train$biomass) -
            min(dat_train$biomass))

        y <- (dat_test$biomass -
          min(dat_train$biomass)) /
          (max(dat_train$biomass) -
            min(dat_train$biomass))

        data.frame(
          y_pred = y_pred,
          y = y,
          species_id = families_sp[i],
          site = dat_test$StationNme,
          year = dat_test$DataYear
        )
      }
  }

pred_r <- map_dbl(data_pred_r, function(x) cor(x$y_pred, x$y))
idx_pred <- which(pred_r > 0.5)
data_pred_r <- data_pred_r[idx_pred]
pred_r <- pred_r[idx_pred]
families_sp_sel <- families_sp[idx_pred]

order(pred_r, decreasing = T)[1:11]

for (i in 1:11) {
  print(
    ggplot() +
      geom_point(
        aes(
          x = data_pred_r[[i]]$y_pred,
          y = data_pred_r[[i]]$y,
          group = data_pred_r[[i]]$site
        ),
        color = "pink3",
        alpha = 0.5,
        size = 2
      ) +
      geom_line(
        aes(
          x = data_pred_r[[i]]$y_pred,
          y = data_pred_r[[i]]$y,
          group = data_pred_r[[i]]$site
        ),
        linewidth = 1.25,
        color = "#333D79FF",
        alpha = 0.75,
        se = F,
        method = "lm",
        stat = "smooth"
      ) +
      #scale_y_continuous(limits = c(0, 1)) +
      #scale_x_continuous(limits = c(0, 1)) +
      labs(
        #title = "American Robin",
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
  )
}


g1 <- ggplot() +
  geom_point(
    aes(
      x = data_pred_r[[5]]$y_pred,
      y = data_pred_r[[5]]$y,
      group = data_pred_r[[5]]$site
    ),
    color = "pink3",
    alpha = 0.5,
    size = 2
  ) +
  geom_line(
    aes(
      x = data_pred_r[[5]]$y_pred,
      y = data_pred_r[[5]]$y,
      group = data_pred_r[[5]]$site
    ),
    linewidth = 1.25,
    color = "#333D79FF",
    alpha = 0.75,
    se = F,
    method = "lm",
    stat = "smooth"
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  #scale_x_continuous(limits = c(0, 1)) +
  labs(
    #title = "American Robin",
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
