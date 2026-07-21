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


# Save combined decomposed DBO data for the SVC analysis ------------------
# (prefix above is the data-prep portion of fit_rf_dbo.R, through the join that
# builds the family x station x year table with decomposed environmental
# components + one-year lags.)

saveRDS(data_rf, "data/data_dbo.rds")
cat("Saved data/data_dbo.rds:", nrow(data_rf), "rows,",
    dplyr::n_distinct(data_rf$family), "families,",
    dplyr::n_distinct(data_rf$StationNme), "stations,",
    "years", paste(range(data_rf$DataYear), collapse = "-"), "\n")
