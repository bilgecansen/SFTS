library(tidyverse)
library(bbsBayes2)
library(sp)
library(prism)
library(geodata)
library(foreach)
library(terra)


# Load main bbs data ------------------------------------------------------

##' Use this if bbs data is not installed
##' fetch_bbs_data()

bbs_all <- load_bbs_data()

# Filter routes for weather and quality
routes <- bbs_all$routes %>%
  filter(rpid == 101, run_type == 1)

# main bbs file
bbs <- left_join(
  bbs_all$birds,
  routes,
  by = c(
    "route_data_id",
    "country_num",
    "state_num",
    "route",
    "bcr",
    "year",
    "rpid"
  )
) %>%
  mutate(site_id = state_num * 1000 + route) %>%
  select(aou, site_id, latitude, longitude, year, species_total) %>%
  rename(
    species_id = "aou",
    lat = "latitude",
    long = "longitude",
    abundance = "species_total"
  )

# remove counts with no route info
idx_lat <- which(!is.na(bbs$lat))
bbs <- bbs[idx_lat, ]


# filter BBS data ---------------------------------------------------------

# Combine subspecies into their common species
combine_subspecies <- function(df, species_table = bbs_all$species) {
  species_table <- species_table %>%
    mutate(species_full = paste(genus, species, sep = " "))

  # Subspecies have two spaces separated by non-spaces
  subspecies_names <- species_table %>%
    filter(aou %in% unique(df$species_id)) %>%
    pull(species_full) %>%
    grep(" [^ ]+ ", ., value = TRUE)

  subspecies_ids <- species_table %>%
    filter(species_full %in% subspecies_names) %>%
    pull(aou)

  # Drop all but the first two words to get the root species name,
  # then find the AOU code
  new_subspecies_ids <- species_table %>%
    slice(match(word(subspecies_names, 1, 2), species_table$species_full)) %>%
    pull(aou)

  # replace the full subspecies names with species-level names
  for (i in seq_along(subspecies_ids)) {
    df$species_id[df$species_id == subspecies_ids[i]] <- new_subspecies_ids[i]
  }

  df %>%
    group_by(site_id, year, species_id, lat, long) %>%
    summarise(abundance = sum(abundance)) %>%
    ungroup()
}

#' Filter poorly sampled BBS species
#'
#' Removes waterbirds, shorebirds, owls, kingfishers, knightjars,
#' dippers. These species are poorly sampled due to their aquatic or
#' noctural nature. Also removes taxa that were either partially unidentified
#' (e.g. "sp.") or were considered hybrids (e.g. "A x B") or were listed as more
#' than one species (e.g. "A / B")
filter_species <- function(df, species_table = bbs_all$species) {
  is_unidentified <- function(names) {
    #Before filtering, account for this one hybrid of 2 subspecies so it's kept
    names[names == 'auratus auratus x auratus cafer'] = 'auratus auratus'
    grepl('sp\\.| x |\\/', names)
  }

  valid_taxa <- species_table %>%
    filter(!is_unidentified(species)) %>%
    filter(aou > 2880) %>%
    filter(aou < 3650 | aou > 3810) %>%
    filter(aou < 3900 | aou > 3910) %>%
    filter(aou < 4160 | aou > 4210) %>%
    filter(aou != 7010)

  filter(df, species_id %in% valid_taxa$aou)
}

bbs2 <- bbs %>%
  filter_species() %>%
  combine_subspecies()

# Add 0 abundance for years when a species is not observed
species <- bbs2$species_id %>%
  unique()

bbs3 <- foreach(i = 1:length(species), .combine = "rbind") %do%
  {
    z <- filter(bbs2, species_id == species[i]) %>%
      pivot_wider(names_from = year, values_from = abundance) %>%
      pivot_longer(
        cols = c(-site_id, -species_id, -lat, -long),
        names_to = "year",
        values_to = "abundance"
      )

    z$abundance[is.na(z$abundance)] <- 0

    z
  }
bbs3$year <- as.numeric(bbs3$year)


# Download and process prism data -----------------------------------------

# Downloads the raw prism rasters into the folder specified above
download_prism <- function(years_to_use = 1981:2023) {
  months <- c(1:12)
  clim_vars <- c("ppt", "tmin", "tmean", "tmax")
  for (month in months) {
    for (clim_var in clim_vars) {
      get_prism_monthlys(
        type = clim_var,
        year = years_to_use,
        mon = month,
        keepZip = F
      )
    }
  }
}

get_prism_data <- function(bbs_data) {
  locations <- dplyr::select(bbs_data, site_id, long, lat) %>%
    distinct()
  coordinates(locations) <- c("long", "lat")
  raster::crs(locations) <-
    '+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'

  # Load the prism data and extract using the bbs locations.
  prism_stacked <- pd_stack(prism_archive_ls())

  # Aggregate to 40km. Original cells average ~4.5 across N. America,
  # so aggregation of 9x9 creates mostly 40km cells.
  prism_stacked <- raster::aggregate(prism_stacked, fact = 9)
  extracted <- raster::extract(prism_stacked, locations)
  prism_bbs_data <- data.frame(
    site_id = locations$site_id,
    coordinates(locations),
    extracted
  )
  prism_bbs_data <- prism_bbs_data %>%
    pivot_longer(
      4:ncol(prism_bbs_data),
      names_to = "date",
      values_to = "value",
    ) %>%
    tidyr::extract(
      date,
      c("clim_var", "year", "month"),
      "PRISM_([[:alpha:]]*)_stable_[[:alnum:]]*_([[:digit:]]{4})([[:digit:]]{2})_"
    )

  # Format the data a little
  prism_bbs_data$year <- as.numeric(prism_bbs_data$year)
  prism_bbs_data$month <- as.numeric(prism_bbs_data$month)

  # Now return the data as asked for
  return(prism_bbs_data)
}

#' Helper function to calculate some of the bioclim variables,
#' like "precip in coldest month"
max_min_combo <- function(vec1, vec2, max = TRUE) {
  # Return the value in vec1 in the position where
  # vec2 is either highest or lowest. But 1st check for na
  # values.
  if (any(is.na(vec1)) | any(is.na(vec2))) {
    return(NA)
  } else if (max) {
    return(vec1[which.max(vec2)])
  } else {
    return(vec1[which.min(vec2)])
  }
}

#' From raw monthly climate values calculate all the bioclim variables.
#' You should not call this directly to load bioclim vars. Instead call
#' get_bioclim_data(). Columns present must be
#' c('year','month','value','clim_var','site_id')

process_bioclim_data <- function(monthly_climate_data) {
  #' Offset the year by 6 months so that the window for calculating
  #' bioclim variables will be July 1 - June 30.
  #' See https://github.com/weecology/bbs-forecasting/issues/114
  monthly_climate_data$year <- with(
    monthly_climate_data,
    ifelse(month %in% 7:12, year + 1, year)
  )

  # Spread out the climate variables ppt, tmean, etc into columns
  monthly_climate_data <- monthly_climate_data %>%
    pivot_wider(names_from = clim_var, values_from = value)

  # Process the quarter ones first.
  quarter_info <- data.frame(
    month = 1:12,
    quarter = c(3, 3, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2)
  )
  bioclim_quarter_data <- monthly_climate_data %>%
    left_join(quarter_info, by = 'month') %>%
    group_by(site_id, year, quarter) %>%
    summarise(precip = sum(ppt), temp = mean(tmean)) %>%
    ungroup() %>%
    group_by(site_id, year) %>%
    summarise(
      bio8 = max_min_combo(temp, precip, max = TRUE),
      bio9 = max_min_combo(temp, precip, max = FALSE),
      bio10 = max(temp),
      bio11 = min(temp),
      bio16 = max(precip),
      bio17 = min(precip),
      bio18 = max_min_combo(precip, temp, max = TRUE),
      bio19 = max_min_combo(precip, temp, max = FALSE)
    ) %>%
    ungroup()

  # Next the yearly ones, joining the quartely ones  back in at the end.
  bioclim_data <- monthly_climate_data %>%
    group_by(site_id, year) %>%
    mutate(monthly_temp_diff = tmax - tmin) %>%
    summarise(
      bio1 = mean(tmean),
      bio2 = mean(monthly_temp_diff),
      bio4 = sd(tmean) * 100,
      bio5 = max_min_combo(tmax, tmean, max = TRUE),
      bio6 = max_min_combo(tmin, tmean, max = FALSE),
      bio12 = sum(ppt),
      bio13 = max(ppt),
      bio14 = min(ppt),
      bio15 = 100 * sd(ppt) / mean(ppt)
    ) %>%
    ungroup() %>%
    mutate(bio7 = bio5 - bio6, bio3 = (bio2 / bio7) * 100) %>%
    full_join(bioclim_quarter_data, by = c('site_id', 'year'))

  bioclim_data_lag <- bioclim_data %>%
    mutate(year = year + 1) %>%
    rename_with(~ paste0(.x, "_lag", recycle0 = TRUE), starts_with("bio"))

  bioclim_data <- left_join(
    bioclim_data,
    bioclim_data_lag,
    by = c('site_id', 'year')
  )

  return(bioclim_data)
}

##' Run if prism data is not downloaded
##' download_prism()

options(prism.path = "./data/prismdata")
prism_data <- get_prism_data(bbs2) %>%
  process_bioclim_data()

saveRDS(prism_data, "data/data_prism.rds")

# decompose data for RF
## these covariates vary across space and time and we should decompose them
decomp_vars <- colnames(prism_data)[c(-1, -2)]

prism_data_dec <- prism_data %>%
  # center without scaling
  dplyr::mutate(across(all_of(decomp_vars), ~ scale(., scale = FALSE)[, 1])) %>%

  # compute the spatial component
  # (mean by pixel across all years of centered variables)
  dplyr::group_by(site_id) %>%
  dplyr::mutate(across(
    all_of(decomp_vars),
    ~ mean(., na.rm = TRUE),
    .names = "{.col}_spatial"
  )) %>%
  dplyr::ungroup() %>%
  # compute the temporal component
  # (mean by year across all pixels of centered variables)
  dplyr::group_by(year) %>%
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
  dplyr::arrange(site_id, year)


# Download and process elevation data -------------------------------------

#' Use the raster package to download the average elevation in a 40km radius
#' of each BBS route
get_elev_data <- function(bbs_data) {
  dir.create("data/alt", showWarnings = FALSE)
  elevation <- elevation_30s(country = "US", path = "data/alt")[[1]]
  route_locations <- bbs_data %>%
    distinct(site_id, long, lat)
  route_locations$elevs <-
    terra::extract(
      elevation,
      route_locations[, 2:3],
      buffer = 40000,
      fun = mean
    )[, 2] #buffers are in meters
  #route_locations <- as.data.frame(route_locations)
  elev_data <- dplyr::select(route_locations, site_id, elevs)

  return(elev_data)
}

elev_data <- get_elev_data(bbs2)

# Combine elevation and prism data
env_data <- full_join(prism_data_dec, elev_data, by = c('site_id'))
saveRDS(env_data, "data/data_env.rds")


# Calculate party effort hours --------------------------------------------

# Number of hours on the route multiplied by number of observers
get_ph <- function(start_hr, end_hr) {
  x <- (end_hr - start_hr) / 100
  ((x - round(x)) * 100 + 60 * round(x))
}

bbs_eff <- bbs_all$routes %>%
  rowwise() %>%
  mutate(
    hours = get_ph(start_time, end_time),
    site_id = state_num * 1000 + route
  ) %>%
  group_by(site_id, year) %>%
  summarize(party_hours = sum(hours))

# add strata as bcr x state combination
bbs_str <- bbs_all$routes %>%
  rowwise() %>%
  mutate(
    site_id = state_num * 1000 + route,
    strata = paste(state, bcr, sep = "-")
  ) %>%
  dplyr::select(site_id, year, strata)


# Process global climate data ---------------------------------------------

# Atlantic multidecadal oscillation
amo <- read.table("data/amo.txt", sep = "", skip = 1, header = F) %>%
  filter(V1 >= 1981 & V1 <= 2023)
colnames(amo) <- c(
  "year",
  "J",
  "F",
  "M",
  "A",
  "M",
  "J",
  "J",
  "A",
  "S",
  "O",
  "N",
  "D"
)
amo <- as.matrix(amo)
amo_sp <- foreach(i = 1:(nrow(amo) - 1), .combine = "rbind") %do%
  {
    data.frame(
      year = amo[i + 1, 1],
      amo = round(mean(amo[i + 1, 4:6]), 3),
      amo_lag = round(mean(amo[i, 4:6]), 3)
    )
  }

# El-Nino southern oscillation
enso <- read.csv("data/censo.csv")
enso$year <- as.numeric(str_split_i(enso$Date, "-", 1))
enso$month <- as.numeric(str_split_i(enso$Date, "-", 2))
enso <- enso[, -1]
colnames(enso)[1] <- "censo"
enso <- pivot_wider(enso, names_from = "month", values_from = "censo") %>%
  filter(year >= 1981 & year <= 2023)
colnames(enso) <- c(
  "year",
  "J",
  "F",
  "M",
  "A",
  "M",
  "J",
  "J",
  "A",
  "S",
  "O",
  "N",
  "D"
)

enso <- as.matrix(enso)
enso_sp <- foreach(i = 1:(nrow(enso) - 1), .combine = "rbind") %do%
  {
    data.frame(
      year = enso[i + 1, 1],
      enso = round(mean(enso[i + 1, 4:6]), 3),
      enso_lag = round(mean(enso[i, 4:6]), 3)
    )
  }

# Pacific decadal oscillation
pdo <- read.table("data/pdo.txt", sep = "", skip = 1, header = T) %>%
  filter(Year >= 1981 & Year <= 2023)

pdo <- as.matrix(pdo)
pdo_sp <- foreach(i = 1:(nrow(pdo) - 1), .combine = "rbind") %do%
  {
    data.frame(
      year = pdo[i + 1, 1],
      pdo = round(mean(pdo[i + 1, 4:6]), 3),
      pdo_lag = round(mean(pdo[i, 4:6]), 3)
    )
  }

# North Atlantic oscillation
nao <- read.table("data/nao.txt", header = T)
nao$year <- as.numeric(rownames(nao))
nao <- filter(nao, year >= 1981 & year <= 2023)

nao <- as.matrix(nao)
nao_sp <- foreach(i = 1:(nrow(nao) - 1), .combine = "rbind") %do%
  {
    data.frame(
      year = nao[i + 1, 13],
      nao = round(mean(nao[i + 1, 3:5]), 3),
      nao_lag = round(mean(nao[i, 3:5]), 3)
    )
  }


# Combine all data --------------------------------------------------------

# remove large objects to ease memory
rm(bbs_all)
rm(bbs)
rm(prism_data_dec)
gc()

bbs_final <-
  left_join(bbs3, bbs_eff, by = c("site_id", "year")) %>%
  left_join(bbs_str, by = c("site_id", "year")) %>%
  filter(year >= 1982 & year <= 2023) %>%
  left_join(env_data, by = c("site_id", "year"), copy = T) %>%
  filter(!is.na(bio1_lag), !is.na(elevs)) %>%
  left_join(amo_sp, by = "year") %>%
  left_join(enso_sp, by = "year") %>%
  left_join(pdo_sp, by = "year") %>%
  left_join(nao_sp, by = "year")

bbs_final_nozero <-
  left_join(bbs2, bbs_eff, by = c("site_id", "year")) %>%
  left_join(bbs_str, by = c("site_id", "year")) %>%
  filter(year >= 1982 & year <= 2023) %>%
  left_join(env_data, by = c("site_id", "year"), copy = T) %>%
  filter(!is.na(bio1_lag), !is.na(elevs)) %>%
  left_join(amo_sp, by = "year") %>%
  left_join(enso_sp, by = "year") %>%
  left_join(pdo_sp, by = "year") %>%
  left_join(nao_sp, by = "year")


saveRDS(bbs_final, "data/data_bbs.rds")
saveRDS(bbs_final_nozero, "data/data_bbs_nozero.rds")
