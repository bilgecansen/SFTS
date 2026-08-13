# Extract the daily NSIDC sea ice concentration series at each DBO station, once.
#
# Rewritten from the version that produced data/data_sic_per.rds. That script counted
# ice days over the CALENDAR year, which for a July cruise merges the tail of the winter
# preceding the sample with the freeze-up of the following winter, in November and
# December, after the biology was collected. It also wrote data/data_sic.rds from a
# monthly block that referenced an undefined object (`file_ice`) and pasted paths without
# a separator, which is why that file is an unlabelled vector with no join key.
#
# This script does one thing: pull the daily value at every station for every day of the
# record and save it. All derived ice variables are built from the saved series by
# build_ice_vars.R, so the rasters are only ever read once.
#
# DATA. data/Sea_ice/daily/{year}/{MM_Mon}/N_YYYYMMDD_concentration_v3.0.tif
#   25 km NSIDC polar stereographic north, EPSG:3411, values scaled by 10 so 0-1000 is
#   0-100%. Flags above 1000 (2510 pole hole, 2530 coast, 2540 land, 2550 missing) are set
#   to NA -- none occur at these 23 stations, but a bare `> 150` test would have counted
#   them as ice. Station coordinates are transformed into the raster's own CRS; the earlier
#   use of EPSG:3413 displaces them by at most 65 m, which is immaterial on a 25 km grid.
#
# Output: data/sic_daily.rds   long, one row per station-day: StationNme, date, sic

suppressMessages({ library(tidyverse); library(readxl); library(terra); library(sf) })

ROOT <- "data/Sea_ice/daily"

# ---- stations ---------------------------------------------------------------------------
e1 <- read_excel("data/EnvtlVariables_UniqueID_0-5nm_DBO_063023_v3_abwk122023.xlsx",
                 sheet = 4)
names(e1)[3] <- "StationNme"
e2 <- read_excel("data/EnvironmentalVariables_UTN_2000-2019.xlsx", sheet = 2)
names(e2)[3] <- "StationNme"

keep <- unique(readRDS("data/data_sic_per.rds")$StationNme)   # the 23 already in use
coord <- bind_rows(select(e1, StationNme, Longitude, Latitude),
                   select(e2, StationNme, Longitude, Latitude)) %>%
  filter(StationNme %in% keep) %>%
  group_by(StationNme) %>%
  summarise(lon = mean(as.numeric(Longitude)), lat = mean(as.numeric(Latitude)),
            .groups = "drop")
stopifnot(setequal(coord$StationNme, keep))
cat(sprintf("%d stations\n", nrow(coord)))

# ---- extract ----------------------------------------------------------------------------
years <- sort(as.integer(list.files(ROOT)))
cat(sprintf("years %d-%d\n", min(years), max(years)))

xy <- NULL                                   # built once, from the first raster's CRS
out <- map_dfr(years, function(yr) {
  map_dfr(list.files(file.path(ROOT, yr)), function(mo) {
    f <- list.files(file.path(ROOT, yr, mo), pattern = "concentration",
                    full.names = TRUE)
    if (!length(f)) return(NULL)
    r <- rast(f)
    if (is.null(xy)) {
      xy <<- st_as_sf(coord, coords = c("lon", "lat"), crs = 4326) %>%
        st_transform(crs(r)) %>% st_coordinates()
    }
    v <- terra::extract(r, xy)
    v <- v[, setdiff(names(v), "ID"), drop = FALSE]
    dates <- as.Date(str_extract(basename(f), "\\d{8}"), format = "%Y%m%d")
    tibble(StationNme = rep(coord$StationNme, times = length(dates)),
           date = rep(dates, each = nrow(coord)),
           sic = as.numeric(as.matrix(v)))
  })
})

# flags: anything above full concentration is not a measurement
out <- out %>% mutate(sic = ifelse(sic > 1000, NA_real_, sic))

cat(sprintf("\n%d station-days, %d stations, %s to %s\n", nrow(out),
            n_distinct(out$StationNme), min(out$date), max(out$date)))
cat(sprintf("missing values: %d (%.2f%%)\n", sum(is.na(out$sic)),
            100 * mean(is.na(out$sic))))
gaps <- setdiff(seq(min(out$date), max(out$date), by = "day"), unique(out$date))
cat(sprintf("calendar days with no file at all: %d\n", length(gaps)))
if (length(gaps)) cat("  ", paste(head(as.Date(gaps, origin = "1970-01-01"), 12),
                                  collapse = ", "), "\n")

saveRDS(out, "data/sic_daily.rds")
cat("saved -> data/sic_daily.rds\n")
