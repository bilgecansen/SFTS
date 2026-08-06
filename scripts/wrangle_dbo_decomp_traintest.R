# Honest train/test decomposition for DBO population-level (mirrors the BBS
# wrangle_decomp_blocks.R). Train 2001-2010, test 2011-2019 are decomposed SEPARATELY
# within their own years, so the spatial climatology and year-marginal never cross the
# train/test boundary (the global decomposition in wrangle_dbo_data.R let test years
# into the climatology). Same math: center by grand mean, spatial = station-marginal
# mean, temporal = year-marginal mean, residual = the rest. RAW env is additionally
# centered by the TRAIN grand mean so the static test projection (raw fed into the
# _spatial slots) is on the train-centered scale.
#
# Env in data/data_dbo.rds is already global-centered; that is irrelevant here because
# every component is re-centered per set (centering-invariant). Output ->
# data/decomp_dbo_traintest.rds (StationNme, DataYear, set, Depth, Lat/Long, raw base,
# base_spatial/temporal/residual).
library(tidyverse)
base <- c("Temp","Salinity","integchla","sedchla","Ammonia","Phosphate",
          "NiTriTra","Silicate","phigte5","TOC","cn")

env <- readRDS("data/data_dbo.rds") %>%
  select(StationNme, DataYear, Depth, Latitude, Longitude, all_of(base)) %>%
  distinct(StationNme, DataYear, .keep_all = TRUE)

decompose_years <- function(clim, yrs) {
  x <- clim %>% filter(DataYear %in% yrs)
  cc <- paste0(base, "__c")
  x <- x %>% mutate(across(all_of(base), ~ .x - mean(.x, na.rm = TRUE), .names = "{.col}__c"))
  x <- x %>% group_by(StationNme) %>% mutate(across(all_of(cc), ~ mean(.x, na.rm = TRUE), .names = "{.col}__sp")) %>% ungroup()
  x <- x %>% group_by(DataYear)   %>% mutate(across(all_of(cc), ~ mean(.x, na.rm = TRUE), .names = "{.col}__tm")) %>% ungroup()
  out <- x %>% select(StationNme, DataYear, Depth, Latitude, Longitude, all_of(base))
  for (v in base) {
    ct <- x[[paste0(v,"__c")]]; sp <- x[[paste0(v,"__c__sp")]]; tm <- x[[paste0(v,"__c__tm")]]
    out[[paste0(v,"_spatial")]]  <- sp
    out[[paste0(v,"_temporal")]] <- tm
    out[[paste0(v,"_residual")]] <- ct - sp - tm
  }
  out
}

build <- function(train_yrs, test_yrs) {
  tr <- decompose_years(env, train_yrs) %>% mutate(set = "train")
  te <- decompose_years(env, test_yrs)  %>% mutate(set = "test")
  g  <- env %>% filter(DataYear %in% train_yrs) %>% summarise(across(all_of(base), ~ mean(.x, na.rm = TRUE)))
  for (v in base) { tr[[v]] <- tr[[v]] - g[[v]]; te[[v]] <- te[[v]] - g[[v]] }
  bind_rows(tr, te)
}

saveRDS(build(2001:2010, 2011:2019), "data/decomp_dbo_traintest.rds")
cat("saved -> data/decomp_dbo_traintest.rds\n")
