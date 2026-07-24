# Per-block spatial/temporal/residual decomposition (CANONICAL).
#
# Replaces the single global decomposition that used to live in wrangle_bbs_data.R
# (which decomposed once over the whole 1982-2023 record, test years included).
# Here each temporal block is decomposed WITHIN its own years, and TRAIN and TEST
# are decomposed SEPARATELY, so no information crosses the train/test boundary:
#   * train rows: components fitted on the block's training years only
#   * test  rows: components fitted on the test years (2011-2020) only
# Same math as before -- center by grand mean, spatial = site-marginal mean,
# temporal = year-marginal mean, residual = centered - spatial - temporal -- just
# applied per set of years. The RAW columns are additionally centered by the TRAIN
# grand mean (both sets) so the static test projection (raw fed into the _spatial
# slots) is on the same centered scale as the train spatial component.
#
# Reads data/data_prism.rds (raw climate, site x year); the base BBS wrangling
# (wrangle_bbs_data.R -> data/data_bbs_nozero.rds) is unchanged and still supplies
# the non-climate columns to the fit scripts.
#
# Blocks:  Temporal train 2001-2010 ; Buffer train 1991-2000 ; test 2011-2020.
# (Spatiotemporal reuses the Temporal block decomposition + spatial blocking.)

library(tidyverse)

# 8 model vars (RF + GAM dynamic/decomp/svc) plus bio1,bio12 (used by the 2-var SI)
dvars <- c("bio1","bio12","bio2","bio3","bio5","bio8","bio9","bio15","bio16","bio18")
prism <- readRDS("data/data_prism.rds") %>% select(site_id, year, all_of(dvars))

decompose_years <- function(clim, yrs) {
  d <- clim %>% filter(year %in% yrs)
  cc <- paste0(dvars, "__c")
  d <- d %>% mutate(across(all_of(dvars), ~ .x - mean(.x, na.rm = TRUE), .names = "{.col}__c"))
  d <- d %>% group_by(site_id) %>% mutate(across(all_of(cc), ~ mean(.x, na.rm = TRUE), .names = "{.col}__sp")) %>% ungroup()
  d <- d %>% group_by(year)    %>% mutate(across(all_of(cc), ~ mean(.x, na.rm = TRUE), .names = "{.col}__tm")) %>% ungroup()
  out <- d %>% select(site_id, year, all_of(dvars))
  for (v in dvars) {
    ct <- d[[paste0(v,"__c")]]; sp <- d[[paste0(v,"__c__sp")]]; tm <- d[[paste0(v,"__c__tm")]]
    out[[paste0(v,"_spatial")]]  <- sp
    out[[paste0(v,"_temporal")]] <- tm
    out[[paste0(v,"_residual")]] <- ct - sp - tm
  }
  out
}

build_block <- function(train_yrs, test_yrs = 2011:2020) {
  tr <- decompose_years(prism, train_yrs) %>% mutate(set = "train")
  te <- decompose_years(prism, test_yrs)  %>% mutate(set = "test")
  gtrain <- prism %>% filter(year %in% train_yrs) %>% summarise(across(all_of(dvars), ~ mean(.x, na.rm = TRUE)))
  for (v in dvars) { tr[[v]] <- tr[[v]] - gtrain[[v]]; te[[v]] <- te[[v]] - gtrain[[v]] }
  bind_rows(tr, te)
}

saveRDS(build_block(2001:2010), "data/decomp_temporal.rds")
saveRDS(build_block(1991:2000), "data/decomp_buffer.rds")
cat("saved -> data/decomp_temporal.rds, data/decomp_buffer.rds\n")
