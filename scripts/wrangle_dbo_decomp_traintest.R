# Train/test decomposition of the DBO environment, each block decomposed within its
# own years so nothing crosses the boundary.
#
# SPLIT BY YEAR COUNT, NOT YEAR RANGE. The 15 sampled years are 2001, 2003-2007 and
# 2011-2019; a calendar-range split at 2010 therefore gives 6 train years against 9
# test years. Taking the first 8 sampled years and the last 7 gives
#     train  2001, 2003-2007, 2011, 2012      (115 station-years)
#     test   2013-2019                        (111 station-years)
# which is nearly balanced, costs no families at the gates used downstream, and
# leaves ~78 residual df for the four-way decomposition in each block (the 6/9
# split left ~53 in train). Test years now begin immediately after the last training
# year, so prediction lead time is 1-7 years rather than 4-12.
#
# FOUR COMPONENTS, by sequential orthogonal projection within each block:
#     station  ->  year  ->  region:year  ->  residual
#   -> _spatial / _global / _regional / _residual, summing exactly to the
# block-centered value. Same construction as scripts/wrangle_dbo_data.R (which does
# it globally, for LOSO) and the same term order as the error decomposition in
# scripts/skill_dbo_4rt.R. Station is nested in region, so station is entered first
# and the time-invariant latitudinal gradient lands in _spatial.
#
# RAW env is additionally centered by the TRAIN grand mean, so the static model's
# test projection (raw values fed into the _spatial slots) is on the train-centered
# scale.
#
# Env in data/data_dbo.rds is already globally centered; irrelevant here, since every
# component is re-centered per block.
#
# Output -> data/decomp_dbo_traintest.rds
#   StationNme, DataYear, DBOreg, set, Depth, Latitude, Longitude,
#   raw base, base_spatial / _global / _regional / _residual
library(tidyverse)

N_TRAIN <- 8                       # number of sampled years in the training block

base <- c("Temp","Salinity","integchla","sedchla","Ammonia","Phosphate",
          "NiTriTra","Silicate","phigte5","TOC","cn")

env <- readRDS("data/data_dbo.rds") %>%
  select(StationNme, DataYear, DBOreg, Depth, Latitude, Longitude, all_of(base)) %>%
  distinct(StationNme, DataYear, .keep_all = TRUE)

yrs <- sort(unique(env$DataYear))
stopifnot(N_TRAIN < length(yrs))
train_yrs <- yrs[seq_len(N_TRAIN)]
test_yrs  <- yrs[-seq_len(N_TRAIN)]

decompose_block <- function(x) {
  stn <- factor(x$StationNme)
  yr  <- factor(x$DataYear)
  ry  <- factor(paste(x$DBOreg, x$DataYear))
  out <- x
  for (v in base) {
    cx <- x[[v]] - mean(x[[v]], na.rm = TRUE)
    f1 <- fitted(lm(cx ~ stn, na.action = na.exclude))
    f2 <- fitted(lm(cx ~ stn + yr, na.action = na.exclude))
    f3 <- fitted(lm(cx ~ stn + yr + ry, na.action = na.exclude))
    out[[paste0(v, "_spatial")]]  <- as.numeric(f1)
    out[[paste0(v, "_global")]]   <- as.numeric(f2 - f1)
    out[[paste0(v, "_regional")]] <- as.numeric(f3 - f2)
    out[[paste0(v, "_residual")]] <- as.numeric(cx - f3)
  }
  out
}

tr <- decompose_block(filter(env, DataYear %in% train_yrs)) %>% mutate(set = "train")
te <- decompose_block(filter(env, DataYear %in% test_yrs))  %>% mutate(set = "test")

g <- env %>% filter(DataYear %in% train_yrs) %>%
  summarise(across(all_of(base), ~ mean(.x, na.rm = TRUE)))
for (v in base) { tr[[v]] <- tr[[v]] - g[[v]]; te[[v]] <- te[[v]] - g[[v]] }

out <- bind_rows(tr, te)
saveRDS(out, "data/decomp_dbo_traintest.rds")

cat(sprintf("train %d years (%s): %d station-years\n", length(train_yrs),
            paste(range(train_yrs), collapse = "-"), nrow(tr)))
cat(sprintf("test  %d years (%s): %d station-years\n", length(test_yrs),
            paste(range(test_yrs), collapse = "-"), nrow(te)))
cat(sprintf("max additivity error: %.1e\n",
    max(map_dbl(base, function(v) {
      s <- out[[paste0(v,"_spatial")]] + out[[paste0(v,"_global")]] +
           out[[paste0(v,"_regional")]] + out[[paste0(v,"_residual")]]
      cx <- ave(out[[v]], out$set, FUN = function(z) z - mean(z))
      max(abs(cx - s))
    }))))
cat("saved -> data/decomp_dbo_traintest.rds\n")
