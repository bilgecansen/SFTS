# DBO RF on RAW (untransformed) biomass WITH ZEROS RETAINED, canonical train/test split.
# Companion to fit_rf_dbo_traintest.R, which models log(biomass) with zeros dropped.
#
#   data   : data/data_dbo_zeros.rds  (complete panel over every sampled station-year;
#            zeros are observations of absence and never counted toward sample size)
#   split  : train 2001-2010 / test 2011-2019, honest per-set decomposition
#            (data/decomp_dbo_traintest.rds -- train and test decomposed independently)
#   models : static / dynamic / decomposed / SVC, Depth as control
#   response: biomass, RAW scale, zeros included
#
# Model definitions (response always spatiotemporal; only the ENV representation differs):
#   static     : TRAIN on per-station means (*_spatial); PREDICT with spatiotemporal env fed
#                into those slots. Abundance-distribution level uses per-station MEAN env.
#   dynamic    : TRAIN and TEST both raw spatiotemporal env.
#   decomposed : TRAIN and TEST both decomposed (spatial + temporal + residual).
#   SVC        : decomposed + Latitude + Longitude.
#
# Filters match fit_rf_dbo_traintest.R, with the rule that ZEROS DO NOT COUNT toward
# sample size: min_years is applied to POSITIVE test records per station, and a station
# counts toward min_stn only if the family was actually observed there in that period.
#
# Prediction levels (new naming):
#   abundance distribution (_ad) : station means over the test years, correlated across stations
#   regional trajectory    (_rt) : stations averaged within a DBO region per year, correlated over years
#   station trajectory     (_st) : within-station correlation over years, median across stations
#
# THE REGIONAL AVERAGE IS COMPUTED TWO WAYS (the open question):
#   _incl : average over ALL stations sampled in that region-year, absences entering as 0.
#           This is regional biomass density -- it tracks occupancy AND local abundance, and
#           is what the model actually predicts (an unconditional expected biomass).
#   _excl : average over only the stations where the family was observed. This is density
#           GIVEN presence. Predictions are averaged over the SAME stations, so any
#           sampling-composition shift is shared by observed and predicted and largely
#           cancels in the correlation -- but the quantity still changes composition
#           year to year, which the _incl version does not.
#
# Output: data/rf_dbo_traintest_zeros_results.rds

library(tidyverse)
library(ranger)

set.seed(1)
n_trees <- 2000
min_stn <- 10; min_test_stn <- 5; min_years <- 5
min_reg_stn <- 3; min_reg_years <- 5

base <- c("Temp", "Salinity", "integchla", "sedchla", "Ammonia", "Phosphate",
          "NiTriTra", "Silicate", "phigte5", "TOC", "cn")
comp_terms <- as.vector(t(outer(base, c("spatial", "temporal", "residual"), paste, sep = "_")))
ctrl <- "Depth"

dec <- readRDS("data/decomp_dbo_traintest.rds")
bio <- readRDS("data/data_dbo_zeros.rds") %>%
  distinct(family, StationNme, DataYear, DBOreg, biomass)
dat <- inner_join(bio, dec, by = c("StationNme", "DataYear"))
train <- filter(dat, set == "train"); test <- filter(dat, set == "test")

f_static  <- as.formula(paste("biomass ~", paste(c(ctrl, paste0(base, "_spatial")), collapse = " + ")))
f_dynamic <- as.formula(paste("biomass ~", paste(c(ctrl, base), collapse = " + ")))
f_decomp  <- as.formula(paste("biomass ~", paste(c(ctrl, comp_terms), collapse = " + ")))
f_svc     <- as.formula(paste("biomass ~", paste(c(ctrl, comp_terms, "Latitude", "Longitude"), collapse = " + ")))

rf <- function(f, d) ranger(f, data = d, num.trees = n_trees, seed = 1, num.threads = 0)
pr <- function(m, newd) as.numeric(predict(m, data = as.data.frame(newd))$predictions)

adcor <- function(stn, y, pred) {                       # abundance distribution
  q <- tibble(stn, y, pred) %>% group_by(stn) %>%
    summarise(o = mean(y), p = mean(pred), .groups = "drop")
  suppressWarnings(cor(q$o, q$p))
}
stcor <- function(stn, y, pred) tibble(stn, y, pred) %>%   # station trajectory
  group_by(stn) %>% filter(sum(y > 0) >= min_years) %>%    # zeros don't count toward n
  summarise(r = suppressWarnings(cor(y, pred)), .groups = "drop") %>%
  pull(r) %>% median(na.rm = TRUE)
rtcor <- function(reg, year, y, pred, include_zeros) {     # regional trajectory
  d <- tibble(reg, year, y, pred)
  if (!include_zeros) d <- filter(d, y > 0)                # same stations for obs AND pred
  d %>% group_by(reg, year) %>% filter(n() >= 2) %>%
    summarise(o = mean(y), p = mean(pred), .groups = "drop") %>%
    group_by(reg) %>% filter(n() >= min_reg_years) %>%
    summarise(r = suppressWarnings(cor(o, p)), .groups = "drop") %>%
    pull(r) %>% median(na.rm = TRUE)
}

families <- unique(train$family)
out <- map(families, function(fm) {
  tr <- filter(train, family == fm)
  te <- filter(test, family == fm) %>%
    filter(StationNme %in% unique(tr$StationNme[tr$biomass > 0])) %>%
    group_by(StationNme) %>% filter(sum(biomass > 0) >= min_years) %>% ungroup()
  if (n_distinct(tr$StationNme[tr$biomass > 0]) < min_stn ||
      n_distinct(te$StationNme) < min_test_stn) return(NULL)
  y <- te$biomass

  te_static <- te; te_static[paste0(base, "_spatial")] <- te[base]
  tem <- te %>% group_by(StationNme) %>%
    summarise(o = mean(biomass), across(all_of(c(base, ctrl)), mean), .groups = "drop")
  te_avg <- tem; te_avg[paste0(base, "_spatial")] <- tem[base]

  tryCatch({
    m_s <- rf(f_static, tr)
    p <- list(static = pr(m_s, te_static), dynamic = pr(rf(f_dynamic, tr), te),
              decomp = pr(rf(f_decomp, tr), te), svc = pr(rf(f_svc, tr), te))
    p_s_ad <- pr(m_s, te_avg)
    o <- tibble(family = fm, n_stn = n_distinct(te$StationNme), n_test = nrow(te),
                pct_zero = mean(te$biomass == 0))
    o$static_ad <- suppressWarnings(cor(tem$o, p_s_ad))
    for (m in names(p)) {
      if (m != "static") o[[paste0(m, "_ad")]] <- adcor(te$StationNme, y, p[[m]])
      o[[paste0(m, "_rt_incl")]] <- rtcor(te$DBOreg, te$DataYear, y, p[[m]], TRUE)
      o[[paste0(m, "_rt_excl")]] <- rtcor(te$DBOreg, te$DataYear, y, p[[m]], FALSE)
      o[[paste0(m, "_st")]]      <- stcor(te$StationNme, y, p[[m]])
    }
    # station-year predictions, retained so the error decomposition can be run
    # without refitting (see scripts/skill_dbo_zeros.R)
    pd <- tibble(family = fm, DBOreg = te$DBOreg, stn = te$StationNme,
                 year = te$DataYear, y = y,
                 p_static = p$static, p_dynamic = p$dynamic,
                 p_decomp = p$decomp, p_svc = p$svc)
    list(summary = o, preds = pd)
  }, error = function(e) NULL)
})
out <- out[!map_lgl(out, is.null)]
res <- map_dfr(out, "summary")
saveRDS(map_dfr(out, "preds"), "data/rf_dbo_traintest_zeros_preds.rds")
saveRDS(res, "data/rf_dbo_traintest_zeros_results.rds")

q <- function(x) sprintf("%+.2f (%+.2f,%+.2f)", median(x, na.rm = TRUE),
  quantile(x, .05, na.rm = TRUE), quantile(x, .95, na.rm = TRUE))
cat(sprintf("\n=== DBO RF, RAW biomass + zeros, train 2001-2010 / test 2011-2019 (%d families) ===\n", nrow(res)))
cat(sprintf("median %% zeros in test rows: %.0f%%\n\n", 100 * median(res$pct_zero)))
cat(sprintf("%-9s %-22s %-22s %-22s %-22s\n", "model",
  "abundance distrib.", "regional (incl 0s)", "regional (excl 0s)", "station trajectory"))
for (m in c("static", "dynamic", "decomp", "svc"))
  cat(sprintf("%-9s %-22s %-22s %-22s %-22s\n", m,
    q(res[[paste0(m, "_ad")]]), q(res[[paste0(m, "_rt_incl")]]),
    q(res[[paste0(m, "_rt_excl")]]), q(res[[paste0(m, "_st")]])))
cat("\nsaved -> data/rf_dbo_traintest_zeros_results.rds\n")
