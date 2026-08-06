# DBO LINEAR (OLS) train/test analysis -- simple-model companion to fit_rf_dbo_traintest.R.
# Same split, same folds, same metrics, same family filters; only the learner and the
# predictor set change, so the linear and RF results are directly comparable.
#
#   split      : train 2001-2010, test 2011-2019 (canonical), honest per-set decomposition
#                (data/decomp_dbo_traintest.rds -- train and test decomposed independently
#                within their own years)
#   predictors : Temp, Salinity, integchla, phigte5 ONLY. No Depth control, no coordinates.
#   models     : static / dynamic / decomposed  (no SVC)
#
# Model definitions (response = log biomass, always spatiotemporal; only the ENV
# representation differs):
#   static     : TRAIN on per-station means (*_spatial). PREDICT with spatiotemporal env
#                fed into those same slots -> a spatial relationship used to predict across
#                time (the space-for-time test). Species-wide uses per-station MEAN env,
#                one prediction per station, matching the RF script.
#   dynamic    : TRAIN and TEST both raw spatiotemporal env.
#   decomposed : TRAIN and TEST both decomposed (spatial + temporal + residual) = 12 terms.
#
# Prediction levels:
#   species-wide (_sw) : across-station, station-mean observed vs predicted
#   regional     (_rg) : average obs & pred across stations within a DBO region per year,
#                        correlate the regional trajectory over years, median across the
#                        family's regions
#   population   (_pl) : within-station correlation over years, median across stations
#
# Output: data/lin_dbo_traintest_results.rds

library(tidyverse)

min_stn <- 10; min_test_stn <- 5; min_years <- 5
base <- c("Temp", "Salinity", "integchla", "phigte5")
comp_terms <- as.vector(t(outer(base, c("spatial", "temporal", "residual"), paste, sep = "_")))

dec <- readRDS("data/decomp_dbo_traintest.rds")
bio <- readRDS("data/data_dbo.rds") %>%
  filter(biomass > 0) %>%
  distinct(family, StationNme, DataYear, DBOreg, biomass)
dat <- inner_join(bio, dec, by = c("StationNme", "DataYear"))
train <- filter(dat, set == "train"); test <- filter(dat, set == "test")

f_static  <- as.formula(paste("log(biomass) ~", paste(paste0(base, "_spatial"), collapse = " + ")))
f_dynamic <- as.formula(paste("log(biomass) ~", paste(base, collapse = " + ")))
f_decomp  <- as.formula(paste("log(biomass) ~", paste(comp_terms, collapse = " + ")))

lmfit <- function(f, d) lm(f, data = d)
pr <- function(m, newd) as.numeric(predict(m, newdata = as.data.frame(newd)))

plcor <- function(stn, y, pred) tibble(stn, y, pred) %>%
  group_by(stn) %>% filter(n() >= min_years) %>%
  summarise(r = suppressWarnings(cor(y, pred)), .groups = "drop") %>%
  pull(r) %>% median(na.rm = TRUE)

swcor <- function(stn, y, pred) {
  q <- tibble(stn, y, pred) %>% group_by(stn) %>%
    summarise(o = mean(y), p = mean(pred), .groups = "drop")
  suppressWarnings(cor(q$o, q$p))
}

rgcor <- function(reg, year, y, pred) tibble(reg, year, y, pred) %>%
  group_by(reg, year) %>% filter(n() >= 2) %>%              # >= 2 stations per regional mean
  summarise(o = mean(y), p = mean(pred), .groups = "drop") %>%
  group_by(reg) %>% filter(n() >= min_years) %>%            # >= 5 region-years to correlate
  summarise(r = suppressWarnings(cor(o, p)), .groups = "drop") %>%
  pull(r) %>% median(na.rm = TRUE)

families <- unique(train$family)
res <- map_dfr(families, function(fm) {
  tr <- filter(train, family == fm)
  te <- filter(test, family == fm) %>%
    filter(StationNme %in% unique(tr$StationNme)) %>%
    group_by(StationNme) %>% filter(n() >= min_years) %>% ungroup()
  if (n_distinct(tr$StationNme) < min_stn || n_distinct(te$StationNme) < min_test_stn) return(NULL)
  y <- log(te$biomass)

  te_static <- te; te_static[paste0(base, "_spatial")] <- te[base]   # raw yearly -> pop-level
  # static species-wide: per-station MEAN env -> one prediction per station
  tem <- te %>% group_by(StationNme) %>%
    summarise(o = mean(log(biomass)), across(all_of(base), mean), .groups = "drop")
  te_avg <- tem; te_avg[paste0(base, "_spatial")] <- tem[base]

  tryCatch({
    m_s <- lmfit(f_static, tr)
    p_s <- pr(m_s, te_static); p_s_sw <- pr(m_s, te_avg)
    p_d <- pr(lmfit(f_dynamic, tr), te)
    p_c <- pr(lmfit(f_decomp, tr), te)
    tibble(family = fm, n_stn = n_distinct(te$StationNme),
      static_sw  = suppressWarnings(cor(tem$o, p_s_sw)),
      static_rg  = rgcor(te$DBOreg, te$DataYear, y, p_s),
      static_pl  = plcor(te$StationNme, y, p_s),
      dynamic_sw = swcor(te$StationNme, y, p_d),
      dynamic_rg = rgcor(te$DBOreg, te$DataYear, y, p_d),
      dynamic_pl = plcor(te$StationNme, y, p_d),
      decomp_sw  = swcor(te$StationNme, y, p_c),
      decomp_rg  = rgcor(te$DBOreg, te$DataYear, y, p_c),
      decomp_pl  = plcor(te$StationNme, y, p_c))
  }, error = function(e) NULL)
})
saveRDS(res, "data/lin_dbo_traintest_results.rds")

q <- function(x) sprintf("%+.2f (%+.2f, %+.2f)", median(x, na.rm = TRUE),
  quantile(x, .05, na.rm = TRUE), quantile(x, .95, na.rm = TRUE))
cat(sprintf("\n=== DBO LINEAR, train 2001-2010 / test 2011-2019 (%d families) ===\n", nrow(res)))
cat(sprintf("predictors: %s | no Depth, no SVC\n\n", paste(base, collapse = ", ")))
cat(sprintf("%-9s %-24s %-24s %-24s\n", "model", "species-wide", "regional", "population-level"))
for (m in c("static", "dynamic", "decomp"))
  cat(sprintf("%-9s %-24s %-24s %-24s\n", m,
    q(res[[paste0(m, "_sw")]]), q(res[[paste0(m, "_rg")]]), q(res[[paste0(m, "_pl")]])))
cat("\nsaved -> data/lin_dbo_traintest_results.rds\n")
