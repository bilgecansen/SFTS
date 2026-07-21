# DBO random forests, ADAPTIVE PER-STATION TEMPORAL SPLIT (species-wide only).
#
# Third species-wide holdout for DBO, completing the 2x2 of what is held out:
#            temporal in-sample        temporal out-of-sample
#   spat in  LOYO  (fit_rf_dbo_3models) THIS  (train-early / test-late)
#   spat out LOSO  (fit_rf_dbo_loso)    (joint station x year -- not run)
# LOYO interpolates on both axes (species-wide near ceiling); LOSO removes the
# spatial crutch; this removes the TEMPORAL crutch while keeping every station
# in-sample -- "can a model fit on early years reproduce the spatial gradient in
# a later period?" The DBO analog of the BBS Temporal treatment.
#
# Split is adaptive per station because DBO series are short and unevenly spread
# across calendar time: a single fixed cutoff year would drop stations sampled
# mostly on one side. For each (family, station) with n positive-biomass years we
# take the LAST ceiling(n/3) years (>= 2) as test and the earlier years (>= 3) as
# train. This keeps the same 45 families (>= 10 stations) as LOYO/LOSO, so the
# three species-wide numbers are directly comparable.
#
# Species-wide ONLY: the per-station test window is too short (median 3 yrs) for a
# within-station year-to-year correlation, and LOYO already gives the honest
# population-level. Three models (static/dynamic/decomposed), Depth control.
#
# NOTE: _spatial is a precomputed full-period station climatology, so the static
# model's spatial predictor carries slight test-period information -- identical to
# how LOYO, LOSO and the BBS scripts treat it, so the comparison stays fair.

library(tidyverse)
library(ranger)

set.seed(1)
NTREE <- 2000
min_stn <- 10       # family needs >= this many qualifying stations
min_train_yr <- 3   # station needs >= this many training (early) years
min_test_yr <- 2    # station needs >= this many test (late) years

d <- readRDS("data/data_dbo.rds") %>% filter(biomass > 0)

base <- c(
  "Temp", "Salinity", "integchla", "sedchla", "Ammonia", "Phosphate",
  "NiTriTra", "Silicate", "phigte5", "TOC", "cn"
)
comp_terms <- as.vector(t(outer(base,
  c("spatial", "temporal", "residual"), paste, sep = "_")))
ctrl <- "Depth"

swcor <- function(stn, y, pred) {
  s <- tibble(stn, y, pred) %>%
    group_by(stn) %>%
    summarise(o = mean(y), p = mean(pred), .groups = "drop")
  suppressWarnings(cor(s$o, s$p))
}

fit_pred <- function(train, newd) {
  m <- ranger(log(biomass) ~ ., data = train, num.trees = NTREE)
  as.numeric(predict(m, data = newd)$predictions)
}

families <- unique(d$family)
cat(sprintf("DBO RF: %d families, adaptive per-station temporal split (species-wide), %d trees\n",
  length(families), NTREE))
flush.console()
t0 <- Sys.time()

res <- map_dfr(seq_along(families), function(i) {
  fm <- families[i]

  # adaptive per-station split: last ceiling(n/3) years -> test, earlier -> train
  df <- d %>%
    filter(family == fm) %>%
    arrange(StationNme, DataYear) %>%
    group_by(StationNme) %>%
    mutate(ny = n(), te_n = pmax(min_test_yr, ceiling(ny / 3)),
      is_test = row_number() > (ny - te_n)) %>%
    ungroup()

  # keep only stations with enough train AND test years
  qual <- df %>% group_by(StationNme) %>%
    summarise(ny = first(ny), te_n = first(te_n), .groups = "drop") %>%
    filter(ny - te_n >= min_train_yr, te_n >= min_test_yr)
  df <- filter(df, StationNme %in% qual$StationNme)
  if (n_distinct(df$StationNme) < min_stn) return(NULL)

  tr <- filter(df, !is_test)
  te <- filter(df, is_test)

  # static: train on climatological (_spatial) means, project onto the test
  # rows' actual environment (raw values renamed to _spatial)
  tr_st <- select(tr, biomass, all_of(ctrl), all_of(paste0(base, "_spatial")))
  te_st <- te %>% select(all_of(ctrl), all_of(base)) %>%
    rename_with(~ paste0(., "_spatial"), all_of(base))

  tr_dy <- select(tr, biomass, all_of(ctrl), all_of(base))
  te_dy <- select(te, all_of(ctrl), all_of(base))

  tr_dc <- select(tr, biomass, all_of(ctrl), all_of(comp_terms))
  te_dc <- select(te, all_of(ctrl), all_of(comp_terms))

  p_static  <- fit_pred(tr_st, te_st)
  p_dynamic <- fit_pred(tr_dy, te_dy)
  p_decomp  <- fit_pred(tr_dc, te_dc)
  y <- log(te$biomass)

  if (i %% 10 == 0 || i == length(families)) {
    cat(sprintf("  %d/%d (%.0fs)\n", i, length(families),
      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    flush.console()
  }

  tibble(
    family = fm, n_stn = n_distinct(df$StationNme),
    train_rows = nrow(tr), test_rows = nrow(te),
    static_sw = swcor(te$StationNme, y, p_static),
    dynamic_sw = swcor(te$StationNme, y, p_dynamic),
    decomp_sw = swcor(te$StationNme, y, p_decomp)
  )
})
saveRDS(res, "data/rf_dbo_tempsplit_results.rds")

q <- function(x) sprintf("%.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
  quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))

# species-wide comparison across the three holdouts (same 45 families) ----------
loyo <- tryCatch(readRDS("data/rf_dbo_results.rds"), error = function(e) NULL)
loso <- tryCatch(readRDS("data/rf_dbo_loso_results.rds"), error = function(e) NULL)
med <- function(r, m) if (is.null(r)) "-" else sprintf("%.2f", median(r[[paste0(m, "_sw")]], na.rm = TRUE))

cat(sprintf("\n\n=== DBO RF species-wide, adaptive temporal split (%d families) ===\n", nrow(res)))
cat(sprintf("%-14s %-22s | %-14s %-14s\n", "Model", "Temporal split (SW)", "LOYO (SW)", "LOSO (SW)"))
for (m in c("static", "dynamic", "decomp")) {
  cat(sprintf("%-14s %-22s | %-14s %-14s\n", m, q(res[[paste0(m, "_sw")]]),
    med(loyo, m), med(loso, m)))
}
cat("\nsaved -> data/rf_dbo_tempsplit_results.rds\n")
