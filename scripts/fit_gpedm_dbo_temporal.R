# DBO GP-EDM, leave-one-year-out, TEMPORAL prediction levels (regional + population-level).
# Matched to scripts/fit_rf_dbo_temporal.R: same 45 families (>= 10 stations), same LOYO
# folds, same two metrics computed by identical code. The one thing that is NOT matched --
# deliberately -- is the information set: the RF models predict from ENVIRONMENT, GP-EDM
# predicts from the population's OWN PAST. That contrast is the point of the comparison.
#
# Model: hierarchical GP-EDM (GPEDM::fitGP), E = 1, tau = 1, one-step-ahead.
#   response  : log biomass at year t          (RAW -- same response as the RF models)
#   predictor : log biomass at year t-1        (the station's own previous year)
#   pop       : StationNme (stations pooled within a family, hierarchical dynamics)
#   scaling   : "global" -- ONE standardization over the whole family dataset, NOT per
#               station. This is the choice that keeps the comparison honest: any
#               per-station centring of the response would either (a) recompute a station
#               mean per fold, injecting the -1/(T-1) leave-one-out anti-correlation that
#               the RF pipeline was rebuilt to remove, or (b) use a full-period station
#               mean, letting the held-out year contribute to its own centring -- response
#               leakage the environment never has. Raw log biomass avoids both.
#
# LOYO fold handling (lenient): training drops rows whose TARGET year is the held-out
# year. Rows predicting year y+1 keep their observed N(y) as the lag input, so the
# held-out observation still appears in training as a predictor. The strict alternative
# (also dropping rows whose lag comes from year y) costs roughly twice the rows per fold.
#
# Coverage note: only consecutive-year pairs are usable, and a family absent in a year
# breaks the chain (zero-biomass rows are dropped upstream), so ~5.5k of the RF's 7.7k
# scored rows survive. All 45 families still qualify.
#
# Input : data/data_dbo.rds
# Output: data/gpedm_dbo_temporal_results.rds (edm_rg / edm_pl per family)
#         data/gpedm_dbo_temporal_preds.rds   (station-level LOYO predictions)

suppressMessages({
  library(GPEDM)
  library(tidyverse)
})

set.seed(1)
min_years <- 5 # station (or region) needs >= this many years for a correlation
min_stn <- 10  # family needs >= this many stations, matching the RF scripts

d <- readRDS("data/data_dbo.rds") %>%
  filter(biomass > 0) %>%
  mutate(logN = log(biomass))

families <- d %>% group_by(family) %>%
  filter(n_distinct(StationNme) >= min_stn) %>% ungroup() %>%
  pull(family) %>% unique() %>% sort()

cat(sprintf("DBO GP-EDM: %d families, LOYO, E=1, scaling=global\n", length(families)))
flush.console()
t0 <- Sys.time()

edm_family <- function(fm) {
  sub <- d %>% filter(family == fm) %>% arrange(StationNme, DataYear) %>%
    select(family, DBOreg, StationNme, DataYear, logN)
  lg <- tryCatch(makelags(as.data.frame(sub), y = "logN", pop = "StationNme",
    time = "DataYear", E = 1, tau = 1), error = function(e) NULL)
  if (is.null(lg)) return(NULL)
  df <- cbind(as.data.frame(sub), lg)
  df <- df[is.finite(df$logN_1), ]
  if (nrow(df) < 20 || n_distinct(df$StationNme) < 3) return(NULL)

  map_dfr(sort(unique(df$DataYear)), function(y) {
    tr <- df[df$DataYear != y, ]                       # lenient: drop targets only
    te <- df[df$DataYear == y, ]
    te <- te[te$StationNme %in% unique(tr$StationNme), ]
    if (nrow(tr) < 15 || nrow(te) < 1 || n_distinct(tr$StationNme) < 3) return(NULL)
    te <- te %>% arrange(StationNme) %>% group_by(StationNme) %>%
      mutate(timestep = row_number()) %>% ungroup()

    fit <- tryCatch(suppressWarnings(fitGP(data = tr, y = "logN", x = "logN_1",
      pop = "StationNme", scaling = "global")), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    p2 <- tryCatch(suppressWarnings(predict(fit, newdata = as.data.frame(te))),
      error = function(e) NULL)
    if (is.null(p2)) return(NULL)

    p2$outsampresults %>%
      left_join(distinct(te, StationNme, timestep, DataYear, DBOreg, family),
        by = c("pop" = "StationNme", "timestep")) %>%
      transmute(family, DBOreg, stn = pop, year = DataYear, y = obs, p_edm = predmean)
  })
}

preds <- map_dfr(seq_along(families), function(i) {
  if (i %% 10 == 0) {
    cat(sprintf("  %d/%d (%.0fs)\n", i, length(families),
      as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
  }
  tryCatch(edm_family(families[i]), error = function(e) NULL)
}) %>% filter(is.finite(y), is.finite(p_edm))
saveRDS(preds, "data/gpedm_dbo_temporal_preds.rds")

# --- metrics: identical code to scripts/fit_rf_dbo_temporal.R -------------------------
plcor <- function(dd, col) dd %>%
  group_by(stn) %>% filter(n() >= min_years) %>%
  summarise(r = suppressWarnings(cor(y, .data[[col]])), .groups = "drop") %>%
  pull(r) %>% median(na.rm = TRUE)

rgcor <- function(dd, col) dd %>%
  group_by(DBOreg, year) %>% filter(n() >= 2) %>%          # >= 2 stations per regional mean
  summarise(o = mean(y), p = mean(.data[[col]]), .groups = "drop") %>%
  group_by(DBOreg) %>% filter(n() >= min_years) %>%        # >= 5 region-years to correlate
  summarise(r = suppressWarnings(cor(o, p)), .groups = "drop") %>%
  pull(r) %>% median(na.rm = TRUE)

res <- preds %>% group_split(family) %>% map_dfr(function(dd) tibble(
  family = dd$family[1], n_stn = n_distinct(dd$stn),
  edm_rg = rgcor(dd, "p_edm"), edm_pl = plcor(dd, "p_edm")))
saveRDS(res, "data/gpedm_dbo_temporal_results.rds")

wsd <- preds %>% group_by(family, stn) %>% filter(n() >= min_years) %>%
  summarise(x = sd(p_edm) / sd(y), .groups = "drop") %>% pull(x) %>% median(na.rm = TRUE)

q <- function(x) sprintf("%+.2f (%.2f - %.2f)", median(x, na.rm = TRUE),
  quantile(x, 0.05, na.rm = TRUE), quantile(x, 0.95, na.rm = TRUE))
cat(sprintf("\n\n=== DBO GP-EDM, LOYO, temporal levels (%d families, %d rows scored) ===\n",
  nrow(res), nrow(preds)))
cat(sprintf("%-14s %-22s %-22s %s\n", "Model", "Regional", "Population-level", "SD(pred)/SD(obs)"))
cat(sprintf("%-14s %-22s %-22s %.2f\n", "GP-EDM", q(res$edm_rg), q(res$edm_pl), wsd))
cat("\nsaved -> data/gpedm_dbo_temporal_results.rds, data/gpedm_dbo_temporal_preds.rds\n")
