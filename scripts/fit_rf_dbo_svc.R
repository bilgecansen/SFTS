# DBO RF SVC-CANDIDATE (exploratory): the RF analog of the GAM SVC = decomposed
# environment PLUS Latitude/Longitude, so RF can split on location (spatial intercept)
# and interact location with the components (spatially varying slopes). Runs BOTH
# holdouts the DBO figure uses, fitting ONLY the SVC model (static/dynamic/decomposed
# already exist in data/rf_dbo_loso_results.rds and data/rf_dbo_results.rds and are NOT
# re-run):
#   LOSO (leave-one-station-out) -> species-wide panel   -> data/rf_dbo_svc_loso_results.rds
#   LOYO (leave-one-year-out)    -> population-level panel -> data/rf_dbo_svc_results.rds
# Everything else copied from fit_rf_dbo_loso.R / fit_rf_dbo_3models.R.
library(tidyverse); library(ranger)

set.seed(1)
NTREE <- 2000; min_years <- 5; min_stn <- 10
d <- readRDS("data/data_dbo.rds") %>% filter(biomass > 0)
base <- c("Temp","Salinity","integchla","sedchla","Ammonia","Phosphate",
          "NiTriTra","Silicate","phigte5","TOC","cn")
comp_terms <- as.vector(t(outer(base, c("spatial","temporal","residual"), paste, sep = "_")))
ctrl <- "Depth"

plcor <- function(stn, y, pred) tibble(stn, y, pred) %>% group_by(stn) %>% filter(n() >= min_years) %>%
  summarise(r = suppressWarnings(cor(y, pred)), .groups = "drop") %>% pull(r) %>% median(na.rm = TRUE)
swcor <- function(stn, y, pred){ s <- tibble(stn, y, pred) %>% group_by(stn) %>%
  summarise(o = mean(y), p = mean(pred), .groups = "drop"); suppressWarnings(cor(s$o, s$p)) }
fit_pred <- function(train, newd){ m <- ranger(log(biomass) ~ ., data = train, num.trees = NTREE)
  as.numeric(predict(m, data = newd)$predictions) }

families <- unique(d$family)

run_svc <- function(holdout, out_path) {
  t0 <- Sys.time()
  res <- map_dfr(seq_along(families), function(i) {
    fm <- families[i]; df <- filter(d, family == fm)
    if (n_distinct(df$StationNme) < min_stn) return(NULL)
    keys <- if (holdout == "station") sort(unique(df$StationNme)) else sort(unique(df$DataYear))
    oos <- map_dfr(keys, function(k) {
      if (holdout == "station") { tr <- filter(df, StationNme != k); te <- filter(df, StationNme == k) }
      else                      { tr <- filter(df, DataYear   != k); te <- filter(df, DataYear   == k) }
      if (nrow(tr) < 10 || nrow(te) < 1) return(NULL)
      # SVC candidate = decomposed components + coordinates
      tr_sv <- select(tr, biomass, all_of(ctrl), all_of(comp_terms), Latitude, Longitude)
      te_sv <- select(te, all_of(ctrl), all_of(comp_terms), Latitude, Longitude)
      tibble(stn = te$StationNme, year = te$DataYear, y = log(te$biomass), p_svc = fit_pred(tr_sv, te_sv))
    })
    if (nrow(oos) < min_years) return(NULL)
    tibble(family = fm, n_stn = n_distinct(df$StationNme), rows = nrow(df),
      svc_sw = swcor(oos$stn, oos$y, oos$p_svc), svc_pl = plcor(oos$stn, oos$y, oos$p_svc))
  })
  saveRDS(res, out_path)
  q <- function(x) sprintf("%.2f (%.2f-%.2f)", median(x, na.rm = TRUE),
    quantile(x, .05, na.rm = TRUE), quantile(x, .95, na.rm = TRUE))
  cat(sprintf("saved -> %s (%d families, %.0fs)  svc_sw %s | svc_pl %s\n",
    out_path, nrow(res), as.numeric(difftime(Sys.time(), t0, units = "secs")),
    q(res$svc_sw), q(res$svc_pl)))
}

run_svc("station", "data/rf_dbo_svc_loso_results.rds")  # species-wide (honest spatial extrapolation)
run_svc("year",    "data/rf_dbo_svc_results.rds")        # population-level (leave-one-year-out)
