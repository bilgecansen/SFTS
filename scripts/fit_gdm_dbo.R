# DBO community turnover with Generalized Dissimilarity Modelling (GDM).
# Complementary to the single-family SDM analysis: predicts COMPOSITIONAL TURNOVER between
# pairs of station-years rather than one taxon's biomass.
#
# TWO MODELS, both ordinary GDMs, each fitted and tested on its own class of site pairs:
#   SPATIAL TURNOVER  : pairs = same year, different station. Predictors are the 11 raw
#                       environmental variables + Depth, with geographic distance (geo=TRUE).
#                       The standard spatial-turnover GDM, pooled across years.
#   TEMPORAL TURNOVER : pairs = same station, different year. Predictors are the 11 raw
#                       environmental variables only. Depth and geographic distance are
#                       CONSTANT within a station, so their pairwise difference is
#                       identically zero and they carry no information here -- they are
#                       dropped rather than fitted as inert splines (geo = FALSE).
#
# Response: Bray-Curtis dissimilarity on FOURTH-ROOT biomass across 57 families -- the same
# scale as the single-family analysis, and the benthic convention.
#
# TWO HOLDOUTS, both applied at the UNIT level (every pair touching a held-out unit is
# removed; each station-year sits in ~225 pairs, so random pair holdout would leak badly):
#   LOSO       : hold out a station. Global full-period decomposition of the environment is
#                irrelevant here since only raw values are used.
#   train/test : fit on pairs with both members in 2001-2007, test on pairs with both
#                members in 2011-2019.
#
# LOSO DOUBLE-COUNTING. A pair has two members, so either being held out puts the pair in a
# test set: a pair between stations A and B is scored twice, once by the model trained
# without A and once by the model trained without B. Cross-station (spatial) pairs are
# therefore doubled; temporal pairs are not, since both members are the same station. The
# two predictions are also correlated, the models sharing 14 of 16 stations. Here each pair
# is scored ONCE, by averaging its predictions, and n reports true pair counts.
#
# Metric: skill = 1 - SSE / SS_obs, with SS_obs taken about the MEAN OBSERVED DISSIMILARITY
# of the test pairs. That null matters: sampling noise puts an unknown floor under temporal
# dissimilarity, so scoring against the mean asks whether the model explains VARIATION in
# turnover rather than its absolute level. Correlation and RMSE reported alongside.
#
# Output: data/gdm_dbo_results.rds (summary), data/gdm_dbo_preds.rds (pair-level predictions
#         with pair class and region-pair labels attached, for scripts/plot_gdm_dbo.R)

suppressMessages({
  library(gdm)
  library(tidyverse)
  library(vegan)
})

set.seed(1)
base <- c("Temp", "Salinity", "integchla", "sedchla", "Ammonia", "Phosphate",
          "NiTriTra", "Silicate", "phigte5", "TOC", "cn")
TRAIN_YRS <- 2001:2007; TEST_YRS <- 2011:2019

z <- readRDS("data/data_dbo_zeros.rds")

# ---- community matrix + Bray-Curtis on fourth-root biomass ---------------------------
comm <- z %>% mutate(v = biomass^0.25) %>%
  select(StationNme, DataYear, family, v) %>%
  pivot_wider(names_from = family, values_from = v, values_fill = 0) %>%
  arrange(StationNme, DataYear)
key <- comm %>% select(StationNme, DataYear)
D <- as.matrix(vegdist(as.matrix(comm[, -(1:2)]), method = "bray"))

sites <- z %>% distinct(StationNme, DataYear, .keep_all = TRUE) %>%
  select(StationNme, DataYear, DBOreg, Latitude, Longitude, Depth, all_of(base))
sord <- key %>% inner_join(sites, by = c("StationNme", "DataYear")) %>% as.data.frame()
stopifnot(identical(paste(sord$StationNme, sord$DataYear), paste(key$StationNme, key$DataYear)))

# ---- all pairs, with class and a stable pair id --------------------------------------
ij <- t(combn(nrow(sord), 2))
pairs_all <- tibble(ki = ij[, 1], kj = ij[, 2],
    distance = D[cbind(ij[, 1], ij[, 2])],
    s1 = sord$StationNme[ki], s2 = sord$StationNme[kj],
    y1 = sord$DataYear[ki],   y2 = sord$DataYear[kj]) %>%
  mutate(cls = case_when(y1 == y2 & s1 != s2 ~ "spatial",
                         s1 == s2 & y1 != y2 ~ "temporal",
                         TRUE ~ "spatiotemporal"),
         pair_id = paste(ki, kj),
         # region labels, carried through to the saved predictions purely for plotting --
         # DBOreg is never a GDM predictor (gdm_table selects only the requested vars)
         r1 = sord$DBOreg[ki], r2 = sord$DBOreg[kj],
         regpair = ifelse(r1 == r2, paste0("Within DBO ", r1),
                          paste0("DBO ", pmin(r1, r2), " - DBO ", pmax(r1, r2))),
         kind = ifelse(r1 == r2, "Within region", "Between regions")) %>%
  filter(is.finite(distance))

cat(sprintf("%d station-years, %d pairs (spatial %d, temporal %d, spatiotemporal %d)\n\n",
  nrow(sord), nrow(pairs_all), sum(pairs_all$cls == "spatial"),
  sum(pairs_all$cls == "temporal"), sum(pairs_all$cls == "spatiotemporal")))

# NOTE: rownames MUST be reset. Indexing sord by repeated row indices yields mangled
# rownames ("15", "15.1", ...), and gdm() SEGFAULTS on those in its C code rather than
# raising a catchable R error -- the values are otherwise identical.
gdm_table <- function(pr, vars) {
  v1 <- sord[pr$ki, vars, drop = FALSE]; names(v1) <- paste0("s1.", vars)
  v2 <- sord[pr$kj, vars, drop = FALSE]; names(v2) <- paste0("s2.", vars)
  rownames(v1) <- NULL; rownames(v2) <- NULL
  out <- data.frame(distance = pr$distance, weights = 1,
    s1.xCoord = sord$Longitude[pr$ki], s1.yCoord = sord$Latitude[pr$ki],
    s2.xCoord = sord$Longitude[pr$kj], s2.yCoord = sord$Latitude[pr$kj],
    v1, v2, check.names = FALSE)
  rownames(out) <- NULL
  class(out) <- c("gdmData", "data.frame")
  out
}

MODELS <- list(
  `Spatial turnover`  = list(cls = "spatial",  vars = c(base, "Depth"), geo = TRUE),
  `Temporal turnover` = list(cls = "temporal", vars = base,             geo = FALSE))

fit_predict <- function(spec, trp, tep) {
  if (nrow(trp) < 60 || nrow(tep) < 5) return(NULL)
  vars <- spec$vars[map_lgl(spec$vars, function(v)
    sd(c(sord[[v]][trp$ki], sord[[v]][trp$kj]), na.rm = TRUE) > 1e-10)]
  if (!length(vars)) return(NULL)
  m <- tryCatch(gdm(gdm_table(trp, vars), geo = spec$geo), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  p <- tryCatch(predict(m, gdm_table(tep, vars)), error = function(e) NULL)
  if (is.null(p)) return(NULL)
  tibble(pair_id = tep$pair_id, obs = tep$distance, pred = as.numeric(p),
         n_train = nrow(trp), dev = m$explained)
}

# each pair scored once: average duplicate predictions (see LOSO note above)
dedup <- function(d) d %>% group_by(pair_id) %>%
  summarise(obs = first(obs), pred = mean(pred), n_rep = n(), .groups = "drop")
score <- function(u, n_scored) tibble(n_pairs = nrow(u), n_scored = n_scored,
  skill = 1 - sum((u$obs - u$pred)^2) / sum((u$obs - mean(u$obs))^2),
  r = suppressWarnings(cor(u$obs, u$pred)),
  rmse = sqrt(mean((u$obs - u$pred)^2)))

# ======================= LOSO ==========================================================
loso <- map(names(MODELS), function(nm) {
  spec <- MODELS[[nm]]
  cls_pairs <- filter(pairs_all, cls == spec$cls)
  out <- map_dfr(sort(unique(sord$StationNme)), function(st) {
    trp <- filter(cls_pairs, s1 != st, s2 != st)
    tep <- filter(cls_pairs, s1 == st | s2 == st)
    r <- fit_predict(spec, trp, tep); if (is.null(r)) return(NULL)
    mutate(r, fold = st)
  })
  if (!nrow(out)) return(NULL)
  u <- dedup(out)
  list(summary = bind_cols(tibble(model = nm, holdout = "LOSO"), score(u, nrow(out)),
                           tibble(mean_dev = round(mean(out$dev), 1))),
       preds = mutate(u, model = nm, holdout = "LOSO"))
})

# ======================= TRAIN / TEST ==================================================
tt <- map(names(MODELS), function(nm) {
  spec <- MODELS[[nm]]
  cls_pairs <- filter(pairs_all, cls == spec$cls)
  trp <- filter(cls_pairs, y1 %in% TRAIN_YRS, y2 %in% TRAIN_YRS)
  tep <- filter(cls_pairs, y1 %in% TEST_YRS,  y2 %in% TEST_YRS)
  r <- fit_predict(spec, trp, tep); if (is.null(r)) return(NULL)
  u <- dedup(r)
  list(summary = bind_cols(tibble(model = nm, holdout = "Train/test"), score(u, nrow(r)),
                           tibble(mean_dev = round(r$dev[1], 1))),
       preds = mutate(u, model = nm, holdout = "Train/test"))
})

all <- c(loso, tt); all <- all[!map_lgl(all, is.null)]
res <- map_dfr(all, "summary")
saveRDS(res, "data/gdm_dbo_results.rds")
saveRDS(map_dfr(all, "preds") %>%
          left_join(distinct(pairs_all, pair_id, cls, regpair, kind), by = "pair_id"),
        "data/gdm_dbo_preds.rds")

cat("=== GDM community turnover: predictive performance ===\n")
cat(sprintf("%-18s %-12s %8s %9s %8s %8s %8s\n",
    "model", "holdout", "n pairs", "n scored", "skill", "r", "rmse"))
pwalk(res, function(model, holdout, n_pairs, n_scored, skill, r, rmse, mean_dev)
  cat(sprintf("%-18s %-12s %8d %9d %+8.3f %+8.3f %8.3f\n",
      model, holdout, n_pairs, n_scored, skill, r, rmse)))
cat("\n(n scored > n pairs under LOSO = cross-station pairs tested in two folds;\n",
    " predictions averaged so each pair contributes once)\n")
cat("\nsaved -> data/gdm_dbo_results.rds, data/gdm_dbo_preds.rds\n")
