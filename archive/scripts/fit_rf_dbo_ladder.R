# Biological-aggregation ladder, random-forest rungs.
#
# Question: does prediction skill rise as the response is aggregated over taxa, and if
# so for which components? Three rungs, identical machinery throughout:
#
#   family   57 families     data/dbo_resp_family.rds
#   class     8 classes      data/dbo_resp_class.rds
#   total     1 series       data/dbo_resp_total.rds   (TotGC)
#
# ONE model per rung: the decomposed environment plus coordinates (the SVC model).
# Predictors are the four components -- _spatial, _global, _regional, _residual -- of
# each of 11 environmental variables and their one-year lags, with Depth, Latitude and
# Longitude. Static / dynamic / decomposed are not fitted here.
#
# Three holdouts:
#   LOSO        hold out one station at a time, 16 folds; environment decomposed over
#               all years (data/data_dbo_zeros.rds).
#   LORO        hold out one DBO region at a time, 3 folds -- every station in the
#               region leaves training together. Regions are uneven (5 / 4 / 7
#               stations), which is fine: the taxon-level gate is applied to the whole
#               data set, not per fold, so holding out DBO 3 does not disqualify taxa.
#               Latitude and Longitude are predictors, and a random forest cannot
#               extrapolate beyond the coordinate range it was trained on, so the
#               held-out region is predicted from the nearest training region.
#   train/test  8 training years (2001, 2003-2007, 2011-2012) against 7 test years
#               (2013-2019), each block decomposed within its own years
#               (data/decomp_dbo_traintest.rds).
#
# Under LOSO and LORO the folds are concatenated per taxon before the error is
# decomposed, so all 16 stations and all 3 regions are present and the region:year term
# is estimable. Decomposing a single LORO fold on its own would collapse region:year
# onto year.
#
# Scale: fourth root, applied AFTER aggregation -- the response tables are raw grams
# carbon and were summed there.
#
# Gates are identical at every rung so the comparison is like for like: a taxon needs
# >=10 stations with positive biomass in training, >=5 test stations under train/test,
# and a test station needs >=5 positive years. This costs classes that are patchily
# distributed; the script reports which taxa enter at each rung.
#
# ERROR DECOMPOSITION. Per taxon, the error e = y - yhat and the observed y are both
# split by sequential orthogonal projection onto station -> year -> region:year ->
# residual, giving abundance distribution / global trajectory / regional trajectory /
# station trajectory. Station is nested in region, so station is entered first and the
# time-invariant latitudinal gradient lands in the abundance component; region:year
# then carries only the interaction, which makes the split order-free.
#
# Skill is POOLED across taxa: skill_k = 1 - sum(SSE_k) / sum(SS_obs_k), summing over
# taxa within a rung. At the total rung there is one taxon and this is just its skill.
# SS_obs and SSE are reported alongside so the change in n up the ladder is visible.
#
# Output: data/rf_dbo_ladder_skill.rds, data/rf_dbo_ladder_preds.rds

suppressMessages({
  library(tidyverse)
  library(ranger)
})

set.seed(1)
NTREE <- 2000
MIN_STN <- 10; MIN_TEST_STN <- 5; MIN_YEARS <- 5; MIN_RESID_DF <- 15

RESP <- function(b) b^0.25

base <- c("Temp", "Salinity", "integchla", "sedchla", "Ammonia", "Phosphate",
          "NiTriTra", "Silicate", "phigte5", "TOC", "cn")
comps <- c("spatial", "global", "regional", "residual")

RUNGS <- c(family = "data/dbo_resp_family.rds",
           class  = "data/dbo_resp_class.rds",
           total  = "data/dbo_resp_total.rds")
CLABS <- c(o1 = "abundance", o2 = "global", o3 = "regional", o4 = "station")

# ---- environment ----------------------------------------------------------------------
# LOSO/LORO use the globally decomposed table, which carries the lagged variables too;
# train/test uses the per-block table, which carries the 11 base variables only. The
# variable list is therefore read off whatever each table actually holds.
env_cv <- readRDS("data/data_dbo_zeros.rds") %>%
  distinct(StationNme, DataYear, .keep_all = TRUE) %>%
  select(-family, -biomass, -DBOreg)

env_tt <- readRDS("data/decomp_dbo_traintest.rds") %>% select(-DBOreg)

# ---- the three model types -------------------------------------------------------------
# The RESPONSE is always spatiotemporal; only the environmental representation differs.
#   static   trained on the per-station MEANS (the _spatial columns), then predicted with
#            the spatiotemporal values fed into those same slots -- a relationship fitted
#            across space and applied through time.
#   dynamic  raw spatiotemporal values in training and in prediction.
#   svc      all four components in both, plus Latitude and Longitude.
# Depth is the control in all three.
model_terms <- function(e) {
  sp <- grep("_spatial$", names(e), value = TRUE)
  v  <- sub("_spatial$", "", sp)
  stopifnot(all(v %in% names(e)))
  list(vars = v, sp = sp,
       static  = c("Depth", sp),
       dynamic = c("Depth", v),
       svc     = c("Depth", "Latitude", "Longitude",
                   as.vector(t(outer(v, comps, paste, sep = "_")))))
}
TM_CV <- model_terms(env_cv); TM_TT <- model_terms(env_tt)
MODELS <- c("static", "dynamic", "svc")

fml <- function(rhs) as.formula(paste("resp ~", paste(rhs, collapse = " + ")))
rf <- function(f, d) ranger(f, data = d, num.trees = NTREE, seed = 1, num.threads = 0)
pr <- function(m, nd) as.numeric(predict(m, data = as.data.frame(nd))$predictions)

fit_models <- function(tr, te, tm) {
  te_st <- te; te_st[tm$sp] <- te[tm$vars]      # raw yearly env -> the _spatial slots
  list(static  = pr(rf(fml(tm$static),  tr), te_st),
       dynamic = pr(rf(fml(tm$dynamic), tr), te),
       svc     = pr(rf(fml(tm$svc),     tr), te))
}

cat(sprintf("predictors  static %d / %d, dynamic %d / %d, svc %d / %d  (LOSO-LORO / train-test)\n\n",
            length(TM_CV$static), length(TM_TT$static),
            length(TM_CV$dynamic), length(TM_TT$dynamic),
            length(TM_CV$svc), length(TM_TT$svc)))

# ---- the three holdouts ---------------------------------------------------------------
# LOSO and LORO differ only in the column that defines a fold.
fit_cv <- function(rsp, fold_by) {
  d <- inner_join(rsp, env_cv, by = c("StationNme", "DataYear")) %>%
    mutate(resp = RESP(biomass))
  map_dfr(sort(unique(d$taxon)), function(tx) {
    df <- filter(d, taxon == tx)
    if (n_distinct(df$StationNme[df$biomass > 0]) < MIN_STN) return(NULL)
    map_dfr(sort(unique(df[[fold_by]])), function(g) {
      tr <- filter(df, .data[[fold_by]] != g); te <- filter(df, .data[[fold_by]] == g)
      if (nrow(tr) < 10 || nrow(te) < 1) return(NULL)
      p <- fit_models(tr, te, TM_CV)
      obs <- tibble(taxon = tx, DBOreg = te$DBOreg, stn = te$StationNme,
                    year = te$DataYear, y = te$resp)
      map_dfr(MODELS, function(m) mutate(obs, model = m, p = p[[m]]))
    })
  })
}

fit_tt <- function(rsp) {
  d <- inner_join(rsp, env_tt, by = c("StationNme", "DataYear")) %>%
    mutate(resp = RESP(biomass))
  tr_all <- filter(d, set == "train"); te_all <- filter(d, set == "test")
  map_dfr(sort(unique(tr_all$taxon)), function(tx) {
    tr <- filter(tr_all, taxon == tx)
    te <- filter(te_all, taxon == tx) %>%
      filter(StationNme %in% unique(tr$StationNme[tr$biomass > 0])) %>%
      group_by(StationNme) %>% filter(sum(biomass > 0) >= MIN_YEARS) %>% ungroup()
    if (n_distinct(tr$StationNme[tr$biomass > 0]) < MIN_STN ||
        n_distinct(te$StationNme) < MIN_TEST_STN) return(NULL)
    p <- fit_models(tr, te, TM_TT)
    obs <- tibble(taxon = tx, DBOreg = te$DBOreg, stn = te$StationNme,
                  year = te$DataYear, y = te$resp)
    map_dfr(MODELS, function(m) mutate(obs, model = m, p = p[[m]]))
  })
}

# ---- four-component error decomposition -----------------------------------------------
seq4 <- function(v, stn, yr, ry) {
  a <- anova(lm(v ~ factor(stn) + factor(yr) + factor(ry)))[["Sum Sq"]]
  if (length(a) < 4)
    a <- c(a[seq_len(length(a) - 1)], rep(0, 4 - length(a)), a[length(a)])[1:4]
  a
}

per_taxon <- function(p) {
  p %>% mutate(ry = paste(DBOreg, year)) %>% group_split(taxon) %>%
    map_dfr(function(d) {
      if (nrow(d) - (n_distinct(d$stn) + n_distinct(d$ry)) < MIN_RESID_DF) return(NULL)
      e <- d$y - d$p
      se <- seq4(e, d$stn, d$year, d$ry); so <- seq4(d$y, d$stn, d$year, d$ry)
      tibble(taxon = d$taxon[1], n = nrow(d),
             e1 = se[1], e2 = se[2], e3 = se[3], e4 = se[4], e_tot = sum(e^2),
             o1 = so[1], o2 = so[2], o3 = so[3], o4 = so[4],
             o_tot = sum((d$y - mean(d$y))^2), bias = nrow(d) * mean(e)^2)
    })
}

skill_tab <- function(p, rung, holdout) {
  map_dfr(MODELS, function(m) {
    r <- per_taxon(filter(p, model == m))
    if (is.null(r) || !nrow(r)) return(NULL)
    bind_rows(
      map_dfr(names(CLABS), function(k) {
        j <- sub("^o", "e", k)
        tibble(component = CLABS[[k]], SS_obs = sum(r[[k]]), SSE = sum(r[[j]]))
      }),
      tibble(component = "pooled", SS_obs = sum(r$o_tot), SSE = sum(r$e_tot))) %>%
      mutate(rung = rung, holdout = holdout, model = m, skill = 1 - SSE / SS_obs,
             n_taxa = nrow(r), n_rows = sum(r$n),
             bias_pct = 100 * sum(r$bias) / sum(r$e_tot),
             addit = max(abs(r$e_tot - (r$bias + r$e1 + r$e2 + r$e3 + r$e4)))) %>%
      select(rung, holdout, model, component, SS_obs, SSE, skill, n_taxa, n_rows,
             bias_pct, addit)
  })
}

# ---- run ------------------------------------------------------------------------------
HOLDOUTS <- c("LOSO", "LORO", "Train/test")

all_preds <- list(); all_skill <- list()
for (rung in names(RUNGS)) {
  rsp <- readRDS(RUNGS[[rung]])
  for (hold in HOLDOUTS) {
    t0 <- Sys.time()
    p <- switch(hold,
                LOSO = fit_cv(rsp, "StationNme"),
                LORO = fit_cv(rsp, "DBOreg"),
                fit_tt(rsp))
    if (!nrow(p)) { cat(sprintf("%-7s %-11s no taxa passed\n", rung, hold)); next }
    cat(sprintf("%-7s %-11s %2d taxa, %5d station-years x %d models (%.0fs): %s\n",
                rung, hold, n_distinct(p$taxon), nrow(p) / length(MODELS), length(MODELS),
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                if (rung == "family") "" else paste(sort(unique(p$taxon)), collapse = ", ")))
    all_preds[[paste(rung, hold)]] <- mutate(p, rung = rung, holdout = hold)
    all_skill[[paste(rung, hold)]] <- skill_tab(p, rung, hold)
  }
}

preds <- bind_rows(all_preds); skill <- bind_rows(all_skill)
saveRDS(preds, "data/rf_dbo_ladder_preds.rds")
saveRDS(skill, "data/rf_dbo_ladder_skill.rds")

# ---- report ---------------------------------------------------------------------------
for (hold in HOLDOUTS) {
  cat(sprintf("\n===== %s =====\n", hold))
  cat(sprintf("%-8s %-7s %5s %6s | %10s %10s %10s %10s %10s\n", "model", "rung",
              "taxa", "rows", "abundance", "global", "regional", "station", "POOLED"))
  for (m in MODELS) {
    for (rung in names(RUNGS)) {
      s <- filter(skill, rung == !!rung, holdout == hold, model == m)
      if (!nrow(s)) next
      v <- setNames(s$skill, s$component)
      cat(sprintf("%-8s %-7s %5d %6d | %+10.3f %+10.3f %+10.3f %+10.3f %+10.3f\n",
                  m, rung, s$n_taxa[1], s$n_rows[1], v[["abundance"]], v[["global"]],
                  v[["regional"]], v[["station"]], v[["pooled"]]))
    }
    cat("\n")
  }
}
cat("\nSS_obs / SSE, all holdouts:\n")
print(as.data.frame(skill %>% transmute(holdout, model, rung, component,
                                        SS_obs = round(SS_obs, 1), SSE = round(SSE, 1),
                                        skill = round(skill, 3))), row.names = FALSE)
cat(sprintf("\nmax additivity error %.1e | bias %.1f-%.1f%% of SSE\n",
            max(skill$addit), min(skill$bias_pct), max(skill$bias_pct)))
cat("saved -> data/rf_dbo_ladder_{skill,preds}.rds\n")
