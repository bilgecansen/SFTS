# Ladder predictions for the main figures, with breakup in the decomposed model only.
#
# Identical to fit_rf_dbo_ladder.R except for one thing: the decomposed + Lat/Long model
# additionally carries breakup_global, breakup_regional and breakup_residual. Static and
# dynamic are unchanged and remain ice-free.
#
#   static    Depth + the _spatial columns, predicted with the year's values in those slots
#   dynamic   Depth + the raw spatiotemporal values
#   svc       Depth, Latitude, Longitude, all four components of every in-situ variable,
#             PLUS the three temporal components of breakup
#
# breakup only, temporal components only. Between stations the ice variables correlate
# 0.77-0.97 with each other and 0.99 with latitude, so a spatial ice term re-labels position
# that Lat/Long already carry; and adding all four ice variables dilutes rather than helps
# (train/test global at the total rung: breakup alone +0.446, all four +0.183). breakup is
# also the variable Frey and Grebmeier identify as the mechanism -- early retreat, longer
# open-water production, export to the bed.
#
# SELECTION CAVEAT. breakup was chosen after scoring four candidates on the test block, and
# a variable correlated 0.95 with persistence once flipped that variable's global-trajectory
# skill from +0.71 to -0.88. The definition has not yet been perturbed.
#
# LOSO and train/test only; LORO is not used by the main figures.
#
# Output: data/rf_ladder_bu_preds.rds, data/rf_ladder_bu_skill.rds

suppressMessages({ library(tidyverse); library(ranger) })

set.seed(1)
NTREE <- 2000
MIN_STN <- 10; MIN_TEST_STN <- 5; MIN_YEARS <- 5; MIN_RESID_DF <- 15
RESP <- function(b) b^0.25

BU <- "breakup"
TCOMP <- c("global", "regional", "residual")
comps <- c("spatial", TCOMP)
MODELS <- c("static", "dynamic", "svc")
CLABS <- c(o1 = "abundance", o2 = "global", o3 = "regional", o4 = "station")
RUNGS <- c(family = "data/dbo_resp_family.rds", class = "data/dbo_resp_class.rds",
           total  = "data/dbo_resp_total.rds")

ice <- readRDS("data/data_sic_vars.rds") %>% select(StationNme, DataYear, all_of(BU)) %>%
  group_by(StationNme) %>%
  mutate(breakup = ifelse(is.na(breakup), median(breakup, na.rm = TRUE), breakup)) %>%
  ungroup()

decompose <- function(x, v, reg) {
  stn <- factor(x$StationNme); yr <- factor(x$DataYear); ry <- factor(paste(reg, x$DataYear))
  cx <- x[[v]] - mean(x[[v]], na.rm = TRUE)
  f1 <- fitted(lm(cx ~ stn)); f2 <- fitted(lm(cx ~ stn + yr))
  f3 <- fitted(lm(cx ~ stn + yr + ry))
  x[[paste0(v, "_global")]]   <- as.numeric(f2 - f1)
  x[[paste0(v, "_regional")]] <- as.numeric(f3 - f2)
  x[[paste0(v, "_residual")]] <- as.numeric(cx - f3)
  x
}

env_cv <- readRDS("data/data_dbo_zeros.rds") %>%
  distinct(StationNme, DataYear, .keep_all = TRUE) %>% select(-family, -biomass) %>%
  inner_join(ice, by = c("StationNme", "DataYear"))
env_cv <- decompose(env_cv, BU, env_cv$DBOreg)

env_tt <- readRDS("data/decomp_dbo_traintest.rds") %>%
  inner_join(ice, by = c("StationNme", "DataYear"))
env_tt <- bind_rows(
  decompose(filter(env_tt, set == "train"), BU, filter(env_tt, set == "train")$DBOreg),
  decompose(filter(env_tt, set == "test"),  BU, filter(env_tt, set == "test")$DBOreg))

# static and dynamic use the in-situ variables only; svc adds the breakup temporal terms
model_terms <- function(e) {
  v <- setdiff(sub("_spatial$", "", grep("_spatial$", names(e), value = TRUE)), BU)
  list(vars = v, sp = paste0(v, "_spatial"),
       static  = c("Depth", paste0(v, "_spatial")),
       dynamic = c("Depth", v),
       svc     = c("Depth", "Latitude", "Longitude",
                   as.vector(t(outer(v, comps, paste, sep = "_"))),
                   paste0(BU, "_", TCOMP)))
}
TM_CV <- model_terms(env_cv); TM_TT <- model_terms(env_tt)
cat(sprintf("LOSO  static %d / dynamic %d / svc %d\n", length(TM_CV$static),
            length(TM_CV$dynamic), length(TM_CV$svc)))
cat(sprintf("Train/test  static %d / dynamic %d / svc %d\n", length(TM_TT$static),
            length(TM_TT$dynamic), length(TM_TT$svc)))

fml <- function(r) as.formula(paste("resp ~", paste(r, collapse = " + ")))
rf  <- function(f, d) ranger(f, data = d, num.trees = NTREE, seed = 1, num.threads = 0)
pr  <- function(m, nd) as.numeric(predict(m, data = as.data.frame(nd))$predictions)

fit_models <- function(tr, te, tm) {
  te_st <- te; te_st[tm$sp] <- te[tm$vars]
  list(static  = pr(rf(fml(tm$static),  tr), te_st),
       dynamic = pr(rf(fml(tm$dynamic), tr), te),
       svc     = pr(rf(fml(tm$svc),     tr), te))
}

fit_loso <- function(rsp) {
  d <- inner_join(rsp, select(env_cv, -DBOreg), by = c("StationNme", "DataYear")) %>%
    mutate(resp = RESP(biomass))
  map_dfr(sort(unique(d$taxon)), function(tx) {
    df <- filter(d, taxon == tx)
    if (n_distinct(df$StationNme[df$biomass > 0]) < MIN_STN) return(NULL)
    map_dfr(sort(unique(df$StationNme)), function(g) {
      tr <- filter(df, StationNme != g); te <- filter(df, StationNme == g)
      if (nrow(tr) < 10 || nrow(te) < 1) return(NULL)
      p <- fit_models(tr, te, TM_CV)
      obs <- tibble(taxon = tx, DBOreg = te$DBOreg, stn = te$StationNme,
                    year = te$DataYear, y = te$resp)
      map_dfr(MODELS, function(m) mutate(obs, model = m, p = p[[m]]))
    })
  })
}

fit_tt <- function(rsp) {
  d <- inner_join(rsp, select(env_tt, -DBOreg), by = c("StationNme", "DataYear")) %>%
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

preds <- map_dfr(names(RUNGS), function(rg) {
  cat(sprintf("  %s\n", rg))
  rsp <- readRDS(RUNGS[[rg]])
  bind_rows(mutate(fit_loso(rsp), holdout = "LOSO"),
            mutate(fit_tt(rsp),   holdout = "Train/test")) %>% mutate(rung = rg)
})
saveRDS(preds, "data/rf_ladder_bu_preds.rds")

# canonical four-component skill from the same predictions, so the figure scripts have a
# reference to assert their repartitioning against
seq4 <- function(v, stn, yr, ry) {
  a <- anova(lm(v ~ factor(stn) + factor(yr) + factor(ry)))[["Sum Sq"]]
  if (length(a) < 4) a <- c(a[seq_len(length(a)-1)], rep(0, 4-length(a)), a[length(a)])[1:4]
  a
}
skill <- preds %>% group_split(rung, holdout, model) %>% map_dfr(function(x) {
  r <- x %>% mutate(ry = paste(DBOreg, year)) %>% group_split(taxon) %>%
    map_dfr(function(t) {
      if (nrow(t) - (n_distinct(t$stn) + n_distinct(t$ry)) < MIN_RESID_DF) return(NULL)
      e <- t$y - t$p
      se <- seq4(e, t$stn, t$year, t$ry); so <- seq4(t$y, t$stn, t$year, t$ry)
      tibble(e1=se[1],e2=se[2],e3=se[3],e4=se[4],o1=so[1],o2=so[2],o3=so[3],o4=so[4]) })
  map_dfr(names(CLABS), function(k) {
    j <- sub("^o","e",k)
    tibble(rung = x$rung[1], holdout = x$holdout[1], model = x$model[1],
           component = CLABS[[k]], skill = 1 - sum(r[[j]]) / sum(r[[k]]))
  })
})
saveRDS(skill, "data/rf_ladder_bu_skill.rds")

cat(sprintf("\n%d prediction rows -> data/rf_ladder_bu_preds.rds\n", nrow(preds)))
print(as.data.frame(skill %>% mutate(skill = round(skill, 3)) %>%
  pivot_wider(names_from = component, values_from = skill) %>%
  arrange(holdout, rung, model)), row.names = FALSE)
