# DBO RF on FOURTH-ROOT biomass with ZEROS RETAINED, both holdouts in one place.
#
# Scale. biomass^0.25. Zeros map to exactly 0, so no offset constant has to be invented
# (log1p barely transforms these data -- the positive median is 0.03, far below 1 -- and
# log(x + c) promotes near-detection-limit values to the most influential points). Fourth
# root takes skew from 17.6 to 2.0 and p99/p50 from 511 to 4.8, and leaves the among-station
# variance share at 47%, close to the 49% of the canonical log analysis, so results stay
# comparable. It is also the conventional transform for benthic biomass, though that
# precedent (Clarke & Warwick / PRIMER) is for multivariate Bray-Curtis work rather than
# univariate regression -- the empirical justification above is the load-bearing one.
#
# Balance. Uses the rebalanced grid from wrangle_dbo_data.R: the seven late-starting SEC
# stations and the sparse years 2008/2010 are dropped, giving 16 stations x 15 years at
# 94% fill (train 88%, test 99%) and stations per region of 5/4/7 instead of 5/4/14.
#
# Holdouts:
#   LOSO       : hold out an entire station (all years); global full-period decomposition
#                (data/data_dbo_zeros.rds), legitimate because the environment is exogenous.
#   train/test : train 2001-2007, test 2011-2019; honest per-set decomposition
#                (data/decomp_dbo_traintest.rds), each block decomposed within its own years.
#
# Models (response always spatiotemporal; only the ENV representation differs):
#   static / dynamic / decomposed / SVC, Depth as control. Static trains on per-station
#   means and predicts with spatiotemporal env fed into those slots.
#
# Levels: abundance distribution (_ad), regional trajectory (_rt), station trajectory (_st).
# Reported BOTH as correlations and as component-wise SKILL (1 - SSE/SS_obs) from the
# sequential orthogonal-projection error decomposition (station, then region-year | station,
# then residual), which partitions exactly under unbalance. Skill < 0 means worse than
# predicting no variation on that axis.
#
# Output: data/rf_dbo_4rt_{loso,traintest}_{results,preds}.rds

library(tidyverse)
library(ranger)

set.seed(1)
NTREE <- 2000
min_stn <- 10; min_test_stn <- 5; min_years <- 5
min_reg_stn <- 3; min_reg_years <- 5
MIN_RESID_DF <- 15

RESP <- function(b) b^0.25          # fourth root; 0 -> 0

base <- c("Temp", "Salinity", "integchla", "sedchla", "Ammonia", "Phosphate",
          "NiTriTra", "Silicate", "phigte5", "TOC", "cn")
comp_terms <- as.vector(t(outer(base, c("spatial", "temporal", "residual"), paste, sep = "_")))
ctrl <- "Depth"
sp_terms <- paste0(base, "_spatial")

rf <- function(f, d) ranger(f, data = d, num.trees = NTREE, seed = 1, num.threads = 0)
pr <- function(m, nd) as.numeric(predict(m, data = as.data.frame(nd))$predictions)
fml <- function(rhs) as.formula(paste("resp ~", paste(rhs, collapse = " + ")))
F_STATIC <- fml(c(ctrl, sp_terms)); F_DYN <- fml(c(ctrl, base))
F_DEC <- fml(c(ctrl, comp_terms));  F_SVC <- fml(c(ctrl, comp_terms, "Latitude", "Longitude"))

# --- the four models, given a train and test frame -----------------------------------
fit_all <- function(tr, te) {
  te_st <- te; te_st[sp_terms] <- te[base]              # raw yearly env -> _spatial slots
  m_s <- rf(F_STATIC, tr)
  list(static  = pr(m_s, te_st),
       dynamic = pr(rf(F_DYN, tr), te),
       decomp  = pr(rf(F_DEC, tr), te),
       svc     = pr(rf(F_SVC, tr), te),
       m_static = m_s)
}

# --- level correlations ---------------------------------------------------------------
adcor <- function(stn, y, p) { q <- tibble(stn, y, p) %>% group_by(stn) %>%
    summarise(o = mean(y), pp = mean(p), .groups = "drop"); suppressWarnings(cor(q$o, q$pp)) }
stcor <- function(stn, y, p) tibble(stn, y, p) %>% group_by(stn) %>%
  filter(sum(y > 0) >= min_years) %>%
  summarise(r = suppressWarnings(cor(y, p)), .groups = "drop") %>% pull(r) %>% median(na.rm = TRUE)
rtcor <- function(reg, yr, y, p) tibble(reg, yr, y, p) %>%
  group_by(reg, yr) %>% filter(n() >= 2) %>%
  summarise(o = mean(y), pp = mean(p), .groups = "drop") %>%
  group_by(reg) %>% filter(n() >= min_reg_years) %>%
  summarise(r = suppressWarnings(cor(o, pp)), .groups = "drop") %>% pull(r) %>% median(na.rm = TRUE)

# --- sequential error decomposition ----------------------------------------------------
seq_ss <- function(v, a, b) anova(lm(v ~ factor(a) + factor(b)))[["Sum Sq"]]
skill_tab <- function(preds, pc) {
  r <- preds %>% mutate(ry = paste(DBOreg, year)) %>% group_split(family) %>%
    map_dfr(function(d) {
      npar <- n_distinct(d$stn) + n_distinct(d$ry)
      if (nrow(d) - npar < MIN_RESID_DF) return(NULL)
      e <- d$y - d[[pc]]
      se <- seq_ss(e, d$stn, d$ry); so <- seq_ss(d$y, d$stn, d$ry)
      tibble(e_ad = se[1], e_rt = se[2], e_st = se[3], e_tot = sum(e^2), bias = nrow(d)*mean(e)^2,
             o_ad = so[1], o_rt = so[2], o_st = so[3], o_tot = sum((d$y - mean(d$y))^2))
    })
  if (!nrow(r)) return(NULL)
  tibble(ad = 1 - sum(r$e_ad)/sum(r$o_ad), rt = 1 - sum(r$e_rt)/sum(r$o_rt),
         st = 1 - sum(r$e_st)/sum(r$o_st), pooled = 1 - sum(r$e_tot)/sum(r$o_tot),
         bias_pct = 100*sum(r$bias)/sum(r$e_tot), nfam = nrow(r),
         addit = max(abs(r$e_tot - (r$bias + r$e_ad + r$e_rt + r$e_st))))
}

# ======================= LOSO =========================================================
d <- readRDS("data/data_dbo_zeros.rds") %>% mutate(resp = RESP(biomass))
cat("LOSO: fitting...\n"); flush.console(); t0 <- Sys.time()
loso <- map(unique(d$family), function(fm) {
  df <- filter(d, family == fm)
  if (n_distinct(df$StationNme[df$biomass > 0]) < min_stn) return(NULL)
  per <- map(sort(unique(df$StationNme)), function(s) {
    tr <- filter(df, StationNme != s); te <- filter(df, StationNme == s)
    if (nrow(tr) < 10 || nrow(te) < 1) return(NULL)
    p <- fit_all(tr, te)
    te_avg <- te %>% summarise(across(all_of(c(base, ctrl)), mean))
    te_avg[sp_terms] <- te_avg[base]
    list(yr = tibble(family = fm, DBOreg = te$DBOreg, stn = te$StationNme, year = te$DataYear,
                     y = te$resp, p_static = p$static, p_dynamic = p$dynamic,
                     p_decomp = p$decomp, p_svc = p$svc),
         ad = tibble(stn = s, o = mean(te$resp), p = pr(p$m_static, te_avg)))
  })
  per <- per[!map_lgl(per, is.null)]
  if (length(per) < 2) return(NULL)
  list(yr = bind_rows(map(per, "yr")), ad = bind_rows(map(per, "ad")))
})
loso <- loso[!map_lgl(loso, is.null)]
pl <- map_dfr(loso, "yr"); ad_static_l <- map(loso, "ad")
saveRDS(pl, "data/rf_dbo_4rt_loso_preds.rds")
cat(sprintf("  done (%.0fs), %d families\n", as.numeric(difftime(Sys.time(), t0, units="secs")), length(loso)))

# ======================= TRAIN / TEST ==================================================
dec <- readRDS("data/decomp_dbo_traintest.rds")
bio <- readRDS("data/data_dbo_zeros.rds") %>% distinct(family, StationNme, DataYear, DBOreg, biomass)
dat <- inner_join(bio, dec, by = c("StationNme", "DataYear")) %>% mutate(resp = RESP(biomass))
train <- filter(dat, set == "train"); test <- filter(dat, set == "test")
cat("train/test: fitting...\n"); flush.console()
tt <- map(unique(train$family), function(fm) {
  tr <- filter(train, family == fm)
  te <- filter(test, family == fm) %>%
    filter(StationNme %in% unique(tr$StationNme[tr$biomass > 0])) %>%
    group_by(StationNme) %>% filter(sum(biomass > 0) >= min_years) %>% ungroup()
  if (n_distinct(tr$StationNme[tr$biomass > 0]) < min_stn ||
      n_distinct(te$StationNme) < min_test_stn) return(NULL)
  p <- tryCatch(fit_all(tr, te), error = function(e) NULL); if (is.null(p)) return(NULL)
  tem <- te %>% group_by(StationNme) %>%
    summarise(o = mean(resp), across(all_of(c(base, ctrl)), mean), .groups = "drop")
  te_avg <- tem; te_avg[sp_terms] <- tem[base]
  list(yr = tibble(family = fm, DBOreg = te$DBOreg, stn = te$StationNme, year = te$DataYear,
                   y = te$resp, p_static = p$static, p_dynamic = p$dynamic,
                   p_decomp = p$decomp, p_svc = p$svc),
       ad = tibble(stn = tem$StationNme, o = tem$o, p = pr(p$m_static, te_avg)))
})
tt <- tt[!map_lgl(tt, is.null)]
pt <- map_dfr(tt, "yr"); ad_static_t <- map(tt, "ad")
saveRDS(pt, "data/rf_dbo_4rt_traintest_preds.rds")

# ======================= REPORT ========================================================
mods <- c(static="p_static", dynamic="p_dynamic", decomp="p_decomp", svc="p_svc")
report <- function(preds, ad_static, lab) {
  cat(sprintf("\n===== %s | %d families, %d station-years =====\n", lab,
              n_distinct(preds$family), nrow(preds)))
  cors <- map_dfr(names(mods), function(m) {
    pc <- mods[[m]]
    per <- preds %>% group_split(family) %>% imap_dfr(function(d, i) tibble(
      ad = if (m == "static") suppressWarnings(cor(ad_static[[i]]$o, ad_static[[i]]$p))
           else adcor(d$stn, d$y, d[[pc]]),
      rt = rtcor(d$DBOreg, d$year, d$y, d[[pc]]),
      st = stcor(d$stn, d$y, d[[pc]])))
    tibble(model = m, ad = median(per$ad, na.rm=TRUE), rt = median(per$rt, na.rm=TRUE),
           st = median(per$st, na.rm=TRUE))
  })
  sks <- map_dfr(names(mods), function(m) bind_cols(tibble(model=m), skill_tab(preds, mods[[m]])))
  cat(sprintf("%-9s | %-27s | %-35s\n", "", "CORRELATION (median)", "SKILL (pooled)"))
  cat(sprintf("%-9s | %8s %8s %8s | %8s %8s %8s %8s\n", "model",
      "abund", "region", "station", "abund", "region", "station", "POOLED"))
  for (i in seq_len(nrow(cors)))
    cat(sprintf("%-9s | %+8.3f %+8.3f %+8.3f | %+8.3f %+8.3f %+8.3f %+8.3f\n",
        cors$model[i], cors$ad[i], cors$rt[i], cors$st[i],
        sks$ad[i], sks$rt[i], sks$st[i], sks$pooled[i]))
  cat(sprintf("  bias %.1f%% of SSE, additivity %.1e\n", sks$bias_pct[1], max(sks$addit)))
  cat("  (variance partition and the 4-component split: scripts/skill_dbo_4rt.R)\n")
  list(cors = cors, skill = sks)
}
r_l <- report(pl, ad_static_l, "LOSO, fourth root + zeros")
r_t <- report(pt, ad_static_t, "TRAIN/TEST 2001-2007 -> 2011-2019, fourth root + zeros")
saveRDS(r_l, "data/rf_dbo_4rt_loso_results.rds")
saveRDS(r_t, "data/rf_dbo_4rt_traintest_results.rds")
cat("\nsaved -> data/rf_dbo_4rt_{loso,traintest}_{results,preds}.rds\n")
