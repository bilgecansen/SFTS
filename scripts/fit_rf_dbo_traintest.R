# DBO population-level under a TRAIN/TEST split (mirrors BBS fit_rf_bbs.R). Replaces
# leave-one-year-out for the population-level panel: LOYO injects a regression-to-mean
# artifact into within-station plcor (permutation null sat at the observed -0.34; a
# pure LOO station mean = -1.00; spatial-only RF = -0.73 -- all reproduce the negative
# with NO real signal). A single train-early/test-late split removes the per-year
# exclusion: held-out years share one trained model, so within-station predictions
# move only with the environment.
#   train 2001-2010, test 2011-2019, honest per-set decomposition (decomp_dbo_traintest.rds)
#   four models: static / dynamic / decomposed / SVC(=decomposed + lat/long)
#   test stations must appear in training (transfer, not extrapolation to new sites)
# Output data/rf_dbo_traintest_results.rds (static/dynamic/decomp/svc _sw and _pl).
# Species-wide stays LOSO (rf_dbo_loso_*), which has no such artifact; only the DBO
# figure's population-level panel is rebuilt from this file.
library(tidyverse); library(ranger)
min_stn <- 10; min_test_stn <- 5; min_years <- 5; n_trees <- 2000
base <- c("Temp","Salinity","integchla","sedchla","Ammonia","Phosphate",
          "NiTriTra","Silicate","phigte5","TOC","cn")
comp_terms <- as.vector(t(outer(base, c("spatial","temporal","residual"), paste, sep = "_")))
ctrl <- "Depth"

dec <- readRDS("data/decomp_dbo_traintest.rds")
bio <- readRDS("data/data_dbo.rds") %>% distinct(family, StationNme, DataYear, biomass)
dat <- inner_join(bio, dec, by = c("StationNme", "DataYear"))
train <- filter(dat, set == "train"); test <- filter(dat, set == "test")

f_static  <- as.formula(paste("log(biomass) ~", paste(c(ctrl, paste0(base,"_spatial")), collapse = " + ")))
f_dynamic <- as.formula(paste("log(biomass) ~", paste(c(ctrl, base), collapse = " + ")))
f_decomp  <- as.formula(paste("log(biomass) ~", paste(c(ctrl, comp_terms), collapse = " + ")))
f_svc     <- as.formula(paste("log(biomass) ~", paste(c(ctrl, comp_terms, "Latitude","Longitude"), collapse = " + ")))

plcor <- function(stn,y,pred) tibble(stn,y,pred) %>% group_by(stn) %>% filter(n()>=min_years) %>%
  summarise(r=suppressWarnings(cor(y,pred)),.groups="drop") %>% pull(r) %>% median(na.rm=TRUE)
swcor <- function(stn,y,pred){ q <- tibble(stn,y,pred) %>% group_by(stn) %>% summarise(o=mean(y),p=mean(pred),.groups="drop"); suppressWarnings(cor(q$o,q$p)) }
rf <- function(f,d) ranger(f, data=d, num.trees=n_trees, seed=1, num.threads=0)
pr <- function(m,newd) predict(m, data=as.data.frame(newd))$predictions

families <- unique(train$family)
res <- map_dfr(families, function(fm) {
  tr <- filter(train, family==fm)
  te <- filter(test, family==fm) %>% filter(StationNme %in% unique(tr$StationNme)) %>%
    group_by(StationNme) %>% filter(n() >= min_years) %>% ungroup()
  if (n_distinct(tr$StationNme) < min_stn || n_distinct(te$StationNme) < min_test_stn) return(NULL)
  y <- log(te$biomass)
  te_static <- te; te_static[paste0(base,"_spatial")] <- te[base]   # raw yearly -> pop-level
  # static species-wide: per-station MEAN env -> one prediction per station
  tem <- te %>% group_by(StationNme) %>% summarise(o=mean(log(biomass)),
           across(all_of(c(base, ctrl)), mean), .groups="drop")
  te_avg <- tem; te_avg[paste0(base,"_spatial")] <- tem[base]
  tryCatch({
    m_s <- rf(f_static, tr)
    p_s <- pr(m_s, te_static); p_s_sw <- pr(m_s, te_avg)
    p_d <- pr(rf(f_dynamic, tr), te)
    p_c <- pr(rf(f_decomp, tr), te); p_v <- pr(rf(f_svc, tr), te)
    tibble(family=fm, n_stn=n_distinct(te$StationNme),
      static_sw=suppressWarnings(cor(tem$o, p_s_sw)), static_pl=plcor(te$StationNme,y,p_s),
      dynamic_sw=swcor(te$StationNme,y,p_d), dynamic_pl=plcor(te$StationNme,y,p_d),
      decomp_sw=swcor(te$StationNme,y,p_c), decomp_pl=plcor(te$StationNme,y,p_c),
      svc_sw=swcor(te$StationNme,y,p_v),    svc_pl=plcor(te$StationNme,y,p_v))
  }, error=function(e) NULL)
})
saveRDS(res, "data/rf_dbo_traintest_results.rds")
q <- function(x) sprintf("%+.2f (%+.2f,%+.2f)", median(x,na.rm=T), quantile(x,.05,na.rm=T), quantile(x,.95,na.rm=T))
cat(sprintf("saved -> data/rf_dbo_traintest_results.rds  (%d families)\n", nrow(res)))
cat(sprintf("%-9s %-22s %-22s\n","model","population-level","species-wide(ref)"))
for (m in c("static","dynamic","decomp","svc"))
  cat(sprintf("%-9s %-22s %-22s\n", m, q(res[[paste0(m,"_pl")]]), q(res[[paste0(m,"_sw")]])))
