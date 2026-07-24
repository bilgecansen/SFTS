# CANONICAL (see wrangle_decomp_blocks.R) -- GAM strata, new structure (per-block decomp, new windows), all 3
# treatments. Mirrors scripts/fit_gam_bbs.R + fit_gam_bbs_spatialcv.R, EXCEPT
# static now uses the 8 vars' spatial components (like RF), not bio1+bio12.
# 4 models: static(8 spatial) / dynamic(8) / decomposed(8x3) / svc(smooth).
library(tidyverse); library(mgcv); library(foreach)
SP <- "/private/tmp/claude-504/-Users-bsen3-Library-Mobile-Documents-com-apple-CloudDocs-Documents-SFTS/5062e77b-19cf-49df-8568-6314674f5c71/scratchpad"
min_strata <- 15; min_test_years <- 5; min_strata_cv <- 25; min_pool_strata <- 15; K <- 8

vars <- c("bio2","bio3","bio5","bio8","bio9","bio15","bio16","bio18")
static_comp <- paste0(vars, "_spatial")            # 8-var static
comp_terms  <- as.vector(t(outer(vars, c("spatial","temporal","residual"), paste, sep="_")))
resid_terms <- paste0(vars, "_residual")
keep <- c("abundance","elevs","party_hours","lat","long", vars, comp_terms)

base <- readRDS("data/data_bbs_nozero.rds") %>%
  select(species_id, site_id, year, strata, lat, long, abundance, elevs, party_hours)

sm <- function(x) paste(sprintf("s(%s)", x), collapse=" + ")
ctrl <- "elevs + party_hours"
rhs_static  <- paste(sm(static_comp), "+", ctrl)
rhs_dynamic <- paste(sm(vars), "+", ctrl)
rhs_decomp  <- paste(sm(comp_terms), "+", ctrl)
rhs_svc_fix <- paste(sm(c(paste0(vars,"_spatial"), paste0(vars,"_temporal"))), "+", ctrl)

plcor <- function(s,y,pred) tibble(s,y,pred) %>% group_by(s) %>% summarise(r=suppressWarnings(cor(y,pred)),.groups="drop") %>% pull(r) %>% median(na.rm=TRUE)
swcor <- function(s,y,pred){ q <- tibble(s,y,pred) %>% group_by(s) %>% summarise(o=mean(y),p=mean(pred),.groups="drop"); suppressWarnings(cor(q$o,q$p)) }
gfit <- function(f,d) bam(f, data=d, discrete=TRUE, method="fREML")
gp   <- function(m,newd) as.numeric(predict(m, newd))
svc_forms <- function(nst) {
  k_sp <- min(30, nst-1); k_by <- min(10, nst-1)
  list(int=sprintf("s(long, lat, k=%d)", k_sp),
       by =paste(sprintf("s(long, lat, by=%s, k=%d)", resid_terms, k_by), collapse=" + "))
}
mods <- function(tr, te) {
  s <- svc_forms(n_distinct(tr$strata))
  f_static  <- as.formula(paste("log(abundance) ~", rhs_static))
  f_dynamic <- as.formula(paste("log(abundance) ~", rhs_dynamic))
  f_decomp  <- as.formula(paste("log(abundance) ~", rhs_decomp))
  f_svc     <- as.formula(paste("log(abundance) ~", rhs_svc_fix, "+", s$int, "+", s$by))
  te_static <- te; te_static[static_comp] <- te[vars]
  list(p_static=gp(gfit(f_static,tr),te_static), p_dynamic=gp(gfit(f_dynamic,tr),te),
       p_decomp=gp(gfit(f_decomp,tr),te),        p_svc=gp(gfit(f_svc,tr),te))
}

run <- function(decomp_path, train_yrs, blocking, out_path, tag) {
  clim <- readRDS(decomp_path) %>% select(-set)
  dat  <- inner_join(base, clim, by=c("site_id","year"))
  train <- filter(dat, year %in% train_yrs); test <- filter(dat, year >= 2011 & year <= 2020)
  year_n <- train %>% group_by(species_id, site_id) %>% summarise(n=sum(abundance>0), .groups="drop") %>% filter(n>4)
  train2 <- semi_join(train, year_n, by=c("species_id","site_id"))
  site_n <- train2 %>% group_by(species_id) %>% summarize(n=length(unique(site_id)), .groups="drop") %>% filter(n>4)
  train3 <- filter(train2, species_id %in% site_n$species_id)
  species <- unique(train3$species_id)
  agg <- function(df, sp, sk=NULL){ d <- filter(df, species_id==sp); if(!is.null(sk)) d <- filter(d, strata %in% sk)
    d %>% group_by(strata,year) %>% summarise(across(all_of(keep),mean),.groups="drop") }
  cat(sprintf("\n[GAM strata %s] %d species, blocking=%s\n", tag, length(species), blocking)); flush.console()

  if (!blocking) {
    res <- foreach(i=seq_along(species)) %do% {
      sp <- species[i]
      tr <- agg(train3, sp) %>% mutate(strata=factor(strata))
      te <- agg(test, sp) %>% filter(strata %in% levels(tr$strata)) %>%
        group_by(strata) %>% filter(n()>=min_test_years) %>% ungroup() %>% mutate(strata=factor(strata, levels=levels(tr$strata)))
      if (i %% 25 == 0) { cat(sprintf("  %d/%d\n", i, length(species))); flush.console() }
      if (n_distinct(tr$strata) < min_strata || nrow(te) < min_test_years) return(NULL)
      y <- log(te$abundance)
      tryCatch({ p <- mods(tr, te)
        tibble(species_id=sp, n_strata=n_distinct(tr$strata),
          static_sw=swcor(te$strata,y,p$p_static),  static_pl=plcor(te$strata,y,p$p_static),
          dynamic_sw=swcor(te$strata,y,p$p_dynamic),dynamic_pl=plcor(te$strata,y,p$p_dynamic),
          decomp_sw=swcor(te$strata,y,p$p_decomp),  decomp_pl=plcor(te$strata,y,p$p_decomp),
          svc_sw=swcor(te$strata,y,p$p_svc),        svc_pl=plcor(te$strata,y,p$p_svc)) }, error=function(e) NULL)
    }
  } else {
    coord <- train3 %>% group_by(species_id, strata) %>% summarise(lat=mean(lat), long=mean(long), .groups="drop")
    res <- foreach(i=seq_along(species)) %do% {
      sp <- species[i]; co <- filter(coord, species_id==sp)
      if (i %% 25 == 0) { cat(sprintf("  %d/%d\n", i, length(species))); flush.console() }
      if (nrow(co) < min_strata_cv) return(NULL)
      xy <- cbind(x=co$long*cos(mean(co$lat)*pi/180), y=co$lat)
      km <- tryCatch({ set.seed(1); kmeans(xy, centers=K, nstart=10) }, error=function(e) NULL); if (is.null(km)) return(NULL)
      co$block <- km$cluster
      pooled <- foreach(k=1:K, .combine=bind_rows) %do% {
        tr <- agg(train3, sp, co$strata[co$block!=k]) %>% mutate(strata=factor(strata))
        te <- agg(test, sp, co$strata[co$block==k]) %>% group_by(strata) %>% filter(n()>=min_test_years) %>% ungroup()
        if (n_distinct(tr$strata) < 5 || nrow(te) < min_test_years) return(NULL)
        tryCatch({ p <- mods(tr, te)
          tibble(strata=te$strata, y=log(te$abundance), p_static=p$p_static, p_dynamic=p$p_dynamic, p_decomp=p$p_decomp, p_svc=p$p_svc) },
          error=function(e) NULL)
      }
      if (is.null(pooled) || n_distinct(pooled$strata) < min_pool_strata) return(NULL)
      tibble(species_id=sp, n_strata=nrow(co), n_pool=n_distinct(pooled$strata),
        static_sw=swcor(pooled$strata,pooled$y,pooled$p_static),  static_pl=plcor(pooled$strata,pooled$y,pooled$p_static),
        dynamic_sw=swcor(pooled$strata,pooled$y,pooled$p_dynamic),dynamic_pl=plcor(pooled$strata,pooled$y,pooled$p_dynamic),
        decomp_sw=swcor(pooled$strata,pooled$y,pooled$p_decomp),  decomp_pl=plcor(pooled$strata,pooled$y,pooled$p_decomp),
        svc_sw=swcor(pooled$strata,pooled$y,pooled$p_svc),        svc_pl=plcor(pooled$strata,pooled$y,pooled$p_svc))
    }
  }
  res <- bind_rows(res); saveRDS(res, out_path)
  q <- function(x) sprintf("%.2f", median(x, na.rm=T))
  cat(sprintf("=== GAM strata [%s] (%d sp) sw sta/dyn/dec/svc: %s/%s/%s/%s\n", tag, nrow(res),
    q(res$static_sw), q(res$dynamic_sw), q(res$decomp_sw), q(res$svc_sw)))
  cat(sprintf("saved -> %s\n", out_path)); flush.console()
}

run("data/decomp_temporal.rds", 2001:2010, FALSE, "data/gam_bbs_noclamp_results.rds", "Temporal")
run("data/decomp_buffer.rds",   1991:2000, FALSE, "data/gam_bbs_noclamp_buffer_results.rds",   "Buffer")
run("data/decomp_temporal.rds", 2001:2010, TRUE,  "data/gam_bbs_noclamp_spatialcv_k8_results.rds", "Spatiotemporal")
