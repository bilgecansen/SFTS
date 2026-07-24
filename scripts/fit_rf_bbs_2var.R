# CANONICAL (see wrangle_decomp_blocks.R) -- RF strata, per-block decomp, all 3 treatments, 2 VARIABLES only
# (bio1 = annual mean temp, bio12 = annual precip). Same structure as the 8-var
# block run. Output -> *_2var_*.rds.
library(tidyverse); library(ranger); library(foreach)
SP <- "/private/tmp/claude-504/-Users-bsen3-Library-Mobile-Documents-com-apple-CloudDocs-Documents-SFTS/5062e77b-19cf-49df-8568-6314674f5c71/scratchpad"
min_strata <- 15; min_test_years <- 5; min_strata_cv <- 25; min_pool_strata <- 15; K <- 8; n_trees <- 1000

vars <- c("bio1","bio12")
comp_terms <- as.vector(t(outer(vars, c("spatial","temporal","residual"), paste, sep="_")))
keep <- c("abundance","elevs","party_hours", vars, comp_terms)

base <- readRDS("data/data_bbs_nozero.rds") %>%
  select(species_id, site_id, year, strata, lat, long, abundance, elevs, party_hours)

f_static  <- as.formula(paste("log(abundance) ~", paste(c("elevs","party_hours", paste0(vars,"_spatial")), collapse=" + ")))
f_dynamic <- as.formula(paste("log(abundance) ~", paste(c("elevs","party_hours", vars), collapse=" + ")))
f_decomp  <- as.formula(paste("log(abundance) ~", paste(c("elevs","party_hours", comp_terms), collapse=" + ")))
plcor <- function(s,y,pred) tibble(s,y,pred) %>% group_by(s) %>% summarise(r=suppressWarnings(cor(y,pred)),.groups="drop") %>% pull(r) %>% median(na.rm=TRUE)
swcor <- function(s,y,pred){ q <- tibble(s,y,pred) %>% group_by(s) %>% summarise(o=mean(y),p=mean(pred),.groups="drop"); suppressWarnings(cor(q$o,q$p)) }
rf <- function(f,d) ranger(f, data=d, num.trees=n_trees, seed=1, num.threads=0)
pr <- function(m,newd) predict(m, data=as.data.frame(newd))$predictions
preds <- function(tr, te){ te_static <- te; te_static[paste0(vars,"_spatial")] <- te[vars]
  list(s=pr(rf(f_static,tr),te_static), d=pr(rf(f_dynamic,tr),te), c=pr(rf(f_decomp,tr),te)) }

run <- function(decomp_path, train_yrs, blocking, out_path, tag) {
  clim <- readRDS(decomp_path) %>% select(-set)
  dat  <- inner_join(base, clim, by=c("site_id","year"))
  train <- filter(dat, year %in% train_yrs); test <- filter(dat, year >= 2011 & year <= 2020)
  year_n <- train %>% group_by(species_id, site_id) %>% summarise(n=length(abundance[abundance>0]), .groups="drop") %>% filter(n>4)
  train2 <- semi_join(train, year_n, by=c("species_id","site_id"))
  site_n <- train2 %>% group_by(species_id) %>% summarize(n=length(unique(site_id)), .groups="drop") %>% filter(n>4)
  train3 <- filter(train2, species_id %in% site_n$species_id)
  species <- unique(train3$species_id)
  agg <- function(df, sp, sk=NULL){ d <- filter(df, species_id==sp); if(!is.null(sk)) d <- filter(d, strata %in% sk)
    d %>% group_by(strata,year) %>% summarise(across(all_of(keep),mean),.groups="drop") }
  cat(sprintf("\n[RF strata 2var %s] %d species, blocking=%s\n", tag, length(species), blocking)); flush.console()
  if (!blocking) {
    res <- foreach(i=seq_along(species)) %do% {
      sp <- species[i]; tr <- agg(train3, sp)
      te <- agg(test, sp) %>% filter(strata %in% unique(tr$strata)) %>% group_by(strata) %>% filter(n()>=min_test_years) %>% ungroup()
      if (n_distinct(tr$strata) < min_strata || nrow(te) < min_test_years) return(NULL)
      y <- log(te$abundance)
      tryCatch({ p <- preds(tr, te)
        tibble(species_id=sp, n_strata=n_distinct(tr$strata),
          static_sw=swcor(te$strata,y,p$s), static_pl=plcor(te$strata,y,p$s),
          dynamic_sw=swcor(te$strata,y,p$d),dynamic_pl=plcor(te$strata,y,p$d),
          decomp_sw=swcor(te$strata,y,p$c), decomp_pl=plcor(te$strata,y,p$c)) }, error=function(e) NULL)
    }
  } else {
    coord <- train3 %>% group_by(species_id, strata) %>% summarise(lat=mean(lat), long=mean(long), .groups="drop")
    res <- foreach(i=seq_along(species)) %do% {
      sp <- species[i]; co <- filter(coord, species_id==sp)
      if (nrow(co) < min_strata_cv) return(NULL)
      xy <- cbind(x=co$long*cos(mean(co$lat)*pi/180), y=co$lat)
      km <- tryCatch({ set.seed(1); kmeans(xy, centers=K, nstart=10) }, error=function(e) NULL); if (is.null(km)) return(NULL)
      co$block <- km$cluster
      pooled <- foreach(k=1:K, .combine=bind_rows) %do% {
        tr <- agg(train3, sp, co$strata[co$block!=k]); te <- agg(test, sp, co$strata[co$block==k]) %>% group_by(strata) %>% filter(n()>=min_test_years) %>% ungroup()
        if (n_distinct(tr$strata) < 5 || nrow(te) < min_test_years) return(NULL)
        tryCatch({ p <- preds(tr, te); tibble(strata=te$strata, y=log(te$abundance), p_s=p$s, p_d=p$d, p_c=p$c) }, error=function(e) NULL)
      }
      if (is.null(pooled) || n_distinct(pooled$strata) < min_pool_strata) return(NULL)
      tibble(species_id=sp, n_strata=nrow(co), n_pool=n_distinct(pooled$strata),
        static_sw=swcor(pooled$strata,pooled$y,pooled$p_s), static_pl=plcor(pooled$strata,pooled$y,pooled$p_s),
        dynamic_sw=swcor(pooled$strata,pooled$y,pooled$p_d),dynamic_pl=plcor(pooled$strata,pooled$y,pooled$p_d),
        decomp_sw=swcor(pooled$strata,pooled$y,pooled$p_c), decomp_pl=plcor(pooled$strata,pooled$y,pooled$p_c))
    }
  }
  res <- bind_rows(res); saveRDS(res, out_path)
  q <- function(x) sprintf("%.2f", median(x, na.rm=T))
  cat(sprintf("=== RF strata 2var [%s] (%d sp) sw %s/%s/%s\n", tag, nrow(res), q(res$static_sw), q(res$dynamic_sw), q(res$decomp_sw)))
  cat(sprintf("saved -> %s\n", out_path)); flush.console()
}
run("data/decomp_temporal.rds", 2001:2010, FALSE, "data/rf_bbs_2var_results.rds", "Temporal")
run("data/decomp_buffer.rds",   1991:2000, FALSE, "data/rf_bbs_2var_buffer_results.rds",   "Buffer")
run("data/decomp_temporal.rds", 2001:2010, TRUE,  "data/rf_bbs_2var_spatialcv_k8_results.rds", "Spatiotemporal")
