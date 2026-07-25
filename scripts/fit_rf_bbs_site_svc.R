# RF site/route SVC-CANDIDATE (exploratory). Site-level counterpart to
# fit_rf_bbs_svc.R: the RF analog of the GAM SVC = decomposed climate components PLUS
# lat/long, so RF can split on location (spatial intercept) and interact location with
# the components (spatially varying slopes). Fits ONLY the SVC model; the existing
# static/dynamic/decomposed live in data/rf_bbs_site_*.rds and are NOT re-run.
# Everything else is copied verbatim from fit_rf_bbs_site.R. Outputs
# data/rf_bbs_site_svc_results.rds (Temporal), _buffer_results.rds,
# _spatialcv_k8_results.rds, each with svc_sw / svc_pl per species.
library(tidyverse); library(ranger); library(foreach)
min_sites <- 15; min_test_years <- 5; min_sites_cv <- 25; min_pool_sites <- 15; K <- 8; n_trees <- 1000

vars <- c("bio2","bio3","bio5","bio8","bio9","bio15","bio16","bio18")
comp_terms <- as.vector(t(outer(vars, c("spatial","temporal","residual"), paste, sep="_")))
keep <- c("abundance","elevs","party_hours","lat","long", vars, comp_terms)  # + lat/long vs fit_rf_bbs_site.R

base <- readRDS("data/data_bbs_nozero.rds") %>%
  select(species_id, site_id, year, lat, long, abundance, elevs, party_hours)

f_svc <- as.formula(paste("log(abundance) ~", paste(c("elevs","party_hours","lat","long", comp_terms), collapse=" + ")))
plcor <- function(s,y,pred) tibble(s,y,pred) %>% group_by(s) %>% summarise(r=suppressWarnings(cor(y,pred)),.groups="drop") %>% pull(r) %>% median(na.rm=TRUE)
swcor <- function(s,y,pred){ q <- tibble(s,y,pred) %>% group_by(s) %>% summarise(o=mean(y),p=mean(pred),.groups="drop"); suppressWarnings(cor(q$o,q$p)) }
rf <- function(f,d) ranger(f, data=d, num.trees=n_trees, seed=1, num.threads=0)
pr <- function(m,newd) predict(m, data=as.data.frame(newd))$predictions
preds <- function(tr, te) list(v=pr(rf(f_svc,tr),te))

run <- function(decomp_path, train_yrs, blocking, out_path, tag) {
  clim <- readRDS(decomp_path) %>% select(-set)
  dat  <- inner_join(base, clim, by=c("site_id","year"))
  train <- filter(dat, year %in% train_yrs); test <- filter(dat, year >= 2011 & year <= 2020)
  year_n <- train %>% group_by(species_id, site_id) %>% summarise(n=length(abundance[abundance>0]), .groups="drop") %>% filter(n>4)
  train2 <- semi_join(train, year_n, by=c("species_id","site_id"))
  site_n <- train2 %>% group_by(species_id) %>% summarize(n=length(unique(site_id)), .groups="drop") %>% filter(n>4)
  train3 <- filter(train2, species_id %in% site_n$species_id)
  species <- unique(train3$species_id)
  agg <- function(df, sp, sk=NULL){ d <- filter(df, species_id==sp); if(!is.null(sk)) d <- filter(d, site_id %in% sk)
    d %>% group_by(site_id,year) %>% summarise(across(all_of(keep),mean),.groups="drop") }
  cat(sprintf("\n[RF site SVC %s] %d species, blocking=%s\n", tag, length(species), blocking)); flush.console()
  if (!blocking) {
    res <- foreach(i=seq_along(species)) %do% {
      sp <- species[i]; tr <- agg(train3, sp)
      te <- agg(test, sp) %>% filter(site_id %in% unique(tr$site_id)) %>% group_by(site_id) %>% filter(n()>=min_test_years) %>% ungroup()
      if (n_distinct(tr$site_id) < min_sites || nrow(te) < min_test_years) return(NULL)
      y <- log(te$abundance)
      tryCatch({ p <- preds(tr, te)
        tibble(species_id=sp, n_sites=n_distinct(tr$site_id),
          svc_sw=swcor(te$site_id,y,p$v), svc_pl=plcor(te$site_id,y,p$v)) }, error=function(e) NULL)
    }
  } else {
    coord <- train3 %>% group_by(species_id, site_id) %>% summarise(lat=mean(lat), long=mean(long), .groups="drop")
    res <- foreach(i=seq_along(species)) %do% {
      sp <- species[i]; co <- filter(coord, species_id==sp)
      if (nrow(co) < min_sites_cv) return(NULL)
      xy <- cbind(x=co$long*cos(mean(co$lat)*pi/180), y=co$lat)
      km <- tryCatch({ set.seed(1); kmeans(xy, centers=K, nstart=10) }, error=function(e) NULL); if (is.null(km)) return(NULL)
      co$block <- km$cluster
      pooled <- foreach(k=1:K, .combine=bind_rows) %do% {
        tr <- agg(train3, sp, co$site_id[co$block!=k]); te <- agg(test, sp, co$site_id[co$block==k]) %>% group_by(site_id) %>% filter(n()>=min_test_years) %>% ungroup()
        if (n_distinct(tr$site_id) < 5 || nrow(te) < min_test_years) return(NULL)
        tryCatch({ p <- preds(tr, te); tibble(site_id=te$site_id, y=log(te$abundance), p_v=p$v) }, error=function(e) NULL)
      }
      if (is.null(pooled) || n_distinct(pooled$site_id) < min_pool_sites) return(NULL)
      tibble(species_id=sp, n_sites=nrow(co), n_pool=n_distinct(pooled$site_id),
        svc_sw=swcor(pooled$site_id,pooled$y,pooled$p_v), svc_pl=plcor(pooled$site_id,pooled$y,pooled$p_v))
    }
  }
  res <- bind_rows(res); saveRDS(res, out_path)
  cat(sprintf("saved -> %s (%d species)\n", out_path, nrow(res))); flush.console()
}
run("data/decomp_temporal.rds", 2001:2010, FALSE, "data/rf_bbs_site_svc_results.rds",             "Temporal")
run("data/decomp_buffer.rds",   1991:2000, FALSE, "data/rf_bbs_site_svc_buffer_results.rds",       "Buffered")
run("data/decomp_temporal.rds", 2001:2010, TRUE,  "data/rf_bbs_site_svc_spatialcv_k8_results.rds", "Spatiotemporal")
