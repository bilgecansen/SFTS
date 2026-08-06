# Four-component error decomposition and component-wise SKILL for the DBO fourth-root
# analysis. Replaces the two-term (station + region-year) structure, which was
# ORDER-DEPENDENT for a structural reason that balancing cannot fix.
#
# The problem. Station is NESTED within region: every station belongs to exactly one DBO
# band, so any between-region contrast is also a between-station contrast. The region main
# effect (3 regions -> 2 dims) therefore lies inside BOTH the station subspace and the
# region-year subspace; their intersection is non-empty and sequential SS depends on term
# order. This persists under a perfectly balanced grid (verified by rank: overlap = 2 either
# way), which is why the SEC/sparse-year fix did not remove it. Those 2 shared dimensions are
# the latitudinal gradient -- the largest spatial signal present -- so the swing was large:
# abundance-distribution skill moved +0.404 (station first) to -0.152 (region-year first).
#
# The fix. Split region-year into a common year effect plus a region-specific deviation:
#     station  ->  year  ->  region:year  ->  residual
# Station and year are FULLY CROSSED with no nesting, so their subspaces meet only at the
# intercept (rank overlap 0, as in the balanced simulation). Entering station first assigns
# the latitudinal gradient to the spatial component -- which is where it belongs, since it is
# static in time -- and region:year then contributes only the region x year interaction, i.e.
# purely time-varying structure that cannot absorb the static gradient. Any residual order
# sensitivity is then due to unbalance alone, and is reported.
#
# Components -> prediction levels:
#   station      -> abundance distribution   (spatial, incl. the between-region gradient)
#   year         -> synchronous trajectory   (the temporal signal common to all three bands)
#   region:year  -> regional divergence      (how the bands' trajectories differ)
#   residual     -> station trajectory       (station-specific deviation)
#
# Skill_k = 1 - SSE_k / SS_observed_k;  < 0 means worse than predicting no variation there.
#
# Input : data/rf_dbo_4rt_{loso,traintest}_preds.rds
# Output: data/rf_dbo_4rt_skill4.rds

library(tidyverse)

MIN_RESID_DF <- 15
mods <- c(static = "p_static", dynamic = "p_dynamic", decomp = "p_decomp", svc = "p_svc")

# sequential SS for v ~ station + year + region:year ; returns 4 components
seq4 <- function(v, stn, yr, ry) {
  a <- anova(lm(v ~ factor(stn) + factor(yr) + factor(ry)))[["Sum Sq"]]
  if (length(a) < 4) a <- c(a[1:(length(a) - 1)], rep(0, 4 - length(a)), a[length(a)])[1:4]
  a
}
seq4_rev <- function(v, stn, yr, ry) {          # year before station, for sensitivity
  a <- anova(lm(v ~ factor(yr) + factor(stn) + factor(ry)))[["Sum Sq"]]
  c(a[2], a[1], a[3], a[4])                      # reorder to station, year, ry, resid
}

decompose <- function(preds, pc, rev = FALSE) {
  f <- if (rev) seq4_rev else seq4
  preds %>% mutate(ry = paste(DBOreg, year)) %>% group_split(family) %>%
    map_dfr(function(d) {
      npar <- n_distinct(d$stn) + n_distinct(d$ry)
      if (nrow(d) - npar < MIN_RESID_DF) return(NULL)
      e <- d$y - d[[pc]]
      se <- f(e, d$stn, d$year, d$ry); so <- f(d$y, d$stn, d$year, d$ry)
      tibble(e_ad = se[1], e_yr = se[2], e_rd = se[3], e_st = se[4],
             o_ad = so[1], o_yr = so[2], o_rd = so[3], o_st = so[4],
             bias = nrow(d) * mean(e)^2, e_tot = sum(e^2),
             o_tot = sum((d$y - mean(d$y))^2))
    })
}
sk <- function(r, a, b) 1 - sum(r[[a]]) / sum(r[[b]])

run <- function(path, lab) {
  p <- readRDS(path)
  cat(sprintf("\n===== %s | %d families, %d station-years =====\n", lab,
              n_distinct(p$family), nrow(p)))
  base_r <- decompose(p, mods[[1]])
  tot <- sum(base_r$o_ad + base_r$o_yr + base_r$o_rd + base_r$o_st)
  cat(sprintf("observed variance: abundance %.0f%% | synchronous %.0f%% | regional div %.0f%% | station %.0f%%\n",
      100*sum(base_r$o_ad)/tot, 100*sum(base_r$o_yr)/tot,
      100*sum(base_r$o_rd)/tot, 100*sum(base_r$o_st)/tot))
  cat(sprintf("additivity: %.1e | families used: %d\n\n",
      max(abs(base_r$e_tot - (base_r$bias + base_r$e_ad + base_r$e_yr + base_r$e_rd + base_r$e_st))),
      nrow(base_r)))
  cat(sprintf("%-9s %10s %12s %13s %10s %10s %8s\n", "model",
      "abundance", "synchronous", "regional div", "station", "POOLED", "bias%"))
  out <- map_dfr(names(mods), function(m) {
    r <- decompose(p, mods[[m]]); rr <- decompose(p, mods[[m]], rev = TRUE)
    o <- tibble(dataset = lab, model = m,
      ad = sk(r, "e_ad", "o_ad"), yr = sk(r, "e_yr", "o_yr"),
      rd = sk(r, "e_rd", "o_rd"), st = sk(r, "e_st", "o_st"),
      pooled = sk(r, "e_tot", "o_tot"), bias_pct = 100*sum(r$bias)/sum(r$e_tot),
      ad_rev = sk(rr, "e_ad", "o_ad"), yr_rev = sk(rr, "e_yr", "o_yr"))
    cat(sprintf("%-9s %+10.3f %+12.3f %+13.3f %+10.3f %+10.3f %7.1f%%\n",
        m, o$ad, o$yr, o$rd, o$st, o$pooled, o$bias_pct))
    o
  })
  cat("\n  order check (year entered before station):\n")
  pwalk(out, function(dataset, model, ad, yr, rd, st, pooled, bias_pct, ad_rev, yr_rev)
    cat(sprintf("    %-9s abundance %+.3f -> %+.3f   synchronous %+.3f -> %+.3f\n",
        model, ad, ad_rev, yr, yr_rev)))
  out
}

res <- bind_rows(
  run("data/rf_dbo_4rt_loso_preds.rds",      "LOSO, fourth root + zeros"),
  run("data/rf_dbo_4rt_traintest_preds.rds", "TRAIN/TEST, fourth root + zeros"))
saveRDS(res, "data/rf_dbo_4rt_skill4.rds")
cat("\nsaved -> data/rf_dbo_4rt_skill4.rds\n")
