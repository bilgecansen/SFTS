# ANOVA-style ERROR decomposition and component-wise SKILL for the DBO raw-scale
# (zeros retained) train/test analysis. Companion to the correlations in
# fit_rf_dbo_traintest_zeros.R -- same predictions, a metric that is not blind to bias
# and amplitude.
#
# Why skill and not correlation. Correlation is invariant to bias and to any rescaling of
# the predictions, so a model that gets the pattern right but systematically compresses
# its predictions still scores well. Squared error is not, and it is the only common
# metric that partitions additively.
#
# Decomposition. Marginal-mean effects are NOT orthogonal in an unbalanced panel, so the
# sums of squares do not add (on the zeros-dropped data the cross-term reached -48% of the
# observed SS). Instead use SEQUENTIAL orthogonal projections onto nested subspaces:
#     intercept  <  + station  <  + region-year  <  everything
# Each component is the projection onto the part of the next subspace orthogonal to all
# previous ones. Those projections are mutually orthogonal and sum to the identity, so
#     ||e||^2 = bias + SSE_station + SSE_regionyear + SSE_residual
# holds EXACTLY whether or not the design is balanced. anova(lm(...)) returns precisely
# this (Type I / sequential SS). The design matrix depends only on which station-years
# exist, which is identical for observed and predicted, so the same projections apply to
# both and the error vector can be decomposed directly.
#
# Components map onto the three prediction levels:
#   station main effect            -> abundance distribution
#   region-year effect | station   -> regional trajectory
#   residual (station x year)      -> station trajectory
#
# Skill_k = 1 - SSE_k / SS_observed_k.  Skill < 0 means the model is WORSE than predicting
# no variation at all on that axis.
#
# ORDER. Sequential SS depends on term order. Station is entered first: it is the dominant
# axis, and the question is what remains for time after space is accounted for -- the
# conservative ordering for a temporal claim. The reversed order is reported as a
# sensitivity check.
#
# Input : data/rf_dbo_traintest_zeros_preds.rds
# Output: data/rf_dbo_traintest_zeros_skill.rds

library(tidyverse)

MIN_RESID_DF <- 15

p <- readRDS("data/rf_dbo_traintest_zeros_preds.rds") %>% mutate(ry = paste(DBOreg, year))
mods <- c(static = "p_static", dynamic = "p_dynamic", decomp = "p_decomp", svc = "p_svc")

seq_ss <- function(v, a, b) {          # sequential SS for v ~ a + b, in that order
  s <- anova(lm(v ~ factor(a) + factor(b)))[["Sum Sq"]]
  c(first = s[1], second = s[2], resid = s[3])
}

decomp_family <- function(d, pc, reverse = FALSE) {
  n <- nrow(d)
  npar <- n_distinct(d$stn) + n_distinct(d$ry)
  if (n - npar < MIN_RESID_DF) return(NULL)
  e <- d$y - d[[pc]]
  if (!reverse) { se <- seq_ss(e, d$stn, d$ry); so <- seq_ss(d$y, d$stn, d$ry) }
  else          { se <- seq_ss(e, d$ry, d$stn); so <- seq_ss(d$y, d$ry, d$stn) }
  tibble(family = d$family[1], n = n, resid_df = n - npar,
    bias = n * mean(e)^2,
    e_ad = if (reverse) se["second"] else se["first"],
    e_rt = if (reverse) se["first"]  else se["second"],
    e_st = se["resid"], e_tot = sum(e^2),
    o_ad = if (reverse) so["second"] else so["first"],
    o_rt = if (reverse) so["first"]  else so["second"],
    o_st = so["resid"], o_tot = sum((d$y - mean(d$y))^2))
}

run <- function(pc, reverse = FALSE)
  p %>% group_split(family) %>% map_dfr(~ decomp_family(.x, pc, reverse))

cat("=== DBO RAW biomass + zeros, train 2001-2010 / test 2011-2019 ===\n")
cat("Component-wise SKILL = 1 - SSE/SS_obs  (pooled: sum SSE / sum SS_obs across families)\n")

chk <- run(mods[[1]])
cat(sprintf("\nfamilies decomposed: %d (median residual df %.0f)\n", nrow(chk), median(chk$resid_df)))
cat(sprintf("ADDITIVITY: max|SSE_total - sum(components)| = %.2e\n",
  with(chk, max(abs(e_tot - (bias + e_ad + e_rt + e_st))))))
cat(sprintf("\nobserved variation partition (pooled): abundance distribution %.0f%% | regional %.0f%% | station %.0f%%\n",
  100*sum(chk$o_ad)/sum(chk$o_tot), 100*sum(chk$o_rt)/sum(chk$o_tot), 100*sum(chk$o_st)/sum(chk$o_tot)))

sk <- function(r, a, b) 1 - sum(r[[a]]) / sum(r[[b]])
all_res <- map_dfr(names(mods), function(m) {
  r <- run(mods[[m]]); rr <- run(mods[[m]], reverse = TRUE)
  tibble(model = m,
    ad = sk(r, "e_ad", "o_ad"), rt = sk(r, "e_rt", "o_rt"),
    st = sk(r, "e_st", "o_st"), pooled = sk(r, "e_tot", "o_tot"),
    bias_pct = 100 * sum(r$bias) / sum(r$e_tot),
    ad_rev = sk(rr, "e_ad", "o_ad"), rt_rev = sk(rr, "e_rt", "o_rt"))
})
saveRDS(all_res, "data/rf_dbo_traintest_zeros_skill.rds")

cat(sprintf("\n%-9s %11s %11s %11s %11s %8s\n", "model",
  "abund.dist", "regional", "station", "POOLED", "bias%"))
pwalk(all_res, function(model, ad, rt, st, pooled, bias_pct, ad_rev, rt_rev)
  cat(sprintf("%-9s %+11.3f %+11.3f %+11.3f %+11.3f %7.1f%%\n",
    model, ad, rt, st, pooled, bias_pct)))

cat("\n-- sensitivity: region-year entered BEFORE station --\n")
cat(sprintf("%-9s %11s %11s\n", "model", "abund.dist", "regional"))
pwalk(all_res, function(model, ad, rt, st, pooled, bias_pct, ad_rev, rt_rev)
  cat(sprintf("%-9s %+11.3f %+11.3f\n", model, ad_rev, rt_rev)))
cat("\nsaved -> data/rf_dbo_traintest_zeros_skill.rds\n")
