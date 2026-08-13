# SHAP attribution of the global and regional trajectory skill, dynamic model, total biomass.
#
# An exact alternative to the permutation version in plot_dyn_importance.R. TreeSHAP splits
# every prediction into a baseline plus one contribution per variable,
#
#     p_i = phi_0 + sum_j  phi_ij
#
# and the four-way decomposition is a set of linear projections. Both operations are linear
# and every projection reproduces a constant, so the constant baseline drops out and the
# component of the PREDICTION decomposes exactly into a sum over variables:
#
#     g(p) = sum_j g(phi_j)          g = the global (or a band's regional) projection
#
# That gives an exact per-variable split of the component SKILL. With
# SS = ||g(y)||^2 and SSE = ||g(y) - g(p)||^2,
#
#     skill = 1 - SSE/SS = sum_j [ 2<g(y), g(phi_j)> - <g(phi_j), g(p)> ] / SS
#
# so each variable's contribution is that bracket over SS, POSITIVE if the variable's
# contribution moves the component toward the observations and NEGATIVE if it moves it away.
# The contributions sum to the component skill, asserted below. No permutation, no error
# bars, no Monte Carlo noise -- unlike the permutation version, where the leading effects
# had standard deviations as large as the effects themselves.
#
# FOUR PANELS matching plot_dyn_trajectory.R: the network-wide global trajectory, then each
# band's regional trajectory scored on its own rows.
#
# What this does NOT do: SHAP attributes what the model USES. A variable can be given a
# large contribution and still be scored negatively here, which is the informative case --
# the model leans on it and is worse off for it out of sample.
#
# Predictors are strongly correlated (salinity with the nutrients, the sediment variables
# with each other), and SHAP divides credit among correlated variables in ways the data do
# not uniquely determine. Treat the ordering as indicative.
#
# Output: figures/dyn_shap.{png,pdf}, data/dyn_shap.rds

suppressMessages({ library(tidyverse); library(ranger); library(treeshap) })

set.seed(1)
NTREE <- 2000
base <- c("Temp","Salinity","integchla","sedchla","Ammonia","Phosphate",
          "NiTriTra","Silicate","phigte5","TOC","cn")
TERMS <- c("Depth", base)
NICE <- c(Depth = "Depth", Temp = "Bottom temperature", Salinity = "Bottom salinity",
          integchla = "Integrated chl a", sedchla = "Sediment chl a",
          Ammonia = "Ammonia", Phosphate = "Phosphate", NiTriTra = "Nitrite + nitrate",
          Silicate = "Silicate", phigte5 = "Silt and clay (>=5 phi)",
          TOC = "Sediment TOC", cn = "Sediment C:N")

env <- readRDS("data/decomp_dbo_traintest.rds") %>% select(-DBOreg)
tot <- readRDS("data/dbo_resp_total.rds")
d <- inner_join(tot, env, by = c("StationNme", "DataYear")) %>%
  transmute(stn = StationNme, reg = as.integer(as.character(DBOreg)), year = DataYear, set,
            y = biomass^0.25, across(all_of(TERMS))) %>%
  mutate(ry = paste(reg, year))
stopifnot(nrow(d) == nrow(tot), !anyNA(d))

tr <- filter(d, set == "train"); te <- filter(d, set == "test")
m <- ranger(as.formula(paste("y ~", paste(TERMS, collapse = " + "))), data = tr,
            num.trees = NTREE, seed = 1, num.threads = 0)
p  <- as.numeric(predict(m, data = as.data.frame(te[TERMS]))$predictions)

sh <- treeshap(ranger.unify(m, as.data.frame(tr[TERMS])),
               as.data.frame(te[TERMS]), verbose = FALSE)$shaps
sh <- as.matrix(sh)[, TERMS, drop = FALSE]
stopifnot(max(abs(p - (mean(p) - mean(rowSums(sh)) + rowSums(sh)))) < 1e-8)

# ---- projections --------------------------------------------------------------------------
S <- factor(te$stn); Y <- factor(te$year); R <- factor(te$ry)
H <- function(f) { M <- model.matrix(f); M %*% MASS::ginv(M) }
H1 <- H(~ S); H2 <- H(~ S + Y); H3 <- H(~ S + Y + R)
g_global   <- function(v) as.numeric((H2 - H1) %*% v)
g_regional <- function(v) as.numeric((H3 - H2) %*% v)

regs <- sort(unique(te$reg))
PANELS <- c("Global", paste("DBO", regs))

# contribution of each variable to a component's skill, on a given set of rows
attribute <- function(gfun, keep) {
  gy <- gfun(te$y)[keep]; gp <- gfun(p)[keep]
  SS <- sum(gy^2); SSE <- sum((gy - gp)^2)
  gj <- map(TERMS, ~ gfun(sh[, .x])[keep])
  contrib <- map_dbl(gj, ~ (2 * sum(gy * .x) - sum(.x * gp)) / SS)
  stopifnot(abs(sum(contrib) - (1 - SSE / SS)) < 1e-8,          # contributions sum to skill
            max(abs(gp - reduce(gj, `+`))) < 1e-8)              # g(p) = sum_j g(phi_j)
  list(skill = 1 - SSE / SS, contrib = set_names(contrib, TERMS))
}

res <- c(list(Global = attribute(g_global, rep(TRUE, nrow(te)))),
         set_names(map(regs, ~ attribute(g_regional, te$reg == .x)), paste("DBO", regs)))

sk  <- map_dbl(res, "skill")
imp <- imap_dfr(res, ~ tibble(component = .y, var = TERMS, contrib = .x$contrib))
saveRDS(list(contrib = imp, skill = sk), "data/dyn_shap.rds")

cat("component skill (the contributions sum to these):\n")
print(round(sk, 3))

# ---- figure ---------------------------------------------------------------------------------
ord <- imp %>% group_by(var) %>% summarise(o = mean(contrib), .groups = "drop") %>% arrange(o)
lab <- tibble(component = names(sk),
              t = sprintf("component skill %+.3f", sk)) %>%
  mutate(component = factor(component, PANELS))
pd <- imp %>% mutate(var = factor(NICE[var], NICE[ord$var]),
                     component = factor(component, PANELS),
                     sign = contrib > 0)

g <- ggplot(pd, aes(contrib, var, fill = sign)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(width = 0.62) +
  geom_text(data = lab, aes(x = Inf, y = Inf, label = t), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.6, size = 2.6, colour = "grey25") +
  facet_wrap(~ component, nrow = 2, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = "#4E79A7", `FALSE` = "#E15759"), guide = "none") +
  labs(x = "Contribution to component skill  (sums to the component skill)", y = NULL,
       title = "Dynamic model, total biomass, train/test: SHAP attribution") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.25),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        plot.title = element_text(size = 10))

ggsave("figures/dyn_shap.png", g, width = 180, height = 150, units = "mm", dpi = 300)
ggsave("figures/dyn_shap.pdf", g, width = 180, height = 150, units = "mm")
cat("\nsaved -> figures/dyn_shap.png\n\n")

print(as.data.frame(imp %>% mutate(variable = NICE[var]) %>%
  arrange(factor(component, PANELS), desc(contrib)) %>%
  transmute(component, variable, contribution = round(contrib, 3))), row.names = FALSE)
