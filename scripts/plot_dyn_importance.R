# Which variables carry the global and regional trajectories in the dynamic model?
#
# Standard permutation importance measures the increase in overall prediction error, which
# is dominated by the abundance distribution and says little about the trajectory terms.
# This script targets the components directly: a variable is permuted in the HELD-OUT rows,
# the predictions are recomputed, the error is decomposed into the usual four components,
# and importance is the DROP IN THAT COMPONENT'S SKILL.
#
# FOUR PANELS, matching plot_dyn_trajectory.R: the network-wide GLOBAL trajectory, then each
# band's REGIONAL trajectory scored on its own rows. The regional component is a per-row
# quantity, so its sum of squares partitions exactly by region for the observations and the
# error alike, and each band gets its own skill and its own importance ranking.
#
# Dynamic model, total biomass, TRAIN/TEST only: Depth plus the eleven in-situ variables
# as raw values.
# Models are fitted once on unpermuted data; permutation happens only at prediction time,
# so this measures what the fitted model USES, not what it would learn without the variable.
#
# NREP permutations per variable; the figure shows the mean and +/- 1 sd across them.
#
# Output: figures/dyn_importance.{png,pdf}, data/dyn_importance.rds

suppressMessages({ library(tidyverse); library(ranger) })

set.seed(1)
NTREE <- 2000; NREP <- 40
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

FRM <- as.formula(paste("y ~", paste(TERMS, collapse = " + ")))
regs <- sort(unique(d$reg))
CLABS <- c("Global", paste("DBO", regs))

# global trajectory (network) and the regional trajectory split by band
comp_skill <- function(x, p) {
  g <- function(v) {
    f1 <- fitted(lm(v ~ factor(x$stn)))
    f2 <- fitted(lm(v ~ factor(x$stn) + factor(x$year)))
    f3 <- fitted(lm(v ~ factor(x$stn) + factor(x$year) + factor(x$ry)))
    by_reg <- function(z) as.numeric(tapply(z, x$reg, sum)[as.character(regs)])
    c(sum((f2 - f1)^2), by_reg((f3 - f2)^2))
  }
  set_names(1 - g(x$y - p) / g(x$y), CLABS)
}

# ---- train/test --------------------------------------------------------------------------
tr <- filter(d, set == "train"); te <- filter(d, set == "test")
m_tt <- ranger(FRM, data = tr, num.trees = NTREE, seed = 1, num.threads = 0)
pred_tt <- function(x) as.numeric(predict(m_tt, data = as.data.frame(x))$predictions)
base_tt <- comp_skill(te, pred_tt(te))

imp_tt <- map_dfr(TERMS, function(v) map_dfr(seq_len(NREP), function(i) {
  x <- te; x[[v]] <- sample(x[[v]])
  tibble(var = v, rep = i, component = CLABS,
         drop = base_tt - comp_skill(te, pred_tt(x)))
}))

imp <- imp_tt
saveRDS(list(importance = imp, baseline = base_tt), "data/dyn_importance.rds")

cat("baseline skill: global trajectory and each band's regional trajectory\n")
print(round(base_tt, 3))

# ---- figure ---------------------------------------------------------------------------------
S <- imp %>% group_by(component, var) %>%
  summarise(m = mean(drop), s = sd(drop), .groups = "drop")
ord <- S %>% group_by(var) %>% summarise(o = mean(m), .groups = "drop") %>% arrange(o)
S <- S %>% mutate(var = factor(NICE[var], NICE[ord$var]),
                  component = factor(component, CLABS))

g <- ggplot(S, aes(m, var)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_linerange(aes(xmin = m - s, xmax = m + s), colour = "#E15759", linewidth = 0.5) +
  geom_point(size = 2.3, colour = "#E15759") +
  facet_wrap(~ component, nrow = 2, scales = "free_x") +
  labs(x = "Drop in component skill when the variable is permuted", y = NULL,
       title = "Dynamic model, total biomass, train/test") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(linewidth = 0.25),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        plot.title = element_text(size = 10))

ggsave("figures/dyn_importance.png", g, width = 180, height = 150, units = "mm", dpi = 300)
ggsave("figures/dyn_importance.pdf", g, width = 180, height = 150, units = "mm")
cat("\nsaved -> figures/dyn_importance.png\n\n")

print(as.data.frame(S %>% arrange(component, desc(m)) %>%
  transmute(component, variable = as.character(var),
            drop = round(m, 3), sd = round(s, 3))), row.names = FALSE)
