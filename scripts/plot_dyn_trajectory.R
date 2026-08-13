# Global and regional trajectory components of total biomass, observed and predicted by the
# DYNAMIC model, across the training and test blocks.
#
# The dynamic model is Depth plus the eleven in-situ variables as raw spatiotemporal values
# -- no decomposition, no coordinates, no ice. It is the model that carries the trajectory
# terms at the total rung under train/test (global +0.179, regional +0.055), which is why
# it is the one plotted here.
#
#   Network panel   the global trajectory: the year effect after station means
#   DBO panels      the regional trajectory: that band's departure from the year effect
#
# THE TWO HALVES ARE NOT THE SAME QUANTITY, and the dashed rule marks the boundary.
#   2001-2012   out-of-bag predictions, so each tree predicts only rows it did not see
#   2013-2019   ordinary predictions for a block the model never saw
#
# The decomposition is applied WITHIN EACH BLOCK, as the canonical analysis does, so the
# test-block numbers reproduce the ladder. A consequence is that each half is centred on its
# own means, so the two halves are not on a common baseline and the step at the boundary is
# not readable from this figure.
#
# 2002 and 2008-2010 were not sampled, so the lines break there.
#
# A second figure, dyn_series.png, shows the SAME model and the same panels without any
# demeaning: the raw network and region mean series, which retain the level and the 2013
# step that the component view removes.
#
# Output: figures/dyn_trajectory.{png,pdf}, figures/dyn_series.{png,pdf}

suppressMessages({ library(tidyverse); library(ranger) })

set.seed(1)
NTREE <- 2000
OBS <- "grey15"; PRD <- "#DD9100"
base <- c("Temp","Salinity","integchla","sedchla","Ammonia","Phosphate",
          "NiTriTra","Silicate","phigte5","TOC","cn")
TERMS <- c("Depth", base)

env <- readRDS("data/decomp_dbo_traintest.rds") %>% select(-DBOreg)
tot <- readRDS("data/dbo_resp_total.rds")

d <- inner_join(tot, env, by = c("StationNme", "DataYear")) %>%
  transmute(stn = StationNme, reg = as.integer(as.character(DBOreg)), year = DataYear, set,
            y = biomass^0.25, across(all_of(TERMS)))
stopifnot(nrow(d) == nrow(tot), !anyNA(d))

tr <- filter(d, set == "train"); te <- filter(d, set == "test")
m <- ranger(as.formula(paste("y ~", paste(TERMS, collapse = " + "))), data = tr,
            num.trees = NTREE, seed = 1, num.threads = 0)
p <- bind_rows(mutate(tr, p = m$predictions),                       # out-of-bag
               mutate(te, p = as.numeric(predict(m, data = as.data.frame(te))$predictions))) %>%
  mutate(ry = paste(reg, year))
BOUND <- (max(tr$year) + min(te$year)) / 2

# global and regional components, computed within each block
comp <- function(v, stn, yr, ry) {
  f1 <- fitted(lm(v ~ factor(stn)))
  f2 <- fitted(lm(v ~ factor(stn) + factor(yr)))
  f3 <- fitted(lm(v ~ factor(stn) + factor(yr) + factor(ry)))
  list(global = as.numeric(f2 - f1), regional = as.numeric(f3 - f2))
}
cc <- p %>% group_split(set) %>% map_dfr(function(b) {
  a <- comp(b$y, b$stn, b$year, b$ry); e <- comp(b$p, b$stn, b$year, b$ry)
  tibble(reg = b$reg, year = b$year, set = b$set,
         og = a$global, pg = e$global, orr = a$regional, prr = e$regional)
})

PAN <- c("Global", "DBO 1", "DBO 2", "DBO 3")
ser <- bind_rows(
  cc %>% group_by(year, set) %>%
    summarise(Observed = mean(og), Predicted = mean(pg), .groups = "drop") %>%
    mutate(panel = "Global"),
  cc %>% group_by(panel = paste("DBO", reg), year, set) %>%
    summarise(Observed = mean(orr), Predicted = mean(prr), .groups = "drop")) %>%
  mutate(panel = factor(panel, PAN))

make_fig <- function(ser, PAN, ylab, ttl, stem) {
  ser <- mutate(ser, panel = factor(panel, PAN))
  S <- ser %>% group_by(panel, set) %>%
    summarise(r = cor(Observed, Predicted), RMSE = sqrt(mean((Observed - Predicted)^2)),
              .groups = "drop")
  lab <- S %>%
    mutate(t = sprintf("%-8s r %+.2f   RMSE %.3f",
                       ifelse(set == "train", "fit", "forecast"), r, RMSE)) %>%
    group_by(panel) %>%
    summarise(lab = paste(t[set == "train"], t[set == "test"], sep = "\n"), .groups = "drop")

  plt <- ser %>% complete(panel, year = min(d$year):max(d$year)) %>%
    pivot_longer(c(Observed, Predicted), names_to = "series", values_to = "v")

  g <- ggplot(plt, aes(year, v, colour = series)) +
    geom_vline(xintercept = BOUND, linetype = "22", colour = "grey45", linewidth = 0.4) +
    geom_line(linewidth = 0.85, na.rm = TRUE) + geom_point(size = 1.9, na.rm = TRUE) +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = 0, vjust = 1.25, size = 2.2, colour = "grey25", lineheight = 1.2,
              family = "mono", nudge_x = 0.4) +
    facet_wrap(~ panel, nrow = 2, scales = "free_y") +
    scale_colour_manual(values = c(Observed = OBS, Predicted = PRD), name = NULL) +
    scale_x_continuous(breaks = seq(2002, 2018, 4)) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.26))) +
    labs(x = NULL, y = ylab, title = ttl) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 10),
          strip.background = element_rect(fill = "grey95", colour = "grey80"),
          legend.position = "bottom")
  ggsave(sprintf("figures/%s.png", stem), g, width = 190, height = 145, units = "mm", dpi = 300)
  ggsave(sprintf("figures/%s.pdf", stem), g, width = 190, height = 145, units = "mm")
  cat(sprintf("saved -> figures/%s.png\n", stem))
  S
}

S1 <- make_fig(ser, PAN, "Trajectory component, fourth-root biomass",
               "Dynamic model, total biomass: trajectory components",
               "dyn_trajectory")

# ---- the same panels without demeaning: raw network and region means ---------------------
raw <- bind_rows(
  p %>% group_by(year, set) %>%
    summarise(Observed = mean(y), Predicted = mean(p), .groups = "drop") %>%
    mutate(panel = "Network"),
  p %>% group_by(panel = paste("DBO", reg), year, set) %>%
    summarise(Observed = mean(y), Predicted = mean(p), .groups = "drop"))
S2 <- make_fig(raw, c("Network", "DBO 1", "DBO 2", "DBO 3"), "Fourth-root total biomass",
               "Dynamic model, total biomass: full series, no demeaning", "dyn_series")

cat("\n=== trajectory components ===\n")
print(as.data.frame(S1 %>% transmute(panel, block = ifelse(set == "train", "fit", "forecast"),
                                     r = round(r, 3), RMSE = round(RMSE, 3))), row.names = FALSE)
cat("\n=== full series, no demeaning ===\n")
print(as.data.frame(S2 %>% transmute(panel, block = ifelse(set == "train", "fit", "forecast"),
                                     r = round(r, 3), RMSE = round(RMSE, 3),
                                     .groups = NULL)), row.names = FALSE)
cat("\n=== the raw series ===\n")
print(as.data.frame(raw %>% mutate(panel = factor(panel, c("Network","DBO 1","DBO 2","DBO 3"))) %>%
  arrange(panel, year) %>%
  transmute(panel, year, block = ifelse(set == "train", "fit", "forecast"),
            observed = round(Observed, 3), predicted = round(Predicted, 3),
            error = round(Predicted - Observed, 3))), row.names = FALSE)
