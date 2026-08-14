# SHAP dependence: how the dynamic model's biomass response varies with the three variables
# that carry the most component skill in figures/dyn_shap.png.
#
# For each held-out observation TreeSHAP gives phi_ij, the signed amount variable j moves that
# prediction away from the baseline. Plotting phi_ij against x_ij traces the model's fitted
# response to variable j, on the fourth-root biomass scale, in the units of the prediction.
#
# WHY THIS AND NOT A PARTIAL DEPENDENCE PLOT. A PDP evaluates the model on the full crossing of
# variable j with the observed values of everything else, which includes predictor combinations
# that never occur. SHAP dependence is evaluated only AT the data, so it never leaves the
# observed envelope. The vertical spread at a given x is real information -- it is the part of
# variable j's effect that depends on the other predictors -- where a PDP would average it away.
#
# Same fit as plot_dyn_shap.R: dynamic model, total biomass, trained on the training block and
# explained on the test block, seed 1. The contributions here are the same phi that the bars in
# dyn_shap.png aggregate, so the two figures are views of one object.
#
# Env values are centred on the TRAINING grand mean (as stored in decomp_dbo_traintest.rds),
# so x = 0 is the training-period average, not a physical zero. Depth is uncentred but is not
# one of the three shown.
#
# Output: figures/dyn_shap_dependence.{png,pdf}, data/dyn_shap_dependence.rds

suppressMessages({ library(tidyverse); library(ranger); library(treeshap) })

set.seed(1)
NTREE <- 2000
base <- c("Temp","Salinity","integchla","sedchla","Ammonia","Phosphate",
          "NiTriTra","Silicate","phigte5","TOC","cn")
TERMS <- c("Depth", base)

# the four that carry component skill most consistently in dyn_shap.png: salinity dominates the
# global trajectory, ammonia leads DBO 1 and 2, sediment chl a leads DBO 3, and nitrite+nitrate
# is positive in all four components without ever being the largest. Salinity, ammonia and
# nitrite+nitrate are the only variables positive in every component.
FOCUS <- c("Salinity", "Ammonia", "sedchla", "NiTriTra")
NICE  <- c(Salinity = "Bottom salinity (train-centred)",
           Ammonia  = "Ammonia (train-centred)",
           sedchla  = "Sediment chl a (train-centred)",
           NiTriTra = "Nitrite + nitrate (train-centred)")
COLR  <- c("DBO 1" = "#4E79A7", "DBO 2" = "#F28E2B", "DBO 3" = "#59A14F")

env <- readRDS("data/decomp_dbo_traintest.rds") %>% select(-DBOreg)
tot <- readRDS("data/dbo_resp_total.rds")
d <- inner_join(tot, env, by = c("StationNme", "DataYear")) %>%
  transmute(stn = StationNme, reg = as.integer(as.character(DBOreg)), year = DataYear, set,
            y = biomass^0.25, across(all_of(TERMS)))
stopifnot(nrow(d) == nrow(tot), !anyNA(d))

tr <- filter(d, set == "train"); te <- filter(d, set == "test")
m <- ranger(as.formula(paste("y ~", paste(TERMS, collapse = " + "))), data = tr,
            num.trees = NTREE, seed = 1, num.threads = 0)
p <- as.numeric(predict(m, data = as.data.frame(te[TERMS]))$predictions)

sh <- treeshap(ranger.unify(m, as.data.frame(tr[TERMS])),
               as.data.frame(te[TERMS]), verbose = FALSE)$shaps
sh <- as.matrix(sh)[, TERMS, drop = FALSE]
stopifnot(max(abs(p - (mean(p) - mean(rowSums(sh)) + rowSums(sh)))) < 1e-8)

pd <- map_dfr(FOCUS, ~ tibble(var = .x, x = te[[.x]], phi = sh[, .x],
                              region = paste("DBO", te$reg), stn = te$stn, year = te$year)) %>%
  mutate(var = factor(NICE[var], NICE[FOCUS]))
saveRDS(pd, "data/dyn_shap_dependence.rds")

g <- ggplot(pd, aes(x, phi)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.9, se = FALSE,
              colour = "grey30", linewidth = 1.3, linetype = "22") +
  geom_point(aes(colour = region), size = 2.6, alpha = 0.75) +
  facet_wrap(~ var, nrow = 2, scales = "free_x") +
  scale_colour_manual(values = COLR, name = NULL) +
  labs(x = NULL, y = "SHAP contribution to fourth-root biomass") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        axis.text = element_text(size = 11), axis.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        strip.text = element_text(size = 11))

ggsave("figures/dyn_shap_dependence.png", g, width = 180, height = 165, units = "mm", dpi = 300)
ggsave("figures/dyn_shap_dependence.pdf", g, width = 180, height = 165, units = "mm")
cat("saved -> figures/dyn_shap_dependence.png\n\n")

cat("monotonicity of the fitted response (Spearman of phi on x), n =", nrow(te), "test rows\n")
print(as.data.frame(pd %>% group_by(variable = var) %>%
  summarise(all = cor(x, phi, method = "spearman"),
            `DBO 1` = cor(x[region == "DBO 1"], phi[region == "DBO 1"], method = "spearman"),
            `DBO 2` = cor(x[region == "DBO 2"], phi[region == "DBO 2"], method = "spearman"),
            `DBO 3` = cor(x[region == "DBO 3"], phi[region == "DBO 3"], method = "spearman"),
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

cat("\nrange of the fitted response (max phi - min phi), fourth-root biomass units\n")
print(as.data.frame(pd %>% group_by(variable = var) %>%
  summarise(span = round(diff(range(phi)), 3), sd = round(sd(phi), 3), .groups = "drop")),
  row.names = FALSE)
