# What drives SPATIAL community turnover, and how -- the I-splines of the spatial-turnover GDM.
#
# GDM predicts compositional dissimilarity between a pair of sites from the sum of per-variable
# differences, each passed through a MONOTONE INCREASING function f_j fitted as an I-spline:
#
#     predicted dissimilarity = link( sum_j [ f_j(x2_j) - f_j(x1_j) ] )
#
# The two things asked of a GDM come straight out of those functions:
#
#   IMPORTANCE   the MAXIMUM HEIGHT of f_j is the total amount of compositional turnover
#                associated with variable j across its observed range, holding the others
#                fixed. This is the standard GDM importance measure, and unlike the RF
#                measures it is on an interpretable scale -- units of ecological distance.
#
#   RELATIONSHIP the SHAPE of f_j says WHERE along the gradient turnover happens. A steep
#                stretch is a part of the gradient across which communities change fast; a
#                flat stretch is one across which they do not change at all. Because f_j is
#                constrained monotone, the shape carries all the information -- there is no
#                sign to read, only rate.
#
# THE MODEL. One GDM fitted to ALL spatial pairs (same year, different station), the same
# class and the same 13 predictors as the `Spatial turnover` model in fit_gdm_dbo.R, plus
# geographic distance. This is the descriptive full-data fit, not a holdout fit -- it is what
# the fitted relationship IS, and the held-out skill for it lives in gdm_dbo_region_*.png.
#
# STABILITY. The 16 leave-one-station-out refits are drawn behind each spline as thin lines.
# With 16 stations a single influential station can carry a spline, and the spread shows
# directly which fitted relationships survive dropping one and which do not.
#
# Env variables are globally centred (as stored in data_dbo_zeros.rds); Depth (m), breakup
# (day of year) and geographic distance are in natural units.
#
# Output: figures/gdm_dbo_splines.{png,pdf}, data/gdm_dbo_splines.rds

suppressMessages({ library(gdm); library(tidyverse); library(vegan); library(patchwork) })

set.seed(1)
base <- c("Temp", "Salinity", "integchla", "sedchla", "Ammonia", "Phosphate",
          "NiTriTra", "Silicate", "phigte5", "TOC", "cn")
ICE <- "breakup"
ENV <- c(base, ICE)
VARS <- c(ENV, "Depth")
NTOP <- 4                                    # splines shown in the right-hand grid
# Geographic distance and Depth stay in the importance ranking but are not shown as curves:
# neither is an environmental gradient one could act on or project forward, and geographic
# distance in particular is a stand-in for everything unmeasured that varies with separation.
NOSPLINE <- c("Geographic", "Depth")

NICE <- c(Geographic = "Geographic distance", Depth = "Depth",
          Temp = "Bottom temperature", Salinity = "Bottom salinity",
          integchla = "Integrated chl a", sedchla = "Sediment chl a",
          Ammonia = "Ammonia", Phosphate = "Phosphate", NiTriTra = "Nitrite + nitrate",
          Silicate = "Silicate", phigte5 = "Silt and clay (>=5 phi)",
          TOC = "Sediment TOC", cn = "Sediment C:N", breakup = "Ice breakup (day of year)")

# ---- site table and pairs, identical to fit_gdm_dbo.R -------------------------------------
z <- readRDS("data/data_dbo_zeros.rds")
ice <- readRDS("data/data_sic_vars.rds") %>%
  select(StationNme, DataYear, all_of(ICE)) %>% group_by(StationNme) %>%
  mutate(breakup = ifelse(is.na(breakup), median(breakup, na.rm = TRUE), breakup)) %>%
  ungroup()

key <- z %>% distinct(StationNme, DataYear) %>% arrange(StationNme, DataYear)
sord <- key %>%
  inner_join(z %>% distinct(StationNme, DataYear, .keep_all = TRUE) %>%
               select(StationNme, DataYear, DBOreg, Latitude, Longitude, Depth, all_of(base)),
             by = c("StationNme", "DataYear")) %>%
  inner_join(ice, by = c("StationNme", "DataYear")) %>% as.data.frame()
stopifnot(nrow(sord) == nrow(key), !anyNA(sord[VARS]))

comm <- readRDS("data/dbo_resp_family.rds") %>% mutate(v = biomass^0.25) %>%
  select(StationNme, DataYear, taxon, v) %>%
  pivot_wider(names_from = taxon, values_from = v, values_fill = 0) %>%
  right_join(key, by = c("StationNme", "DataYear")) %>%
  arrange(StationNme, DataYear) %>% mutate(across(-(1:2), ~ replace_na(.x, 0)))
stopifnot(identical(paste(comm$StationNme, comm$DataYear),
                    paste(sord$StationNme, sord$DataYear)))
D <- as.matrix(vegdist(as.matrix(comm[, -(1:2)]), method = "bray"))

ij <- t(combn(nrow(sord), 2))
spat <- tibble(ki = ij[, 1], kj = ij[, 2], distance = D[cbind(ij[, 1], ij[, 2])],
               s1 = sord$StationNme[ki], s2 = sord$StationNme[kj],
               y1 = sord$DataYear[ki],   y2 = sord$DataYear[kj]) %>%
  filter(is.finite(distance), y1 == y2, s1 != s2)

# rownames MUST be reset -- gdm() segfaults on the mangled ones repeated indexing produces
gdm_table <- function(pr, vars) {
  v1 <- sord[pr$ki, vars, drop = FALSE]; names(v1) <- paste0("s1.", vars)
  v2 <- sord[pr$kj, vars, drop = FALSE]; names(v2) <- paste0("s2.", vars)
  rownames(v1) <- NULL; rownames(v2) <- NULL
  out <- data.frame(distance = pr$distance, weights = 1,
    s1.xCoord = sord$Longitude[pr$ki], s1.yCoord = sord$Latitude[pr$ki],
    s2.xCoord = sord$Longitude[pr$kj], s2.yCoord = sord$Latitude[pr$kj],
    v1, v2, check.names = FALSE)
  rownames(out) <- NULL
  class(out) <- c("gdmData", "data.frame")
  out
}

# ---- full fit, then the 16 LOSO refits ----------------------------------------------------
m <- gdm(gdm_table(spat, VARS), geo = TRUE)
cat(sprintf("spatial pairs %d | deviance explained %.1f%% | intercept %.3f\n\n",
            nrow(spat), m$explained, m$intercept))

splines_of <- function(mod) {
  s <- isplineExtract(mod)
  map_dfr(colnames(s$x), ~ tibble(var = .x, x = s$x[, .x], f = s$y[, .x]))
}

full <- splines_of(m) %>% mutate(fit = "full")
loso <- map_dfr(sort(unique(sord$StationNme)), function(st) {
  pr <- filter(spat, s1 != st, s2 != st)
  mm <- tryCatch(gdm(gdm_table(pr, VARS), geo = TRUE), error = function(e) NULL)
  if (is.null(mm)) return(NULL)
  splines_of(mm) %>% mutate(fit = st)
})

imp <- full %>% group_by(var) %>% summarise(height = max(f), .groups = "drop") %>%
  arrange(desc(height)) %>%
  mutate(share = 100 * height / sum(height), label = NICE[var])
loso_imp <- loso %>% group_by(fit, var) %>% summarise(height = max(f), .groups = "drop")
saveRDS(list(splines = full, loso = loso, importance = imp,
             loso_importance = loso_imp, explained = m$explained),
        "data/gdm_dbo_splines.rds")

# ---- figure ------------------------------------------------------------------------------
BAR <- "#254F9A"
ord <- imp$label

pa <- ggplot(imp %>% mutate(label = factor(label, rev(ord))), aes(height, label)) +
  geom_col(fill = BAR, alpha = 0.85, width = 0.7) +
  geom_point(data = loso_imp %>% mutate(label = factor(NICE[var], rev(ord))),
             aes(height, label), colour = "grey25", size = 0.55, alpha = 0.5) +
  labs(x = "Total turnover (I-spline height)", y = NULL,
       title = "a  Variable importance") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.25),
        plot.title = element_text(size = 10, face = "bold"))

top <- setdiff(imp$var, NOSPLINE)[seq_len(NTOP)]
pb <- ggplot(loso %>% filter(var %in% top) %>%
               mutate(label = factor(NICE[var], NICE[top])),
             aes(x, f)) +
  geom_line(aes(group = fit), colour = "grey65", linewidth = 0.3, alpha = 0.65) +
  geom_line(data = full %>% filter(var %in% top) %>%
              mutate(label = factor(NICE[var], NICE[top])),
            colour = BAR, linewidth = 1.1) +
  facet_wrap(~ label, nrow = 2, scales = "free_x") +
  # free_x puts each panel's outermost tick label hard against the panel edge, where it
  # collides with its neighbour's -- fewer breaks plus panel spacing keeps them apart
  scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.04)) +
  labs(x = NULL, y = "Compositional turnover  f(x)",
       title = sprintf("b  Fitted relationships, %d strongest environmental variables", NTOP),
       subtitle = "grey = the 16 leave-one-station-out refits") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        strip.text = element_text(size = 9),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        panel.spacing.x = unit(5, "mm"))

g <- pa + pb + plot_layout(widths = c(1, 1.5))
ggsave("figures/gdm_dbo_splines.png", g, width = 215, height = 115, units = "mm", dpi = 300)
ggsave("figures/gdm_dbo_splines.pdf", g, width = 215, height = 115, units = "mm")
cat("saved -> figures/gdm_dbo_splines.png\n\n")

cat("I-spline heights: total turnover attributable to each variable\n")
print(as.data.frame(imp %>%
  left_join(loso_imp %>% group_by(var) %>%
              summarise(loso_min = min(height), loso_max = max(height),
                        n_zero = sum(height < 1e-8), .groups = "drop"), by = "var") %>%
  transmute(variable = label, height = round(height, 3), `share_%` = round(share, 1),
            loso_min = round(loso_min, 3), loso_max = round(loso_max, 3),
            loso_folds_zero = n_zero)), row.names = FALSE)
