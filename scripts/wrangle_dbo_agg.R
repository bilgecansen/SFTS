# Response tables for the biological-aggregation ladder.
#
#   family   57 families            data/dbo_resp_family.rds
#   class     8 classes             data/dbo_resp_class.rds
#   total     1 series (TotGC)      data/dbo_resp_total.rds
#
# All three sit on the SAME 16 station x 15 year grid (226 station-years) produced by
# scripts/wrangle_dbo_data.R, so the only thing that changes up the ladder is the
# taxonomic grain.
#
# JOIN KEY. Class biomass comes from the NOAA CLASS workbook, joined on uniqueID via
# data/dbo_grabkey.rds rather than on StationNme + DataYear: 21 station-years have two
# grabs in the raw files and only one survives the filters in wrangle_dbo_data.R.
# Joining on station-year would silently duplicate those rows.
#
# GATES. The class level uses the same gates as the family level: a taxon x station
# pair is kept if the taxon has >4 years with positive biomass there, and a taxon is
# kept if >4 stations survive that. Zeros are retained as observations of absence and
# never count toward sample size. Of the ten CLASS columns, gc_Others is excluded (a
# residual bin, not a class) and Echinoidea fails the gate (positive at one station,
# 91% zeros), leaving EIGHT classes.
#
# TOTAL is TotGC, the workbook's own column, which equals the sum of all ten class
# columns exactly -- so it includes gc_Others and Echinoidea. On the analysis grid the
# 57 retained families sum to a median 94% of TotGC, so total biomass is close to, but
# not identical with, the assemblage the family models see.
#
# Carbon is concentrated: median per-station-year shares are Bivalvia 51%, Polychaeta
# 19%, Crustacea 3.5%, and under 0.3% for each of the remaining five classes. Pooled
# skill at the class level is therefore effectively a two-taxon result.
#
# SCALE. These tables are raw grams carbon. Summing happens here, on the raw scale;
# the fourth root is applied downstream, after aggregation.

suppressMessages({
  library(tidyverse)
  library(readxl)
})

CLASS_FILE <- "data/AllStation_BenthicInfauna_CLASS_ALL3_2000-2019_NOAA_abwk122023.xlsx"
DROP_CLASS <- "Others"          # residual bin, excluded from the class level
MIN_YEARS  <- 5                 # positive years required per taxon x station
MIN_STN    <- 5                 # stations required per taxon

key <- readRDS("data/dbo_grabkey.rds")          # uniqueID -> StationNme, DataYear, DBOreg

# gate a long taxon x station x year table, exactly as wrangle_dbo_data.R gates families
gate <- function(d) {
  ok_pair <- d %>% group_by(taxon, StationNme) %>%
    summarise(n = sum(biomass > 0), .groups = "drop") %>%
    filter(n >= MIN_YEARS) %>% mutate(ts = paste(taxon, StationNme))
  d <- filter(d, paste(taxon, StationNme) %in% ok_pair$ts)
  ok_tax <- d %>% group_by(taxon) %>%
    summarise(n = n_distinct(StationNme), .groups = "drop") %>% filter(n >= MIN_STN)
  filter(d, taxon %in% ok_tax$taxon)
}

# ---- family: reuse the existing gated panel ------------------------------------------
fam <- readRDS("data/data_dbo_zeros.rds") %>%
  transmute(StationNme, DataYear, DBOreg, taxon = family, biomass)
saveRDS(fam, "data/dbo_resp_family.rds")

# ---- class and total ------------------------------------------------------------------
cw <- suppressWarnings(read_excel(CLASS_FILE)) %>%
  select(uniqueID = uniqueid_gc, starts_with("gc_"), TotGC) %>%
  inner_join(key, by = "uniqueID")
stopifnot(nrow(cw) == nrow(key), !anyDuplicated(cw$uniqueID))

# TotGC is the sum of the ten class columns
stopifnot(max(abs(cw$TotGC - rowSums(select(cw, starts_with("gc_"))))) < 1e-8)

cls_all <- cw %>%
  pivot_longer(starts_with("gc_"), names_to = "taxon", values_to = "biomass") %>%
  mutate(taxon = sub("^gc_", "", taxon)) %>%
  filter(!taxon %in% DROP_CLASS) %>%
  select(StationNme, DataYear, DBOreg, taxon, biomass)

cls <- gate(cls_all)
saveRDS(cls, "data/dbo_resp_class.rds")

tot <- cw %>% transmute(StationNme, DataYear, DBOreg, taxon = "Total", biomass = TotGC)
stopifnot(nrow(gate(tot)) == nrow(tot))
saveRDS(tot, "data/dbo_resp_total.rds")

# ---- report ---------------------------------------------------------------------------
report <- function(d, lab) cat(sprintf(
  "%-7s %3d taxa | %5d rows | %d station-years | %5.1f%% zeros\n",
  lab, n_distinct(d$taxon), nrow(d), n_distinct(paste(d$StationNme, d$DataYear)),
  100 * mean(d$biomass == 0)))
report(fam, "family"); report(cls, "class"); report(tot, "total")

cat("\ndropped classes:",
    paste(setdiff(unique(cls_all$taxon), unique(cls$taxon)), collapse = ", "),
    "(gate) +", DROP_CLASS, "(excluded)\n")
cat("\nclass shares of TotGC (median %):\n")
cls_all %>% group_by(StationNme, DataYear) %>%
  mutate(sh = 100 * biomass / sum(biomass)) %>% ungroup() %>%
  group_by(taxon) %>% summarise(median_pct = round(median(sh), 2),
                                kept = first(taxon) %in% unique(cls$taxon)) %>%
  arrange(desc(median_pct)) %>% as.data.frame() %>% print(row.names = FALSE)

cat(sprintf("\nfamily sum / TotGC: median %.3f\n",
    median((fam %>% group_by(StationNme, DataYear) %>% summarise(s = sum(biomass), .groups="drop") %>%
            inner_join(tot, by = c("StationNme","DataYear")) %>% mutate(r = s/biomass))$r)))
cat("\nsaved -> data/dbo_resp_{family,class,total}.rds\n")
