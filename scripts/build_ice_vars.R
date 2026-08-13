# Build the ice variables from the saved daily series.
#
# THE ICE YEAR. For a station sampled in July of year t, the relevant ice season is the
# winter BEFORE that cruise. The previous persistence variable counted ice days over the
# calendar year, so roughly a third of the count came from November and December of year t
# -- the freeze-up of the following winter, after the biology had been collected. Every
# variable here is defined on the ice year
#
#     1 September (t-1)  ->  31 August (t)
#
# so the whole covariate precedes the July cruise. per_jun is the same count truncated at
# 30 June, which drops the two months that follow the cruise; the two differ only if ice is
# present in July or August, which is reported below.
#
# THE 15% THRESHOLD is the standard presence/absence criterion. Events are taken from the
# LONGEST CONTINUOUS RUN of ice days in the ice year rather than from the first crossing
# after a fixed calendar date. Fixed search dates fail here for two reasons: DBO 1 is the
# St Lawrence Island Polynya and is already open on 1 March in 12% of station-years, and
# midwinter leads open briefly elsewhere (UTBS1 falls to 18% in late February 2019 and
# returns to 78% in March). Both would be read as the seasonal retreat. The longest run is
# the main ice season and is unaffected by short excursions either side of it.
#
#   freezeup   first day of the longest ice run,
#              reported as days after 1 September -- larger means later freeze-up
#   breakup    day after the last day of that run,
#              reported as day of year in t -- larger means later breakup
#   ice_season breakup date minus freezeup date, in days
#   per_ice    days >= 15% over the whole ice year
#   per_jun    days >= 15% from 1 Sep (t-1) to 30 Jun (t)
#   per_cal    the OLD calendar-year count, kept so the change can be seen
#
# The ice record starts 1 Jan 2000, so ice year 2000 is incomplete and dropped; 2001-2019
# are complete, which covers every sampled year.
#
# Output: data/data_sic_vars.rds   StationNme, DataYear, and the six variables above

suppressMessages({ library(tidyverse) })

THR <- 150                                        # 15%, values scaled by 10

d <- readRDS("data/sic_daily.rds") %>% arrange(StationNme, date) %>%
  mutate(ice = sic >= THR,
         yr  = as.integer(format(date, "%Y")),
         mo  = as.integer(format(date, "%m")),
         doy = as.integer(format(date, "%j")),
         ice_year = ifelse(mo >= 9, yr + 1L, yr))  # Sep-Dec belong to the NEXT ice year

# start and end of the longest continuous run of TRUE
longest_run <- function(flag) {
  r <- rle(flag)
  k <- which(r$values)
  if (!length(k)) return(c(NA_integer_, NA_integer_))
  j <- k[which.max(r$lengths[k])]
  end <- cumsum(r$lengths)[j]
  c(end - r$lengths[j] + 1L, end)
}

LAST <- 2019L                     # ice year 2020 would need Jan-Aug 2020
vars <- d %>% filter(ice_year >= 2001, ice_year <= LAST) %>%
  group_split(StationNme, ice_year) %>%
  map_dfr(function(x) {
    t   <- x$ice_year[1]
    sep <- as.Date(sprintf("%d-09-01", t - 1))
    jun <- as.Date(sprintf("%d-06-30", t))
    rr  <- longest_run(x$ice)                      # the main ice season
    fu  <- x$date[rr[1]]                           # first day of ice
    bu  <- x$date[rr[2]] + 1                       # first day of open water after it
    tibble(StationNme = x$StationNme[1], DataYear = t,
           freezeup   = as.numeric(fu - sep),
           breakup    = as.integer(format(bu, "%j")),
           ice_season = as.numeric(bu - fu),
           per_ice    = sum(x$ice),
           per_jun    = sum(x$ice[x$date <= jun]))
  }) %>%
  # a station that never froze has no freeze-up or breakup date, and an ice season of zero
  mutate(no_ice     = per_ice == 0,
         freezeup   = ifelse(no_ice, NA_real_, freezeup),
         breakup    = ifelse(no_ice, NA_integer_, breakup),
         ice_season = ifelse(no_ice, 0, ice_season))

# the old calendar-year count, for comparison
cal <- d %>% group_by(StationNme, DataYear = yr) %>%
  summarise(per_cal = sum(ice), .groups = "drop")
vars <- left_join(vars, cal, by = c("StationNme", "DataYear"))

cat(sprintf("%d station-years, %d stations, %d-%d\n", nrow(vars),
            n_distinct(vars$StationNme), min(vars$DataYear), max(vars$DataYear)))
cat(sprintf("station-years with NO ice at all: %d  (%s)\n", sum(vars$no_ice),
            paste(sprintf("%s %d", vars$StationNme[vars$no_ice], vars$DataYear[vars$no_ice]),
                  collapse = ", ")))
cat(sprintf("missing freezeup %d, breakup %d\n",
            sum(is.na(vars$freezeup)), sum(is.na(vars$breakup))))
cat(sprintf("ice days in July-August (the months per_ice includes and per_jun does not): %d of %d station-years have any; median %d\n",
            sum(vars$per_ice - vars$per_jun > 0), nrow(vars),
            median(vars$per_ice - vars$per_jun)))

# agreement with the file this replaces
old <- readRDS("data/data_sic_per.rds")
chk <- inner_join(vars, rename(old, per_old = per), by = c("StationNme", "DataYear"))
cat(sprintf("\nvs data_sic_per.rds: per_cal matches the old value in %d of %d rows (max diff %d)\n",
            sum(chk$per_cal == chk$per_old), nrow(chk), max(abs(chk$per_cal - chk$per_old))))
cat(sprintf("correlation of the fixed count with the old one: %.3f | mean shift %+.1f days\n",
            cor(chk$per_ice, chk$per_old), mean(chk$per_ice - chk$per_old)))

saveRDS(vars, "data/data_sic_vars.rds")
cat("saved -> data/data_sic_vars.rds\n")

# ---- what the variables look like -------------------------------------------------------
tot <- readRDS("data/dbo_resp_total.rds") %>%
  distinct(StationNme, DataYear, DBOreg) %>%
  mutate(reg = as.integer(as.character(DBOreg)))
v <- inner_join(vars, tot, by = c("StationNme", "DataYear"))

cat("\n=== region means by year ===\n")
for (nm in c("per_ice", "per_cal", "breakup", "freezeup", "ice_season")) {
  cat(sprintf("\n%s\n", nm))
  print(as.data.frame(v %>% group_by(reg, DataYear) %>%
    summarise(m = round(mean(.data[[nm]], na.rm = TRUE)), .groups = "drop") %>%
    pivot_wider(names_from = reg, values_from = m, names_prefix = "DBO")),
    row.names = FALSE)
}

# ---- figure ------------------------------------------------------------------------------
suppressMessages(library(patchwork))
COLR <- c("DBO 1" = "#4E79A7", "DBO 2" = "#F28E2B", "DBO 3" = "#59A14F")
LAB <- c(per_ice = "Ice days in the ice year", ice_season = "Ice season length (days)",
         freezeup = "Freeze-up (days after 1 Sep)", breakup = "Breakup (day of year)")
pd <- v %>% mutate(region = paste("DBO", reg)) %>%
  select(region, DataYear, all_of(names(LAB))) %>%
  pivot_longer(-c(region, DataYear), names_to = "var", values_to = "val") %>%
  group_by(region, DataYear, var) %>%
  summarise(val = mean(val, na.rm = TRUE), .groups = "drop") %>%
  mutate(var = factor(LAB[var], LAB))

g <- ggplot(pd, aes(DataYear, val, colour = region)) +
  geom_vline(xintercept = 2012.5, linetype = "22", colour = "grey45", linewidth = 0.4) +
  geom_line(linewidth = 0.8, na.rm = TRUE) + geom_point(size = 1.7, na.rm = TRUE) +
  facet_wrap(~ var, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = COLR, name = NULL) +
  scale_x_continuous(breaks = seq(2002, 2018, 4)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "grey95", colour = "grey80"))
ggsave("figures/ice_vars.png", g, width = 185, height = 130, units = "mm", dpi = 300)
cat("\nsaved -> figures/ice_vars.png\n")
