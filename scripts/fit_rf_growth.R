# Growth RF models --------------------------------------------------------

growth <- bbs_train3 %>%
  group_by(species_id, site_id) %>%
  reframe(
    growth = abundance[2:length(abundance)] /
      abundance[1:(length(abundance) - 1)],
    year_diff = year[2:length(year)] -
      year[1:(length(year) - 1)],
    effort_diff = party_hours[2:length(year)] /
      party_hours[1:(length(year) - 1)],
    year = year[1:(length(year) - 1)]
  ) %>%
  filter(year_diff == 1) %>%
  select(-year_diff)

bbs_growth <- left_join(
  growth,
  select(bbs_train3, -abundance),
  by = c("species_id", "site_id", "year")
)

# Route level
var_imp_growth <- list()
R2_growth <- c()
pb <- txtProgressBar(max = length(species), style = 3)
for (i in 1:length(species)) {
  bbs_species <- filter(bbs_growth, species_id == species[i])

  data_rf <- select(
    bbs_species,
    growth,
    elevs,
    lat,
    long,
    effort_diff,
    contains(c("spatial", "temporal", "residual")),
    -contains("lag_spatial")
  )

  res_rf <- ranger(
    log(growth) ~ .,
    data = data_rf,
    num.trees = 500,
    importance = "permutation"
  )

  R2_growth[i] <- res_rf$r.squared
  var_imp_growth[[i]] <- res_rf$variable.importance

  setTxtProgressBar(pb, i)
}
