# CANONICAL main BBS figures: the four "two-axis" figures (metric x scale).
# x = model type (relaxing STE ->), rows = treatment (harder to transfer, down),
# GAM vs RF by colour. See make_2axis() in plot_helpers.R. Writes to figures/.
source("scripts/plot_helpers.R")

make_2axis("sw", "strata", "species_wide_strata")
make_2axis("sw", "site",   "species_wide_site")
make_2axis("pl", "strata", "population_level_strata")
make_2axis("pl", "site",   "population_level_site")
