# SFTS
R scripts and data for the manuscript titled "Violations of space-time equivalence limit the near-term forecasting ability of species distribution models".

**scripts/**

*fit_rf_bbs.R*: Fits random forests models to BBS data and calculate predictions performance of static, dynamic and decomposed SDMS.

*fit_rf_dbo.R*: Wrangles DBO data, then fits random forests models and calculate predictions performance of static, dynamic and decomposed SDMS.

*fit_rf_simulations.R*: Simulates data under different space-time equivalence scenarios and fits random forest models.

*plots.R*: Script to replicate plots in the manuscript.

*wrangle_bbs_data.R*: Wrangles BBS data from scract. It can install BSS, PRISM and elevation data. It will wrangle and combine these to make the data ready for *fit_rf_bbs.R*. This will take a while to run. You don't have to run this to replicate the analysis. The output of this code is already in the **data** folder.

**data/**: Contains the necessary data to run *fit_rf_bbs.R* and *fit_rf_dbo.R*.

**figures/**: Figures made by *plots.R* will be saved here.
