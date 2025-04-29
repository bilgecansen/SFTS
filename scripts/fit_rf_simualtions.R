library(tidyverse)
library(ranger)
library(foreach)
library(patchwork)
library(truncnorm)


# Univariate Relationships ------------------------------------------------

simulate_ab <- function(response, t_x = 0:19) {
  s_x <- 1:20
  m_x <- foreach(i = 1:20, .combine = "rbind") %do%
    {
      s_x + t_x[i]
    }
  m_x <- t(apply(m_x, c(1, 2), function(x) x + rnorm(1, 0, 1)))

  dat_x <- data.frame(m_x)
  colnames(dat_x) <- 1:20
  dat_x$site <- 1:20

  dat_x <- pivot_longer(dat_x, !site, names_to = "time", values_to = "env")
  dat_x$env_rand1 <- rnorm(nrow(dat_x), 0, 10)
  dat_x$env_rand2 <- rnorm(nrow(dat_x), 0, 10)
  dat_x$env_rand3 <- rnorm(nrow(dat_x), 0, 10)

  if (response == "spatial") {
    a <- 1000
    b <- 10
    mu_x <- apply(m_x, 2, mean)
    mu_ab <- a + b * mu_x

    m_ab <- lapply(mu_ab, function(x) rnorm(20, x, 10))
    m_ab <- do.call(rbind, m_ab)

    dat_ab <- data.frame(m_ab)
    colnames(dat_ab) <- 1:20
    dat_ab$site <- 1:20

    dat_ab <- pivot_longer(
      dat_ab,
      !site,
      names_to = "time",
      values_to = "abundance"
    )

    dat_raw <- left_join(dat_ab, dat_x, by = c("site", "time"))
  }

  if (response == "temporal-fixed") {
    a <- 1000
    b <- 10
    mu_x <- apply(m_x, 1, mean)
    site_effects <- rnorm(20, 0, 10)
    m_ab <- foreach(i = 1:length(mu_x), .combine = "cbind") %do%
      {
        a + b * mu_x[i] + site_effects
      }

    dat_ab <- data.frame(m_ab)
    colnames(dat_ab) <- 1:20
    dat_ab$site <- 1:20

    dat_ab <- pivot_longer(
      dat_ab,
      !site,
      names_to = "time",
      values_to = "abundance"
    )

    dat_raw <- left_join(dat_ab, dat_x, by = c("site", "time"))
  }

  if (response == "temporal-variable") {
    a <- 1000
    #b <- rtruncnorm(20, a = 0, mean = 10, sd = 10)
    b <- rnorm(20, 10, 10)
    mu_x <- apply(m_x, 1, mean)
    site_effects <- rnorm(20, 0, 10)
    m_ab <- foreach(i = 1:length(mu_x), .combine = "cbind") %do%
      {
        a + b * mu_x[i] + site_effects
      }

    dat_ab <- data.frame(m_ab)
    colnames(dat_ab) <- 1:20
    dat_ab$site <- 1:20

    dat_ab <- pivot_longer(
      dat_ab,
      !site,
      names_to = "time",
      values_to = "abundance"
    )

    dat_raw <- left_join(dat_ab, dat_x, by = c("site", "time"))
  }

  if (response == "spatio-temporal") {
    a <- 1000
    b <- 10
    mu_ab <- a + b * m_x

    m_ab <- apply(mu_ab, c(1, 2), function(x) rnorm(1, x, 10))

    dat_ab <- data.frame(m_ab)
    colnames(dat_ab) <- 1:20
    dat_ab$site <- 1:20

    dat_ab <- pivot_longer(
      dat_ab,
      !site,
      names_to = "time",
      values_to = "abundance"
    )

    dat_raw <- left_join(dat_ab, dat_x, by = c("site", "time"))
  }

  if (response == "spatial + temporal-variable") {
    a <- 1000
    b <- 10
    mu_sx <- apply(m_x, 2, mean)
    mu_ab <- a + b * mu_x

    b2 <- rnorm(20, 10, 10)
    mu_tx <- apply(m_x, 1, mean)
    #site_effects <- rnorm(20, 0, 10)
    m_ab <- foreach(i = 1:length(mu_tx), .combine = "cbind") %do%
      {
        mu_ab + b2 * mu_tx[i] #+ site_effects
      }

    dat_ab <- data.frame(m_ab)
    colnames(dat_ab) <- 1:20
    dat_ab$site <- 1:20

    dat_ab <- pivot_longer(
      dat_ab,
      !site,
      names_to = "time",
      values_to = "abundance"
    )

    dat_raw <- left_join(dat_ab, dat_x, by = c("site", "time"))
  }

  # decompose simulated env data
  decomp_vars <- c("env", "env_rand1", "env_rand2", "env_rand3")
  dat_dec <- dat_raw %>%
    # center without scaling
    dplyr::mutate(across(
      all_of(decomp_vars),
      ~ scale(., scale = FALSE)[, 1]
    )) %>%

    # compute the spatial component
    # (mean by pixel across all years of centered variables)
    dplyr::group_by(site) %>%
    dplyr::mutate(across(
      all_of(decomp_vars),
      ~ mean(., na.rm = TRUE),
      .names = "{.col}_spatial"
    )) %>%
    dplyr::ungroup() %>%
    # compute the temporal component
    # (mean by year across all pixels of centered variables)
    dplyr::group_by(time) %>%
    dplyr::mutate(across(
      all_of(decomp_vars),
      ~ mean(., na.rm = TRUE),
      .names = "{.col}_temporal"
    )) %>%
    dplyr::ungroup() %>%

    # compute residual for each site i and year j as centered variable value
    # i,j - spatial mean i - temporal mean j
    dplyr::mutate(across(
      all_of(decomp_vars),
      ~ . -
        get(paste0(cur_column(), "_spatial")) -
        get(paste0(cur_column(), "_temporal")),
      .names = "{.col}_residual"
    )) %>%
    dplyr::arrange(site, time)

  dat_raw <- dat_raw %>%
    # center without scaling
    dplyr::mutate(across(
      all_of(decomp_vars),
      ~ scale(., scale = FALSE)[, 1]
    ))

  list(raw = dat_raw, dec = dat_dec)
}

dat_sp1 <- simulate_ab(response = "spatial")
dat_temp1 <- simulate_ab(response = "temporal-fixed")
dat_sptemp1 <- simulate_ab(response = "spatio-temporal")
dat_tempvar1 <- simulate_ab(response = "temporal-variable")
dat_sptempvar1 <- simulate_ab(response = "spatial + temporal-variable")

dat_sp2 <- simulate_ab(response = "spatial", t_x = seq(0, 9.5, 0.5))
dat_temp2 <- simulate_ab(response = "temporal-fixed", t_x = seq(0, 9.5, 0.5))
dat_sptemp2 <- simulate_ab(response = "spatio-temporal", t_x = seq(0, 9.5, 0.5))
dat_tempvar2 <- simulate_ab(
  response = "temporal-variable",
  t_x = seq(0, 9.5, 0.5)
)
dat_sptempvar2 <- simulate_ab(
  response = "spatial + temporal-variable",
  t_x = seq(0, 9.5, 0.5)
)

plot_univ <- function(dat) {
  ggplot(data = dat) +
    geom_point(
      aes(x = env, y = abundance),
      color = "pink3",
      alpha = 0.5,
      size = 2
    ) +
    geom_line(
      aes(
        x = env,
        y = abundance,
        group = site
      ),
      linewidth = 1.25,
      color = "#333D79FF",
      alpha = 0.75,
      se = F,
      method = "lm",
      stat = "smooth"
    ) +
    labs(y = "Simulated Abundance", x = "Environmental Variable") +
    #scale_y_continuous(breaks = seq(0, 125, 25), limits = c(0, 135)) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      #axis.title.x = element_text(margin = margin(t = 10)),
      #axis.title.y = element_text(margin = margin(r = 10)),
      axis.title = element_text(size = 10)
    )
}

gm1 <- plot_univ(dat_sp1$raw)
gm2 <- plot_univ(dat_temp1$raw)
gm3 <- plot_univ(dat_sptemp1$raw)
gm4 <- plot_univ(dat_tempvar1$raw)
gm5 <- plot_univ(dat_sptempvar1$raw)

gm1 +
  gm2 +
  gm3 +
  gm4 +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = "a")

#ggsave("fig1.pdf", width = 180, height = 90, units = "mm", dpi = 600)

gm5 <- plot_univ(dat_sp2$raw)
gm6 <- plot_univ(dat_temp2$raw)
gm7 <- plot_univ(dat_sptemp2$raw)
gm8 <- plot_univ(dat_tempvar2$raw)
gm9 <- plot_univ(dat_sptempvar2$raw)

gm5 +
  gm6 +
  gm7 +
  gm8 +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = "a")

#ggsave("fig1.pdf", width = 180, height = 90, units = "mm", dpi = 600)

# Variable Importance -----------------------------------------------------

res_sim_sp1 <- ranger(
  abundance ~ .,
  data = select(
    dat_sp1$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10,
  num.trees = 2000
)

res_sim_sp2 <- ranger(
  abundance ~ .,
  data = select(
    dat_sp2$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10,
  num.trees = 2000
)

res_sim_temp1 <- ranger(
  abundance ~ .,
  data = select(
    dat_temp1$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10,
  num.trees = 2000
)

res_sim_temp2 <- ranger(
  abundance ~ .,
  data = select(
    dat_temp2$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10,
  num.trees = 2000
)

res_sim_sptemp1 <- ranger(
  abundance ~ .,
  data = select(
    dat_sptemp1$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10,
  num.trees = 2000
)

res_sim_sptemp2 <- ranger(
  abundance ~ .,
  data = select(
    dat_sptemp$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10,
  num.trees = 2000
)

res_sim_tempvar1 <- ranger(
  abundance ~ .,
  data = select(
    dat_tempvar2$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10,
  num.trees = 2000
)

res_sim_tempvar2 <- ranger(
  abundance ~ .,
  data = select(
    dat_tempvar2$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10,
  num.trees = 2000
)

res_sim_tempvar1 <- ranger(
  abundance ~ .,
  data = select(
    dat_tempvar2$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10,
  num.trees = 2000
)

res_sim_tempvar2 <- ranger(
  abundance ~ .,
  data = select(
    dat_tempvar2$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10,
  num.trees = 2000
)
plot_imp <- function(res, title) {
  var_imp <- res$variable.importance
  var_imp <- var_imp[order(var_imp, decreasing = T)][1:12]
  names_imp <- factor(
    names(var_imp)[12:1],
    levels = names(var_imp)[12:1]
  )

  ggplot() +
    geom_segment(
      aes(y = names_imp, x = 0, xend = var_imp[12:1]),
      color = "#333D79FF",
      linewidth = 1.1
    ) +
    labs(
      title = title,
      y = "Variables",
      x = "Permutation Importance"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      #axis.title.x = element_text(margin = margin(t = 10)),
      #axis.title.y = element_text(margin = margin(r = 10)),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
}

plot_imp(res_sim_sp1, title = "Spatial Response")
plot_imp(res_sim_temp1, title = "Fixed Temporal Response")
plot_imp(res_sim_sptemp1, title = "Spatio-Temporal Response")
plot_imp(res_sim_tempvar1, title = "Variable Temporal Response")

plot_imp(res_sim_sp2, title = "Spatial Response")
plot_imp(res_sim_temp2, title = "Fixed Temporal Response")
plot_imp(res_sim_sptemp2, title = "Spatio-Temporal Response")
plot_imp(res_sim_tempvar2, title = "Variable Temporal Response")


# Test set predictions ----------------------------------------------------

# Spatial
res_sim_sp1 <- ranger(
  abundance ~
    env_spatial + env_rand1_spatial + env_rand2_spatial + env_rand3_spatial,
  data = dat_sp1$dec,
  mtry = 3,
  num.trees = 2000
)

dat_sp_test1 <- simulate_ab(response = "spatial")
colnames(dat_sp_test1$raw)[-(1:3)] <- paste(
  colnames(dat_sp_test1$raw)[-(1:3)],
  "spatial",
  sep = "_"
)
pred_sim1 <- predict(res_sim_sp1, data = dat_sp_test1$raw)
dat_sp_test1$raw$pred <- pred_sim1$predictions

res_sim_sp2 <- ranger(
  abundance ~
    env_spatial + env_rand1_spatial + env_rand2_spatial + env_rand3_spatial,
  data = dat_sp2$dec,
  mtry = 3,
  num.trees = 2000
)

dat_sp_test2 <- simulate_ab(response = "spatial", t_x = seq(0, 9.5, 0.5))
colnames(dat_sp_test2$raw)[-(1:3)] <- paste(
  colnames(dat_sp_test2$raw)[-(1:3)],
  "spatial",
  sep = "_"
)
pred_sim2 <- predict(res_sim_sp2, data = dat_sp_test2$raw)
dat_sp_test2$raw$pred <- pred_sim2$predictions

# Fixed Temporal
res_sim_temp1 <- ranger(
  abundance ~
    env_temporal + env_rand1_temporal + env_rand2_temporal + env_rand3_temporal,
  data = dat_temp1$dec,
  mtry = 3,
  num.trees = 2000
)

dat_temp_test1 <- simulate_ab(response = "temporal-fixed")
colnames(dat_temp_test1$raw)[-(1:3)] <- paste(
  colnames(dat_temp_test1$raw)[-(1:3)],
  "temporal",
  sep = "_"
)
pred_sim1 <- predict(res_sim_temp1, data = dat_temp_test1$raw)
dat_temp_test1$raw$pred <- pred_sim1$predictions

res_sim_temp2 <- ranger(
  abundance ~
    env_temporal + env_rand1_temporal + env_rand2_temporal + env_rand3_temporal,
  data = dat_temp2$dec,
  mtry = 3,
  num.trees = 2000
)

dat_temp_test2 <- simulate_ab(
  response = "temporal-fixed",
  t_x = seq(0, 9.5, 0.5)
)
colnames(dat_temp_test2$raw)[-(1:3)] <- paste(
  colnames(dat_temp_test2$raw)[-(1:3)],
  "temporal",
  sep = "_"
)
pred_sim2 <- predict(res_sim_temp2, data = dat_temp_test2$raw)
dat_temp_test2$raw$pred <- pred_sim2$predictions

# Spatio-temporal
res_sim_sptemp1 <- ranger(
  abundance ~
    env_temporal +
      env_rand1_temporal +
      env_rand2_temporal +
      env_rand3_temporal +
      env_spatial +
      env_rand1_spatial +
      env_rand2_spatial +
      env_rand3_spatial,
  data = dat_sptemp1$dec,
  mtry = 3,
  num.trees = 2000
)

dat_sptemp_test1 <- simulate_ab(response = "spatio-temporal")
pred_sim1 <- predict(res_sim_sptemp1, data = dat_sptemp_test1$dec)
dat_sptemp_test1$dec$pred <- pred_sim1$predictions

res_sim_sptemp2 <- ranger(
  abundance ~
    env_temporal +
      env_rand1_temporal +
      env_rand2_temporal +
      env_rand3_temporal +
      env_spatial +
      env_rand1_spatial +
      env_rand2_spatial +
      env_rand3_spatial,
  data = dat_sptemp2$dec,
  mtry = 3,
  num.trees = 2000
)

dat_sptemp_test2 <- simulate_ab(
  response = "spatio-temporal",
  t_x = seq(0, 9.5, 0.5)
)
pred_sim2 <- predict(res_sim_sptemp2, data = dat_sptemp_test2$dec)
dat_sptemp_test2$dec$pred <- pred_sim2$predictions

# Variable temporal
res_sim_tempvar1 <- ranger(
  abundance ~
    env_temporal +
      env_rand1_temporal +
      env_rand2_temporal +
      env_rand3_temporal +
      env_spatial +
      env_rand1_spatial +
      env_rand2_spatial +
      env_rand3_spatial,
  data = dat_tempvar1$dec,
  mtry = 6,
  num.trees = 2000
)

dat_tempvar_test1 <- simulate_ab(response = "temporal-variable")
pred_sim1 <- predict(res_sim_tempvar1, data = dat_tempvar_test1$dec)
dat_tempvar_test1$dec$pred <- pred_sim1$predictions

res_sim_tempvar2 <- ranger(
  abundance ~
    env_temporal +
      env_rand1_temporal +
      env_rand2_temporal +
      env_rand3_temporal +
      env_spatial +
      env_rand1_spatial +
      env_rand2_spatial +
      env_rand3_spatial,
  data = dat_tempvar2$dec,
  mtry = 6,
  num.trees = 2000
)

dat_tempvar_test2 <- simulate_ab(
  response = "temporal-variable",
  t_x = seq(0, 9.5, 0.5)
)
pred_sim2 <- predict(res_sim_tempvar2, data = dat_tempvar_test2$dec)
dat_tempvar_test2$dec$pred <- pred_sim2$predictions

# Plots
plot_test <- function(dat) {
  ggplot(data = dat) +
    geom_point(
      aes(x = pred, y = abundance),
      color = "pink3",
      alpha = 0.5,
      size = 2
    ) +
    geom_line(
      aes(
        x = pred,
        y = abundance,
        group = site
      ),
      linewidth = 1.25,
      color = "#333D79FF",
      alpha = 0.75,
      se = F,
      method = "lm",
      stat = "smooth"
    ) +
    labs(y = "Simulated Abundance", x = "Predicted Abundance") +
    #scale_y_continuous(breaks = seq(0, 125, 25), limits = c(0, 135)) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      #axis.title.x = element_text(margin = margin(t = 10)),
      #axis.title.y = element_text(margin = margin(r = 10)),
      axis.title = element_text(size = 10)
    )
}

plot_test(dat_sp_test1$raw)
plot_test(dat_sp_test2$raw)

plot_test(dat_temp_test1$raw)
plot_test(dat_temp_test2$raw)

plot_test(dat_sptemp_test1$dec)
plot_test(dat_sptemp_test2$dec)

plot_test(dat_tempvar_test1$dec)
plot_test(dat_tempvar_test2$dec)
