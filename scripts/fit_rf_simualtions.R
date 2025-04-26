library(tidyverse)
library(ranger)
library(foreach)
library(patchwork)
library(truncnorm)


# Univariate Relationships ------------------------------------------------

simulate_ab <- function(response) {
  x <- 1:20
  m_x <- foreach(i = 0:19, .combine = "rbind") %do%
    {
      x + i
    }
  m_x <- apply(m_x, c(1, 2), function(x) x + rnorm(1, 0, 1))

  dat_x <- data.frame(m_x)
  colnames(dat_x) <- 1:20
  dat_x$site <- 1:20

  dat_x <- pivot_longer(dat_x, !site, names_to = "time", values_to = "env")
  dat_x$env_rand1 <- rnorm(nrow(dat_x), 0, 10)
  dat_x$env_rand2 <- rnorm(nrow(dat_x), 0, 10)
  dat_x$env_rand3 <- rnorm(nrow(dat_x), 0, 10)

  if (response == "spatial") {
    a <- 10
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
    a <- 10
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
    a <- 10
    b <- rtruncnorm(20, a = 0, mean = 10, sd = 10)
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
    a <- 10
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

  list(raw = dat_raw, dec = dat_dec)
}

dat_sp <- simulate_ab(response = "spatial")
dat_temp <- simulate_ab(response = "temporal-fixed")
dat_sptemp <- simulate_ab(response = "spatio-temporal")
dat_tempvar <- simulate_ab(response = "temporal-variable")

gm1 <- ggplot(data = dat_sp$raw) +
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

gm2 <- ggplot(data = dat_temp$raw) +
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

gm3 <- ggplot(data = dat_sptemp$raw) +
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

gm4 <- ggplot(data = dat_tempvar$raw) +
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

gm1 +
  gm2 +
  gm3 +
  gm4 +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = "a")

#ggsave("fig1.pdf", width = 180, height = 90, units = "mm", dpi = 600)

# Variable Importance -----------------------------------------------------

res_sim_sp <- ranger(
  abundance ~ .,
  data = select(
    dat_sp$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10
)

res_sim_temp <- ranger(
  abundance ~ .,
  data = select(
    dat_temp$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10
)

res_sim_sptemp <- ranger(
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
  mtry = 10
)

res_sim_tempvar <- ranger(
  abundance ~ .,
  data = select(
    dat_tempvar$dec,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10
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

plot_imp(res_sim_sp, title = "Spatial Response")
plot_imp(res_sim_temp, title = "Fixed Temporal Response")
plot_imp(res_sim_sptemp, title = "Spatio-Temporal Response")
plot_imp(res_sim_tempvar, title = "Variable Temporal Response")


# Test set predictions ----------------------------------------------------

res_sim_sp <- ranger(
  abundance ~
    env_spatial + env_rand1_spatial + env_rand2_spatial + env_rand3_spatial,
  data = m_train,
  mtry = 3,
  importance = "permutation"
)
pred_sim <- predict(res_sim_sp, data = m_test)

m_test$pred <- pred_sim$predictions

gm2 <- ggplot(data = m_test) +
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


# Temporal effects --------------------------------------------------------

a <- 10 + rnorm(20, 0, 10)
b <- 10
mu_x <- 1:20
m_ab <- foreach(i = 1:length(a), .combine = "rbind") %do%
  {
    a[i] + b * mu_x
  }
m_ab <- data.frame(m_ab)
colnames(m_ab) <- 1:20
m_ab$site <- 1:20

m_x <- matrix(rep(1:20, each = 20), ncol = 20, nrow = 20)
m_x <- t(apply(m_x, 1, function(x) x + rnorm(1, 0, 2)))
m_x <- data.frame(m_x)
colnames(m_x) <- 1:20
m_x$site <- 1:20

m_ab <- pivot_longer(m_ab, !site, names_to = "time", values_to = "abundance")

m_x <- pivot_longer(m_x, !site, names_to = "time", values_to = "env")
m_x$env_rand1 <- rnorm(nrow(m_x), 0, 10)
m_x$env_rand2 <- rnorm(nrow(m_x), 0, 10)
m_x$env_rand3 <- rnorm(nrow(m_x), 0, 10)

m <- left_join(m_ab, m_x, by = c("site", "time"))

gm3 <- ggplot(data = m) +
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

# Create a test set
m_ab_test <- lapply(mu_ab, function(x) rnorm(20, x, 10))
m_ab_test <- do.call(rbind, m_ab_test)
m_ab_test <- data.frame(m_ab_test)
colnames(m_ab_test) <- 1:20
m_ab_test$site <- 1:10

m_x_test <- lapply(mu_x, function(x) rnorm(20, x, 2))
m_x_test <- do.call(rbind, m_x_test)
m_x_test <- data.frame(m_x_test)
colnames(m_x_test) <- 1:20
m_x_test$site <- 1:10

m_ab_test <- pivot_longer(
  m_ab_test,
  !site,
  names_to = "time",
  values_to = "abundance"
)

m_x_test <- pivot_longer(
  m_x_test,
  !site,
  names_to = "time",
  values_to = "env_spatial"
)
m_x_test$env_rand1_spatial <- rnorm(nrow(m_x_test), 0, 10)
m_x_test$env_rand2_spatial <- rnorm(nrow(m_x_test), 0, 10)
m_x_test$env_rand3_spatial <- rnorm(nrow(m_x_test), 0, 10)

m_test <- left_join(m_ab_test, m_x_test, by = c("site", "time"))
m_test <- m_test %>%
  dplyr::mutate(across(
    all_of(
      c(
        "env_spatial",
        "env_rand1_spatial",
        "env_rand2_spatial",
        "env_rand3_spatial"
      )
    ),
    ~ scale(., scale = FALSE)[, 1]
  ))

# decompose simulated env data
decomp_vars <- c("env", "env_rand1", "env_rand2", "env_rand3")
m_train <- m %>%
  # center without scaling
  dplyr::mutate(across(all_of(decomp_vars), ~ scale(., scale = FALSE)[, 1])) %>%

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

res_sim <- ranger(
  abundance ~ .,
  data = select(
    m_train,
    -site,
    -time,
    -env,
    -env_rand1,
    -env_rand2,
    -env_rand3
  ),
  importance = "permutation",
  mtry = 10
)

res_sim_sp <- ranger(
  abundance ~
    env_spatial + env_rand1_spatial + env_rand2_spatial + env_rand3_spatial,
  data = m_train,
  mtry = 3,
  importance = "permutation"
)
pred_sim <- predict(res_sim_sp, data = m_test)

m_test$pred <- pred_sim$predictions
gm1 + gm2 + plot_layout(axes = "collect") + plot_annotation(tag_levels = "a")

ggsave("fig1.pdf", width = 180, height = 90, units = "mm", dpi = 600)
