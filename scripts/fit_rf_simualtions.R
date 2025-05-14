library(tidyverse)
library(ranger)
library(foreach)
library(patchwork)
library(truncnorm)

theme_set(theme_bw())


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
    mu_ab <- a + b * mu_sx

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

  if (response == "quadratic") {
    a <- 1000
    b <- 40
    c <- -1
    mu_ab <- a + b * m_x + c * m_x^2

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

  dat_raw <- dat_raw %>%
    # center without scaling
    dplyr::mutate(across(
      all_of(decomp_vars),
      ~ scale(., scale = FALSE)[, 1]
    ))

  list(raw = dat_raw, dec = dat_dec)
}

scenarios <- c(
  "spatial",
  "temporal-fixed",
  "spatio-temporal",
  "temporal-variable",
  "spatial + temporal-variable",
  "quadratic"
)

dat_equal <- foreach(i = 1:length(scenarios)) %do%
  {
    simulate_ab(response = scenarios[i])
  }
names(dat_equal) <- scenarios

dat_uneven <- foreach(i = 1:length(scenarios)) %do%
  {
    simulate_ab(response = scenarios[i], t_x = seq(0, 9.5, 0.5))
  }
names(dat_uneven) <- scenarios

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

gm_equal <- foreach(i = 1:length(dat_equal)) %do%
  {
    plot_univ(dat_equal[[i]]$raw)
  }

gm_uneven <- foreach(i = 1:length(dat_equal)) %do%
  {
    plot_univ(dat_uneven[[i]]$raw)
  }


# Variable Importance -----------------------------------------------------

res_equal <- foreach(i = 1:length(dat_equal)) %do%
  {
    ranger(
      abundance ~ .,
      data = select(
        dat_equal[[i]]$dec,
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
  }
names(res_equal) <- scenarios

res_uneven <- foreach(i = 1:length(dat_uneven)) %do%
  {
    ranger(
      abundance ~ .,
      data = select(
        dat_uneven[[i]]$dec,
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
  }
names(res_uneven) <- scenarios

# Plots
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

titles <- c(
  "Spatial Response",
  "Fixed Temporal Response",
  "Spatio-Temporal Response",
  "Variable Temporal Response",
  "Spatial + Variable Temporal Response",
  "Quadratic Response"
)

g_imp_equal <- foreach(i = 1:length(res_equal)) %do%
  {
    plot_imp(res_equal[[i]], title = titles[i])
  }

g_imp_uneven <- foreach(i = 1:length(res_uneven)) %do%
  {
    plot_imp(res_uneven[[i]], title = titles[i])
  }


# Test set predictions ----------------------------------------------------

# Models with only spatial components
res_test_equal <- foreach(i = 1:length(dat_equal)) %do%
  {
    ranger(
      abundance ~
        env_spatial + env_rand1_spatial + env_rand2_spatial + env_rand3_spatial,
      data = dat_equal[[i]]$dec,
      mtry = 3,
      num.trees = 2000
    )
  }

res_test_uneven <- foreach(i = 1:length(dat_equal)) %do%
  {
    ranger(
      abundance ~
        env_spatial + env_rand1_spatial + env_rand2_spatial + env_rand3_spatial,
      data = dat_uneven[[i]]$dec,
      mtry = 3,
      num.trees = 2000
    )
  }

# Simulate new test data
dat_test_equal <- foreach(i = 1:length(scenarios)) %do%
  {
    z <- simulate_ab(response = scenarios[i])
    colnames(z$raw)[-(1:3)] <- paste(
      colnames(z$raw)[-(1:3)],
      "spatial",
      sep = "_"
    )

    z
  }
names(dat_test_equal) <- scenarios

dat_test_uneven <- foreach(i = 1:length(scenarios)) %do%
  {
    z <- simulate_ab(response = scenarios[i], t_x = seq(0, 9.5, 0.5))
    colnames(z$raw)[-(1:3)] <- paste(
      colnames(z$raw)[-(1:3)],
      "spatial",
      sep = "_"
    )

    z
  }
names(dat_test_uneven) <- scenarios

# predict test data
pred_equal <- foreach(i = 1:length(res_test_equal)) %do%
  {
    predict(res_test_equal[[i]], data = dat_test_equal[[i]]$raw)
  }
for (i in 1:length(dat_test_equal)) {
  dat_test_equal[[i]]$raw$pred <- pred_equal[[i]]$predictions
}

pred_uneven <- foreach(i = 1:length(res_test_uneven)) %do%
  {
    predict(res_test_uneven[[i]], data = dat_test_uneven[[i]]$raw)
  }
for (i in 1:length(dat_test_uneven)) {
  dat_test_uneven[[i]]$raw$pred <- pred_uneven[[i]]$predictions
}

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

g_pred_equal <- foreach(i = 1:length(dat_test_equal)) %do%
  {
    plot_test(dat_test_equal[[i]]$raw)
  }
names(g_pred_equal) <- scenarios

g_pred_uneven <- foreach(i = 1:length(dat_test_uneven)) %do%
  {
    plot_test(dat_test_uneven[[i]]$raw)
  }
names(g_pred_uneven) <- scenarios


# Plots -------------------------------------------------------------------

# Simulated environment
## Equal trends in space and time
s_x <- 1:20
m_x <- foreach(i = 0:19, .combine = "rbind") %do%
  {
    s_x + i
  }
m_x2 <- t(apply(m_x, c(1, 2), function(x) x + rnorm(1, 0, 1)))

m_x <- as.data.frame(m_x)
colnames(m_x) <- 1:20
m_x$site <- 1:20

m_x <- m_x %>%
  pivot_longer(-site, names_to = "time", values_to = "x")
m_x$time <- as.numeric(m_x$time)
m_x$x <- as.numeric(m_x$x)

ggplot(data = m_x) +
  geom_tile(
    aes(y = site, x = time, fill = x),
    #hjust = 0,
    #vjust = 0,
    color = "black",
    linewidth = 0.5
  ) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_viridis_c(limits = c(-1, 40))
ggsave("fig11.pdf", width = 180, height = 180, units = "mm", dpi = 600)

m_x2 <- as.data.frame(m_x2)
colnames(m_x2) <- 1:20
m_x2$site <- 1:20

m_x2 <- m_x2 %>%
  pivot_longer(-site, names_to = "time", values_to = "x")
m_x2$time <- as.numeric(m_x2$time)
m_x2$x <- as.numeric(m_x2$x)

ggplot(data = m_x2) +
  geom_tile(
    aes(y = site, x = time, fill = x),
    #hjust = 0,
    #vjust = 0,
    color = "black",
    linewidth = 0.5
  ) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_viridis_c(limits = c(-1, 40))

ggsave("fig12.pdf", width = 180, height = 180, units = "mm", dpi = 600)

## Weaker trends across time
t_x <- seq(0, 9.5, 0.5)
m_x3 <- foreach(i = 1:20, .combine = "rbind") %do%
  {
    s_x + t_x[i]
  }
m_x3 <- t(m_x3)
m_x4 <- apply(m_x3, c(1, 2), function(x) x + rnorm(1, 0, 1))

m_x3 <- as.data.frame(m_x3)
colnames(m_x3) <- 1:20
m_x3$site <- 1:20

m_x3 <- m_x3 %>%
  pivot_longer(-site, names_to = "time", values_to = "x")
m_x3$time <- as.numeric(m_x3$time)
m_x3$x <- as.numeric(m_x3$x)

ggplot(data = m_x3) +
  geom_tile(
    aes(y = site, x = time, fill = x),
    #hjust = 0,
    #vjust = 0,
    color = "black",
    linewidth = 0.5
  ) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_viridis_c(limits = c(-1, 40))

ggsave("fig13.pdf", width = 180, height = 180, units = "mm", dpi = 600)

m_x4 <- as.data.frame(m_x4)
colnames(m_x4) <- 1:20
m_x4$site <- 1:20

m_x4 <- m_x4 %>%
  pivot_longer(-site, names_to = "time", values_to = "x")
m_x4$time <- as.numeric(m_x4$time)
m_x4$x <- as.numeric(m_x4$x)

ggplot(data = m_x4) +
  geom_tile(
    aes(y = site, x = time, fill = x),
    #hjust = 0,
    #vjust = 0,
    color = "black",
    linewidth = 0.5
  ) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_viridis_c(limits = c(-1, 40))

ggsave("fig14.pdf", width = 180, height = 180, units = "mm", dpi = 600)
