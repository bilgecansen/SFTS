library(tidyverse)
library(ranger)
library(foreach)
library(patchwork)
library(truncnorm)
library(ggbeeswarm)

theme_set(theme_bw())


# Univariate Relationships ------------------------------------------------

simulate_ab <- function(b_mean, t_x = 0:19) {
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

  a <- 1000
  b <- 50
  mu_sx <- apply(m_x, 2, mean)
  site_effects <- rnorm(20, 0, 10)
  mu_ab <- a + b * mu_sx + site_effects
  b2 <- rnorm(20, b_mean, 1)
  m_ab <- foreach(i = 1:length(mu_ab), .combine = "rbind") %do%
    {
      mu_ab[i] + b2[i] * m_x[i, ] + rnorm(20, 0, 5)
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

dat_uneven <- list()
dat_uneven$pos <- simulate_ab(b_mean = 10, t_x = seq(0, 9.5, 0.5))
dat_uneven$neut <- simulate_ab(b_mean = 0, t_x = seq(0, 9.5, 0.5))
dat_uneven$neg <- simulate_ab(b_mean = -10, t_x = seq(0, 9.5, 0.5))

plot_univ <- function(dat) {
  ggplot(data = dat) +
    geom_line(
      aes(
        x = env,
        y = abundance,
        group = site
      ),
      linewidth = 1.25,
      color = "#90AFC5",
      alpha = 0.75,
      se = F,
      method = "lm",
      stat = "smooth"
    ) +
    geom_point(
      aes(x = env, y = abundance),
      color = "#763626",
      alpha = 0.5,
      size = 2
    ) +
    labs(y = "Simulated Abundance", x = "Environmental Variable") +
    #scale_y_continuous(breaks = seq(0, 125, 25), limits = c(0, 135)) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(size = 10)
    )
}

gm_uneven <- list()
gm_uneven$pos <- plot_univ(dat_uneven$pos$raw)
gm_uneven$neut <- plot_univ(dat_uneven$neut$raw)
gm_uneven$neg <- plot_univ(dat_uneven$neg$raw)


# Test set predictions ----------------------------------------------------

# Models with only spatial components
res_test_uneven <- foreach(i = 1:length(dat_uneven)) %do%
  {
    ranger(
      abundance ~
        env_spatial + env_rand1_spatial + env_rand2_spatial + env_rand3_spatial,
      data = dat_uneven[[i]]$dec,
      mtry = 3,
      num.trees = 2000
    )
  }
names(res_test_uneven) <- c("pos", "neut", "neg")

# Simulate new test data
b_mean <- c(20, 0, -20)
dat_test_uneven <- foreach(i = 1:length(res_test_uneven)) %do%
  {
    z <- simulate_ab(b_mean = b_mean[i], t_x = seq(0, 9.5, 0.5))
    colnames(z$raw)[-(1:3)] <- paste(
      colnames(z$raw)[-(1:3)],
      "spatial",
      sep = "_"
    )

    z
  }
names(dat_test_uneven) <- c("pos", "neut", "neg")

# predict test data
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
    geom_line(
      aes(
        x = pred,
        y = abundance,
        group = site
      ),
      linewidth = 1.25,
      color = "#90AFC5",
      alpha = 0.75,
      se = F,
      method = "lm",
      stat = "smooth"
    ) +
    geom_point(
      data = group_by(dat, site) %>%
        summarise(pred = mean(pred), abundance = mean(abundance)),
      aes(x = pred, y = abundance),
      color = "#763626",
      alpha = 0.5,
      size = 2
    ) +
    labs(y = "Simulated Abundance", x = "Predicted Abundance") +
    #scale_y_continuous(breaks = seq(0, 125, 25), limits = c(0, 135)) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(size = 10)
    )
}

g_pred_uneven <- foreach(i = 1:length(dat_test_uneven)) %do%
  {
    plot_test(dat_test_uneven[[i]]$raw)
  }
names(g_pred_uneven) <- c("pos", "neut", "neg")

plots_sim <- list(
  gm_uneven = gm_uneven,
  g_pred_uneven = g_pred_uneven
)

plots_sim

saveRDS(plots_sim, "figures/plots_sim.rds")
