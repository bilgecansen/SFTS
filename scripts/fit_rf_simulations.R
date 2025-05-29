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

  if (response == "spatial + temporal-opposite") {
    a <- 1000
    b <- 10
    mu_sx <- apply(m_x, 2, mean)
    mu_ab <- a + b * mu_sx

    b2 <- -2
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
  "spatial + temporal-opposite",
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
      color = "#763626",
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
      color = "#90AFC5",
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
names(gm_equal) <- scenarios

gm_uneven <- foreach(i = 1:length(dat_equal)) %do%
  {
    plot_univ(dat_uneven[[i]]$raw)
  }
names(gm_uneven) <- scenarios


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
      color = "#90AFC5",
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
  "Simulation: Spatial(+), Temporal(0)",
  "Simulation: Spatial(0), Temporal(+)",
  "Simulation: Spatial(+), Temporal(+)",
  "Simulation: Spatial(0), Temporal(random)",
  "Simulation: Spatial(+), Temporal(-)",
  "Simulation: Quadratic Response"
)

g_imp_equal <- foreach(i = 1:length(res_equal)) %do%
  {
    plot_imp(res_equal[[i]], title = titles[i])
  }
names(g_imp_equal) <- scenarios

g_imp_uneven <- foreach(i = 1:length(res_uneven)) %do%
  {
    plot_imp(res_uneven[[i]], title = titles[i])
  }
names(g_imp_uneven) <- scenarios


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
plot_test <- function(dat, title) {
  ggplot(data = dat) +
    geom_point(
      aes(x = pred, y = abundance),
      color = "#763626",
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
      color = "#90AFC5",
      alpha = 0.75,
      se = F,
      method = "lm",
      stat = "smooth"
    ) +
    labs(y = "Simulated Abundance", x = "Predicted Abundance", title = title) +
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
    plot_test(dat_test_equal[[i]]$raw, titles[i])
  }
names(g_pred_equal) <- scenarios

g_pred_uneven <- foreach(i = 1:length(dat_test_uneven)) %do%
  {
    plot_test(dat_test_uneven[[i]]$raw, titles[i])
  }
names(g_pred_uneven) <- scenarios


# Simulations with multiple iterations ------------------------------------

# generate data
dat_equal_multi <- foreach(i = 1:length(scenarios)) %:%
  foreach(h = 1:100) %do%
  {
    simulate_ab(response = scenarios[i])
  }
names(dat_equal_multi) <- scenarios

dat_uneven_multi <- foreach(i = 1:length(scenarios)) %:%
  foreach(h = 1:100) %do%
  {
    simulate_ab(response = scenarios[i], t_x = seq(0, 9.5, 0.5))
  }
names(dat_equal_multi) <- scenarios

# run RF
varimp_equal <- foreach(i = 1:length(dat_equal_multi)) %:%
  foreach(h = 1:100) %do%
  {
    res_rf <- ranger(
      abundance ~ .,
      data = select(
        dat_equal_multi[[i]][[h]]$dec,
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

    res_rf$variable.importance
  }
names(varimp_equal) <- scenarios

varimp_uneven <- foreach(i = 1:length(dat_uneven_multi)) %:%
  foreach(h = 1:100) %do%
  {
    res_rf <- ranger(
      abundance ~ .,
      data = select(
        dat_uneven_multi[[i]][[h]]$dec,
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

    res_rf$variable.importance
  }
names(varimp_uneven) <- scenarios

med_imp_equal <- foreach(i = 1:length(varimp_equal)) %:%
  foreach(h = 1:100, .combine = "rbind") %do%
  {
    z <- varimp_equal[[i]][[h]][order(varimp_equal[[i]][[h]], decreasing = T)]

    data.frame(
      best_imp = c(
        min(str_which(names(z), "spatial")),
        min(str_which(names(z), "temporal")),
        min(str_which(names(z), "residual"))
      ),
      type = c("Spatial", "Temporal", "Residual")
    )
  }

med_imp_uneven <- foreach(i = 1:length(varimp_uneven)) %:%
  foreach(h = 1:100, .combine = "rbind") %do%
  {
    z <- varimp_uneven[[i]][[h]][order(varimp_uneven[[i]][[h]], decreasing = T)]

    data.frame(
      best_imp = c(
        min(str_which(names(z), "spatial")),
        min(str_which(names(z), "temporal")),
        min(str_which(names(z), "residual"))
      ),
      type = c("Spatial", "Temporal", "Residual")
    )
  }

g_imp_equal_multi <- foreach(i = 1:length(scenarios)) %do%
  {
    ggplot() +
      geom_jitter(
        data = med_imp_equal[[i]],
        aes(x = type, y = best_imp),
        color = "#90AFC5",
        alpha = 0.5
      ) +
      labs(
        y = "Rank of Top Variable",
        x = "Variable Composition Type",
        title = titles[i]
      ) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10)),
        axis.title = element_text(size = 10)
      ) +
      scale_y_continuous(breaks = seq(1, 15, 2))
  }
names(g_imp_equal_multi) <- scenarios

g_imp_uneven_multi <- foreach(i = 1:length(scenarios)) %do%
  {
    ggplot() +
      geom_jitter(
        data = med_imp_uneven[[i]],
        aes(x = type, y = best_imp),
        color = "#90AFC5",
        alpha = 0.5
      ) +
      labs(
        y = "Rank of Top Variable",
        x = "Variable Composition Type",
        title = titles[i]
      ) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10)),
        axis.title = element_text(size = 10)
      ) +
      scale_y_continuous(breaks = seq(1, 15, 2))
  }
names(g_imp_uneven_multi) <- scenarios

# Test set predictions
pred_equal_multi <- foreach(i = 1:length(dat_equal)) %:%
  foreach(h = 1:100, .combine = "rbind") %do%
  {
    res <- ranger(
      abundance ~
        env_spatial + env_rand1_spatial + env_rand2_spatial + env_rand3_spatial,
      data = dat_equal_multi[[i]][[h]]$dec,
      mtry = 3,
      num.trees = 2000
    )

    z <- simulate_ab(response = scenarios[i])
    colnames(z$raw)[-(1:3)] <- paste(
      colnames(z$raw)[-(1:3)],
      "spatial",
      sep = "_"
    )

    z_pred <- predict(res, data = z$raw)
    z$raw$pred <- z_pred$predictions
    z$raw$iter <- h

    z$raw
  }
names(pred_equal_multi) <- scenarios

pred_uneven_multi <- foreach(i = 1:length(dat_uneven)) %:%
  foreach(h = 1:100, .combine = "rbind") %do%
  {
    res <- ranger(
      abundance ~
        env_spatial + env_rand1_spatial + env_rand2_spatial + env_rand3_spatial,
      data = dat_uneven_multi[[i]][[h]]$dec,
      mtry = 3,
      num.trees = 2000
    )

    z <- simulate_ab(response = scenarios[i])
    colnames(z$raw)[-(1:3)] <- paste(
      colnames(z$raw)[-(1:3)],
      "spatial",
      sep = "_"
    )

    z_pred <- predict(res, data = z$raw)
    z$raw$pred <- z_pred$predictions
    z$raw$iter <- h

    z$raw
  }
names(pred_uneven_multi) <- scenarios

plot_pred <- function(dat) {
  cor_temp <- dat %>%
    group_by(iter, site) %>%
    summarise(r = cor(abundance, pred)) %>%
    ungroup() %>%
    group_by(iter) %>%
    summarise(r = median(r))

  cor_all <- dat %>%
    group_by(iter) %>%
    summarise(r = cor(abundance, pred))

  dat_cor <- data.frame(
    r = c(cor_all$r, cor_temp$r),
    type = rep(
      c("Spatio-temporal", "Temporal"),
      each = length(cor_all$r)
    )
  )

  ggplot(dat_cor) +
    geom_jitter(
      aes(x = type, y = r),
      color = "#763626",
      alpha = 0.5
    ) +
    labs(
      y = "Prediction Correlation",
      x = "Prediction Type",
      title = titles[i]
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      #axis.title.x = element_text(margin = margin(t = 10)),
      #axis.title.y = element_text(margin = margin(r = 10)),
      axis.title = element_text(size = 10)
    ) #+
  #scale_y_continuous(breaks = seq(0, 1, 0.1))

  #ggplot() +
  #geom_point(
  #aes(x = cor_all$r, y = cor_temp$r),
  #color = "#763626",
  #alpha = 0.75,
  #size = 2
  #) +
  #labs(y = "Median Temporal Correlation", x = "Spatiotemporal Correlation") +
  #theme(
  #panel.grid.minor = element_blank(),
  #panel.border = element_blank(),
  #axis.title.x = element_text(margin = margin(t = 10)),
  #axis.title.y = element_text(margin = margin(r = 10)),
  #axis.title = element_text(size = 10)
  #)
}

g_pred_equal_multi <- foreach(i = 1:length(pred_equal_multi)) %do%
  {
    plot_pred(pred_equal_multi[[i]])
  }
names(g_pred_equal_multi) <- scenarios

g_pred_uneven_multi <- foreach(i = 1:length(pred_uneven_multi)) %do%
  {
    plot_pred(pred_uneven_multi[[i]])
  }
names(g_pred_uneven_multi) <- scenarios

plots_sim <- list(
  gm_equal = gm_equal,
  gm_uneven = gm_uneven,
  g_imp_equal = g_imp_equal,
  g_imp_equal_multi = g_imp_equal_multi,
  g_imp_uneven = g_imp_uneven,
  g_imp_uneven_multi = g_imp_uneven_multi,
  g_pred_equal = g_pred_equal,
  g_pred_equal_multi = g_pred_equal_multi,
  g_pred_uneven = g_pred_uneven,
  g_pred_uneven_multi = g_pred_uneven_multi
)

saveRDS(plots_sim, "plots_sim.rds")


# Plots -------------------------------------------------------------------

(g_imp_uneven$`spatio-temporal` + g_pred_uneven$`spatio-temporal`) /
  (g_imp_uneven_multi$`spatio-temporal` + g_pred_uneven_multi$`spatio-temporal`)
ggsave("fig1.pdf", width = 180, height = 180, units = "mm", dpi = 600)


(g_imp_uneven$spatial + g_pred_uneven$spatial) /
  (g_imp_uneven_multi$spatial + g_pred_uneven_multi$spatial)
ggsave("fig2.pdf", width = 180, height = 180, units = "mm", dpi = 600)

(g_imp_uneven$`spatial + temporal-opposite` +
  g_pred_uneven$`spatial + temporal-opposite`) /
  (g_imp_uneven_multi$`spatial + temporal-opposite` +
    g_pred_uneven_multi$`spatial + temporal-opposite`)
ggsave("fig3.pdf", width = 180, height = 180, units = "mm", dpi = 600)


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
