#-------------------------------------------------------------------------------
# packages
library(tidyverse)
library(ggpubr)

#-------------------------------------------------------------------------------
# data
load("simulation_results.rda")
# load("./simulations/simulation_results.rda")

#-------------------------------------------------------------------------------
# ggplot themes and settings
theme_set(theme_bw(base_size = 10))

## reformat the rejection rates such that all colors are always present
## (artificial rejection rates of 0, 0.001, ..., 0.999, 1 outside of plotting
## region)
data <- data %>%
  group_by(n, K, type, method) %>%
  nest() %>%
  mutate(
    data = map(
      .x = data,
      .f = function(df) {
        rbind(
          df,
          tibble(
            bias = rep(-5, 1001),
            dispersion = rep(-5, 1001),
            rej_0.001 = seq(0, 1, length.out = 1001),
            rej_0.01 = rej_0.001,
            rej_0.05 = rej_0.01,
            rej_0.1 = rej_0.001
          )
        )
      }
    )
  ) %>%
  unnest(cols = data) %>%
  ungroup()

#-------------------------------------------------------------------------------
# illustration of simulated PIT
set.seed(12345)
n <- 10000
mu <- rnorm(n)
y <- rnorm(n, mu)

bias_disp <- expand.grid(
  bias = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5),
  dispersion = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5)
)

hist_dat <- as_tibble(bias_disp) %>%
  mutate(pit = map2(
    .x = bias,
    .y = dispersion,
    ~pnorm(y, mean = mu + .x, sd = sqrt(1 + .y))
  )) %>%
  unnest(cols = pit) %>%
  mutate(
    dispersion = factor(
      dispersion,
      levels = rev(sort(unique(dispersion))),
      labels = rev(sort(unique(dispersion))),
      ordered = TRUE
    )
  )

pit_grid_1 <- ggplot(hist_dat) +
  geom_hline(yintercept = 1, col = "red") +
  geom_histogram(aes(x = pit, y = ..density..), boundary = 0, binwidth = 0.05) +
  facet_grid(rows = vars(dispersion), cols = vars(bias), scales = "free_y") +
  theme_bw(base_size = 8) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12)
  ) +
  labs(x = element_blank(), y = element_blank()) +
  ggtitle("(a)")

pit_grid_2 <- hist_dat %>%
  group_by(bias, dispersion) %>%
  filter(seq_along(pit) < 361) %>%
  ggplot() +
  geom_hline(yintercept = 1, col = "red") +
  geom_histogram(aes(x = pit, y = ..density..), boundary = 0, binwidth = 0.05) +
  facet_grid(rows = vars(dispersion), cols = vars(bias), scales = "free_y") +
  theme_bw(base_size = 8) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12)
  ) +
  labs(x = element_blank(), y = element_blank()) +
  ggtitle("(b)")

pdf(width = 8, height = 4, file = "./tex/pit_illustration.pdf")
ggarrange(pit_grid_1, pit_grid_2, ncol = 2)
dev.off()

#-------------------------------------------------------------------------------
# power

nn <- 360
tp <- "pit"
kk <- 0
pit <- data %>%
  filter(type == tp & n == nn & K == kk) %>%
  select(method, rej_0.05, bias, dispersion) %>%
  mutate(
    rej_0.05 = cut(
      rej_0.05,
      c(0, 0.025, 0.055, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    method = fct_recode(
      method,
      "ks.test" = "ptest",
      "beta-e" = "beta",
      "kernel-e" = "kernel"
    )
  ) %>%
  ggplot() +
  geom_tile(
    aes(x = bias, y = dispersion, fill = rej_0.05),
    height = 0.09,
    width = 0.09
  ) +
  facet_grid(cols = vars(method)) +
  scale_fill_brewer(palette = "PuBuGn", direction = -1) +
  labs(y = "Dispersion error", x = element_blank()) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_path(
    data = tibble(
      x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
      y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
    ),
    aes(x = x, y = y),
    col = 2,
    lwd = 0.5
  ) +
  ggtitle("(a)")

tp <- "rank_histogram"
kk <- 20
rhist <- data %>%
  filter(type == tp & n == nn & K == kk) %>%
  select(method, rej_0.05, bias, dispersion) %>%
  mutate(
    rej_0.05 = cut(
      rej_0.05,
      c(0, 0.025, 0.055, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    method = fct_recode(
      method,
      "chisq.test" = "ptest",
      "empirical-e" = "empirical",
      "betabinom-e" = "betabinom"
    )
  ) %>%
  ggplot() +
  geom_tile(
    aes(x = bias, y = dispersion, fill = rej_0.05),
    height = 0.09,
    width = 0.09
  ) +
  facet_grid(cols = vars(method)) +
  scale_fill_brewer(palette = "PuBuGn", direction = -1) +
  labs(y = "Dispersion error", x = element_blank()) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_path(
    data = tibble(
      x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
      y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
    ),
    aes(x = x, y = y),
    col = 2,
    lwd = 0.5
  ) +
  ggtitle("(b)")

tp <- "stochtest"
kk <- 0
stochtest <- data %>%
  filter(type == tp & n == nn & K == kk) %>%
  select(method, rej_0.05, bias, dispersion) %>%
  mutate(
    rej_0.05 = cut(
      rej_0.05,
      c(0, 0.025, 0.055, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    method = fct_recode(
      method,
      "ks.test" = "ptest",
      "bernstein-e" = "bernstein",
      "grenander-e" = "grenander"
    )
  ) %>%
  ggplot() +
  geom_tile(
    aes(x = bias, y = dispersion, fill = rej_0.05),
    height = 0.09,
    width = 0.09
  ) +
  facet_grid(cols = vars(method)) +
  scale_fill_brewer(palette = "PuBuGn", direction = -1) +
  labs(y = "Dispersion error", x = element_blank()) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_path(
    data = tibble(
      x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
      y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
    ),
    aes(x = x, y = y),
    col = 2,
    lwd = 0.5
  ) +
  ggtitle("(c)")

tp <- "quantile_pit"
kk <- 19
qpit <- data %>%
  filter(type == tp & n == nn & K == kk) %>%
  select(method, rej_0.05, bias, dispersion) %>%
  mutate(
    rej_0.05 = cut(
      rej_0.05,
      c(0, 0.025, 0.055, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    method = fct_recode(
      method,
      "ks.test" = "ptest",
      "bernstein-e" = "bernstein",
      "grenander-e" = "grenander"
    )
  ) %>%
  ggplot() +
  geom_tile(
    aes(x = bias, y = dispersion, fill = rej_0.05),
    height = 0.09,
    width = 0.09
  ) +
  facet_grid(cols = vars(method)) +
  scale_fill_brewer(palette = "PuBuGn", direction = -1) +
  labs(
    x = "Bias", y = "Dispersion error", fill = "Rejection rate") +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.direction = "horizontal",
    legend.key = element_rect(color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(6, "mm"),
    legend.spacing.x = unit(2, "mm")
  ) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_path(
    data = tibble(
      x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
      y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
    ),
    aes(x = x, y = y),
    col = 2,
    lwd = 0.5
  ) +
  ggtitle("(d)")

pdf(width = 8, height = 12, file = "./tex/simulations_power.pdf")
ggarrange(pit, rhist, stochtest, qpit, ncol = 1, heights = c(1, 1, 1, 1.2))
dev.off()


#-------------------------------------------------------------------------------
# supplementary figures
specs <- list(c(180, 0.05), c(720, 0.05), c(360, 0.01), c(360, 0.1))
for (spec in specs) {
  nn <- spec[1]
  alph <- spec[2]
  
  ## n = 180, same plots as in paper
  tp <- "pit"
  kk <- 0
  pit <- data %>%
    gather(key = "alpha", value = "rej", starts_with("rej")) %>%
    mutate(alpha = parse_number(alpha)) %>%
    filter(alpha == alph) %>%
    filter(type == tp & n == nn & K == kk) %>%
    select(method, rej, bias, dispersion) %>%
    mutate(
      rej = cut(
        rej,
        c(0, alph / 2, alph * 1.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
        include.lowest = TRUE
      )
    ) %>%
    mutate(
      method = fct_recode(
        method,
        "ks.test" = "ptest",
        "beta-e" = "beta",
        "kernel-e" = "kernel"
      )
    ) %>%
    ggplot() +
    geom_tile(
      aes(x = bias, y = dispersion, fill = rej),
      height = 0.09,
      width = 0.09
    ) +
    facet_grid(cols = vars(method)) +
    scale_fill_brewer(palette = "PuBuGn", direction = -1) +
    labs(y = "Dispersion error", x = element_blank()) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
    scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
    geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
    geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
    geom_path(
      data = tibble(
        x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
        y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
      ),
      aes(x = x, y = y),
      col = 2,
      lwd = 0.5
    ) +
    ggtitle("(a)")
  
  tp <- "rank_histogram"
  kk <- 20
  rhist <- data %>%
    gather(key = "alpha", value = "rej", starts_with("rej")) %>%
    mutate(alpha = parse_number(alpha)) %>%
    filter(alpha == alph) %>%
    filter(type == tp & n == nn & K == kk) %>%
    select(method, rej, bias, dispersion) %>%
    mutate(
      rej = cut(
        rej,
        c(0, alph / 2, alph * 1.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
        include.lowest = TRUE
      )
    ) %>%
    mutate(
      method = fct_recode(
        method,
        "chisq.test" = "ptest",
        "empirical-e" = "empirical",
        "betabinom-e" = "betabinom"
      )
    ) %>%
    ggplot() +
    geom_tile(
      aes(x = bias, y = dispersion, fill = rej),
      height = 0.09,
      width = 0.09
    ) +
    facet_grid(cols = vars(method)) +
    scale_fill_brewer(palette = "PuBuGn", direction = -1) +
    labs(y = "Dispersion error", x = element_blank()) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
    scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
    geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
    geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
    geom_path(
      data = tibble(
        x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
        y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
      ),
      aes(x = x, y = y),
      col = 2,
      lwd = 0.5
    ) +
    ggtitle("(b)")
  
  tp <- "stochtest"
  kk <- 0
  stochtest <- data %>%
    gather(key = "alpha", value = "rej", starts_with("rej")) %>%
    mutate(alpha = parse_number(alpha)) %>%
    filter(alpha == alph) %>%
    filter(type == tp & n == nn & K == kk) %>%
    select(method, rej, bias, dispersion) %>%
    mutate(
      rej = cut(
        rej,
        c(0, alph / 2, alph * 1.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
        include.lowest = TRUE
      )
    ) %>%
    mutate(
      method = fct_recode(
        method,
        "ks.test" = "ptest",
        "bernstein-e" = "bernstein",
        "grenander-e" = "grenander"
      )
    ) %>%
    ggplot() +
    geom_tile(
      aes(x = bias, y = dispersion, fill = rej),
      height = 0.09,
      width = 0.09
    ) +
    facet_grid(cols = vars(method)) +
    scale_fill_brewer(palette = "PuBuGn", direction = -1) +
    labs(y = "Dispersion error", x = element_blank()) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
    scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
    geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
    geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
    geom_path(
      data = tibble(
        x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
        y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
      ),
      aes(x = x, y = y),
      col = 2,
      lwd = 0.5
    ) +
    ggtitle("(c)")
  
  tp <- "quantile_pit"
  kk <- 19
  qpit <- data %>%
    gather(key = "alpha", value = "rej", starts_with("rej")) %>%
    mutate(alpha = parse_number(alpha)) %>%
    filter(alpha == alph) %>%
    filter(type == tp & n == nn & K == kk) %>%
    select(method, rej, bias, dispersion) %>%
    mutate(
      rej = cut(
        rej,
        c(0, alph / 2, alph * 1.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
        include.lowest = TRUE
      )
    ) %>%
    mutate(
      method = fct_recode(
        method,
        "ks.test" = "ptest",
        "bernstein-e" = "bernstein",
        "grenander-e" = "grenander"
      )
    ) %>%
    ggplot() +
    geom_tile(
      aes(x = bias, y = dispersion, fill = rej),
      height = 0.09,
      width = 0.09
    ) +
    facet_grid(cols = vars(method)) +
    scale_fill_brewer(palette = "PuBuGn", direction = -1) +
    labs(
      x = "Bias", y = "Dispersion error", fill = "Rejection rate"
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.direction = "horizontal",
      legend.key = element_rect(color = "black"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.key.size = unit(6, "mm"),
      legend.spacing.x = unit(2, "mm")
    ) +
    coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
    scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
    geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
    geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
    geom_path(
      data = tibble(
        x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
        y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
      ),
      aes(x = x, y = y),
      col = 2,
      lwd = 0.5
    ) +
    ggtitle("(d)")
  
  plt <- 
    ggarrange(pit, rhist, stochtest, qpit, ncol = 1, heights = c(1, 1, 1, 1.2))
  
  pdf(
    width = 8,
    height = 12,
    file = paste0("./tex/simulations_power_n_", nn, "_alpha_", alph, ".pdf")
  )
  print(plt)
  dev.off()
}
  
## different ensemble sizes for discrete uniform distribution
nn <- 360
alph <- 0.05

tp <- "rank_histogram"
kk <- 10
rhist_10 <- data %>%
  gather(key = "alpha", value = "rej", starts_with("rej")) %>%
  mutate(alpha = parse_number(alpha)) %>%
  filter(alpha == alph) %>%
  filter(type == tp & n == nn & K == kk) %>%
  select(method, rej, bias, dispersion) %>%
  mutate(
    rej = cut(
      rej,
      c(0, 0.025, 0.055, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    method = fct_recode(
      method,
      "chisq.test" = "ptest",
      "empirical-e" = "empirical",
      "betabinom-e" = "betabinom"
    )
  ) %>%
  ggplot() +
  geom_tile(
    aes(x = bias, y = dispersion, fill = rej),
    height = 0.09,
    width = 0.09
  ) +
  facet_grid(cols = vars(method)) +
  scale_fill_brewer(palette = "PuBuGn", direction = -1) +
  labs(y = "Dispersion error", x = element_blank()) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_path(
    data = tibble(
      x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
      y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
    ),
    aes(x = x, y = y),
    col = 2,
    lwd = 0.5
  ) +
  ggtitle("(a)")

kk <- 20
rhist_20 <- data %>%
  gather(key = "alpha", value = "rej", starts_with("rej")) %>%
  mutate(alpha = parse_number(alpha)) %>%
  filter(alpha == alph) %>%
  filter(type == tp & n == nn & K == kk) %>%
  select(method, rej, bias, dispersion) %>%
  mutate(
    rej = cut(
      rej,
      c(0, 0.025, 0.055, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    method = fct_recode(
      method,
      "chisq.test" = "ptest",
      "empirical-e" = "empirical",
      "betabinom-e" = "betabinom"
    )
  ) %>%
  ggplot() +
  geom_tile(
    aes(x = bias, y = dispersion, fill = rej),
    height = 0.09,
    width = 0.09
  ) +
  facet_grid(cols = vars(method)) +
  scale_fill_brewer(palette = "PuBuGn", direction = -1) +
  labs(y = "Dispersion error", x = element_blank()) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_path(
    data = tibble(
      x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
      y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
    ),
    aes(x = x, y = y),
    col = 2,
    lwd = 0.5
  ) +
  ggtitle("(b)")

kk <- 50
rhist_50 <- data %>%
  gather(key = "alpha", value = "rej", starts_with("rej")) %>%
  mutate(alpha = parse_number(alpha)) %>%
  filter(alpha == alph) %>%
  filter(type == tp & n == nn & K == kk) %>%
  select(method, rej, bias, dispersion) %>%
  mutate(
    rej = cut(
      rej,
      c(0, 0.025, 0.055, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    method = fct_recode(
      method,
      "chisq.test" = "ptest",
      "empirical-e" = "empirical",
      "betabinom-e" = "betabinom"
    )
  ) %>%
  ggplot() +
  geom_tile(
    aes(x = bias, y = dispersion, fill = rej),
    height = 0.09,
    width = 0.09
  ) +
  facet_grid(cols = vars(method)) +
  scale_fill_brewer(palette = "PuBuGn", direction = -1) +
  labs(
    x = "Bias", y = "Dispersion error", fill = "Rejection rate"
  ) + 
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.direction = "horizontal",
    legend.key = element_rect(color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(6, "mm"),
    legend.spacing.x = unit(2, "mm")
  ) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_path(
    data = tibble(
      x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
      y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
    ),
    aes(x = x, y = y),
    col = 2,
    lwd = 0.5
  ) +
  ggtitle("(c)")

pdf(width = 8, height = 10, file = "./tex/simulations_power_rhist.pdf")
ggarrange(rhist_10, rhist_20, rhist_50, ncol = 1, heights = c(1, 1, 1.2))
dev.off()

## quantile pit
tp <- "quantile_pit"
kk <- 9
qpit_9 <- data %>%
  gather(key = "alpha", value = "rej", starts_with("rej")) %>%
  mutate(alpha = parse_number(alpha)) %>%
  filter(alpha == alph) %>%
  filter(type == tp & n == nn & K == kk) %>%
  select(method, rej, bias, dispersion) %>%
  mutate(
    rej = cut(
      rej,
      c(0, 0.025, 0.055, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    method = fct_recode(
      method,
      "ks.test" = "ptest",
      "bernstein-e" = "bernstein",
      "grenander-e" = "grenander"
    )
  ) %>%
  ggplot() +
  geom_tile(
    aes(x = bias, y = dispersion, fill = rej),
    height = 0.09,
    width = 0.09
  ) +
  facet_grid(cols = vars(method)) +
  scale_fill_brewer(palette = "PuBuGn", direction = -1) +
  labs(y = "Dispersion error", x = element_blank()) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_path(
    data = tibble(
      x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
      y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
    ),
    aes(x = x, y = y),
    col = 2,
    lwd = 0.5
  ) +
  ggtitle("(a)")

kk <- 19
qpit_19 <- data %>%
  gather(key = "alpha", value = "rej", starts_with("rej")) %>%
  mutate(alpha = parse_number(alpha)) %>%
  filter(alpha == alph) %>%
  filter(type == tp & n == nn & K == kk) %>%
  select(method, rej, bias, dispersion) %>%
  mutate(
    rej = cut(
      rej,
      c(0, 0.025, 0.055, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    method = fct_recode(
      method,
      "ks.test" = "ptest",
      "bernstein-e" = "bernstein",
      "grenander-e" = "grenander"
    )
  ) %>%
  ggplot() +
  geom_tile(
    aes(x = bias, y = dispersion, fill = rej),
    height = 0.09,
    width = 0.09
  ) +
  facet_grid(cols = vars(method)) +
  scale_fill_brewer(palette = "PuBuGn", direction = -1) +
  labs(
    x = "Bias", y = "Dispersion error", fill = "Rejection rate"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.direction = "horizontal",
    legend.key = element_rect(color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(6, "mm"),
    legend.spacing.x = unit(2, "mm")
  ) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  geom_hline(yintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_vline(xintercept = seq(-0.45, 0.45, 0.1), color = "lightgray", lwd = 0.5) +
  geom_path(
    data = tibble(
      x = c(-0.05, 0.05, 0.05, -0.05, -0.05),
      y = c(-0.05, -0.05, 0.05, 0.05, -0.05)
    ),
    aes(x = x, y = y),
    col = 2,
    lwd = 0.5
  ) +
  ggtitle("(b)")

pdf(width = 8, height = 7, file = "./tex/simulations_power_qpit.pdf")
ggarrange(qpit_9, qpit_19, ncol = 1, heights = c(1, 1.2))
dev.off()