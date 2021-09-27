#-------------------------------------------------------------------------------
# packages
library(tidyverse)
library(ggpubr)
library(feather)
library(epit)
library(lubridate)
library(bde)
library(KernSmooth)
library(knitr)

#-------------------------------------------------------------------------------
# data (do not run)
precip <- read_feather("./weather/prec_PIT.feather")
temp <- read_feather("./weather/t2m_PIT.feather")
wind <- read_feather("./weather/ws_PIT.feather")
statinfo <- read_feather("./weather/station_info.feather")

#-------------------------------------------------------------------------------
# ggplot themes
theme_set(theme_bw(base_size = 10))

#-------------------------------------------------------------------------------
# station information
statinfo <- statinfo %>%
  select(-height) %>%
  rename(
    Latitude = lat,
    Longitude = lon,
    `WMO ID` = ID,
    Name = name
  ) %>%
  arrange(`WMO ID`)
statinfo_wide <- cbind(statinfo[1:11, ], statinfo[12:22, ])
kable(statinfo_wide, booktabs = TRUE, format = "latex")

#-------------------------------------------------------------------------------
# format data and pre-select relevant variables (do not run)
data <- bind_rows(
  mutate(precip, variable = "Precipitation"),
  mutate(temp, variable = "Temperature"),
  mutate(wind, variable = "Wind speed")
)

data <- data %>%
  select(date, variable, leadtime, station, obs, rank, pit)

data <- data %>%
  filter(leadtime %in% c(24, 48, 72))

data <- data %>%
  filter(!is.na(pit) & !is.na(rank))

date_ranges <- data %>%
  select(-rank, -obs, -pit) %>%
  group_by(station, leadtime, variable) %>%
  summarise(min_date = min(date), max_date = max(date)) %>%
  mutate(date = map2(
    .x = min_date,
    .y = max_date,
    .f = ~seq(.x, .y, as.difftime(1, units = "days"))
  )) %>%
  select(-min_date, -max_date) %>%
  unnest(cols = date)

data <- full_join(
  data,
  date_ranges,
  by = c("station", "leadtime", "variable", "date")
)

data <- data %>%
  arrange(station, variable, leadtime, date)

# check
data %>%
  group_by(variable, leadtime, station) %>%
  summarise(maxdiff = max(diff(date)), mindiff = min(diff(date))) %>%
  summarise(
    minmin = min(mindiff),
    maxmin = max(mindiff),
    minmax = min(maxdiff),
    maxmax = max(maxdiff)
  )

# check date range
data %>%
  group_by(variable, leadtime, station) %>%
  summarise(min_date = min(date), max_date = max(date)) %>%
  print(n = 200)

data %>%
  group_by(station, leadtime, variable) %>%
  summarise(n = n(), na = sum(is.na(pit)))

# export selection
save(list = "data", file = "./weather/processed_weather_data.rda")

#-------------------------------------------------------------------------------
# analysis
load("weather_results.rda")
# load("./weather/weather_results.rda")

## e-value and deviation of PIT from uniform density
mean_dist_dat <- weather_results %>%
  group_by(variable, leadtime, station) %>%
  summarise(
    abs_diff = 0.05 * sum(abs(
      prop.table(table(cut(pit, seq(0, 1, 0.05), include.lowest = TRUE))) /
        0.05 - 1
    )),
    e = tail(vec_kernel_e, n = 1)
  )

cor_dat <- mean_dist_dat %>%
  group_by(variable, leadtime) %>%
  summarise(
    rho = cor(abs_diff, e, method = "spearman"),
    e = max(e)
  ) %>%
  group_by(variable) %>%
  mutate(ypos = sqrt(max(e))) %>%
  mutate(rho = paste0("rho == ", round(rho, 2)))

dist_e <- ggplot() +
  geom_hline(yintercept = c(1, 100), lty = 3) +
  geom_point(data = mean_dist_dat, aes(x = abs_diff, y = e), cex = 0.5) +
  geom_text(
    data = cor_dat,
    aes(label = rho, x = 0.05, y = ypos),
    hjust = 0,
    vjust = 0,
    parse = TRUE
  ) +
  facet_grid(rows = vars(variable), cols = vars(leadtime), scales = "free_y") +
  scale_y_log10() +
  labs(
    x = expression(L[1]~distance~of~PIT~histogram~from~uniform~density),
    y = "E-value"
  )

pdf(width = 8, height = 4, file = "./tex/distance_e_value.pdf")
dist_e
dev.off()

## station plots
stats <- c(10015, 10162, 10729)
vs <- c("Temperature", "Wind speed", "Precipitation")
lags <- c(24, 24, 48)

### station 10015
i <- 1
pithist <- weather_results %>%
  filter(station == stats[i]) %>%
  ggplot() +
  geom_hline(yintercept = 1, col = "red") +
  geom_histogram(
    mapping = aes(x = pit, y = ..density..),
    boundary = 0,
    binwidth = 0.05
  ) +
  facet_grid(rows = vars(variable), cols = vars(leadtime)) +
  labs(x = "PIT", y = "Density") +
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) +
  ggtitle(paste0("(a) Station: ", stats[i]))

e_process <- weather_results %>%
  filter(station == stats[i]) %>%
  select(date, variable, leadtime, vec_kernel_e) %>%
  ggplot() +
  geom_hline(yintercept = c(1, 100), col = 1, lty = 3) +
  geom_line(aes(x = date, y = vec_kernel_e)) +
  facet_grid(rows = vars(variable), cols = vars(leadtime), scales = "free_y") +
  scale_y_log10() +
  labs(x = "Date", y = "E-value", color = element_blank()) +
  ggtitle("(b)")

kdens_dat <- weather_results %>%
  filter(station == stats[i] & variable == vs[i] & leadtime == lags[i]) %>%
  mutate(year = year(date) + ifelse(month(date) < 7, 0, 0.5))

grd <- seq(0.001, 1 - 0.001, 0.001)
yrs <- unique(kdens_dat$year)
bws <- sapply(
  yrs,
  function(x) {
    z <- kdens_dat$pit[kdens_dat$year == x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    KernSmooth::dpik(
      x = z,
      scalest = "min",
      level = 2L,
      kernel = "epanech",
      gridsize = 401L,
      range.x = c(0, 1)
    )
  }
)
dens <- mapply(
  function(x, bws) {
    z <- kdens_dat$pit[kdens_dat$year == x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    kernel_estim <- bde::jonesCorrectionMuller94BoundaryKernel(
      dataPoints = z,
      b = bws,
      mu = 1,
      dataPointsCache = grd
    )
    kernel_estim@densityCache / bde::distribution(kernel_estim, x = 1, FALSE)
  },
  x = yrs,
  bws = bws,
  SIMPLIFY = FALSE
)

dens_dat_single_year <- tibble(
  year = yrs,
  density = dens,
  grid = list(grd)
)

dens_dat_single_year <- dens_dat_single_year %>%
  unnest(cols = c(density, grid)) %>%
  mutate(type = "only given time period")

bws <- sapply(
  yrs,
  function(x) {
    z <- kdens_dat$pit[kdens_dat$year <= x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    KernSmooth::dpik(
      x = z,
      scalest = "min",
      level = 2L,
      kernel = "epanech",
      gridsize = 401L,
      range.x = c(0, 1)
    )
  }
)
dens <- mapply(
  function(x, bws) {
    z <- kdens_dat$pit[kdens_dat$year <= x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    kernel_estim <- bde::jonesCorrectionMuller94BoundaryKernel(
      dataPoints = z,
      b = bws,
      mu = 1,
      dataPointsCache = grd
    )
    kernel_estim@densityCache / bde::distribution(kernel_estim, x = 1, FALSE)
  },
  x = yrs,
  bws = bws,
  SIMPLIFY = FALSE
)

dens_dat_accum <- tibble(
  year = yrs,
  density = dens,
  grid = list(grd)
)

dens_dat <- dens_dat_accum %>%
  unnest(cols = c(density, grid)) %>%
  mutate(type = "all until given time period") %>%
  bind_rows(dens_dat_single_year) %>%
  filter(
    floor(year) > 2009 & year < 2014
  ) %>%
  mutate(
    year = map_chr(
      year,
      function(yr) {
        if (floor(yr) == yr) return(paste0(yr, " (Jan - Jun)"))
        paste0(floor(yr), " (Jul - Dec)")
      }
    )
  )

dens_plot <- dens_dat %>%
  unnest(cols = c(density, grid)) %>%
  ggplot() +
  geom_hline(yintercept = 1, col = "darkgray") +
  geom_line(aes(x = grid, y = density, linetype = type, group = type)) +
  scale_linetype_manual(values = c(5, 1)) +
  facet_wrap(.~year, nrow = 2, ncol = 4) +
  theme(legend.position = "bottom") +
  labs(
    x = "PIT",
    y = "Estimated density",
    linetype = "Data"
  ) +
  ggtitle("(c)") +
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8))

pdf(
  height = 8,
  width = 8,
  file = paste0("./tex/station_", stats[i], ".pdf")
)
print(
  ggarrange(ggarrange(pithist, e_process, ncol = 2), dens_plot, nrow = 2)
)
dev.off()

### station 10162
i <- 2
pithist <- weather_results %>%
  filter(station == stats[i]) %>%
  ggplot() +
  geom_hline(yintercept = 1, col = "red") +
  geom_histogram(
    mapping = aes(x = pit, y = ..density..),
    boundary = 0,
    binwidth = 0.05
  ) +
  facet_grid(rows = vars(variable), cols = vars(leadtime)) +
  labs(x = "PIT", y = "Density") +
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) +
  ggtitle(paste0("(a) Station: ", stats[i]))

e_process <- weather_results %>%
  filter(station == stats[i]) %>%
  select(date, variable, leadtime, vec_kernel_e) %>%
  ggplot() +
  geom_hline(yintercept = c(1, 100), col = 1, lty = 3) +
  geom_line(aes(x = date, y = vec_kernel_e)) +
  facet_grid(rows = vars(variable), cols = vars(leadtime), scales = "free_y") +
  scale_y_log10() +
  labs(x = "Date", y = "E-value", color = element_blank()) +
  ggtitle("(b)")

kdens_dat <- weather_results %>%
  filter(station == stats[i] & variable == vs[i] & leadtime == lags[i]) %>%
  mutate(year = year(date))

grd <- seq(0.001, 1 - 0.001, 0.001)
yrs <- unique(kdens_dat$year)
bws <- sapply(
  yrs,
  function(x) {
    z <- kdens_dat$pit[kdens_dat$year == x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    KernSmooth::dpik(
      x = z,
      scalest = "min",
      level = 2L,
      kernel = "epanech",
      gridsize = 401L,
      range.x = c(0, 1)
    )
  }
)
dens <- mapply(
  function(x, bws) {
    z <- kdens_dat$pit[kdens_dat$year == x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    kernel_estim <- bde::jonesCorrectionMuller94BoundaryKernel(
      dataPoints = z,
      b = bws,
      mu = 1,
      dataPointsCache = grd
    )
    kernel_estim@densityCache / bde::distribution(kernel_estim, x = 1, FALSE)
  },
  x = yrs,
  bws = bws,
  SIMPLIFY = FALSE
)

dens_dat_single_year <- tibble(
  year = yrs,
  density = dens,
  grid = list(grd)
)

dens_dat_single_year <- dens_dat_single_year %>%
  unnest(cols = c(density, grid)) %>%
  mutate(type = "only given time period")

bws <- sapply(
  yrs,
  function(x) {
    z <- kdens_dat$pit[kdens_dat$year <= x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    KernSmooth::dpik(
      x = z,
      scalest = "min",
      level = 2L,
      kernel = "epanech",
      gridsize = 401L,
      range.x = c(0, 1)
    )
  }
)
dens <- mapply(
  function(x, bws) {
    z <- kdens_dat$pit[kdens_dat$year <= x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    kernel_estim <- bde::jonesCorrectionMuller94BoundaryKernel(
      dataPoints = z,
      b = bws,
      mu = 1,
      dataPointsCache = grd
    )
    kernel_estim@densityCache / bde::distribution(kernel_estim, x = 1, FALSE)
  },
  x = yrs,
  bws = bws,
  SIMPLIFY = FALSE
)

dens_dat_accum <- tibble(
  year = yrs,
  density = dens,
  grid = list(grd)
)

dens_dat <- dens_dat_accum %>%
  unnest(cols = c(density, grid)) %>%
  mutate(type = "all until given time period") %>%
  bind_rows(dens_dat_single_year) %>%
  filter(
    floor(year) > 2009 & year < 2014
  )

dens_plot <- dens_dat %>%
  unnest(cols = c(density, grid)) %>%
  ggplot() +
  geom_hline(yintercept = 1, col = "darkgray") +
  geom_line(aes(x = grid, y = density, linetype = type, group = type)) +
  scale_linetype_manual(values = c(5, 1)) +
  facet_wrap(.~year, nrow = 2, ncol = 4) +
  theme(legend.position = "bottom") +
  labs(
    x = "PIT",
    y = "Estimated density",
    linetype = "Data"
  ) +
  ggtitle("(c)") +
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8))

pdf(
  height = 6,
  width = 8,
  file = paste0("./tex/station_", stats[i], ".pdf")
)
print(
  ggarrange(
    ggarrange(pithist, e_process, ncol = 2),
    dens_plot,
    nrow = 2,
    heights = c(3, 2)
  ))
dev.off()

### station 10729
i <- 3
pithist <- weather_results %>%
  filter(station == stats[i]) %>%
  ggplot() +
  geom_hline(yintercept = 1, col = "red") +
  geom_histogram(
    mapping = aes(x = pit, y = ..density..),
    boundary = 0,
    binwidth = 0.05
  ) +
  facet_grid(rows = vars(variable), cols = vars(leadtime)) +
  labs(x = "PIT", y = "Density") +
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) +
  ggtitle(paste0("(a) Station: ", stats[i]))

e_process <- weather_results %>%
  filter(station == stats[i]) %>%
  select(date, variable, leadtime, vec_kernel_e) %>%
  ggplot() +
  geom_hline(yintercept = c(1, 100), col = 1, lty = 3) +
  geom_line(aes(x = date, y = vec_kernel_e)) +
  facet_grid(rows = vars(variable), cols = vars(leadtime), scales = "free_y") +
  scale_y_log10() +
  labs(x = "Date", y = "E-value", color = element_blank()) +
  ggtitle("(b)")

kdens_dat <- weather_results %>%
  filter(station == stats[i] & variable == vs[i] & leadtime == lags[i]) %>%
  mutate(year = year(date))

grd <- seq(0.001, 1 - 0.001, 0.001)
yrs <- unique(kdens_dat$year)
bws <- sapply(
  yrs,
  function(x) {
    z <- kdens_dat$pit[kdens_dat$year == x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    KernSmooth::dpik(
      x = z,
      scalest = "min",
      level = 2L,
      kernel = "epanech",
      gridsize = 401L,
      range.x = c(0, 1)
    )
  }
)
dens <- mapply(
  function(x, bws) {
    z <- kdens_dat$pit[kdens_dat$year == x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    kernel_estim <- bde::jonesCorrectionMuller94BoundaryKernel(
      dataPoints = z,
      b = bws,
      mu = 1,
      dataPointsCache = grd
    )
    kernel_estim@densityCache / bde::distribution(kernel_estim, x = 1, FALSE)
  },
  x = yrs,
  bws = bws,
  SIMPLIFY = FALSE
)

dens_dat_single_year <- tibble(
  year = yrs,
  density = dens,
  grid = list(grd)
)

dens_dat_single_year <- dens_dat_single_year %>%
  unnest(cols = c(density, grid)) %>%
  mutate(type = "only given time period")

bws <- sapply(
  yrs,
  function(x) {
    z <- kdens_dat$pit[kdens_dat$year <= x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    KernSmooth::dpik(
      x = z,
      scalest = "min",
      level = 2L,
      kernel = "epanech",
      gridsize = 401L,
      range.x = c(0, 1)
    )
  }
)
dens <- mapply(
  function(x, bws) {
    z <- kdens_dat$pit[kdens_dat$year <= x]
    z <- z[!is.na(z) & z > 0 & z < 1]
    kernel_estim <- bde::jonesCorrectionMuller94BoundaryKernel(
      dataPoints = z,
      b = bws,
      mu = 1,
      dataPointsCache = grd
    )
    kernel_estim@densityCache / bde::distribution(kernel_estim, x = 1, FALSE)
  },
  x = yrs,
  bws = bws,
  SIMPLIFY = FALSE
)

dens_dat_accum <- tibble(
  year = yrs,
  density = dens,
  grid = list(grd)
)

dens_dat <- dens_dat_accum %>%
  unnest(cols = c(density, grid)) %>%
  mutate(type = "all until given time period") %>%
  bind_rows(dens_dat_single_year) %>%
  filter(
    floor(year) > 2009 & year < 2014
  )

dens_plot <- dens_dat %>%
  unnest(cols = c(density, grid)) %>%
  ggplot() +
  geom_hline(yintercept = 1, col = "darkgray") +
  geom_line(aes(x = grid, y = density, linetype = type, group = type)) +
  scale_linetype_manual(values = c(5, 1)) +
  facet_wrap(.~year, nrow = 2, ncol = 4) +
  theme(legend.position = "bottom") +
  labs(
    x = "PIT",
    y = "Estimated density",
    linetype = "Data"
  ) +
  ggtitle("(c)") +
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8))

pdf(
  height = 6,
  width = 8,
  file = paste0("./tex/station_", stats[i], ".pdf")
)
print(
  ggarrange(
    ggarrange(pithist, e_process, ncol = 2),
    dens_plot,
    nrow = 2,
    heights = c(3, 2)
  )
)
dev.off()

#-------------------------------------------------------------------------------
# table for supplementary material
weather_table <- weather_results %>%
  group_by(station, leadtime, variable) %>%
  summarise(
    empirical = min(which(vec_empirical_rank_e > 10^8)) - 366,
    betabinom = min(which(vec_betabinom_rank_e > 10^8) - 366)
  ) %>%
  pivot_wider(
    names_from = "leadtime",
    values_from = c("empirical", "betabinom")
  ) %>%
  arrange(station, variable)

nwt <- nrow(weather_table)
weather_table$station[!(seq_len(nwt) %% 3 == 1)] <- ""
weather_table <- cbind(weather_table[1:33, ], weather_table[34:66, ])
kable(
  weather_table,
  format = "latex",
  booktabs = TRUE
)
