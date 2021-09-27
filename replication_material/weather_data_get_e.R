#-------------------------------------------------------------------------------
# packages, data, environment variables, functions
library(epit)
library(dplyr)
library(tidyr)
library(purrr)

load("processed_weather_data.rda")
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

e_vec_merge <- function(e_output) {
  if (e_output$h == 1) return(cumprod(e_output$e))
  es <- lapply(e_output$evalues_h, function(x) x$e)
  epit:::evalue_combine_h(es)
}

#-------------------------------------------------------------------------------
# compute e-values

stat_id <- unique(data$station)[id]
df <- data %>%
  filter(station == stat_id) %>%
  group_by(leadtime, variable) %>%
  nest() %>%
  mutate(
    h = as.integer(leadtime / 24),
    data = map(data, ~arrange(., date)),
    beta_e = map2(
      .x = data,
      .y = h,
      ~e_pit(
        z = .x$pit,
        h = .y,
        strategy = "beta",
        options = list(n0 = as.integer(366 / .y))
      )
    ),
    kernel_e = map2(
      .x = data,
      .y = h,
      ~e_pit(
        z = .x$pit,
        h = .y,
        strategy = "kernel",
        options = list(n0 = as.integer(366 / .y))
      )
    ),
    betabinom_rank_e = map2(
      .x = data,
      .y = h,
      ~e_rank_histogram(
        r = .x$rank,
        h = .y,
        m = 50,
        strategy = "betabinom",
        options = list(n0 = as.integer(366 / .y))
      )
    ),
    empirical_rank_e = map2(
      .x = data,
      .y = h,
      ~e_rank_histogram(
        r = .x$rank,
        h = .y,
        m = 50,
        strategy = "empirical",
        options = list(n0 = as.integer(366 / .y))
      )
    ),
    vec_beta_e = map(beta_e, e_vec_merge),
    vec_kernel_e = map(kernel_e, e_vec_merge),
    vec_betabinom_rank_e = map(betabinom_rank_e, e_vec_merge),
    vec_empirical_rank_e = map(empirical_rank_e, e_vec_merge)
  ) %>%
  select(-h, -beta_e, -betabinom_rank_e, -empirical_rank_e, -kernel_e) %>%
  unnest(cols =
    c(
      vec_beta_e,
      vec_kernel_e,
      vec_betabinom_rank_e,
      vec_empirical_rank_e,
      data
    )
  ) %>%
  ungroup()

#-------------------------------------------------------------------------------
# export
save(list = "df", file = paste0("weather_e_values_", id, ".rda"))
