library(tidyverse)

weather_results <- vector("list", 22)
for (j in seq_len(22)) {
  load(paste0("./weather/weather_e_values_", j, ".rda"))
  df <- df %>%
    as.data.frame()
  weather_results[[j]] <- df
}

weather_results <- do.call(rbind, weather_results)
save(list = "weather_results", file = "./weather/weather_results.rda")
