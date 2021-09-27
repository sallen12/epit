id <- 485:2541 # enter indices here

data <- vector("list", length(id))
for (j in seq_along(id)) {
  load(paste0("pit_sim_", id[j], ".rda"))
  data[[j]] <- out
}

data <- do.call(rbind, data)

save(list = "data", file = "aaa.rda")
