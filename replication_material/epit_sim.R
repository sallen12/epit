# Load packages
library(epit)
source("epit_sim_functions.R")

# Load simulation parameters and array index
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
load("sim_params.rda")

method <- sim_params$method[id]
K <- sim_params$K[id]
type <- sim_params$type[id]
n <- sim_params$n[id]
bias <- sim_params$bias[id]
dispersion <- sim_params$dispersion[id]

# Container for output
out <- data.frame(
  n = n,
  K = K,
  type = type,
  bias = bias,
  dispersion = dispersion,
  method = method,
  rej_0.001 = 0,
  rej_0.01 = 0,
  rej_0.05 = 0,
  rej_0.1 = 0
)

# Run simulations
m <- 5000
pb <- txtProgressBar(max = m)
for (i in seq_len(m)) {
  setTxtProgressBar(pb, value = i)
  pit_sim <- simulate_pit(
    n = n,
    type = ifelse(type == "stochtest", "pit", type),
    bias = bias,
    dispersion = dispersion,
    K = K,
    seed = n^2 + i
  )
  pit_sim2 <- simulate_pit(
    n = n,
    type = ifelse(type == "stochtest", "pit", type),
    bias = -bias,
    dispersion = dispersion,
    K = K,
    seed = n^2 + i
  )

  pval <- apply_method(
    type = type,
    method = method,
    pit_sim = pit_sim,
    K = K
  )

  pval2 <- apply_method(
    type = type,
    method = method,
    pit_sim = pit_sim2,
    K = K
  )
  if (abs(pval - pval2) > 1e-13) break

  out$rej_0.001 <- out$rej_0.001 + (pval <= 0.001)
  out$rej_0.01 <- out$rej_0.01 + (pval <= 0.01)
  out$rej_0.05 <- out$rej_0.05 + (pval <= 0.05)
  out$rej_0.1 <- out$rej_0.1 + (pval <= 0.1)
}
close(pb)

out[, 7:10] <- out[, 7:10] / m
out

# Export
save(list = "out", file = paste0("pit_sim_", id, ".rda"))
