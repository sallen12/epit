# Parametrize all settings considered

## Sample sizes
n <- c(180, 360, 720)

## Forecast parameters
bias_disp <- expand.grid(
  n = n,
  bias = seq(-0.5, 0.5, 0.1),
  dispersion = seq(-0.5, 0.5, 0.1)
)

## Construct data.frame of all parameters

### Usual pit
method <- c("ptest", "beta", "kernel")
K <- 0
type <- "pit"

pit_params <- expand.grid(
  method = method,
  K = K,
  stringsAsFactors = FALSE
)

rep_ind <- rep(seq_len(nrow(pit_params)), each = nrow(bias_disp))
rep_ind_2 <- rep(seq_len(nrow(bias_disp)), nrow(pit_params))

pit_params <- data.frame(
  method = pit_params$method[rep_ind],
  K = pit_params$K[rep_ind],
  type = type,
  bias_disp[rep_ind_2, ]
)

### Quantile pit
method <- c("ptest", "grenander", "bernstein")
K <- c(9, 19)
type <- "quantile_pit"

qpit_params <- expand.grid(
  method = method,
  K = K,
  stringsAsFactors = FALSE
)

rep_ind <- rep(seq_len(nrow(qpit_params)), each = nrow(bias_disp))
rep_ind_2 <- rep(seq_len(nrow(bias_disp)), nrow(qpit_params))

qpit_params <- data.frame(
  method = qpit_params$method[rep_ind],
  K = qpit_params$K[rep_ind],
  type = type,
  bias_disp[rep_ind_2, ]
)

### Rank histograms
method <- c("ptest", "empirical", "betabinom")
K <- c(10, 20, 50)
type <- "rank_histogram"

rhist_params <- expand.grid(
  method = method,
  K = K,
  stringsAsFactors = FALSE
)

rep_ind <- rep(seq_len(nrow(rhist_params)), each = nrow(bias_disp))
rep_ind_2 <- rep(seq_len(nrow(bias_disp)), nrow(rhist_params))

rhist_params <- data.frame(
  method = rhist_params$method[rep_ind],
  K = rhist_params$K[rep_ind],
  type = type,
  bias_disp[rep_ind_2, ]
)

### Stochastic dominance testing
K <- 0
type <- "stochtest"
method <- c("ptest", "bernstein", "grenander")

stochtest_params <- expand.grid(
  K = K,
  type = type,
  method = method,
  stringsAsFactors = FALSE
)

rep_ind <- rep(seq_len(nrow(stochtest_params)), each = nrow(bias_disp))
rep_ind_2 <- rep(seq_len(nrow(bias_disp)), nrow(stochtest_params))

stochtest_params <- data.frame(
  method = stochtest_params$method[rep_ind],
  K = stochtest_params$K[rep_ind],
  type = type,
  bias_disp[rep_ind_2, ]
)

# Export
sim_params <- rbind(pit_params, qpit_params, rhist_params, stochtest_params)

sim_params <- rbind(
  sim_params[sim_params$method %in% c("bernstein", "kernel"), ],
  sim_params[sim_params$method == "betabinom", ],
  sim_params[!(sim_params$method %in% c("bernstein", "kernel", "betabinom")), ]
)
row.names(sim_params) <- NULL

range(which(sim_params$method == "bernstein"))
range(which(sim_params$method %in% c("bernstein", "kernel")))
range(which(!(sim_params$method %in% c("bernstein", "kernel"))))

save(list = "sim_params", file = "sim_params.rda")
