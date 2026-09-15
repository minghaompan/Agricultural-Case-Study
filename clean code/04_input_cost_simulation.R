# ============================================================
# 04_input_cost_simulation.R
# Distribution fitting and input-cost simulation
# ============================================================

# ---- 1. Packages and seed --------------------------------------------------

library(fitdistrplus)
library(actuar)
library(truncdist)

set.seed(123)


# ---- 2. Load data ----------------------------------------------------------

cost_data <- read.csv(
  "data/Crop Input Cost in Crop Alternative.csv"
)

names(cost_data) <- tolower(names(cost_data))

crops <- c("swheat", "oats", "barley", "canola")
crop_labels <- c("Spring wheat", "Oats", "Barley", "Canola")


# ---- 3. Candidate distributions -------------------------------------------

distribution_names <- c(
  "Normal",
  "Lognormal",
  "Log-logistic",
  "Uniform",
  "Gamma",
  "Weibull"
)

distribution_codes <- c(
  "norm",
  "lnorm",
  "llogis",
  "unif",
  "gamma",
  "weibull"
)

fit_candidates <- function(x) {
  setNames(
    lapply(
      distribution_codes,
      function(distribution) {
        fitdist(
          x,
          distribution,
          method = "mle"
        )
      }
    ),
    distribution_names
  )
}

cost_fits <- setNames(
  lapply(
    crops,
    function(crop) {
      fit_candidates(cost_data[[crop]])
    }
  ),
  crops
)

cost_gof <- setNames(
  lapply(
    crops,
    function(crop) {
      gofstat(
        cost_fits[[crop]],
        fitnames = distribution_names
      )
    }
  ),
  crops
)

for (crop in crops) {
  cat(
    "\nGoodness-of-fit:",
    crop,
    "\n"
  )
  print(cost_gof[[crop]])
}

for (i in seq_along(crops)) {
  qqcomp(
    cost_fits[[crops[i]]],
    legendtext = distribution_names,
    main = crop_labels[i]
  )
}

for (i in seq_along(crops)) {
  cdfcomp(
    cost_fits[[crops[i]]],
    legendtext = distribution_names,
    main = crop_labels[i]
  )
}

for (i in seq_along(crops)) {
  plot(
    cost_fits[[crops[i]]][["Log-logistic"]],
    demp = TRUE
  )
}


# ---- 4. Simulation settings ------------------------------------------------

n_years <- 20L
n_iter <- 10000L
n_sim <- n_years * n_iter

cost_parameters <- data.frame(
  crop = crops,
  shape = c(14.629, 11.137, 13.290, 13.433),
  scale = c(333.285, 275.896, 296.279, 403.916),
  lower = c(244.825, 182.871, 205.433, 288.920),
  upper = c(457.355, 422.049, 419.667, 577.890)
)

rtrunc_llogis <- function(
  n,
  shape,
  scale,
  lower,
  upper
) {
  rtrunc(
    n,
    spec = "llogis",
    a = lower,
    b = upper,
    shape = shape,
    scale = scale
  )
}


# ---- 5. Simulate input costs -----------------------------------------------

cost_draws <- lapply(
  seq_len(nrow(cost_parameters)),
  function(i) {
    rtrunc_llogis(
      n = n_sim,
      shape = cost_parameters$shape[i],
      scale = cost_parameters$scale[i],
      lower = cost_parameters$lower[i],
      upper = cost_parameters$upper[i]
    )
  }
)

names(cost_draws) <- cost_parameters$crop

simulated_input_cost <- data.frame(
  year = rep(seq_len(n_years), each = n_iter),
  iteration = rep(seq_len(n_iter), times = n_years),
  cost_draws,
  check.names = FALSE
)


# ---- 6. Save simulation ----------------------------------------------------

write.csv(
  simulated_input_cost,
  "data/simulated_inputcost.csv",
  row.names = FALSE
)


# ---- 7. Results ------------------------------------------------------------

cat("\nSimulated input-cost summary\n")
print(summary(simulated_input_cost[crops]))

cat("\nSimulated input-cost standard deviations\n")
print(
  sapply(
    simulated_input_cost[crops],
    sd,
    na.rm = TRUE
  )
)
