# ============================================================
# 02_price_simulation.R
# AR(2)-SUR simulation of crop prices
# ============================================================

# ---- 1. Packages and seed --------------------------------------------------

library(dplyr)
library(systemfit)
library(tmvtnorm)

set.seed(123)


# ---- 2. Load data ----------------------------------------------------------

price_data <- read.csv("data/Historical crop price 1987-2024.csv")
names(price_data) <- tolower(names(price_data))
price_data <- price_data[order(price_data$years), ]


# ---- 3. Detrend oats and barley -------------------------------------------

oats_trend_model <- lm(oats ~ years, data = price_data)
barley_trend_model <- lm(barley ~ years, data = price_data)

price_data$oats_detrended <- residuals(oats_trend_model)
price_data$barley_detrended <- residuals(barley_trend_model)


# ---- 4. AR(2)-SUR model ----------------------------------------------------

model_data <- price_data %>%
  mutate(
    swheat_lag1 = lag(swheat, 1),
    swheat_lag2 = lag(swheat, 2),
    oats_lag1 = lag(oats_detrended, 1),
    oats_lag2 = lag(oats_detrended, 2),
    barley_lag1 = lag(barley_detrended, 1),
    barley_lag2 = lag(barley_detrended, 2),
    canola_lag1 = lag(canola, 1),
    canola_lag2 = lag(canola, 2)
  ) %>%
  na.omit()

equations <- list(
  swheat = swheat ~ swheat_lag1 + swheat_lag2,
  oats = oats_detrended ~ oats_lag1 + oats_lag2,
  barley = barley_detrended ~ barley_lag1 + barley_lag2,
  canola = canola ~ canola_lag1 + canola_lag2
)

sur_fit <- systemfit(
  equations,
  method = "SUR",
  data = model_data
)

beta <- coef(sur_fit)
sigma <- sur_fit$residCov


# ---- 5. Shock bounds -------------------------------------------------------

ols_fits <- lapply(
  equations,
  lm,
  data = model_data
)

ols_residuals <- do.call(
  cbind,
  lapply(ols_fits, residuals)
)

colnames(ols_residuals) <- c("swheat", "oats", "barley", "canola")

lower_bounds <- apply(ols_residuals, 2, min)
upper_bounds <- apply(ols_residuals, 2, max)

residual_bounds <- data.frame(
  crop = names(lower_bounds),
  lower = as.numeric(lower_bounds),
  upper = as.numeric(upper_bounds)
)


# ---- 6. Simulation settings ------------------------------------------------

n_periods <- 20L
n_sim <- 10000L
n_crops <- length(equations)
n_draws <- n_periods * n_sim

last_year <- max(price_data$years)
future_years <- last_year + seq_len(n_periods)

oats_trend <- predict(
  oats_trend_model,
  newdata = data.frame(years = future_years)
)

barley_trend <- predict(
  barley_trend_model,
  newdata = data.frame(years = future_years)
)

oats_trend <- oats_trend +
  mean(price_data$oats, na.rm = TRUE) -
  mean(oats_trend)

barley_trend <- barley_trend +
  mean(price_data$barley, na.rm = TRUE) -
  mean(barley_trend)


# ---- 7. Correlated price shocks --------------------------------------------

shock_draws <- rtmvnorm(
  n = n_draws,
  mean = rep(0, n_crops),
  sigma = sigma,
  lower = lower_bounds,
  upper = upper_bounds,
  algorithm = "gibbs",
  burn.in.samples = 100,
  thinning = 5,
  start.value = rep(0, n_crops)
)

shock_draws <- shock_draws[
  sample.int(n_draws),
  ,
  drop = FALSE
]

simulated_shocks <- aperm(
  array(
    t(shock_draws),
    dim = c(n_crops, n_periods, n_sim)
  ),
  c(2, 1, 3)
)

dimnames(simulated_shocks) <- list(
  year = future_years,
  crop = c("swheat", "oats", "barley", "canola"),
  iteration = seq_len(n_sim)
)


# ---- 8. Simulate crop prices -----------------------------------------------

simulated_prices <- array(
  NA_real_,
  dim = c(n_periods, n_crops, n_sim),
  dimnames = dimnames(simulated_shocks)
)

prev_swheat_1 <- rep(tail(model_data$swheat, 1), n_sim)
prev_swheat_2 <- rep(tail(model_data$swheat, 2)[1], n_sim)

prev_oats_1 <- rep(tail(model_data$oats_detrended, 1), n_sim)
prev_oats_2 <- rep(tail(model_data$oats_detrended, 2)[1], n_sim)

prev_barley_1 <- rep(tail(model_data$barley_detrended, 1), n_sim)
prev_barley_2 <- rep(tail(model_data$barley_detrended, 2)[1], n_sim)

prev_canola_1 <- rep(tail(model_data$canola, 1), n_sim)
prev_canola_2 <- rep(tail(model_data$canola, 2)[1], n_sim)

for (t in seq_len(n_periods)) {
  mu_swheat <-
    beta["swheat_(Intercept)"] +
    beta["swheat_swheat_lag1"] * prev_swheat_1 +
    beta["swheat_swheat_lag2"] * prev_swheat_2

  mu_oats <-
    beta["oats_(Intercept)"] +
    beta["oats_oats_lag1"] * prev_oats_1 +
    beta["oats_oats_lag2"] * prev_oats_2

  mu_barley <-
    beta["barley_(Intercept)"] +
    beta["barley_barley_lag1"] * prev_barley_1 +
    beta["barley_barley_lag2"] * prev_barley_2

  mu_canola <-
    beta["canola_(Intercept)"] +
    beta["canola_canola_lag1"] * prev_canola_1 +
    beta["canola_canola_lag2"] * prev_canola_2

  eps <- simulated_shocks[t, , ]

  oats_detrended_t <- mu_oats + eps["oats", ]
  barley_detrended_t <- mu_barley + eps["barley", ]

  swheat_t <- pmax(0, mu_swheat + eps["swheat", ])
  oats_t <- pmax(0, oats_trend[t] + oats_detrended_t)
  barley_t <- pmax(0, barley_trend[t] + barley_detrended_t)
  canola_t <- pmax(0, mu_canola + eps["canola", ])

  simulated_prices[t, "swheat", ] <- swheat_t
  simulated_prices[t, "oats", ] <- oats_t
  simulated_prices[t, "barley", ] <- barley_t
  simulated_prices[t, "canola", ] <- canola_t

  prev_swheat_2 <- prev_swheat_1
  prev_swheat_1 <- swheat_t

  prev_oats_2 <- prev_oats_1
  prev_oats_1 <- oats_detrended_t

  prev_barley_2 <- prev_barley_1
  prev_barley_1 <- barley_detrended_t

  prev_canola_2 <- prev_canola_1
  prev_canola_1 <- canola_t
}


# ---- 9. Save simulation ----------------------------------------------------

simulated_data <- data.frame(
  year = rep(future_years, times = n_sim),
  iteration = rep(seq_len(n_sim), each = n_periods),
  swheat = as.vector(simulated_prices[, "swheat", ]),
  oats = as.vector(simulated_prices[, "oats", ]),
  barley = as.vector(simulated_prices[, "barley", ]),
  canola = as.vector(simulated_prices[, "canola", ])
)

write.csv(
  simulated_data,
  "data/simulated_price.csv",
  row.names = FALSE
)


# ---- 10. Results -----------------------------------------------------------

cat("\nSUR model\n")
print(summary(sur_fit))

cat("\nOLS residual bounds\n")
print(residual_bounds, row.names = FALSE)

cat("\nSimulated price summary\n")
print(summary(simulated_data[c("swheat", "oats", "barley", "canola")]))

cat("\nSimulated price correlations\n")
print(cor(simulated_data[c("swheat", "oats", "barley", "canola")]))
