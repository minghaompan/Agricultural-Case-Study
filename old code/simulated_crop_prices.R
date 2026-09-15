# --- Load Required Packages ---
library(systemfit)
library(tmvtnorm)
library(dplyr)
library(readxl)
library(writexl)
library(tibble)
library(tseries)
library(urca)
library(forecast)

# --- Reproducibility ---
set.seed(123)

# --- Load and Prepare Data ---
data <- read.csv("data/Historical crop price 1987-2024.csv")
colnames(data) <- tolower(colnames(data))

# --- Detrend Barley and Oats (trend-stationary crops) ---
barley_trend_model <- lm(barley ~ years, data = data)
oats_trend_model   <- lm(oats   ~ years, data = data)

# Keep model summary output (unchanged behavior)
summary(barley_trend_model)

# Store residuals (detrended series used in tests, lag selection, and SUR)
data$barley_detrended <- residuals(barley_trend_model)
data$oats_detrended   <- residuals(oats_trend_model)

# --- Test each crop for unit roots (ADF/KPSS) ---
crop_list <- c("swheat", "oats", "barley", "canola", "barley_detrended", "oats_detrended")
start_year <- min(data$years)

for (crop in crop_list) {
  cat("\n=========================================\n")
  cat("Crop:", crop, "\n")
  cat("=========================================\n")
  
  ts_data <- ts(data[[crop]], start = start_year, frequency = 1)
  
  # Visual check
  plot(ts_data, main = paste("Time Series -", crop), ylab = "Price", xlab = "Year")
  
  # ----- ADF tests (no constant, drift, trend) -----
  cat("\n--- Augmented Dickey-Fuller (ADF) Tests ---\n")
  
  adf_none  <- ur.df(ts_data, type = "none",  selectlags = "AIC")
  cat("\nADF (no constant):\n");  print(summary(adf_none))
  
  adf_drift <- ur.df(ts_data, type = "drift", selectlags = "AIC")
  cat("\nADF (with drift):\n");   print(summary(adf_drift))
  
  adf_trend <- ur.df(ts_data, type = "trend", selectlags = "AIC")
  cat("\nADF (with trend):\n");   print(summary(adf_trend))
  
  # ----- KPSS tests (level and trend stationarity) -----
  cat("\n--- KPSS Tests ---\n")
  kpss_level <- kpss.test(ts_data, null = "Level")
  cat("\nKPSS (level-stationary - drift):\n"); print(kpss_level)
  
  kpss_trend <- kpss.test(ts_data, null = "Trend")
  cat("\nKPSS (trend-stationary):\n");        print(kpss_trend)
}

# --- Time series objects (explicit starts set to 1987; matches original) ---
swheat_ts   <- ts(data$swheat,           start = 1987, frequency = 1)
oats_ts     <- ts(data$oats,             start = 1987, frequency = 1)
barley_ts   <- ts(data$barley,           start = 1987, frequency = 1)
canola_ts   <- ts(data$canola,           start = 1987, frequency = 1)
debarley_ts <- ts(data$barley_detrended, start = 1987, frequency = 1)
deoats_ts   <- ts(data$oats_detrended,   start = 1987, frequency = 1)

# --- Lag selection helper (AIC/BIC/HQIC) ---
select_lag <- function(ts_data, max_lag = 5) {
  ts_data <- na.omit(ts_data)
  n <- length(ts_data)
  results <- data.frame(Lag = 1:max_lag, AIC = NA, BIC = NA, HQIC = NA)
  
  for (lag in 1:max_lag) {
    fit <- Arima(ts_data, order = c(lag, 0, 0))
    loglik <- as.numeric(logLik(fit))
    k <- length(fit$coef)
    
    # AIC/BIC from forecast::AIC/BIC; HQIC manually computed
    aic  <- AIC(fit)
    bic  <- BIC(fit)
    hqic <- -2 * loglik + 2 * k * log(log(n))
    
    results[lag, ] <- c(lag, aic, bic, hqic)
  }
  results
}

# --- Apply lag selection (unchanged) ---
swheat_lags   <- select_lag(swheat_ts)
oats_lags     <- select_lag(oats_ts)
barley_lags   <- select_lag(barley_ts)
canola_lags   <- select_lag(canola_ts)
debarley_lags <- select_lag(debarley_ts)
deoats_lags   <- select_lag(deoats_ts)

# Print results (preserved)
print("Spring Wheat Lag Selection:"); print(swheat_lags)
print("Oats Lag Selection:");         print(oats_lags)
print("Barley Lag Selection:");       print(barley_lags)
print("Canola Lag Selection:");       print(canola_lags)
print("Debarley Lag Selection:");     print(debarley_lags)
print("Deoats Lag Selection:");       print(deoats_lags)

# Extract best lags by criterion (unchanged)
get_best_lag <- function(results) {
  list(
    AIC  = results$Lag[which.min(results$AIC)],
    BIC  = results$Lag[which.min(results$BIC)],
    HQIC = results$Lag[which.min(results$HQIC)]
  )
}

# Print again (as in original)
cat("Spring Wheat Lag Selection\n"); print(swheat_lags)
cat("Oats Lag Selection\n");         print(oats_lags)
cat("Barley Lag Selection\n");       print(barley_lags)
cat("Canola Lag Selection\n");       print(canola_lags)
cat("Debarley Lag Selection\n");     print(debarley_lags)
cat("Deoats Lag Selection\n");       print(deoats_lags)

print(list(
  Swheat   = get_best_lag(swheat_lags),
  Oats     = get_best_lag(oats_lags),
  Barley   = get_best_lag(barley_lags),
  Canola   = get_best_lag(canola_lags),
  Debarley = get_best_lag(debarley_lags),
  Deoats   = get_best_lag(deoats_lags)
))

# --- Create lagged variables for SUR and drop NAs (unchanged design) ---
data <- data %>%
  mutate(
    swheat_lag1            = lag(swheat, 1),
    swheat_lag2            = lag(swheat, 2),
    oats_detrended_lag1    = lag(oats_detrended, 1),
    oats_detrended_lag2    = lag(oats_detrended, 2),
    barley_detrended_lag1  = lag(barley_detrended, 1),
    barley_detrended_lag2  = lag(barley_detrended, 2),
    canola_lag1            = lag(canola, 1),
    canola_lag2            = lag(canola, 2)
  ) %>%
  na.omit()

# --- SUR Model (AR(2) per equation; oats/barley use detrended series) ---
eq1 <- swheat           ~ swheat_lag1 + swheat_lag2
eq2 <- oats_detrended   ~ oats_detrended_lag1 + oats_detrended_lag2
eq3 <- barley_detrended ~ barley_detrended_lag1 + barley_detrended_lag2
eq4 <- canola           ~ canola_lag1 + canola_lag2

system <- list(swheat = eq1, oats = eq2, barley = eq3, canola = eq4)
fit <- systemfit(system, method = "SUR", data = data)
summary(fit)

b <- coef(fit)          # SUR coefficients for mean recursion
resid_cov <- fit$residCov  # Cross-equation covariance for shocks

s <- summary(fit)
s$system$OLS.R2
s$system$McElroy.R2

# --- Simulation settings (unchanged) ---
n_periods <- 20
n_sim     <- 10000
n_crops   <- 4

simulated_prices <- array(NA, dim = c(n_periods, n_crops, n_sim))
dimnames(simulated_prices) <- list(1:n_periods, c("swheat", "oats", "barley", "canola"), 1:n_sim)

# --- Last observed values for AR recursion (unchanged) ---
last_obs <- data[nrow(data), ]

# NOTE: This produces future_years = 1:20 as in your original code.
future_years <- max(0) + 1:n_periods

# --- Demeaned trend projections for barley and oats (add trend back later) ---
barley_trend <- predict(barley_trend_model, newdata = data.frame(years = future_years))
oats_trend   <- predict(oats_trend_model,   newdata = data.frame(years = future_years))

mean_barley_hist <- mean(data$barley, na.rm = TRUE)
mean_oats_hist   <- mean(data$oats,   na.rm = TRUE)

barley_trend_demeaned <- barley_trend + (mean_barley_hist - mean(barley_trend))
oats_trend_demeaned   <- oats_trend   + (mean_oats_hist   - mean(oats_trend))

# --- Build empirical bounds for shocks from historical one-step changes (unchanged) ---
price_changes <- data %>%
  transmute(
    swheat = swheat - lag(swheat),
    oats   = oats_detrended - lag(oats_detrended),
    barley = barley_detrended - lag(barley_detrended),
    canola = canola - lag(canola)
  ) %>%
  na.omit()

# Using 0% and 100% quantiles reproduces min/max bounds from history
lower_bounds <- apply(price_changes, 2, quantile, probs = 0)
upper_bounds <- apply(price_changes, 2, quantile, probs = 1)

# --- Simulate future prices with SUR-based mean + truncated MVN shocks (unchanged) ---
for (s in 1:n_sim) {
  # Initialize with last two observations per series
  prev_swheat_1 <- last_obs$swheat
  prev_swheat_2 <- data$swheat[nrow(data) - 1]
  
  prev_oats_1   <- last_obs$oats_detrended
  prev_oats_2   <- data$oats_detrended[nrow(data) - 1]
  
  prev_barley_1 <- last_obs$barley_detrended
  prev_barley_2 <- data$barley_detrended[nrow(data) - 1]
  
  prev_canola_1 <- last_obs$canola
  prev_canola_2 <- data$canola[nrow(data) - 1]
  
  for (t in 1:n_periods) {
    # SUR mean recursion (AR(2) in each equation, correlated across crops via resid_cov)
    mean_prices <- c(
      b["swheat_(Intercept)"] +
        b["swheat_swheat_lag1"] * prev_swheat_1 +
        b["swheat_swheat_lag2"] * prev_swheat_2,
      
      b["oats_(Intercept)"] +
        b["oats_oats_detrended_lag1"] * prev_oats_1 +
        b["oats_oats_detrended_lag2"] * prev_oats_2,
      
      b["barley_(Intercept)"] +
        b["barley_barley_detrended_lag1"] * prev_barley_1 +
        b["barley_barley_detrended_lag2"] * prev_barley_2,
      
      b["canola_(Intercept)"] +
        b["canola_canola_lag1"] * prev_canola_1 +
        b["canola_canola_lag2"] * prev_canola_2
    )
    
    # Draw one correlated shock vector within empirical min/max bounds
    simulated_errors <- rtmvnorm(
      n = 1,
      mean  = rep(0, n_crops),
      sigma = resid_cov,
      lower = lower_bounds,
      upper = upper_bounds
    )
    
    # Construct prices (oats/barley: add back trend component; all truncated at zero)
    simulated_prices[t, "swheat", s] <- max(0, mean_prices[1] + simulated_errors[1])
    simulated_prices[t, "oats",   s] <- max(0, oats_trend_demeaned[t]   + mean_prices[2] + simulated_errors[2])
    simulated_prices[t, "barley", s] <- max(0, barley_trend_demeaned[t] + mean_prices[3] + simulated_errors[3])
    simulated_prices[t, "canola", s] <- max(0, mean_prices[4] + simulated_errors[4])
    
    # Roll the AR(2) state forward (note: oats/barley states are on the detrended scale)
    prev_swheat_2 <- prev_swheat_1
    prev_swheat_1 <- simulated_prices[t, "swheat", s]
    
    prev_oats_2 <- prev_oats_1
    prev_oats_1 <- mean_prices[2] + simulated_errors[2]
    
    prev_barley_2 <- prev_barley_1
    prev_barley_1 <- mean_prices[3] + simulated_errors[3]
    
    prev_canola_2 <- prev_canola_1
    prev_canola_1 <- simulated_prices[t, "canola", s]
  }
}

# --- Long-format output for downstream use (unchanged) ---
simulated_data <- data.frame(
  year      = rep(future_years, each = n_sim),
  Iteration = rep(1:n_sim, times = n_periods),
  swheat    = as.vector(simulated_prices[, "swheat", ]),
  oats      = as.vector(simulated_prices[, "oats",   ]),
  barley    = as.vector(simulated_prices[, "barley", ]),
  canola    = as.vector(simulated_prices[, "canola", ])
)

# --- Preview + simple diagnostics (unchanged) ---
head(simulated_data)
summary(simulated_data)
sapply(simulated_data, sd, na.rm = TRUE)

# 10th and 90th percentiles per crop
percentiles <- sapply(simulated_data[, c("swheat", "oats", "barley", "canola")],
                      quantile, probs = c(0.10, 0.90))
percentiles <- t(percentiles)
print(percentiles)

# --- Save output ---
write.csv(
  simulated_data,
  "data/simulated_price.csv",
  row.names = FALSE
)


