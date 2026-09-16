# ============================================================
# 01_price_model_selection.R
# Crop-price stationarity tests and AR lag selection
# ============================================================

# ---- 1. Packages -----------------------------------------------------------

library(urca)
library(forecast)


# ---- 2. Load data ----------------------------------------------------------

price_data <- read.csv("data/Historical crop price 1987-2024.csv")
names(price_data) <- tolower(names(price_data))
price_data <- price_data[order(price_data$years), ]


# ---- 3. Detrend oats and barley -------------------------------------------

oats_trend_model <- lm(oats ~ years, data = price_data)
barley_trend_model <- lm(barley ~ years, data = price_data)

price_data$oats_detrended <- residuals(oats_trend_model)
price_data$barley_detrended <- residuals(barley_trend_model)


# ---- 4. Stationarity tests -------------------------------------------------

test_series <- list(
  swheat = price_data$swheat,
  canola = price_data$canola,
  barley = price_data$barley,
  oats = price_data$oats,
  barley_detrended = price_data$barley_detrended,
  oats_detrended = price_data$oats_detrended
)

run_adf <- function(x, type, lag = 1L) {
  fit <- ur.df(
    x,
    type = type,
    lags = lag,
    selectlags = "Fixed"
  )

  tau_name <- switch(
    type,
    none = "tau1",
    drift = "tau2",
    trend = "tau3"
  )

  unname(fit@teststat[1, tau_name])
}

adf_results <- do.call(
  rbind,
  lapply(names(test_series), function(series_name) {
    x <- test_series[[series_name]]

    data.frame(
      series = series_name,
      lag = 1L,
      random_walk = run_adf(x, "none"),
      random_walk_with_drift = run_adf(x, "drift"),
      trend_stationary = run_adf(x, "trend")
    )
  })
)

run_kpss <- function(x, type, lag = 3L) {
  fit <- ur.kpss(
    x,
    type = type,
    lags = "nil",
    use.lag = lag
  )

  unname(fit@teststat)
}

kpss_results <- do.call(
  rbind,
  lapply(names(test_series), function(series_name) {
    x <- test_series[[series_name]]

    data.frame(
      series = series_name,
      lag = 3L,
      level_stationarity = run_kpss(x, "mu"),
      trend_stationarity = run_kpss(x, "tau")
    )
  })
)


# ---- 5. Final series -------------------------------------------------------

series_treatment <- data.frame(
  crop = c("Spring wheat", "Canola", "Oats", "Barley"),
  series_used = c(
    "Level",
    "Level",
    "Detrended residual",
    "Detrended residual"
  )
)

final_series <- list(
  swheat = price_data$swheat,
  oats = price_data$oats_detrended,
  barley = price_data$barley_detrended,
  canola = price_data$canola
)


# ---- 6. AR lag selection ---------------------------------------------------

select_ar_lag <- function(x, max_lag = 5L) {
  do.call(
    rbind,
    lapply(seq_len(max_lag), function(p) {
      fit <- Arima(
        x,
        order = c(p, 0, 0),
        include.mean = TRUE,
        method = "ML"
      )

      log_likelihood <- logLik(fit)
      n_parameters <- attr(log_likelihood, "df")
      n_observations <- length(x)

      data.frame(
        lag = p,
        AIC = AIC(fit),
        BIC = BIC(fit),
        HQIC = -2 * as.numeric(log_likelihood) +
          2 * n_parameters * log(log(n_observations))
      )
    })
  )
}

lag_results <- do.call(
  rbind,
  lapply(names(final_series), function(series_name) {
    result <- select_ar_lag(final_series[[series_name]])
    result$series <- series_name
    result[, c("series", "lag", "AIC", "BIC", "HQIC")]
  })
)

selected_lags <- do.call(
  rbind,
  lapply(split(lag_results, lag_results$series), function(x) {
    data.frame(
      series = x$series[1],
      AIC_lag = x$lag[which.min(x$AIC)],
      BIC_lag = x$lag[which.min(x$BIC)],
      HQIC_lag = x$lag[which.min(x$HQIC)]
    )
  })
)

row.names(selected_lags) <- NULL


# ---- 7. Results ------------------------------------------------------------

cat("\nADF tests: lag = 1\n")
print(adf_results, row.names = FALSE)

cat("\nKPSS tests: lag = 3\n")
print(kpss_results, row.names = FALSE)

cat("\nFinal series treatment\n")
print(series_treatment, row.names = FALSE)

cat("\nAR lag selection\n")
print(selected_lags, row.names = FALSE)
