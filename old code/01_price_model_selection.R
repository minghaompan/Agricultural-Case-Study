# ============================================================
# 01_model_selection.R
# ============================================================

library(urca)
library(forecast)

# ---- 1. Load data ----------------------------------------------------------
data <- read.csv("data/Historical crop price 1987-2024.csv")
names(data) <- tolower(names(data))
data <- data[order(data$years), ]

required_vars <- c("years", "swheat", "oats", "barley", "canola")

# ---- 2. Detrend oats and barley -------------------------------------------
oats_trend_model   <- lm(oats ~ years, data = data)
barley_trend_model <- lm(barley ~ years, data = data)

data$oats_detrended   <- residuals(oats_trend_model)
data$barley_detrended <- residuals(barley_trend_model)

# ---- 3. Stationarity tests -------------------------------------------------
test_series <- list(
  swheat            = data$swheat,
  canola            = data$canola,
  barley            = data$barley,
  oats              = data$oats,
  barley_detrended  = data$barley_detrended,
  oats_detrended    = data$oats_detrended
)

run_adf <- function(x, type, lag = 1L) {
  fit <- ur.df(
    x,
    type = type,
    lags = lag,
    selectlags = "Fixed"
  )

  tau <- switch(
    type,
    none  = "tau1",
    drift = "tau2",
    trend = "tau3"
  )

  unname(fit@teststat[1, tau])
}

adf_results <- do.call(
  rbind,
  lapply(names(test_series), function(series_name) {
    x <- test_series[[series_name]]

    data.frame(
      Series = series_name,
      Lag = 1L,
      Random_Walk = run_adf(x, "none"),
      Random_Walk_with_Drift = run_adf(x, "drift"),
      Trend_Stationary = run_adf(x, "trend")
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
      Series = series_name,
      Lag = 3L,
      Level_Stationarity = run_kpss(x, "mu"),
      Trend_Stationarity = run_kpss(x, "tau")
    )
  })
)

# ---- 4. Final series treatment ---------------------------------------------
final_treatment <- data.frame(
  Crop = c("Spring wheat", "Canola", "Oats", "Barley"),
  Series_Used = c(
    "Level",
    "Level",
    "Detrended residual",
    "Detrended residual"
  )
)

# ---- 5. AR lag-order selection ---------------------------------------------
final_series <- list(
  swheat = data$swheat,
  oats   = data$oats_detrended,
  barley = data$barley_detrended,
  canola = data$canola
)

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

      ll <- logLik(fit)
      k  <- attr(ll, "df")
      n  <- length(x)

      data.frame(
        Lag = p,
        AIC = AIC(fit),
        BIC = BIC(fit),
        HQIC = -2 * as.numeric(ll) + 2 * k * log(log(n))
      )
    })
  )
}

lag_results <- do.call(
  rbind,
  lapply(names(final_series), function(series_name) {
    result <- select_ar_lag(final_series[[series_name]])
    result$Series <- series_name
    result[, c("Series", "Lag", "AIC", "BIC", "HQIC")]
  })
)

selected_lags <- do.call(
  rbind,
  lapply(split(lag_results, lag_results$Series), function(x) {
    data.frame(
      Series = x$Series[1],
      AIC_Lag = x$Lag[which.min(x$AIC)],
      BIC_Lag = x$Lag[which.min(x$BIC)],
      HQIC_Lag = x$Lag[which.min(x$HQIC)]
    )
  })
)

row.names(selected_lags) <- NULL

# ---- 6. Console summary ----------------------------------------------------
cat("\n================ ADF Tests: Lag = 1 ================\n")
print(adf_results, row.names = FALSE)

cat("\n================ KPSS Tests: Lag = 3 ===============\n")
print(kpss_results, row.names = FALSE)

cat("\n================ Final Series Treatment ============\n")
print(final_treatment, row.names = FALSE)

cat("\n================ AR Lag Selection ==================\n")
print(selected_lags, row.names = FALSE)
