# -----------------------------
# Packages
# -----------------------------
library(vars)
library(tseries)
library(readxl)
library(urca)
library(writexl)
library(fitdistrplus)
library(pastecs)
library(stats)
library(actuar)
library(MASS)

# =========================================
# PART A: De-trend yields and save to CSV
# =========================================

# Import the raw CSV
data <- read.csv("data/Crop Yield in AthbascaAB.csv", check.names = FALSE)

# Make all column names lower-case
names(data) <- tolower(names(data))

# Quick peek
head(data)

# -----------------------------
# Helper: detrend one crop series
# -----------------------------
detrend_to_base_year <- function(df, ycol, year_col = "cropyear", base_year = 2024) {
  ok <- is.finite(df[[ycol]]) & is.finite(df[[year_col]])
  df_fit <- df[ok, , drop = FALSE]
  
  if (nrow(df_fit) < 2) {
    stop(sprintf("Not enough valid observations to fit OLS for %s.", ycol))
  }
  
  m <- lm(as.formula(paste0(ycol, " ~ ", year_col)), data = df_fit)
  
  # predict trend for all rows
  trend_all <- predict(m, newdata = df)
  
  # base-year trend (first match if multiple)
  idx_base <- which(df[[year_col]] == base_year & is.finite(trend_all))
  if (length(idx_base) == 0) {
    stop(sprintf("Base year %s not found (or trend NA) for %s.", base_year, ycol))
  }
  base_trend <- trend_all[idx_base[1]]
  
  detrended <- df[[ycol]] - trend_all + base_trend
  detrended[!is.finite(df[[ycol]])] <- NA_real_
  
  list(model = m, trend = trend_all, detrended = detrended)
}

# -----------------------------
# Spring wheat
# -----------------------------
res_sw <- detrend_to_base_year(data, "swheat", base_year = 2024)
data$trend_swheat <- res_sw$trend
data$detrended_swheat <- res_sw$detrended
summary(res_sw$model)
par(mfrow = c(1, 1))
plot(data$cropyear, data$swheat, type = "l",
     main = "Spring Wheat: Original & De-trended (Base Year 2024)",
     xlab = "Year", ylab = "Yield")
lines(data$cropyear, data$trend_swheat, lty = 2)
lines(data$cropyear, data$detrended_swheat)
abline(v = 2024, lty = 2)
legend("bottomright", legend = c("Original", "Trend", "De-trended"),
       lty = c(1, 2, 1), bty = "n")

# -----------------------------
# Oats
# -----------------------------
res_oats <- detrend_to_base_year(data, "oats", base_year = 2024)
data$trend_oats <- res_oats$trend
data$detrended_oats <- res_oats$detrended
summary(res_oats$model)

plot(data$cropyear, data$oats, type = "l",
     main = "Oats: Original & De-trended (Base Year 2024)",
     xlab = "Year", ylab = "Yield")
lines(data$cropyear, data$trend_oats, lty = 2)
lines(data$cropyear, data$detrended_oats)
abline(v = 2024, lty = 2)
legend("bottomright", legend = c("Original", "Trend", "De-trended"),
       lty = c(1, 2, 1), bty = "n")

# -----------------------------
# Barley
# -----------------------------
res_barley <- detrend_to_base_year(data, "barley", base_year = 2024)
data$trend_barley <- res_barley$trend
data$detrended_barley <- res_barley$detrended
summary(res_barley$model)

plot(data$cropyear, data$barley, type = "l",
     main = "Barley: Original & De-trended (Base Year 2024)",
     xlab = "Year", ylab = "Yield")
lines(data$cropyear, data$trend_barley, lty = 2)
lines(data$cropyear, data$detrended_barley)
abline(v = 2024, lty = 2)
legend("bottomright", legend = c("Original", "Trend", "De-trended"),
       lty = c(1, 2, 1), bty = "n")

# -----------------------------
# Canola
# -----------------------------
res_canola <- detrend_to_base_year(data, "canola", base_year = 2024)
data$trend_canola <- res_canola$trend
data$detrended_canola <- res_canola$detrended
summary(res_canola$model)

plot(data$cropyear, data$canola, type = "l",
     main = "Canola: Original & De-trended (Base Year 2024)",
     xlab = "Year", ylab = "Yield")
lines(data$cropyear, data$trend_canola, lty = 2)
lines(data$cropyear, data$detrended_canola)
abline(v = 2024, lty = 2)
legend("bottomright", legend = c("Original", "Trend", "De-trended"),
       lty = c(1, 2, 1), bty = "n")

# Summary checks
summary(data)
sapply(data[, c("swheat","oats","barley","canola",
                "detrended_swheat","detrended_oats","detrended_barley","detrended_canola")],
       sd, na.rm = TRUE)

# Save
write.csv(
  data,
  "data/de-trended_crop_yield_in_athabasca.csv",
  row.names = FALSE
)
# ==========================================================
# PART B: Fit distributions to de-trended series (as given)
# ==========================================================

# Re-import the saved 
data <- read.csv("data/de-trended_crop_yield_in_athabasca.csv")
head(data)

# ---------- Spring Wheat ----------
fn.sw  <- fitdist(data$detrended_swheat, "norm",  method = "mle"); summary(fn.sw)
fln.sw <- fitdist(data$detrended_swheat, "lnorm", method = "mle"); summary(fln.sw)
fll.sw <- fitdist(data$detrended_swheat, "llogis", method = "mle"); summary(fll.sw)
fu.sw  <- fitdist(data$detrended_swheat, "unif",  method = "mle"); summary(fu.sw)
fg.sw  <- fitdist(data$detrended_swheat, "gamma", method = "mle"); summary(fg.sw)
fw.sw  <- fitdist(data$detrended_swheat, "weibull", method = "mle"); summary(fw.sw)

# Weibull parameter adjustment
alpha_old <- 9.88
beta_old  <- 1849.79
mean_original <- beta_old * gamma(1 + 1/alpha_old)
sd_original   <- sqrt(beta_old^2 * (gamma(1 + 2/alpha_old) - (gamma(1 + 1/alpha_old))^2))
sd_new <- 263.86

weibull_sd <- function(alpha_new, beta_new) {
  gamma1 <- gamma(1 + 1/alpha_new)
  gamma2 <- gamma(1 + 2/alpha_new)
  sqrt(beta_new^2 * (gamma2 - gamma1^2))
}
adjust_weibull <- function(alpha_old, beta_old, mean_original, sd_new) {
  alpha_new <- uniroot(function(alpha_new) {
    beta_new <- beta_old * gamma(1 + 1/alpha_old) / gamma(1 + 1/alpha_new)
    weibull_sd(alpha_new, beta_new) - sd_new
  }, lower = 0.1, upper = 10)$root
  beta_new <- beta_old * gamma(1 + 1/alpha_old) / gamma(1 + 1/alpha_new)
  list(alpha_new = alpha_new, beta_new = beta_new)
}
adjusted_params <- adjust_weibull(alpha_old, beta_old, mean_original, sd_new)
adjusted_params

quantile(fln.sw, probs = 0.05)
quantile(fln.sw, probs = 0.01)

gofstat(list(fn.sw, fln.sw, fll.sw, fu.sw, fg.sw, fw.sw),
        fitnames = c("Normal", "Lognormal", "Log-logistic", "uniform", "Gamma", "Weibull"))

# ---------- Oats ----------
fn.o  <- fitdist(data$detrended_oats, "norm",  method = "mle"); summary(fn.o)
fln.o <- fitdist(data$detrended_oats, "lnorm", method = "mle"); summary(fln.o)
fll.o <- fitdist(data$detrended_oats, "llogis", method = "mle"); summary(fll.o)
fu.o  <- fitdist(data$detrended_oats, "unif",  method = "mle"); summary(fu.o)
fg.o  <- fitdist(data$detrended_oats, "gamma", method = "mle"); summary(fg.o)
fw.o  <- fitdist(data$detrended_oats, "weibull", method = "mle"); summary(fw.o)

alpha_old <- 7.37
beta_old  <- 1965.27
mean_original <- beta_old * gamma(1 + 1/alpha_old)
sd_original   <- sqrt(beta_old^2 * (gamma(1 + 2/alpha_old) - (gamma(1 + 1/alpha_old))^2))
sd_new <- 353.27

weibull_sd <- function(alpha_new, beta_new) {
  gamma1 <- gamma(1 + 1/alpha_new)
  gamma2 <- gamma(1 + 2/alpha_new)
  sqrt(beta_new^2 * (gamma2 - gamma1^2))
}
adjust_weibull <- function(alpha_old, beta_old, mean_original, sd_new) {
  alpha_new <- uniroot(function(alpha_new) {
    beta_new <- beta_old * gamma(1 + 1/alpha_old) / gamma(1 + 1/alpha_new)
    weibull_sd(alpha_new, beta_new) - sd_new
  }, lower = 0.1, upper = 10)$root
  beta_new <- beta_old * gamma(1 + 1/alpha_old) / gamma(1 + 1/alpha_new)
  list(alpha_new = alpha_new, beta_new = beta_new)
}
adjusted_params <- adjust_weibull(alpha_old, beta_old, mean_original, sd_new)
adjusted_params

quantile(fln.o, probs = 0.05)
quantile(fln.o, probs = 0.01)

gofstat(list(fn.o, fln.o, fll.o, fu.o, fg.o, fw.o),
        fitnames = c("Normal", "Lognormal", "Log-logistic", "uniform", "Gamma", "Weibull"))

# ---------- Barley ----------
fn.b  <- fitdist(data$detrended_barley, "norm",  method = "mle"); summary(fn.b)
fln.b <- fitdist(data$detrended_barley, "lnorm", method = "mle"); summary(fln.b)
fll.b <- fitdist(data$detrended_barley, "llogis", method = "mle"); summary(fll.b)
fu.b  <- fitdist(data$detrended_barley, "unif",  method = "mle"); summary(fu.b)
fg.b  <- fitdist(data$detrended_barley, "gamma", method = "mle"); summary(fg.b)
fw.b  <- fitdist(data$detrended_barley, "weibull", method = "mle"); summary(fw.b)

alpha_old <- 7.75
beta_old  <- 1770.16
mean_original <- beta_old * gamma(1 + 1/alpha_old)
sd_original   <- sqrt(beta_old^2 * (gamma(1 + 2/alpha_old) - (gamma(1 + 1/alpha_old))^2))
sd_new <- 321.38

weibull_sd <- function(alpha_new, beta_new) {
  gamma1 <- gamma(1 + 1/alpha_new)
  gamma2 <- gamma(1 + 2/alpha_new)
  sqrt(beta_new^2 * (gamma2 - gamma1^2))
}
adjust_weibull <- function(alpha_old, beta_old, mean_original, sd_new) {
  alpha_new <- uniroot(function(alpha_new) {
    beta_new <- beta_old * gamma(1 + 1/alpha_old) / gamma(1 + 1/alpha_new)
    weibull_sd(alpha_new, beta_new) - sd_new
  }, lower = 0.1, upper = 10)$root
  beta_new <- beta_old * gamma(1 + 1/alpha_old) / gamma(1 + 1/alpha_new)
  list(alpha_new = alpha_new, beta_new = beta_new)
}
adjusted_params <- adjust_weibull(alpha_old, beta_old, mean_original, sd_new)
adjusted_params

quantile(fln.b, probs = 0.05)
quantile(fln.b, probs = 0.01)

gofstat(list(fn.b, fln.b, fll.b, fu.b, fg.b, fw.b),
        fitnames = c("Normal", "Lognormal", "Log-logistic", "uniform", "Gamma", "Weibull"))

# ---------- Canola ----------
fn.c  <- fitdist(data$detrended_canola, "norm",  method = "mle"); summary(fn.c)
fln.c <- fitdist(data$detrended_canola, "lnorm", method = "mle"); summary(fln.c)
fll.c <- fitdist(data$detrended_canola, "llogis", method = "mle"); summary(fll.c)
fu.c  <- fitdist(data$detrended_canola, "unif",  method = "mle"); summary(fu.c)
fg.c  <- fitdist(data$detrended_canola, "gamma", method = "mle"); summary(fg.c)
fw.c  <- fitdist(data$detrended_canola, "weibull", method = "mle"); summary(fw.c)

alpha_old <- 7.02
beta_old  <- 951.92
mean_original <- beta_old * gamma(1 + 1/alpha_old)
sd_original   <- sqrt(beta_old^2 * (gamma(1 + 2/alpha_old) - (gamma(1 + 1/alpha_old))^2))
sd_new <- 175.05

weibull_sd <- function(alpha_new, beta_new) {
  gamma1 <- gamma(1 + 1/alpha_new)
  gamma2 <- gamma(1 + 2/alpha_new)
  sqrt(beta_new^2 * (gamma2 - gamma1^2))
}
adjust_weibull <- function(alpha_old, beta_old, mean_original, sd_new) {
  alpha_new <- uniroot(function(alpha_new) {
    beta_new <- beta_old * gamma(1 + 1/alpha_old) / gamma(1 + 1/alpha_new)
    weibull_sd(alpha_new, beta_new) - sd_new
  }, lower = 0.1, upper = 10)$root
  beta_new <- beta_old * gamma(1 + 1/alpha_old) / gamma(1 + 1/alpha_new)
  list(alpha_new = alpha_new, beta_new = beta_new)
}
adjusted_params <- adjust_weibull(alpha_old, beta_old, mean_original, sd_new)
adjusted_params

quantile(fln.c, probs = 0.05)
quantile(fln.c, probs = 0.01)

gofstat(list(fn.c, fln.c, fll.c, fu.c, fg.c, fw.c),
        fitnames = c("Normal", "Lognormal", "Log-logistic", "uniform", "Gamma", "Weibull"))

# --------- Diagnostic plots (kept) ----------
plot.legend <- c("Normal", "Lognormal", "Log-logistic", "uniform", "Gamma", "Weibull")

par(mfrow = c(1, 1))
qqcomp(list(fn.sw, fln.sw, fll.sw, fu.sw, fg.sw, fw.sw), legendtext = plot.legend, main = "Spring Wheat")
qqcomp(list(fn.o,  fln.o,  fll.o,  fu.o,  fg.o,  fw.o),  legendtext = plot.legend, main = "Oats")
qqcomp(list(fn.b,  fln.b,  fll.b,  fu.b,  fg.b,  fw.b),  legendtext = plot.legend, main = "Barley")
qqcomp(list(fn.c,  fln.c,  fll.c,  fu.c,  fg.c,  fw.c),  legendtext = plot.legend, main = "Canola")
par(mfrow = c(1, 1))

qqcomp(list(fn.sw, fln.sw, fll.sw, fu.sw, fg.sw, fw.sw),
       legendtext = plot.legend, main = "Spring Wheat", fitcol = 1:6)
qqcomp(list(fn.o, fln.o, fll.o, fu.o, fg.o, fw.o),
       legendtext = plot.legend, main = "Oats", fitcol = 1:6)
qqcomp(list(fn.b, fln.b, fll.b, fu.b, fg.b, fw.b),
       legendtext = plot.legend, main = "Barley", fitcol = 1:6)
qqcomp(list(fn.c, fln.c, fll.c, fu.c, fg.c, fw.c),
       legendtext = plot.legend, main = "Canola", fitcol = 1:6)

cdfcomp(list(fn.sw, fln.sw, fll.sw, fu.sw, fg.sw, fw.sw),
        legendtext = plot.legend, main = "Spring Wheat")
cdfcomp(list(fn.o, fln.o, fll.o, fu.o, fg.o, fw.o),
        legendtext = plot.legend, main = "Oats")
cdfcomp(list(fn.b, fln.b, fll.b, fu.b, fg.b, fw.b),
        legendtext = plot.legend, main = "Barley")
cdfcomp(list(fn.c, fln.c, fll.c, fu.c, fg.c, fw.c),
        legendtext = plot.legend, main = "Canola")

# =====================================================
# PART C: Copula-based simulation with Weibull margins
# =====================================================

# 1) Simulation sizes 
n_years <- 20
n_iter  <- 10000
n_sim   <- n_years * n_iter  # 20,000

# 2) Weibull parameters per crop (shape α, scale β) 
shape_params <- c(7.91, 6.07, 6.02, 5.91)
scale_params <- c(1868.78, 1985.56, 1793.74, 960.81)

# 3) Hard caps for yields 
max_yields <- c(2390.52, 2573.89, 2342.07, 1248.06)

# 4) Target correlation matrix for the Gaussian copula 
cor_matrix <- matrix(c(
  1,      0.653,  0.665, 0.543,
  0.653,  1,      0.654, 0.503,
  0.665, 0.654, 1,      0.547,
  0.543, 0.503, 0.547, 1
), nrow = 4)

# 5) Generate correlated standard normals (seed kept for identical results)
set.seed(123)
norm_copula_input <- mvrnorm(n = n_sim, mu = rep(0, 4), Sigma = cor_matrix)

# 6) Map to uniforms via Φ
uniform_copula_output <- pnorm(norm_copula_input)

# 7) Map uniforms to Weibull via inverse CDF (qweibull), margin-wise
weibull_sim <- mapply(function(u, shape, scale) {
  qweibull(u, shape = shape, scale = scale)
}, as.data.frame(uniform_copula_output), shape_params, scale_params)
weibull_sim <- as.matrix(weibull_sim)

# 8) Apply maximum yield caps (elementwise min)
weibull_sim_restricted <- sweep(weibull_sim, 2, max_yields, pmin)

# 9) Build final dataframe with Year & Iteration
simulated_yield <- as.data.frame(weibull_sim_restricted)
colnames(simulated_yield) <- c("swheat", "oats", "barley", "canola")
simulated_yield$Year <- rep(1:n_years, each = n_iter)
simulated_yield$Iteration <- rep(1:n_iter, times = n_years)
simulated_yield <- simulated_yield[, c("Year", "Iteration", "swheat", "oats", "barley", "canola")]

# 10) Diagnostics 
head(simulated_yield)
cat("Correlation (before restriction):\n"); print(round(cor(weibull_sim), 3))
cat("Correlation (after restriction):\n");  print(round(cor(simulated_yield[, 3:6]), 3))
summary(simulated_yield)
sapply(simulated_yield, sd, na.rm = TRUE)
boxplot(simulated_yield[, 3:6], main = "Boxplots of Restricted Simulated Yields", ylab = "Yield")

# 11) Save simulated yields 
write.csv(
  simulated_yield,
  "data/simulated_yield.csv",
  row.names = FALSE
)
