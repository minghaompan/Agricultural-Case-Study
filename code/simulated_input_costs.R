# ============================
# Packages
# ============================
library(vars)
library(tseries)
library(readxl)
library(urca)
library(fitdistrplus)
library(pastecs)
library(stats)
library(actuar)

# ==========================================
# PART A — Load data & fit marginal models
# ==========================================

# Import Excel
data <- read.csv("data/Crop Input Cost in Crop Alternative.csv")
colnames(data) <- tolower(colnames(data))

# Quick peek
head(data)

# --------------------------
# Spring wheat input costs
# --------------------------
fn.sw  <- fitdist(data$swheat, "norm",    method = "mle"); summary(fn.sw)
fln.sw <- fitdist(data$swheat, "lnorm",   method = "mle"); summary(fln.sw)
fll.sw <- fitdist(data$swheat, "llogis",  method = "mle"); summary(fll.sw)
fu.sw  <- fitdist(data$swheat, "unif",    method = "mle"); summary(fu.sw)
fg.sw  <- fitdist(data$swheat, "gamma",   method = "mle"); summary(fg.sw)
fw.sw  <- fitdist(data$swheat, "weibull", method = "mle"); summary(fw.sw)

gofstat(list(fn.sw, fln.sw, fll.sw, fu.sw, fg.sw, fw.sw),
        fitnames = c("Normal", "Lognormal", "Log-logistic", "uniform", "Gamma", "Weibull"))

# -------------
# Oats
# -------------
fn.o  <- fitdist(data$oats, "norm",    method = "mle"); summary(fn.o)
fln.o <- fitdist(data$oats, "lnorm",   method = "mle"); summary(fln.o)
fll.o <- fitdist(data$oats, "llogis",  method = "mle"); summary(fll.o)
fu.o  <- fitdist(data$oats, "unif",    method = "mle"); summary(fu.o)
fg.o  <- fitdist(data$oats, "gamma",   method = "mle"); summary(fg.o)
fw.o  <- fitdist(data$oats, "weibull", method = "mle"); summary(fw.o)

gofstat(list(fn.o, fln.o, fll.o, fu.o, fg.o, fw.o),
        fitnames = c("Normal", "Lognormal", "Log-logistic", "uniform", "Gamma", "Weibull"))

# -------------
# Barley
# -------------
fn.b  <- fitdist(data$barley, "norm",    method = "mle"); summary(fn.b)
fln.b <- fitdist(data$barley, "lnorm",   method = "mle"); summary(fln.b)
fll.b <- fitdist(data$barley, "llogis",  method = "mle"); summary(fll.b)
fu.b  <- fitdist(data$barley, "unif",    method = "mle"); summary(fu.b)
fg.b  <- fitdist(data$barley, "gamma",   method = "mle"); summary(fg.b)
fw.b  <- fitdist(data$barley, "weibull", method = "mle"); summary(fw.b)

gofstat(list(fn.b, fln.b, fll.b, fu.b, fg.b, fw.b),
        fitnames = c("Normal", "Lognormal", "Log-logistic", "uniform", "Gamma", "Weibull"))

# -------------
# Canola
# -------------
fn.c  <- fitdist(data$canola, "norm",    method = "mle"); summary(fn.c)
fln.c <- fitdist(data$canola, "lnorm",   method = "mle"); summary(fln.c)
fll.c <- fitdist(data$canola, "llogis",  method = "mle"); summary(fll.c)
fu.c  <- fitdist(data$canola, "unif",    method = "mle"); summary(fu.c)
fg.c  <- fitdist(data$canola, "gamma",   method = "mle"); summary(fg.c)
fw.c  <- fitdist(data$canola, "weibull", method = "mle"); summary(fw.c)

gofstat(list(fn.c, fln.c, fll.c, fu.c, fg.c, fw.c),
        fitnames = c("Normal", "Lognormal", "Log-logistic", "uniform", "Gamma", "Weibull"))

# ------------------------
# Diagnostics: QQ & CDFs
# ------------------------
plot.legend <- c("Normal", "Lognormal", "Log-logistic", "uniform", "Gamma", "Weibull")

par(mfrow = c(2, 2))
qqcomp(list(fn.sw, fln.sw, fll.sw, fu.sw, fg.sw, fw.sw), legendtext = plot.legend, main = "Spring Wheat")
qqcomp(list(fn.o,  fln.o,  fll.o,  fu.o,  fg.o,  fw.o), legendtext = plot.legend, main = "Oats")
qqcomp(list(fn.b,  fln.b,  fll.b,  fu.b,  fg.b,  fw.b), legendtext = plot.legend, main = "Barley")
qqcomp(list(fn.c,  fln.c,  fll.c,  fu.c,  fg.c,  fw.c), legendtext = plot.legend, main = "Canola")
par(mfrow = c(1, 1))

# (Kept intentionally; duplicates your original extra QQ calls)
qqcomp(list(fn.sw, fln.sw, fll.sw, fu.sw, fg.sw, fw.sw),
       legendtext = c("Normal", "Lognormal", "Log-logistic", "Uniform", "Gamma", "Weibull"),
       main = "Spring Wheat", fitcol = 1:6)
qqcomp(list(fn.o, fln.o, fll.o, fu.o, fg.o, fw.o),
       legendtext = c("Normal", "Lognormal", "Log-logistic", "Uniform", "Gamma", "Weibull"),
       main = "Oats", fitcol = 1:6)
qqcomp(list(fn.b, fln.b, fll.b, fu.b, fg.b, fw.b),
       legendtext = c("Normal", "Lognormal", "Log-logistic", "Uniform", "Gamma", "Weibull"),
       main = "Barley", fitcol = 1:6)
qqcomp(list(fn.c, fln.c, fll.c, fu.c, fg.c, fw.c),
       legendtext = c("Normal", "Lognormal", "Log-logistic", "Uniform", "Gamma", "Weibull"),
       main = "Canola", fitcol = 1:6)

cdfcomp(list(fn.sw, fln.sw, fll.sw, fu.sw, fg.sw, fw.sw), legendtext = plot.legend, main = "Spring Wheat")
cdfcomp(list(fn.o,  fln.o,  fll.o,  fu.o,  fg.o,  fw.o), legendtext = plot.legend, main = "Oats")
cdfcomp(list(fn.b,  fln.b,  fll.b,  fu.b,  fg.b,  fw.b), legendtext = plot.legend, main = "Barley")
cdfcomp(list(fn.c,  fln.c,  fll.c,  fu.c,  fg.c,  fw.c), legendtext = plot.legend, main = "Canola")


# ===========================================================
# PART B — Truncated log-logistic simulation of input costs
# ===========================================================

# Required for truncated sampling & data handling
library(truncdist)
library(dplyr)
library(tidyr)
library(writexl)

# --- Custom sampler: truncated log-logistic ---
rtrunc_llogis <- function(n, shape, scale, lower, upper) {
  rtrunc(n, spec = "llogis", a = lower, b = upper, shape = shape, scale = scale)
}

# --- Simulation setup ---
set.seed(123)   # critical for result reproducibility
years      <- 20
iterations <- 1000
total_obs  <- years * iterations

crops <- c("swheat", "oats", "barley", "canola")

# Crop-specific truncated log-logistic parameters and bounds
shape_params <- c(14.629, 11.137, 13.290, 13.433)
scale_params <- c(333.285, 275.896, 296.279, 403.916)
lower_bounds <- c(244.825, 182.871, 205.433, 288.92)
upper_bounds <- c(457.355, 422.049, 419.667, 577.89)

# --- Output container  ---
simulated_inputcost <- data.frame(
  Year      = rep(1:years, each = iterations),
  Iteration = rep(1:iterations, times = years)
)

# --- Simulate crop input costs (per crop, same loop logic) ---
for (i in seq_along(crops)) {
  crop_name <- crops[i]
  shape     <- shape_params[i]
  scale     <- scale_params[i]
  lower     <- lower_bounds[i]
  upper     <- upper_bounds[i]
  
  # Draw from truncated log-logistic
  cost_samples <- rtrunc_llogis(total_obs, shape, scale, lower, upper)
  
  # Append to result frame
  simulated_inputcost[[crop_name]] <- cost_samples
}

# Quick checks 
head(simulated_inputcost)
summary(simulated_inputcost)
sapply(simulated_inputcost, sd, na.rm = TRUE)

# Save
# --- Save output ---
write.csv(
  simulated_inputcost,
  "data/simulated_inputcost.csv",
  row.names = FALSE
)