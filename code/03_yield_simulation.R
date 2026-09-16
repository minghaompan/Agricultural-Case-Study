# ============================================================
# 03_yield_simulation.R
# Detrending, distribution fitting, and copula yield simulation
# ============================================================

# ---- 1. Packages and seed --------------------------------------------------

library(fitdistrplus)
library(actuar)
library(MASS)

set.seed(123)


# ---- 2. Load data ----------------------------------------------------------

yield_data <- read.csv(
  "data/Crop Yield in AthbascaAB.csv",
  check.names = FALSE
)

names(yield_data) <- tolower(names(yield_data))

crops <- c("swheat", "oats", "barley", "canola")
crop_labels <- c("Spring wheat", "Oats", "Barley", "Canola")


# ---- 3. Detrend yields to 2024 --------------------------------------------

detrend_to_base_year <- function(
  data,
  crop,
  year_col = "cropyear",
  base_year = 2024
) {
  trend_model <- lm(
    reformulate(year_col, response = crop),
    data = data,
    na.action = na.exclude
  )

  trend <- predict(
    trend_model,
    newdata = data
  )

  base_data <- data.frame(base_year)
  names(base_data) <- year_col

  base_trend <- as.numeric(
    predict(
      trend_model,
      newdata = base_data
    )
  )

  detrended <- data[[crop]] - trend + base_trend

  list(
    model = trend_model,
    trend = trend,
    detrended = detrended
  )
}

detrended_results <- setNames(
  lapply(
    crops,
    function(crop) {
      detrend_to_base_year(
        data = yield_data,
        crop = crop,
        base_year = 2024
      )
    }
  ),
  crops
)

for (crop in crops) {
  yield_data[[paste0("trend_", crop)]] <-
    detrended_results[[crop]]$trend

  yield_data[[paste0("detrended_", crop)]] <-
    detrended_results[[crop]]$detrended
}

plot_yield_trend <- function(crop, crop_label) {
  plot(
    yield_data$cropyear,
    yield_data[[crop]],
    type = "l",
    main = paste0(
      crop_label,
      ": Original and detrended yield"
    ),
    xlab = "Year",
    ylab = "Yield"
  )

  lines(
    yield_data$cropyear,
    yield_data[[paste0("trend_", crop)]],
    lty = 2
  )

  lines(
    yield_data$cropyear,
    yield_data[[paste0("detrended_", crop)]]
  )

  abline(v = 2024, lty = 2)

  legend(
    "bottomright",
    legend = c("Original", "Trend", "Detrended"),
    lty = c(1, 2, 1),
    bty = "n"
  )
}

invisible(
  Map(
    plot_yield_trend,
    crops,
    crop_labels
  )
)

write.csv(
  yield_data,
  "data/de-trended_crop_yield_in_athabasca.csv",
  row.names = FALSE
)


# ---- 4. Candidate distributions -------------------------------------------

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

yield_fits <- setNames(
  lapply(
    crops,
    function(crop) {
      fit_candidates(
        yield_data[[paste0("detrended_", crop)]]
      )
    }
  ),
  crops
)

yield_gof <- setNames(
  lapply(
    crops,
    function(crop) {
      gofstat(
        yield_fits[[crop]],
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
  print(yield_gof[[crop]])
}

for (i in seq_along(crops)) {
  qqcomp(
    yield_fits[[crops[i]]],
    legendtext = distribution_names,
    main = crop_labels[i]
  )
}

for (i in seq_along(crops)) {
  cdfcomp(
    yield_fits[[crops[i]]],
    legendtext = distribution_names,
    main = crop_labels[i]
  )
}


# ---- 5. Weibull parameters -------------------------------------------------

weibull_sd <- function(shape, scale) {
  gamma_1 <- gamma(1 + 1 / shape)
  gamma_2 <- gamma(1 + 2 / shape)

  scale * sqrt(gamma_2 - gamma_1^2)
}

adjust_weibull <- function(
  shape_old,
  scale_old,
  target_sd
) {
  mean_original <-
    scale_old *
    gamma(1 + 1 / shape_old)

  shape_new <- uniroot(
    function(shape) {
      scale_new <-
        mean_original /
        gamma(1 + 1 / shape)

      weibull_sd(
        shape,
        scale_new
      ) -
        target_sd
    },
    lower = 0.1,
    upper = 10
  )$root

  scale_new <-
    mean_original /
    gamma(1 + 1 / shape_new)

  c(
    shape = shape_new,
    scale = scale_new
  )
}

weibull_parameter_data <- data.frame(
  crop = crops,
  shape_old = c(9.88, 7.37, 7.75, 7.02),
  scale_old = c(1849.79, 1965.27, 1770.16, 951.92),
  target_sd = c(263.86, 353.27, 321.38, 175.05)
)

adjusted_weibull <- t(
  mapply(
    adjust_weibull,
    weibull_parameter_data$shape_old,
    weibull_parameter_data$scale_old,
    weibull_parameter_data$target_sd
  )
)

weibull_parameter_data$shape <-
  round(adjusted_weibull[, "shape"], 2)

weibull_parameter_data$scale <-
  round(adjusted_weibull[, "scale"], 2)

print(weibull_parameter_data)


# ---- 6. Gaussian copula simulation ----------------------------------------

n_years <- 20L
n_iter <- 10000L
n_sim <- n_years * n_iter

shape_params <- weibull_parameter_data$shape
scale_params <- weibull_parameter_data$scale

max_yields <- c(
  2390.52,
  2573.89,
  2342.07,
  1248.06
)

cor_matrix <- matrix(
  c(
    1.000, 0.653, 0.665, 0.543,
    0.653, 1.000, 0.654, 0.503,
    0.665, 0.654, 1.000, 0.547,
    0.543, 0.503, 0.547, 1.000
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(crops, crops)
)

normal_draws <- mvrnorm(
  n = n_sim,
  mu = rep(0, length(crops)),
  Sigma = cor_matrix
)

uniform_draws <- pnorm(normal_draws)

weibull_draws <- mapply(
  function(u, shape, scale) {
    qweibull(
      u,
      shape = shape,
      scale = scale
    )
  },
  as.data.frame(uniform_draws),
  shape_params,
  scale_params
)

weibull_draws <- as.matrix(weibull_draws)

restricted_yields <- sweep(
  weibull_draws,
  2,
  max_yields,
  pmin
)

simulated_yield <- data.frame(
  year = rep(seq_len(n_years), each = n_iter),
  iteration = rep(seq_len(n_iter), times = n_years),
  restricted_yields
)

names(simulated_yield)[3:6] <- crops


# ---- 7. Save simulation ----------------------------------------------------

write.csv(
  simulated_yield,
  "data/simulated_yield.csv",
  row.names = FALSE
)


# ---- 8. Results ------------------------------------------------------------

cat("\nCorrelation before yield caps\n")
print(
  round(
    cor(weibull_draws),
    3
  )
)

cat("\nCorrelation after yield caps\n")
print(
  round(
    cor(simulated_yield[crops]),
    3
  )
)

cat("\nSimulated yield summary\n")
print(summary(simulated_yield[crops]))

cat("\nSimulated yield standard deviations\n")
print(
  sapply(
    simulated_yield[crops],
    sd,
    na.rm = TRUE
  )
)
