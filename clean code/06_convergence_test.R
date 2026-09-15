# ============================================================
# 06_convergence_test.R
# Monte Carlo convergence test for farm NPV
# ============================================================

# ---- 1. Packages -----------------------------------------------------------

library(dplyr)
library(tidyr)


# ---- 2. Run farm simulation ------------------------------------------------

source("updated code/05_farm_npv_simulation.R")


# ---- 3. Convergence settings -----------------------------------------------

sample_sizes <- c(
  1000L,
  5000L,
  10000L
)


# ---- 4. NPV convergence summary --------------------------------------------

summarise_convergence <- function(n) {
  npv_all %>%
    filter(iteration <= n) %>%
    group_by(rotation, scenario) %>%
    summarise(
      n_iterations = n,
      mean_npv = mean(npv_per_ha, na.rm = TRUE),
      sd_npv = sd(npv_per_ha, na.rm = TRUE),
      p05_npv = quantile(npv_per_ha, 0.05, na.rm = TRUE),
      median_npv = median(npv_per_ha, na.rm = TRUE),
      p95_npv = quantile(npv_per_ha, 0.95, na.rm = TRUE),
      p_npv_negative = mean(npv_per_ha < 0, na.rm = TRUE),
      mcse_mean_npv = sd_npv / sqrt(n),
      .groups = "drop"
    )
}

convergence_summary <- bind_rows(
  lapply(
    sample_sizes,
    summarise_convergence
  )
) %>%
  arrange(
    rotation,
    scenario,
    n_iterations
  )


# ---- 5. Change relative to 10,000 iterations -------------------------------

reference_10000 <- convergence_summary %>%
  filter(n_iterations == 10000L) %>%
  dplyr::select(
    rotation,
    scenario,
    mean_npv_10000 = mean_npv,
    sd_npv_10000 = sd_npv,
    p_npv_negative_10000 = p_npv_negative
  )

convergence_change <- convergence_summary %>%
  left_join(
    reference_10000,
    by = c("rotation", "scenario")
  ) %>%
  mutate(
    mean_npv_change =
      mean_npv - mean_npv_10000,
    mean_npv_pct_change =
      100 *
      (mean_npv - mean_npv_10000) /
      abs(mean_npv_10000),
    sd_npv_change =
      sd_npv - sd_npv_10000,
    sd_npv_pct_change =
      100 *
      (sd_npv - sd_npv_10000) /
      sd_npv_10000,
    p_negative_change_pp =
      100 *
      (p_npv_negative - p_npv_negative_10000)
  ) %>%
  dplyr::select(
    rotation,
    scenario,
    n_iterations,
    mean_npv,
    mean_npv_change,
    mean_npv_pct_change,
    sd_npv,
    sd_npv_change,
    sd_npv_pct_change,
    p05_npv,
    median_npv,
    p95_npv,
    p_npv_negative,
    p_negative_change_pp,
    mcse_mean_npv
  )


# ---- 6. Stability assessment -----------------------------------------------

npv_tolerance_pct <- 1.0
negative_npv_tolerance_pp <- 0.5

stability_assessment <- convergence_change %>%
  filter(n_iterations %in% c(1000L, 5000L)) %>%
  mutate(
    mean_npv_stable =
      abs(mean_npv_pct_change) <= npv_tolerance_pct,
    negative_npv_stable =
      abs(p_negative_change_pp) <= negative_npv_tolerance_pp,
    convergence =
      case_when(
        mean_npv_stable &
          negative_npv_stable ~ "Stable",
        TRUE ~ "Review"
      )
  ) %>%
  dplyr::select(
    rotation,
    scenario,
    n_iterations,
    mean_npv,
    mean_npv_pct_change,
    sd_npv,
    p05_npv,
    median_npv,
    p95_npv,
    p_npv_negative,
    p_negative_change_pp,
    mcse_mean_npv,
    convergence
  )

convergence_5000_vs_10000 <- stability_assessment %>%
  filter(n_iterations == 5000L)

convergence_overview <- convergence_5000_vs_10000 %>%
  group_by(rotation) %>%
  summarise(
    max_abs_mean_npv_change_pct =
      max(abs(mean_npv_pct_change), na.rm = TRUE),
    max_abs_negative_npv_change_pp =
      max(abs(p_negative_change_pp), na.rm = TRUE),
    scenarios_stable =
      sum(convergence == "Stable"),
    scenarios_total = n(),
    .groups = "drop"
  )


# ---- 7. Console summary ----------------------------------------------------

cat("\n================ NPV Convergence Summary ================\n")
print(convergence_summary)

cat("\n================ Change from 10,000 Iterations ==========\n")
print(convergence_change)

cat("\n================ 5,000 vs. 10,000 Iterations ============\n")
print(convergence_5000_vs_10000)

cat("\n================ Convergence Overview ===================\n")
print(convergence_overview)
