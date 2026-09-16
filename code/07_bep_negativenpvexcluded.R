# ============================================================
# 07_bep_negativenpvexcluded_scenario.R
# ============================================================

# ---- 1. Packages and seed --------------------------------------------------

library(dplyr)
library(tidyr)
library(zoo)

set.seed(123)


# ---- 2. Model settings -----------------------------------------------------

total_acres <- 1920
coverage_level <- 0.80
machinery_cost_per_acre <- 73.04
baseline_discount_rate <- 0.10

carbon_stock_tC_per_ha <- 394
decay_rate <- 0.02


# ---- 3. Load simulation inputs --------------------------------------------

add_model_year <- function(data) {
  data %>%
    arrange(iteration, year) %>%
    group_by(iteration) %>%
    mutate(model_year = row_number()) %>%
    ungroup()
}

sim_price <- read.csv("data/simulated_price.csv") %>%
  rename_with(tolower) %>%
  add_model_year()

sim_yield <- read.csv("data/simulated_yield.csv") %>%
  rename_with(tolower) %>%
  add_model_year()

sim_cost <- read.csv("data/simulated_inputcost.csv") %>%
  rename_with(tolower) %>%
  add_model_year()

iterations <- sort(unique(sim_price$iteration))

drainage_draw <- tibble(
  iteration = iterations,
  drain_cost_per_ha = runif(
    length(iterations),
    min = 559,
    max = 1810
  )
)


# ---- 4. Carbon accounting -------------------------------------------------

calc_pv_avoided_tCO2e <- function(
  carbon_stock_tC_per_ha,
  decay_rate,
  discount_rate
) {
  carbon_stock_tC_per_ha *
    (44 / 12) *
    decay_rate /
    (discount_rate + decay_rate)
}

calc_bep <- function(npv_per_ha, pv_avoided_tCO2e) {
  npv_per_ha / pv_avoided_tCO2e
}


# ---- 5. Baseline farm simulation ------------------------------------------

run_rotation <- function(
  rotation,
  crops,
  acres,
  raay_year1,
  sip_base
) {

  rotation_acres <- sum(acres)

  prices_long <- sim_price %>%
    dplyr::select(iteration, model_year, all_of(crops)) %>%
    pivot_longer(
      cols = all_of(crops),
      names_to = "crop",
      values_to = "price"
    )

  yields_long <- sim_yield %>%
    dplyr::select(iteration, model_year, all_of(crops)) %>%
    pivot_longer(
      cols = all_of(crops),
      names_to = "crop",
      values_to = "yield"
    )

  costs_long <- sim_cost %>%
    dplyr::select(iteration, model_year, all_of(crops)) %>%
    pivot_longer(
      cols = all_of(crops),
      names_to = "crop",
      values_to = "cost"
    )

  # Crop cash flows under baseline assumptions
  farm_data <- yields_long %>%
    left_join(
      prices_long,
      by = c("iteration", "model_year", "crop")
    ) %>%
    left_join(
      costs_long,
      by = c("iteration", "model_year", "crop")
    ) %>%
    mutate(
      cost = cost + machinery_cost_per_acre,
      yield = yield / 1000,
      crop_acres = acres[crop]
    ) %>%
    arrange(iteration, crop, model_year) %>%
    group_by(iteration, crop) %>%
    mutate(
      raay = cummean(yield),
      raay = replace(
        raay,
        model_year == 1,
        unname(raay_year1[crop][model_year == 1])
      ),
      insured_yield = raay * coverage_level,
      spring_insurance_price = sip_base[crop],
      yield_shortfall = pmax(insured_yield - yield, 0),
      production_loss_payment =
        yield_shortfall * spring_insurance_price,
      variable_price_payment =
        yield_shortfall *
        pmax(
          0,
          pmin(
            price - spring_insurance_price,
            0.5 * spring_insurance_price
          )
        ) *
        as.numeric(price >= 1.1 * spring_insurance_price)
    ) %>%
    ungroup() %>%
    mutate(
      crop_revenue = yield * price * crop_acres,
      insurance_payment =
        (production_loss_payment + variable_price_payment) *
        crop_acres,
      allowable_income = crop_revenue + insurance_payment,
      allowable_expense = cost * crop_acres,
      net_crop_cashflow = allowable_income - allowable_expense
    )

  # AgriStability under baseline assumptions
  farm_margin <- farm_data %>%
    group_by(iteration, model_year) %>%
    summarise(
      total_income = sum(allowable_income),
      total_expense = sum(allowable_expense),
      production_margin = total_income - total_expense,
      .groups = "drop"
    ) %>%
    arrange(iteration, model_year) %>%
    group_by(iteration) %>%
    mutate(
      reference_margin_5 = zoo::rollapplyr(
        production_margin,
        width = 5,
        FUN = function(x) mean(sort(x)[2:4]),
        fill = NA,
        partial = FALSE
      ),
      reference_margin_3 = zoo::rollapplyr(
        production_margin,
        width = 3,
        FUN = mean,
        fill = NA,
        partial = FALSE
      ),
      reference_margin_3 = replace(
        reference_margin_3,
        row_number() < 3,
        NA_real_
      ),
      reference_margin = coalesce(
        reference_margin_5,
        reference_margin_3
      ),
      trigger_margin = 0.7 * reference_margin,
      margin_shortfall = pmax(
        coalesce(trigger_margin - production_margin, 0),
        0
      ),
      agristability_payment = 0.8 * margin_shortfall
    ) %>%
    ungroup()

  # Annual farm cash flows
  farm_cashflow <- farm_data %>%
    group_by(iteration, model_year) %>%
    summarise(
      net_crop_cashflow = sum(net_crop_cashflow),
      .groups = "drop"
    ) %>%
    left_join(
      farm_margin %>%
        dplyr::select(
          iteration,
          model_year,
          agristability_payment
        ),
      by = c("iteration", "model_year")
    ) %>%
    mutate(
      total_cashflow =
        net_crop_cashflow + agristability_payment,
      discounted_cashflow =
        total_cashflow /
        (1 + baseline_discount_rate)^model_year
    )

  # Baseline NPV
  farm_cashflow %>%
    group_by(iteration) %>%
    summarise(
      npv_cashflows = sum(discounted_cashflow),
      final_year = max(model_year),
      final_cashflow =
        total_cashflow[which.max(model_year)],
      terminal_value =
        (final_cashflow / baseline_discount_rate) /
        (1 + baseline_discount_rate)^final_year,
      npv_total = npv_cashflows + terminal_value,
      .groups = "drop"
    ) %>%
    left_join(drainage_draw, by = "iteration") %>%
    mutate(
      drainage_total =
        drain_cost_per_ha * (rotation_acres / 2.47),
      npv_total_adj =
        npv_total - drainage_total,
      npv_per_ha =
        (npv_total_adj / rotation_acres) * 2.47,
      rotation = rotation,
      discount_rate = baseline_discount_rate,
      pv_avoided_tCO2e =
        calc_pv_avoided_tCO2e(
          carbon_stock_tC_per_ha,
          decay_rate,
          baseline_discount_rate
        ),
      bep_cad_per_tCO2e =
        calc_bep(
          npv_per_ha,
          pv_avoided_tCO2e
        )
    ) %>%
    dplyr::select(
      iteration,
      rotation,
      discount_rate,
      npv_per_ha,
      bep_cad_per_tCO2e
    )
}


# ---- 6. Crop rotations ----------------------------------------------------

acres_4 <- c(
  swheat = round(total_acres * 0.27),
  canola = round(total_acres * 0.41),
  barley = round(total_acres * 0.18)
)

acres_4 <- c(
  acres_4,
  oats = total_acres - sum(acres_4)
)

acres_2 <- c(
  swheat = round(total_acres * 0.40),
  canola = total_acres - round(total_acres * 0.40)
)

raay <- c(
  swheat = 1.756,
  oats = 1.843,
  barley = 1.661,
  canola = 0.889
)

sip <- c(
  swheat = 285,
  oats = 250,
  barley = 245,
  canola = 575
)


# ---- 7. Run baseline simulations ------------------------------------------

npv_4 <- run_rotation(
  rotation = "4-crop",
  crops = c("swheat", "canola", "barley", "oats"),
  acres = acres_4,
  raay_year1 = raay[c("swheat", "canola", "barley", "oats")],
  sip_base = sip[c("swheat", "canola", "barley", "oats")]
)

npv_2 <- run_rotation(
  rotation = "2-crop",
  crops = c("swheat", "canola"),
  acres = acres_2,
  raay_year1 = raay[c("swheat", "canola")],
  sip_base = sip[c("swheat", "canola")]
)

npv_all <- bind_rows(npv_4, npv_2) %>%
  mutate(
    rotation = factor(
      rotation,
      levels = c("4-crop", "2-crop")
    )
  )


# ---- 8. Baseline NPV summary ----------------------------------------------

npv_summary <- npv_all %>%
  group_by(rotation) %>%
  summarise(
    mean_npv = mean(npv_per_ha, na.rm = TRUE),
    sd_npv = sd(npv_per_ha, na.rm = TRUE),
    p05 = quantile(npv_per_ha, 0.05, na.rm = TRUE),
    median = median(npv_per_ha, na.rm = TRUE),
    p95 = quantile(npv_per_ha, 0.95, na.rm = TRUE),
    p_npv_negative = mean(npv_per_ha < 0, na.rm = TRUE),
    .groups = "drop"
  )

print(npv_summary)


# ---- 9. Baseline break-even carbon price ----------------------------------
# BEP is summarized only for simulations with non-negative NPV.

bep_all <- npv_all %>%
  filter(npv_per_ha >= 0)

bep_summary <- bep_all %>%
  group_by(rotation) %>%
  summarise(
    n_positive_npv = n(),
    mean_bep = mean(bep_cad_per_tCO2e, na.rm = TRUE),
    sd_bep = sd(bep_cad_per_tCO2e, na.rm = TRUE),
    p05 = quantile(bep_cad_per_tCO2e, 0.05, na.rm = TRUE),
    median = median(bep_cad_per_tCO2e, na.rm = TRUE),
    p95 = quantile(bep_cad_per_tCO2e, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

print(bep_summary)

