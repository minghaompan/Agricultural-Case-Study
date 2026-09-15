# ============================================================
# Farm NPV Monte Carlo — BASELINE ONLY (4/3/2-crop rotations)
# Keeps NPV mechanics:
#   - Crop insurance (PL + VPB) exactly as in your script
#   - AgriStability exactly as in your script
#   - Drainage deduction (shared draw across rotations)
# Removes:
#   - All plots / ggplot code
# Assumptions (from your earlier request):
#   - 20-year horizon only (years 1–20)
#   - discount rate = 0 (no discounting)
#   - NO terminal value
# Outputs:
#   - iteration-level NPV (CAD/ha) after drainage deduction
#   - summary table by rotation
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(readr)
})

# --------------------------
# Settings
# --------------------------
HORIZON_YEARS <- 20
DISCOUNT_RATE <- 0.1

price_path <- "data/simulated_price.csv"
yield_path <- "data/simulated_yield.csv"
cost_path  <- "data/simulated_inputcost.csv"

# --------------------------
# Load
# --------------------------
sim_price <- read.csv(price_path) |> rename_with(tolower)
sim_yield <- read.csv(yield_path) |> rename_with(tolower)
sim_cost  <- read.csv(cost_path)  |> rename_with(tolower)

# Common drainage draw (shared across all rotations)
set.seed(123)
all_iterations <- sort(unique(sim_price$iteration))
stopifnot(length(all_iterations) == 10000)

drainage_draw <- tibble(
  iteration = all_iterations,
  drain_cost_per_ha = runif(length(all_iterations), min = 494, max = 1600)
)

# --------------------------
# Core runner: BASELINE NPV only (no plots)
# --------------------------
run_rotation_baseline_npv <- function(
    rot_name,
    crops,                  # crop names in CSV columns
    acres,                  # named acres per crop, sums to total_acres
    raay_year1,             # named t/ha for year 1 (as in your script)
    sip_base,               # named SIP
    coverage_level = 0.80,
    machinery_fix_per_acre = 73.04,
    total_acres = sum(acres),
    include_crop_insurance = TRUE,
    include_agristability  = TRUE
) {
  stopifnot(setequal(names(acres), crops))
  stopifnot(setequal(names(raay_year1), crops))
  stopifnot(setequal(names(sip_base), crops))
  
  prices_long <- sim_price |>
    pivot_longer(-c(year, iteration), names_to = "crop", values_to = "price") |>
    filter(crop %in% crops)
  
  yields_long <- sim_yield |>
    pivot_longer(-c(year, iteration), names_to = "crop", values_to = "yield") |>
    filter(crop %in% crops)
  
  costs_long <- sim_cost |>
    pivot_longer(-c(year, iteration), names_to = "crop", values_to = "cost") |>
    filter(crop %in% crops)
  
  # Join + baseline transforms
  data <- yields_long %>%
    left_join(prices_long, by = c("year", "iteration", "crop")) %>%
    left_join(costs_long,  by = c("year", "iteration", "crop")) %>%
    mutate(
      yield      = yield / 1000,                 # kg -> tonnes
      cost       = cost + machinery_fix_per_acre,
      crop_acres = acres[crop]
    ) %>%
    filter(year >= 1, year <= HORIZON_YEARS)
  
  # Crop insurance (PL + VPB) — identical logic to your script
  if (include_crop_insurance) {
    data <- data %>%
      arrange(iteration, crop, year) %>%
      group_by(iteration, crop) %>%
      mutate(
        raay          = if_else(year == 1, unname(raay_year1[crop]), cummean(yield)),
        insured_yield = raay * coverage_level,
        sip           = sip_base[crop],
        pl_payment    = if_else(yield < insured_yield, (insured_yield - yield) * sip, 0),
        vpb_payment   = case_when(
          yield < insured_yield & price > 1.5 * sip ~ 0.5 * sip * (insured_yield - yield),
          yield < insured_yield & price >= 1.1 * sip & price <= 1.5 * sip ~ (insured_yield - yield) * (price - sip),
          TRUE ~ 0
        )
      ) %>%
      ungroup()
  } else {
    data <- data %>% mutate(pl_payment = 0, vpb_payment = 0)
  }
  
  # Crop-level cash flows
  data <- data %>%
    mutate(
      crop_revenue      = yield * price * crop_acres,
      insurance_total   = (pl_payment + vpb_payment) * crop_acres,
      allowable_income  = if (include_crop_insurance) crop_revenue + insurance_total else crop_revenue,
      allowable_expense = cost * crop_acres,
      net_crop_cf       = allowable_income - allowable_expense
    )
  
  # Whole-farm production margin (for AgriStability)
  farm_margin <- data %>%
    group_by(iteration, year) %>%
    summarise(
      total_income  = sum(allowable_income),
      total_expense = sum(allowable_expense),
      .groups = "drop"
    ) %>%
    mutate(production_margin = total_income - total_expense)
  
  # AgriStability — identical logic to your script
  if (include_agristability) {
    farm_margin <- farm_margin %>%
      group_by(iteration) %>%
      arrange(year) %>%
      mutate(
        ref_margin_5 = zoo::rollapplyr(
          production_margin, 5,
          function(x) mean(x[order(x)][2:4]),
          fill = NA, partial = FALSE
        ),
        ref_margin_3 = zoo::rollapplyr(
          production_margin, 3,
          mean, fill = NA, partial = FALSE
        ),
        ref_margin = ifelse(
          !is.na(ref_margin_5), ref_margin_5,
          ifelse(row_number() >= 3, ref_margin_3, NA)
        )
      ) %>%
      ungroup() %>%
      mutate(
        trigger_margin      = 0.7 * ref_margin,
        agstability_payment = ifelse(
          !is.na(ref_margin) & production_margin < trigger_margin,
          0.8 * (trigger_margin - production_margin), 0
        )
      )
  } else {
    farm_margin <- farm_margin %>% mutate(agstability_payment = 0)
  }
  
  # Farm cash flow (years 1..20) and NPV (r=0, no terminal value)
  farm_cashflow <- data %>%
    group_by(iteration, year) %>%
    summarise(net_crop_cashflow = sum(net_crop_cf), .groups = "drop") %>%
    left_join(
      farm_margin %>% select(iteration, year, agstability_payment),
      by = c("iteration", "year")
    ) %>%
    mutate(total_cashflow = net_crop_cashflow + agstability_payment)
  
  npv_iter <- farm_cashflow %>%
    mutate(discounted_cf = total_cashflow / ((1 + DISCOUNT_RATE) ^ year)) %>%
    group_by(iteration) %>%
    summarise(NPV_total = sum(discounted_cf), .groups = "drop") %>%
    mutate(NPV_per_acre = NPV_total / total_acres) %>%
    left_join(drainage_draw, by = "iteration") %>%
    mutate(
      drainage_total   = drain_cost_per_ha * (total_acres / 2.47), # ha = acres/2.47
      NPV_total_adj    = NPV_total - drainage_total,
      NPV_per_acre_adj = NPV_total_adj / total_acres,
      NPV_per_ha_adj   = NPV_per_acre_adj * 2.47,
      rotation         = rot_name
    ) %>%
    select(rotation, iteration, NPV_per_ha_adj, NPV_total_adj)
  
  npv_sum <- npv_iter %>%
    summarise(
      rotation = first(rotation),
      mean_ha  = mean(NPV_per_ha_adj, na.rm = TRUE),
      sd_ha    = sd(NPV_per_ha_adj, na.rm = TRUE),
      p05      = quantile(NPV_per_ha_adj, 0.05, na.rm = TRUE),
      p50      = quantile(NPV_per_ha_adj, 0.50, na.rm = TRUE),
      p95      = quantile(NPV_per_ha_adj, 0.95, na.rm = TRUE),
      p_neg    = mean(NPV_per_ha_adj < 0, na.rm = TRUE)
    )
  
  list(npv_iter = npv_iter, npv_sum = npv_sum)
}

# --------------------------
# Rotations (same as your script)
# --------------------------
total_acres <- 1920

# 4-crop: swheat 27%, canola 41%, barley 18%, oats remainder
acres_4 <- c(
  swheat = round(total_acres * 0.27),
  canola = round(total_acres * 0.41),
  barley = round(total_acres * 0.18)
)
acres_4 <- c(acres_4, oats = total_acres - sum(acres_4))
raay_4 <- c(swheat = 1.756, oats = 1.843, barley = 1.661, canola = 0.889)
sip_4  <- c(swheat = 285,   oats = 250,   barley = 245,   canola = 575)

# 3-crop: swheat 32%, canola 47%, barley 21%
acres_3 <- round(total_acres * c(swheat = 0.32, canola = 0.47, barley = 0.21))
raay_3  <- c(swheat = 1.756, barley = 1.661, canola = 0.889)
sip_3   <- c(swheat = 285,   barley = 245,   canola = 575)

# 2-crop: swheat 40%, canola 60%
acres_2 <- c(
  swheat = round(total_acres * 0.40),
  canola = total_acres - round(total_acres * 0.40)
)
raay_2 <- c(swheat = 1.756, canola = 0.889)
sip_2  <- c(swheat = 285,   canola = 575)

# --------------------------
# Run baseline NPVs (with crop insurance + AgriStability)
# --------------------------
res_4 <- run_rotation_baseline_npv(
  rot_name = "4-crop (Swheat–Canola–Barley–Oats)",
  crops = c("swheat", "canola", "barley", "oats"),
  acres = acres_4, raay_year1 = raay_4, sip_base = sip_4,
  include_crop_insurance = TRUE,
  include_agristability  = TRUE
)

res_3 <- run_rotation_baseline_npv(
  rot_name = "3-crop (Swheat–Canola–Barley)",
  crops = c("swheat", "barley", "canola"),
  acres = acres_3, raay_year1 = raay_3, sip_base = sip_3,
  include_crop_insurance = TRUE,
  include_agristability  = TRUE
)

res_2 <- run_rotation_baseline_npv(
  rot_name = "2-crop (Swheat–Canola)",
  crops = c("swheat", "canola"),
  acres = acres_2, raay_year1 = raay_2, sip_base = sip_2,
  include_crop_insurance = TRUE,
  include_agristability  = TRUE
)

# --------------------------
# Combine + export
# --------------------------
npv_iter_all <- bind_rows(res_4$npv_iter, res_3$npv_iter, res_2$npv_iter)
npv_sum_all  <- bind_rows(res_4$npv_sum,  res_3$npv_sum,  res_2$npv_sum)

print(npv_sum_all)

write_csv(npv_iter_all, "data/NPV_baseline_20yr_nodiscount_iterations.csv")
write_csv(npv_sum_all,  "data/NPV_baseline_20yr_nodiscount_summary.csv")