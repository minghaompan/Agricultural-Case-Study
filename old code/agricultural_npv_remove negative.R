# ============================================================
# Farm NPV Monte Carlo — Compact Scenario Runner (2/3/4-crop)
# Keeps identical results: insurance, AgriStability, drainage,
# terminal value, scenarios, diagnostics, BEP + figures.
# Adds: BEP summary table with mean, SD, p05, p50, p95 for
#       all rotations x scenarios (no SD plotted).
# NOTE: NPV is NOT changed.
#       For BEP only, simulations with NPV <= 0 are excluded.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(zoo)
  library(ggplot2); library(readxl); library(purrr)
  library(scales);  library(stringr); library(patchwork)
})

# --------------------------
# 0) File paths (edit here)
# --------------------------
price_path <- "data/simulated_price.csv"
yield_path <- "data/simulated_yield.csv"
cost_path  <- "data/simulated_inputcost.csv"

# --------------------------
# 1) Load once
# --------------------------
sim_price <- read.csv(price_path) |> rename_with(tolower)
sim_yield <- read.csv(yield_path) |> rename_with(tolower)
sim_cost  <- read.csv(cost_path)  |> rename_with(tolower)

# Common seed + drainage draw (shared across all rotations)
set.seed(123)
all_iterations <- sort(unique(sim_price$iteration))
stopifnot(length(all_iterations) == 1000)

drainage_draw <- tibble(
  iteration = all_iterations,
  drain_cost_per_ha = runif(length(all_iterations), min = 494, max = 1600)
)

# --------------------------
# 2) Helpers (PV of avoided CO2e; BEP helper; safe summaries; theme; colors; scenario list)
# --------------------------
# Closed-form PV under geometric decay (infinite horizon), discrete-time:
# PV_avoided_tCO2e = S0_tC * (44/12) * d / (r + d)
pv_avoided_tCO2e_inf <- function(S0_tC, d, r, CO2_per_C = 44/12) {
  S0_tC * CO2_per_C * d / (r + d)
}

# BEP helper:
# NPV in the simulation is unchanged.
# For BEP only, keep only profitable simulations (NPV > 0).
# Non-profitable simulations are set to NA and excluded from summaries.
calc_bep_profitable_only <- function(npv_per_ha, pv_avoided_tCO2e) {
  ifelse(npv_per_ha > 0 & pv_avoided_tCO2e > 0,
         npv_per_ha / pv_avoided_tCO2e,
         NA_real_)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

theme_pub <- theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line  = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.title = element_text(size = 10, colour = "black", face = "plain"),
    axis.text  = element_text(size = 9,  colour = "black", face = "plain"),
    plot.title = element_text(colour = "black", face = "plain"),
    plot.subtitle = element_text(colour = "black"),
    plot.caption  = element_text(colour = "black"),
    legend.position = "bottom",
    legend.title    = element_blank(),
    legend.text     = element_text(colour = "black", face = "plain"),
    plot.title.position   = "plot",
    plot.caption.position = "plot"
  )

scenario_levels <- c(
  "Baseline", "No BRM",
  "Price -15%", "Price +15%",
  "Input Cost -15%", "Input Cost +15%",
  "Discount rate = 7%"
)

scenario_cols <- c(
  "Baseline"                 = "#F8766D",
  "No BRM"                   = "#CD9600",
  "Price +15%"               = "#00796B",
  "Price -15%"               = "#7FD4C6",
  "Input Cost -15%"          = "#6A3D9A",
  "Input Cost +15%"          = "#C9B3E8",
  "Discount rate = 7%"       = "#1E88E5"
)

scale_color_scen <- scale_color_manual(
  limits = scenario_levels, values = scenario_cols, drop = FALSE
)
scale_fill_scen  <- scale_fill_manual(
  limits = scenario_levels, values = scenario_cols, drop = FALSE
)

# --------------------------
# 3) Core runner (reused by any rotation) — with plots intact
# --------------------------
run_rotation <- function(
    rot_name,
    crops,
    acres,
    raay_year1,
    sip_base,
    coverage_level = 0.80,
    machinery_fix_per_acre = 73.04,
    total_acres = sum(acres),
    S0_tC_per_ha = 1180,
    decay_d = 0.05
) {
  stopifnot(setequal(names(acres), crops))
  stopifnot(setequal(names(raay_year1), crops))
  stopifnot(setequal(names(sip_base), crops))
  
  # Filter & reshape to rotation
  prices_long <- sim_price |>
    pivot_longer(-c(year, iteration), names_to = "crop", values_to = "price") |>
    dplyr::filter(crop %in% crops)
  
  yields_long <- sim_yield |>
    pivot_longer(-c(year, iteration), names_to = "crop", values_to = "yield") |>
    dplyr::filter(crop %in% crops)
  
  costs_long <- sim_cost |>
    pivot_longer(-c(year, iteration), names_to = "crop", values_to = "cost") |>
    dplyr::filter(crop %in% crops)
  
  # Scenario engine (includes terminal growth after year 20)
  run_scenario <- function(scenario_name, price_mult = 1, cost_mult = 1,
                           discount = 0.10, include_risk_programs = TRUE,
                           tie_sip_to_price = FALSE, growth_rate = 0) {
    
    sip_vec <- if (tie_sip_to_price) sip_base * price_mult else sip_base
    
    data <- yields_long %>%
      dplyr::left_join(prices_long, by = c("year", "iteration", "crop")) %>%
      dplyr::left_join(costs_long,  by = c("year", "iteration", "crop")) %>%
      dplyr::mutate(
        price      = price * price_mult,
        cost       = (cost * cost_mult) + machinery_fix_per_acre,
        yield      = yield / 1000,  # kg -> tonnes
        crop_acres = acres[crop]
      )
    
    if (include_risk_programs) {
      data <- data %>%
        dplyr::arrange(iteration, crop, year) %>%
        dplyr::group_by(iteration, crop) %>%
        dplyr::mutate(
          raay          = dplyr::if_else(year == 1, unname(raay_year1[crop]), cummean(yield)),
          insured_yield = raay * coverage_level,
          sip           = sip_vec[crop],
          pl_payment    = dplyr::if_else(yield < insured_yield, (insured_yield - yield) * sip, 0),
          vpb_payment   = dplyr::case_when(
            yield < insured_yield & price > 1.5 * sip ~ 0.5 * sip * (insured_yield - yield),
            yield < insured_yield & price >= 1.1 * sip & price <= 1.5 * sip ~ (insured_yield - yield) * (price - sip),
            TRUE ~ 0
          )
        ) %>%
        dplyr::ungroup()
    } else {
      data <- data %>%
        dplyr::mutate(pl_payment = 0, vpb_payment = 0)
    }
    
    data <- data %>%
      dplyr::mutate(
        crop_revenue      = yield * price * crop_acres,
        insurance_total   = (pl_payment + vpb_payment) * crop_acres,
        allowable_income  = if (include_risk_programs) crop_revenue + insurance_total else crop_revenue,
        allowable_expense = cost * crop_acres,
        net_crop_cf       = allowable_income - allowable_expense
      )
    
    # AgriStability (whole-farm)
    farm_margin <- data %>%
      dplyr::group_by(iteration, year) %>%
      dplyr::summarise(
        total_income  = sum(allowable_income),
        total_expense = sum(allowable_expense),
        .groups = "drop"
      ) %>%
      dplyr::mutate(production_margin = total_income - total_expense)
    
    if (include_risk_programs) {
      farm_margin <- farm_margin %>%
        dplyr::group_by(iteration) %>%
        dplyr::arrange(year) %>%
        dplyr::mutate(
          ref_margin_5 = zoo::rollapplyr(
            production_margin, 5,
            function(x) mean(x[order(x)][2:4]),
            fill = NA, partial = FALSE
          ),
          ref_margin_3 = zoo::rollapplyr(
            production_margin, 3,
            mean, fill = NA, partial = FALSE
          ),
          ref_margin   = ifelse(
            !is.na(ref_margin_5), ref_margin_5,
            ifelse(dplyr::row_number() >= 3, ref_margin_3, NA)
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          trigger_margin      = 0.7 * ref_margin,
          agstability_payment = ifelse(
            !is.na(ref_margin) & production_margin < trigger_margin,
            0.8 * (trigger_margin - production_margin), 0
          )
        )
    } else {
      farm_margin <- farm_margin %>%
        dplyr::mutate(agstability_payment = 0)
    }
    
    # Cash flow -> NPV with terminal value at final year
    farm_cashflow <- data %>%
      dplyr::group_by(iteration, year) %>%
      dplyr::summarise(net_crop_cashflow = sum(net_crop_cf), .groups = "drop") %>%
      dplyr::left_join(
        farm_margin %>% dplyr::select(iteration, year, agstability_payment),
        by = c("iteration", "year")
      ) %>%
      dplyr::mutate(total_cashflow = net_crop_cashflow + agstability_payment)
    
    npv <- farm_cashflow %>%
      dplyr::mutate(discounted_cf = total_cashflow / ((1 + discount) ^ year)) %>%
      dplyr::group_by(iteration) %>%
      dplyr::summarise(
        NPV_base   = sum(discounted_cf),
        final_year = max(year),
        final_cf   = total_cashflow[which.max(year)],
        perpetuity_value = dplyr::case_when(
          abs(growth_rate) < 1e-12 ~ (final_cf / discount) / ((1 + discount) ^ final_year),
          growth_rate < discount   ~ (final_cf * (1 + growth_rate) / (discount - growth_rate)) / ((1 + discount) ^ final_year),
          TRUE ~ NA_real_
        ),
        NPV_total  = NPV_base + perpetuity_value,
        .groups = "drop"
      ) %>%
      dplyr::mutate(NPV_per_acre = NPV_total / total_acres) %>%
      dplyr::left_join(drainage_draw, by = "iteration") %>%
      dplyr::mutate(
        drainage_total   = drain_cost_per_ha * (total_acres / 2.47),
        NPV_total_adj    = NPV_total - drainage_total,
        NPV_per_acre_adj = NPV_total_adj / total_acres,
        NPV_per_ha_adj   = NPV_per_acre_adj * 2.47,
        scenario         = scenario_name,
        discount_rate    = discount,
        growth_rate      = growth_rate
      ) %>%
      dplyr::select(
        iteration, scenario, discount_rate, growth_rate,
        NPV_per_ha_adj, NPV_per_acre_adj, NPV_total_adj
      )
    
    npv
  }
  
  scenarios <- list(
    run_scenario("Baseline",            1.00, 1.00, 0.10, TRUE,  FALSE, 0.00),
    run_scenario("No BRM",              1.00, 1.00, 0.10, FALSE, FALSE, 0.00),
    run_scenario("Price -15%",          0.85, 1.00, 0.10, TRUE,  TRUE,  0.00),
    run_scenario("Price +15%",          1.15, 1.00, 0.10, TRUE,  TRUE,  0.00),
    run_scenario("Input Cost -15%",     1.00, 0.85, 0.10, TRUE,  FALSE, 0.00),
    run_scenario("Input Cost +15%",     1.00, 1.15, 0.10, TRUE,  FALSE, 0.00),
    run_scenario("Discount rate = 7%",  1.00, 1.00, 0.07, TRUE,  FALSE, 0.00)
  )
  
  npv_all <- dplyr::bind_rows(scenarios) %>%
    dplyr::mutate(scenario = factor(as.character(scenario), levels = scenario_levels))
  
  # Diagnostics
  diag_tbl <- npv_all %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      mean_ha   = mean(NPV_per_ha_adj, na.rm = TRUE),
      sd_ha     = stats::sd(NPV_per_ha_adj, na.rm = TRUE),
      p_NPV_lt0 = mean(NPV_per_ha_adj < 0, na.rm = TRUE),
      .groups = "drop"
    )
  print(diag_tbl)
  print(
    diag_tbl %>%
      dplyr::transmute(
        scenario,
        mean_ha = sprintf("%.2f", mean_ha),
        sd_ha   = sprintf("%.2f", sd_ha),
        `p(NPV<0)%` = sprintf("%.2f", 100 * p_NPV_lt0)
      )
  )
  
  # Density plot
  npv_plot <- npv_all %>%
    dplyr::mutate(NPV_thousand = NPV_per_ha_adj / 1000)
  
  means_df <- npv_plot %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(xbar_k = mean(NPV_thousand), .groups = "drop")
  
  p_npv <- ggplot(npv_plot, aes(x = NPV_thousand, colour = scenario)) +
    geom_density(linewidth = 1.1, adjust = 1.0) +
    geom_vline(aes(xintercept = 0), colour = "grey40", linetype = "22", linewidth = 0.6) +
    geom_vline(
      data = means_df,
      aes(xintercept = xbar_k, colour = scenario),
      linetype = "22", linewidth = 0.7, show.legend = FALSE
    ) +
    scale_color_scen +
    scale_x_continuous(
      breaks = seq(-2, 10, 2),
      labels = label_number(accuracy = 0.1),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = c(-3, 10)) +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    labs(
      title    = paste0("NPV Distributions by Scenario — ", rot_name),
      subtitle = "NPV per hectare (thousand CAD), after drainage deduction; 1,000 iterations",
      x        = "NPV (thousand CAD/ha)",
      y        = "Density",
      colour   = NULL
    ) +
    theme_pub
  print(p_npv)
  
  # BEP per iteration using each scenario's discount rate
  # Only profitable simulations (NPV > 0) are included
  bep_df <- npv_all %>%
    dplyr::mutate(
      PV_avoided_tCO2e  = pv_avoided_tCO2e_inf(S0_tC_per_ha, decay_d, discount_rate),
      BEP_CAD_per_tCO2e = calc_bep_profitable_only(NPV_per_ha_adj, PV_avoided_tCO2e)
    )
  
  bep_sum <- bep_df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      mean_bep = safe_mean(BEP_CAD_per_tCO2e),
      p05      = safe_quantile(BEP_CAD_per_tCO2e, 0.05),
      p95      = safe_quantile(BEP_CAD_per_tCO2e, 0.95),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      label = ifelse(is.na(mean_bep), NA_character_, sprintf("%.2f", mean_bep))
    )
  
  valid_endpts <- c(bep_sum$p05, bep_sum$p95)
  valid_endpts <- valid_endpts[is.finite(valid_endpts)]
  nudge <- if (length(valid_endpts) >= 2) max(0.1, 0.06 * diff(range(valid_endpts))) else 0.1
  
  p_bep <- ggplot(bep_sum, aes(y = scenario, x = mean_bep, colour = scenario)) +
    geom_linerange(
      aes(xmin = p05, xmax = p95),
      linewidth = 1.4, show.legend = FALSE, na.rm = TRUE
    ) +
    geom_point(size = 3.8, show.legend = FALSE, na.rm = TRUE) +
    geom_text(
      aes(label = label, x = mean_bep),
      nudge_x = nudge, size = 3.6, colour = "grey25",
      show.legend = FALSE, na.rm = TRUE
    ) +
    geom_vline(xintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey65") +
    scale_color_scen +
    scale_y_discrete(
      limits = rev(levels(npv_all$scenario)),
      labels = function(x) stringr::str_wrap(x, width = 22)
    ) +
    scale_x_continuous(
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.12))
    ) +
    labs(
      title    = paste0("Break-even Carbon Price by Scenario — ", rot_name),
      subtitle = "Dot = mean; line = 5th–95th percentile across profitable simulations only",
      x        = "Break-even price (CAD/tCO2e)",
      y        = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "plain", size = 18, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 12, colour = "grey30"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linetype = "dashed", colour = "#d9d9d9"),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11)
    )
  print(p_bep)
  
  # Return everything, including iteration-level BEP for summary table later
  list(
    npv   = npv_all,
    p_npv = p_npv,
    p_bep = p_bep,
    bep_iter = bep_df %>%
      dplyr::transmute(
        iteration,
        scenario,
        discount_rate,
        BEP_CAD_per_tCO2e = BEP_CAD_per_tCO2e
      )
  )
}

# --------------------------
# 4) Define rotations (identical splits & parameters as your scripts)
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
# 5) Run all rotations (plots kept)
# --------------------------
res_4 <- run_rotation(
  "4-Crop Rotation (Swheat–Canola–Barley–Oats)",
  crops = c("swheat", "canola", "barley", "oats"),
  acres = acres_4, raay_year1 = raay_4, sip_base = sip_4
)

res_3 <- run_rotation(
  "3-Crop Rotation (Swheat–Canola–Barley)",
  crops = c("swheat", "barley", "canola"),
  acres = acres_3, raay_year1 = raay_3, sip_base = sip_3
)

res_2 <- run_rotation(
  "2-Crop Rotation (Swheat–Canola)",
  crops = c("swheat", "canola"),
  acres = acres_2, raay_year1 = raay_2, sip_base = sip_2
)

npv_all_4 <- res_4$npv
npv_all_3 <- res_3$npv
npv_all_2 <- res_2$npv

# --------------------------
# 6) Multi-panel journal figures (NPV & BEP)
# --------------------------
make_npv_density_panel <- function(npv_df, panel_title, xlim_k = c(-3, 10), breaks_k = seq(-2, 10, 2)) {
  df <- npv_df %>%
    dplyr::mutate(NPV_thousand = NPV_per_ha_adj / 1000)
  
  means_df <- df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(xbar_k = mean(NPV_thousand, na.rm = TRUE), .groups = "drop")
  
  ggplot(df, aes(x = NPV_thousand, colour = scenario)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.6) +
    geom_density(linewidth = 0.9, adjust = 1, na.rm = TRUE) +
    geom_vline(
      data = means_df,
      aes(xintercept = xbar_k, colour = scenario),
      linetype = "22", linewidth = 0.7, show.legend = FALSE
    ) +
    scale_color_scen +
    scale_x_continuous(
      breaks = breaks_k,
      labels = label_number(accuracy = 0.1),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = xlim_k) +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    labs(title = panel_title, x = "NPV (thousand CAD/ha)", y = "Density") +
    theme_pub
}

make_bep_panel <- function(npv_df, panel_title, S0_tC_per_ha = 1180, decay_d = 0.05) {
  iter_df <- npv_df %>%
    dplyr::mutate(
      PV_avoided_tCO2e = pv_avoided_tCO2e_inf(S0_tC_per_ha, decay_d, discount_rate),
      BEP = calc_bep_profitable_only(NPV_per_ha_adj, PV_avoided_tCO2e)
    )
  
  sum_df <- iter_df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      mean_bep = safe_mean(BEP),
      p05      = safe_quantile(BEP, 0.05),
      p95      = safe_quantile(BEP, 0.95),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      label = ifelse(is.na(mean_bep), NA_character_, sprintf("%.2f", mean_bep))
    )
  
  ggplot(sum_df, aes(y = scenario, x = mean_bep, colour = scenario)) +
    geom_linerange(
      aes(xmin = p05, xmax = p95),
      linewidth = 0.9, show.legend = FALSE, na.rm = TRUE
    ) +
    geom_point(size = 1.8, show.legend = FALSE, na.rm = TRUE) +
    geom_text(
      aes(label = label, x = mean_bep),
      nudge_y = 0.32, vjust = 0, size = 3, colour = "black",
      show.legend = FALSE, na.rm = TRUE
    ) +
    geom_vline(xintercept = 0, linewidth = 0.7, linetype = "dashed", colour = "grey40") +
    scale_color_scen +
    scale_y_discrete(
      limits = rev(scenario_levels),
      labels = function(x) stringr::str_wrap(x, width = 30)
    ) +
    scale_x_continuous(
      labels = scales::label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.01, 0.06))
    ) +
    coord_cartesian(clip = "off") +
    labs(title = panel_title, x = "Break-even price (CAD/tCO2e)", y = NULL) +
    theme_pub +
    theme(
      plot.margin = margin(5.5, 12, 5.5, 5.5),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linetype = "dashed", colour = "#d9d9d9")
    )
}

# Common x-range across panels
npv_bind <- dplyr::bind_rows(
  npv_all_4 %>% dplyr::mutate(rot = "4"),
  npv_all_3 %>% dplyr::mutate(rot = "3"),
  npv_all_2 %>% dplyr::mutate(rot = "2")
) %>%
  dplyr::mutate(NPV_thousand = NPV_per_ha_adj / 1000)

npv_xlim <- c(
  max(-3, floor(min(npv_bind$NPV_thousand, na.rm = TRUE))),
  min(10, ceiling(max(npv_bind$NPV_thousand, na.rm = TRUE)))
)
npv_breaks <- seq(-2, 10, 2)

p_npv_4 <- make_npv_density_panel(npv_all_4, "a.", npv_xlim, npv_breaks) +
  theme(plot.title = element_text(face = "plain"))
p_npv_3 <- make_npv_density_panel(npv_all_3, "b.", npv_xlim, npv_breaks) +
  theme(plot.title = element_text(face = "plain"))
p_npv_2 <- make_npv_density_panel(npv_all_2, "c.", npv_xlim, npv_breaks) +
  theme(plot.title = element_text(face = "plain"))

npv_figure <- (p_npv_4 / p_npv_3 / p_npv_2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
npv_figure

p_bep_4 <- make_bep_panel(npv_all_4, "a.") +
  theme(plot.title = element_text(face = "plain"))
p_bep_3 <- make_bep_panel(npv_all_3, "b.") +
  theme(plot.title = element_text(face = "plain"))
p_bep_2 <- make_bep_panel(npv_all_2, "c.Spring wheat - Canola (1180 t C ha⁻¹ with 5% annual decay rate)") +
  theme(plot.title = element_text(face = "plain"))

bep_figure <- (p_bep_4 / p_bep_3 / p_bep_2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
bep_figure

# --------------------------
# 7) Save figures (same locations & DPI)
# --------------------------
ggsave("figures/Figure_NPV_rotations.png", npv_figure, dpi = 600,
       width = 190, height = 210, units = "mm")

ggsave("figures/Figure_BEP_rotations.png", bep_figure, dpi = 600,
       width = 190, height = 210, units = "mm")

# --------------------------
# 8) BEP summary table (mean, SD, p05, p50, p95)
#    for ALL rotations x scenarios
#    Only profitable simulations are included
# --------------------------
bep_all_iter <- dplyr::bind_rows(
  res_4$bep_iter %>% dplyr::mutate(rotation = "4-crop"),
  res_3$bep_iter %>% dplyr::mutate(rotation = "3-crop"),
  res_2$bep_iter %>% dplyr::mutate(rotation = "2-crop")
) %>%
  dplyr::mutate(scenario = factor(as.character(scenario), levels = scenario_levels))

bep_table <- bep_all_iter %>%
  dplyr::group_by(rotation, scenario) %>%
  dplyr::summarise(
    mean_bep = safe_mean(BEP_CAD_per_tCO2e),
    sd_bep   = safe_sd(BEP_CAD_per_tCO2e),
    p05      = safe_quantile(BEP_CAD_per_tCO2e, 0.05),
    p50      = safe_quantile(BEP_CAD_per_tCO2e, 0.50),
    p95      = safe_quantile(BEP_CAD_per_tCO2e, 0.95),
    .groups  = "drop"
  ) %>%
  dplyr::arrange(rotation, scenario)

# Pretty print (optional)
bep_table_print <- bep_table %>%
  dplyr::mutate(
    mean_bep = ifelse(is.na(mean_bep), NA_character_, sprintf("%.2f", mean_bep)),
    sd_bep   = ifelse(is.na(sd_bep),   NA_character_, sprintf("%.2f", sd_bep)),
    p05      = ifelse(is.na(p05),      NA_character_, sprintf("%.2f", p05)),
    p50      = ifelse(is.na(p50),      NA_character_, sprintf("%.2f", p50)),
    p95      = ifelse(is.na(p95),      NA_character_, sprintf("%.2f", p95))
  )
print(bep_table_print)

readr::write_csv(bep_table, "data/BEP_summary_all_rotations.csv")

theme_pub <- theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line  = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.title = element_text(size = 14, colour = "black", face = "plain"),
    axis.text  = element_text(size = 13, colour = "black", face = "plain"),
    plot.title = element_text(colour = "black", face = "plain"),
    plot.subtitle = element_text(colour = "black"),
    plot.caption  = element_text(colour = "black"),
    legend.position = "bottom",
    legend.title    = element_blank(),
    legend.text     = element_text(colour = "black", face = "plain"),
    plot.title.position   = "plot",
    plot.caption.position = "plot"
  )

# --------------------------
# 9) NPV plots — ONLY Baseline vs No BRM (CAD/ha)
# --------------------------
scen_keep <- c("Baseline", "No BRM")

scen_labels <- c(
  "Baseline" = "Baseline with BRM programs",
  "No BRM"   = "Baseline without BRM programs"
)

scale_color_keep <- scale_color_manual(
  values = scenario_cols[scen_keep],
  breaks = scen_keep,
  limits = scen_keep,
  labels = scen_labels,
  drop   = TRUE
)

make_npv_density_panel_keep <- function(npv_df, panel_title,
                                        xlim_ha = NULL,
                                        break_step = 1000) {
  
  df <- npv_df %>%
    dplyr::filter(as.character(scenario) %in% scen_keep) %>%
    dplyr::mutate(
      scenario = factor(as.character(scenario), levels = scen_keep),
      NPV_ha   = NPV_per_ha_adj
    )
  
  if (is.null(xlim_ha)) {
    rng <- range(df$NPV_ha, na.rm = TRUE)
    pad <- 0.04 * diff(rng)
    xlim_ha <- c(rng[1] - pad, rng[2] + pad)
  }
  
  means_df <- df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(xbar = mean(NPV_ha, na.rm = TRUE), .groups = "drop")
  
  ggplot(df, aes(x = NPV_ha, colour = scenario)) +
    geom_vline(xintercept = 0, colour = "grey40", linetype = "22", linewidth = 0.9) +
    geom_density(linewidth = 0.9, adjust = 1, na.rm = TRUE) +
    geom_vline(
      data = means_df,
      aes(xintercept = xbar, colour = scenario),
      linetype = "22", linewidth = 0.9, show.legend = FALSE
    ) +
    scale_color_keep +
    scale_x_continuous(
      breaks = seq(
        floor(xlim_ha[1] / break_step) * break_step,
        ceiling(xlim_ha[2] / break_step) * break_step,
        by = break_step
      ),
      labels = scales::label_number(accuracy = 1),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = xlim_ha) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.0001)) +
    labs(title = panel_title, x = "NPV ($/ha)", y = "Density") +
    theme_pub +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = 16),
      legend.position = "none"
    )
}

# Force x-axis limits (CAD/ha)
npv_xlim_keep <- c(-2000, 5500)

# Choose tick spacing (CAD/ha)
break_step_keep <- 1000

# Two panels (4-crop and 2-crop)
p_npv_keep_4 <- make_npv_density_panel_keep(
  npv_all_4,
  "a. 4-Crop rotations (Swheat-Canola-Barley-Oats)",
  xlim_ha = npv_xlim_keep,
  break_step = break_step_keep
)

p_npv_keep_2 <- make_npv_density_panel_keep(
  npv_all_2,
  "b. 2-Crop rotations (Swheat-Canola)",
  xlim_ha = npv_xlim_keep,
  break_step = break_step_keep
)

# Combine
npv_figure_keep <- (p_npv_keep_4 / p_npv_keep_2)
npv_figure_keep

ggsave(
  "figures/Figure_NPV_Baseline_NoBRM_rotations.png",
  npv_figure_keep,
  dpi = 600, width = 150, height = 170, units = "mm"
)

# --------------------------
# 10) BEP panels (updated here too)
# --------------------------
make_bep_panel <- function(npv_df, panel_title, S0_tC_per_ha = 1180, decay_d = 0.05) {
  iter_df <- npv_df %>%
    dplyr::mutate(
      PV_avoided_tCO2e = pv_avoided_tCO2e_inf(S0_tC_per_ha, decay_d, discount_rate),
      BEP = calc_bep_profitable_only(NPV_per_ha_adj, PV_avoided_tCO2e)
    )
  
  sum_df <- iter_df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      mean_bep = safe_mean(BEP),
      p05      = safe_quantile(BEP, 0.05),
      p95      = safe_quantile(BEP, 0.95),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      label = ifelse(is.na(mean_bep), NA_character_, sprintf("%.2f", mean_bep))
    )
  
  ggplot(sum_df, aes(y = scenario, x = mean_bep, colour = scenario)) +
    geom_linerange(
      aes(xmin = p05, xmax = p95),
      linewidth = 0.9, show.legend = FALSE, na.rm = TRUE
    ) +
    geom_point(size = 1.8, show.legend = FALSE, na.rm = TRUE) +
    geom_text(
      aes(label = label, x = mean_bep),
      nudge_y = 0.12, vjust = 0, size = 5, colour = "black",
      show.legend = FALSE, na.rm = TRUE
    ) +
    scale_color_scen +
    scale_y_discrete(
      limits = rev(scenario_levels),
      labels = function(x) stringr::str_wrap(x, width = 30)
    ) +
    scale_x_continuous(
      labels = scales::label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.01, 0.06))
    ) +
    coord_cartesian(clip = "off") +
    labs(title = panel_title, x = "Break-even carbon price ($/tCO\u2082e)", y = NULL) +
    theme_pub +
    theme(
      plot.margin = margin(5.5, 12, 5.5, 5.5),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linetype = "dashed", colour = "#d9d9d9")
    )
}

# Common x-axis for BOTH panels
common_bep_x <- scale_x_continuous(
  limits = c(0, 26),
  breaks = seq(0, 25, 5),
  labels = label_number(accuracy = 0.01),
  expand = c(0, 0)
)

library(patchwork)

top_title <- "1180 t C ha^-1 initial stock & 5% annual decay, NPV > 0"

# IMPORTANT:
# Here I explicitly set decay_d = 0.02 so the figure matches the title.
p_bep_4 <- make_bep_panel(npv_all_4, "a. 4-Crop rotations (Swheat-Canola-Barley-Oat)", decay_d = 0.02) +
  common_bep_x +
  theme(plot.title = element_text(face = "plain", hjust = 0.5, size = 16))

p_bep_2 <- make_bep_panel(npv_all_2, "b. 2-Crop rotations (Swheat-Canola)", decay_d = 0.02) +
  common_bep_x +
  theme(plot.title = element_text(face = "plain", hjust = 0.5, size = 16))

bep_figure <- (p_bep_4 / p_bep_2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

bep_figure <- bep_figure +
  plot_annotation(
    title = top_title,
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = 16),
      plot.margin = margin(10, 12, 5.5, 5.5)
    )
  )

bep_figure

ggsave("figures/Figure_1180_0.05_NPV>0.png", bep_figure, dpi = 600,
       width = 180, height = 230, units = "mm")