# ============================================================
# 03_farm_npv_simulation.R
# Stochastic farm NPV and break-even carbon price analysis
# Four-crop and spring wheat-canola rotations
# ============================================================

# ---- 1. Packages and seed --------------------------------------------------

library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)
library(scales)
library(patchwork)

set.seed(123)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)


# ---- 2. Model settings -----------------------------------------------------

total_acres <- 1920
coverage_level <- 0.80
machinery_cost_per_acre <- 73.04
carbon_stock_tC_per_ha <- 394
decay_rate <- 0.02

scenario_levels <- c(
  "Baseline",
  "No BRM",
  "Price -15%",
  "Price +15%",
  "Input Cost -15%",
  "Input Cost +15%",
  "Discount rate = 7%"
)

scenarios <- tibble(
  scenario = scenario_levels,
  price_mult = c(1.00, 1.00, 0.85, 1.15, 1.00, 1.00, 1.00),
  cost_mult = c(1.00, 1.00, 1.00, 1.00, 0.85, 1.15, 1.00),
  discount_rate = c(0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.07),
  brm_flag = c(1, 0, 1, 1, 1, 1, 1),
  sip_price_flag = c(0, 0, 1, 1, 0, 0, 0)
)


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
  pmax(npv_per_ha, 0) / pv_avoided_tCO2e
}


# ---- 5. Farm simulation ---------------------------------------------------

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

  run_scenario <- function(
    scenario_name,
    price_mult,
    cost_mult,
    discount_rate,
    brm_flag,
    sip_price_flag
  ) {

    sip_scenario <- sip_base *
      (1 + sip_price_flag * (price_mult - 1))

    # Crop cash flows
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
        price = price * price_mult,
        cost = cost * cost_mult + machinery_cost_per_acre,
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
        spring_insurance_price = sip_scenario[crop],
        yield_shortfall = pmax(insured_yield - yield, 0),
        production_loss_payment =
          brm_flag *
          yield_shortfall *
          spring_insurance_price,
        variable_price_payment =
          brm_flag *
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

    # AgriStability
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
        agristability_payment =
          brm_flag * 0.8 * margin_shortfall
      ) %>%
      ungroup()

    # Farm cash flows
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
          (1 + discount_rate)^model_year
      )

    # NPV and break-even carbon price
    farm_cashflow %>%
      group_by(iteration) %>%
      summarise(
        npv_cashflows = sum(discounted_cashflow),
        final_year = max(model_year),
        final_cashflow =
          total_cashflow[which.max(model_year)],
        terminal_value =
          (final_cashflow / discount_rate) /
          (1 + discount_rate)^final_year,
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
        scenario = scenario_name,
        discount_rate = discount_rate,
        pv_avoided_tCO2e =
          calc_pv_avoided_tCO2e(
            carbon_stock_tC_per_ha,
            decay_rate,
            discount_rate
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
        scenario,
        discount_rate,
        npv_per_ha,
        bep_cad_per_tCO2e
      )
  }

  scenario_results <- lapply(
    seq_len(nrow(scenarios)),
    function(i) {
      run_scenario(
        scenario_name = scenarios$scenario[i],
        price_mult = scenarios$price_mult[i],
        cost_mult = scenarios$cost_mult[i],
        discount_rate = scenarios$discount_rate[i],
        brm_flag = scenarios$brm_flag[i],
        sip_price_flag = scenarios$sip_price_flag[i]
      )
    }
  )

  bind_rows(scenario_results)
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


# ---- 7. Run simulations ---------------------------------------------------

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
    scenario = factor(
      scenario,
      levels = scenario_levels
    ),
    rotation = factor(
      rotation,
      levels = c("4-crop", "2-crop")
    )
  )


# ---- 8. Summary tables ----------------------------------------------------

npv_summary <- npv_all %>%
  group_by(rotation, scenario) %>%
  summarise(
    mean_npv = mean(npv_per_ha, na.rm = TRUE),
    sd_npv = sd(npv_per_ha, na.rm = TRUE),
    p05 = quantile(npv_per_ha, 0.05, na.rm = TRUE),
    median = median(npv_per_ha, na.rm = TRUE),
    p95 = quantile(npv_per_ha, 0.95, na.rm = TRUE),
    p_npv_negative = mean(npv_per_ha < 0, na.rm = TRUE),
    .groups = "drop"
  )

bep_summary <- npv_all %>%
  group_by(rotation, scenario) %>%
  summarise(
    mean_bep = mean(bep_cad_per_tCO2e, na.rm = TRUE),
    sd_bep = sd(bep_cad_per_tCO2e, na.rm = TRUE),
    p05 = quantile(bep_cad_per_tCO2e, 0.05, na.rm = TRUE),
    median = median(bep_cad_per_tCO2e, na.rm = TRUE),
    p95 = quantile(bep_cad_per_tCO2e, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

print(npv_summary)
print(bep_summary)

write.csv(
  npv_summary,
  "results/NPV_summary.csv",
  row.names = FALSE
)

write.csv(
  bep_summary,
  "results/BEP_summary.csv",
  row.names = FALSE
)


# ---- 9. Figure settings ---------------------------------------------------

scenario_colours <- c(
  "Baseline" = "#F8766D",
  "No BRM" = "#D89000",
  "Price -15%" = "#72C9B8",
  "Price +15%" = "#00897B",
  "Input Cost -15%" = "#6F4AA8",
  "Input Cost +15%" = "#C5A5E3",
  "Discount rate = 7%" = "#1E88E5"
)

theme_publication <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "plain",
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    axis.title = element_text(size = 11),
    axis.text = element_text(
      size = 10,
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.55,
      colour = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.45,
      colour = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    plot.margin = margin(
      6, 8, 5, 6,
      unit = "pt"
    )
  )


# ---- 10. NPV distributions ------------------------------------------------

npv_scenarios <- c(
  "Baseline",
  "No BRM"
)

make_npv_density_panel <- function(
  data,
  panel_title,
  legend_position
) {

  plot_data <- data %>%
    filter(scenario %in% npv_scenarios) %>%
    mutate(
      scenario = factor(
        as.character(scenario),
        levels = npv_scenarios
      )
    )

  mean_data <- plot_data %>%
    group_by(scenario) %>%
    summarise(
      mean_npv = mean(npv_per_ha, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(
    plot_data,
    aes(
      x = npv_per_ha,
      colour = scenario
    )
  ) +
    geom_density(
      linewidth = 1.05,
      adjust = 1,
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = 0,
      colour = "grey45",
      linetype = "dashed",
      linewidth = 0.65
    ) +
    geom_vline(
      data = mean_data,
      aes(
        xintercept = mean_npv,
        colour = scenario
      ),
      linetype = "dashed",
      linewidth = 0.70,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(
      values = scenario_colours[npv_scenarios],
      breaks = npv_scenarios
    ) +
    scale_x_continuous(
      breaks = seq(-2000, 5000, by = 1000),
      labels = label_number(
        accuracy = 1,
        big.mark = " "
      ),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = seq(0, 0.0004, by = 0.0001),
      labels = label_number(accuracy = 0.0001),
      expand = expansion(mult = c(0, 0.04))
    ) +
    coord_cartesian(
      xlim = c(-2000, 5500)
    ) +
    labs(
      title = panel_title,
      x = "NPV ($/ha)",
      y = "Density",
      colour = NULL
    ) +
    theme_publication +
    theme(
      panel.grid = element_blank(),
      legend.position = legend_position,
      legend.key.width = grid::unit(1.8, "cm")
    )
}

p_npv_A <- make_npv_density_panel(
  data = filter(npv_all, rotation == "4-crop"),
  panel_title = "Panel A",
  legend_position = "right"
)

p_npv_B <- make_npv_density_panel(
  data = filter(npv_all, rotation == "2-crop"),
  panel_title = "Panel B",
  legend_position = "none"
)

npv_figure <- p_npv_A / p_npv_B +
  plot_layout(heights = c(1, 1))

print(npv_figure)


# ---- 11. Break-even carbon prices -----------------------------------------

make_bep_panel <- function(
  data,
  panel_title
) {

  summary_data <- data %>%
    group_by(scenario) %>%
    summarise(
      mean_bep =
        mean(bep_cad_per_tCO2e, na.rm = TRUE),
      p05 =
        quantile(
          bep_cad_per_tCO2e,
          0.05,
          na.rm = TRUE
        ),
      p95 =
        quantile(
          bep_cad_per_tCO2e,
          0.95,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) %>%
    mutate(
      scenario = factor(
        as.character(scenario),
        levels = scenario_levels
      ),
      mean_label = sprintf("%.2f", mean_bep)
    )

  ggplot(
    summary_data,
    aes(
      y = scenario,
      x = mean_bep,
      colour = scenario
    )
  ) +
    geom_linerange(
      aes(
        xmin = p05,
        xmax = p95
      ),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    geom_point(
      size = 2.2,
      show.legend = FALSE
    ) +
    geom_text(
      aes(label = mean_label),
      nudge_y = 0.24,
      vjust = 0,
      size = 3.3,
      colour = "black",
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = scenario_colours,
      breaks = scenario_levels,
      drop = FALSE
    ) +
    scale_y_discrete(
      limits = rev(scenario_levels)
    ) +
    scale_x_continuous(
      breaks = seq(0, 30, by = 5),
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(
      xlim = c(0, 30),
      clip = "off"
    ) +
    labs(
      title = panel_title,
      x = expression(
        "Break-even carbon price ($/tCO"[2]*"e)"
      ),
      y = NULL
    ) +
    theme_publication +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_line(
        colour = "grey82",
        linewidth = 0.45,
        linetype = "dashed"
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

p_bep_A <- make_bep_panel(
  data = filter(npv_all, rotation == "4-crop"),
  panel_title = "Panel A"
)

p_bep_B <- make_bep_panel(
  data = filter(npv_all, rotation == "2-crop"),
  panel_title = "Panel B"
)

bep_figure <- p_bep_A / p_bep_B +
  plot_layout(heights = c(1, 1))

print(bep_figure)


# ---- 12. Save figures -----------------------------------------------------

ggsave(
  filename = "figures/Figure_NPV_rotations.png",
  plot = npv_figure,
  width = 180,
  height = 185,
  units = "mm",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = "figures/Figure_BEP_rotations.png",
  plot = bep_figure,
  width = 180,
  height = 220,
  units = "mm",
  dpi = 600,
  bg = "white"
)

