# post-simulation secondary plots and tables:
# WIP 4.1
# note: cost and daly outcomes have already been scaled to full pop

## 1. region-level results to combine to give UIs for plotting ----------
# main_val <- read_xlsx(file.path(out_folder, "RES_access_scenarios", "100glp1_0.1sample_5yr_3.28", "regional_result_table_5_100.xlsx"))
# lb_val <- read_xlsx(file.path(out_folder, "RES_TE_SA", "LOWER_100glp1_0.1sample_5yr_3.28", "regional_result_table_5_100.xlsx")) %>%
#   select(-outcome_value_0_base, -outcome_value_end_base) |> 
#   rename_with(~ paste0(.x, "_lb"), -c(who_region, condition))
# ub_val <- read_xlsx(file.path(out_folder, "RES_TE_SA", "UPPER_100glp1_0.1sample_5yr_3.28", "regional_result_table_5_100.xlsx"))%>%
#   select(-outcome_value_0_base, -outcome_value_end_base) |> 
#   rename_with(~ paste0(.x, "_ub"), -c(who_region, condition))
# 
# 
# wip <- main_val |>  full_join(lb_val,  by = c("who_region", "condition")) |> 
#   full_join(ub_val, by = c("who_region", "condition")) |> 
#   select(condition, who_region,  "outcome_value_0_base" ,    "outcome_value_end_base",
#          "outcome_value_end_int", "outcome_value_end_int_lb","outcome_value_end_int_ub",
#          "abs_diff_pct" , "abs_diff_pct_lb", "abs_diff_pct_ub", 
#          "rel_diff_pct" , "rel_diff_pct_lb", "rel_diff_pct_ub" )
# 
# write_xlsx(wip, file.path(out_folder, "region_level_uncertainty_int_res_100_0.1.xlsx")) # done
#' 
#' 
#' 
## 2. country level ICER  ----------
#'  table ofincremental cost | incremental DALY, by country and wb income group (4 cat), and for 3 cost scenarios (0, cost based, and real world)
#'        ICER = (total cost under treatment - total cost under status quo) / (DALYs under status quo - DALYs under treatment)
#'        report exactly the ICER by cost scenario
#'        use to plot panel/faceted plot by country income lelel
#'        * also try making it a heat map for ICER under current costs (sf, rnaturalearth?)
#'
# main_val <- read_xlsx(file.path(out_folder,
#                                 "RES_access_scenarios", "100glp1_0.1sample_5yr_3.28", "sim_result_iso3_100glp1.xlsx")) |>
#   select( "sim", "iso3", "n",
#           "yll_adj", "yld_adj_total", "daly_total_adj", starts_with("cost_total_"), cost_total_tc_zero =  "cost_conditions_total") |>
#   pivot_longer(
#     cols = starts_with("cost_total_"),
#     names_prefix = "cost_total_",
#     names_to = "cost_scenario",
#     values_to = "cost_total"
#   ) |>
#   rename_with(~ paste0(.x, "_val"), -c(sim, iso3, n, cost_scenario))
#
# lb_val <- read_xlsx(file.path(out_folder, "RES_TE_SA", "LOWER_100glp1_0.1sample_5yr_3.28", "sim_result_iso3_100glp1.xlsx")) |>
#   select( "sim", "iso3",
#           "yll_adj", "yld_adj_total", "daly_total_adj", starts_with("cost_total_"),  cost_total_tc_zero = "cost_conditions_total") |>
#   filter(sim!= "5yr_basecase") |>
#   pivot_longer(
#     cols = starts_with("cost_total_"),
#     names_prefix = "cost_total_",
#     names_to = "cost_scenario",
#     values_to = "cost_total"
#   ) |>
#   rename_with(~ paste0(.x, "_lb"), -c(sim, iso3, cost_scenario))
#
# ub_val <- read_xlsx(file.path(out_folder, "RES_TE_SA", "UPPER_100glp1_0.1sample_5yr_3.28", "sim_result_iso3_100glp1.xlsx"))%>%
#   select( "sim", "iso3",
#           "yll_adj", "yld_adj_total", "daly_total_adj", starts_with("cost_total_"), cost_total_tc_zero = "cost_conditions_total") |>
#   filter(sim!= "5yr_basecase") |>
#   pivot_longer(
#     cols = starts_with("cost_total_"),
#     names_prefix = "cost_total_",
#     names_to = "cost_scenario",
#     values_to = "cost_total"
#   ) |>
#   rename_with(~ paste0(.x, "_ub"), -c(sim, iso3, cost_scenario))
#
# country_ce_temp <- main_val |>
#   full_join(lb_val,  by = c("sim", "iso3", "cost_scenario")) |>
#   full_join(ub_val,  by = c("sim", "iso3", "cost_scenario")) |>
#   left_join(country_groups %>% select(iso3, income_group = world_bank_income_group, who_region),
#             by = "iso3") |>
#   mutate(across(matches("_(lb|ub)$"), # base sim doesn't have lb/ub
#                 ~ if_else(is.na(.x),
#                 get(sub("_(lb|ub)$", "_val", cur_column())),
#                 .x))
#       ) |>
#   group_by(iso3, income_group,who_region, cost_scenario, n ) |>
#   mutate(
#   incremental_cost_val = cost_total_val - cost_total_val[sim == "5yr_basecase"],
#   incremental_cost_lb  = cost_total_lb - cost_total_lb[sim == "5yr_basecase"],
#   incremental_cost_ub  = cost_total_ub - cost_total_ub[sim == "5yr_basecase"],
#
#   daly_averted_val = daly_total_adj_val[sim == "5yr_basecase"] - daly_total_adj_val,
#   daly_averted_lb  = daly_total_adj_lb[sim == "5yr_basecase"] - daly_total_adj_lb,
#   daly_averted_ub  = daly_total_adj_ub[sim == "5yr_basecase"] - daly_total_adj_ub,
#
#   icer_val = incremental_cost_val / daly_averted_val,
#   icer_lb  = incremental_cost_lb  / daly_averted_lb,
#   icer_ub  = incremental_cost_ub  / daly_averted_ub,
#   across(ends_with("_lb"),
#            ~ pmin(.x, get(sub("_lb$", "_ub", cur_column())), na.rm = FALSE),
#            .names = "{sub('_lb$', '_lower', .col)}"),
#   across(ends_with("_ub"),
#            ~ pmax(.x, get(sub("_ub$", "_lb", cur_column())), na.rm = FALSE),
#            .names = "{sub('_ub$', '_upper', .col)}")
#   ) |>
#   ungroup() |>
#   filter(sim != "5yr_basecase") |>  select(-sim, -ends_with("lb"), -ends_with("ub")) |>
#   select(iso3, income_group, who_region, n, cost_scenario, starts_with("icer"), starts_with("increm"),
#          starts_with("cost"),
#          starts_with("yll"), starts_with("yld"), starts_with("daly")
#   )
#
#
# write_xlsx(country_ce_temp, file.path(out_folder, "country_level_icer_res_100_0.1.xlsx"))


## 3. plot icer summary ---------
library(ggplot2)
library(forcats)
library(ggrepel)


plot_df <- country_ce_temp |>
  filter(cost_scenario %in% c("intl_rwe", "id_sust_costx1")) |> # "tc_zero"
  filter((icer_val/1e3) >-20) |>
  mutate(n = n*10, # for 10% sample
    icer_k = icer_val / 1e3,
    is_outlier = icer_k > 200,
    label_iso3 = if_else(
      !is.na(iso3) & (icer_k > 25 | n > 5e8),
      iso3,
      NA_character_
    ),
    cost_scenario = recode(
      cost_scenario,
      intl_rwe = "Estimated Current Price",
      id_sust_costx1 = "Cost-Based Price",
      tc_zero = "Zero Treatment Cost"
    ),
    income_group = income_group %>%
      replace_na("not assigned") %>%
      tolower() %>%
      recode("no wb income group" = "not assigned") %>%
      str_to_title(),
    income_group =factor(income_group,
                         levels = c(
                           "Low Income", "Lower Middle Income", "Upper Middle Income",
                           "High Income", "Not Assigned"))
  )

plot_df %>%
  filter(is_outlier) %>%
  group_by(cost_scenario, income_group) %>%
  summarise(
    n_outliers = n(),
    max_icer = max(icer_k, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(plot_df %>% filter(!is_outlier), aes(x = income_group, y = icer_val/1e3, size = n)) +
  coord_cartesian(ylim = c(min(plot_df$icer_val/1e3, na.rm = TRUE), 170)) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.15, height = 0)) +
  geom_hline(yintercept = c(1, 50, 100), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.8) +
  geom_text_repel(
    aes(label = iso3), size = 3,
    max.overlaps = 8, box.padding = 0.4,
    point.padding = 0.4, min.segment.length = 0
  ) +
  facet_wrap(~ cost_scenario, nrow = 1) +
  scale_size_continuous(name = "Population (millions)",
                        labels = scales::label_number(
                          scale = 1e-6, suffix = "M", accuracy = 1)) +
  labs(
    x = NULL,
    y = "ICER\n(1,000 Int$ per DALY averted)"
  ) +
  theme_gray() +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

ggsave(filename =  file.path(out_folder, "FIG_iso3_icer.png"),
       width = 10, height = 5,
       dpi = 600, bg="white")


# v2 grouped by income;
library(scales)
library(ggsci)

plot_df2 <- country_ce_temp |>
  filter(cost_scenario %in% c("tc_zero", "intl_rwe", "id_sust_costx1", "id_sust_costx3")) |>
  mutate(
    n = n * 10,  # for 10% sample
    cost_scenario = recode(
      cost_scenario,
      intl_rwe = "Estimated Current Price",
      id_sust_costx1 = "Cost-Based Price",
      id_sust_costx3 = "Implementation Adjusted Cost-Based Price",
      tc_zero = "Zero Treatment Cost"
    ),
    income_group = income_group %>%
      replace_na("not assigned") %>%
      tolower() %>%
      recode("no wb income group" = "not assigned") %>%
      str_to_title(),
    income_group = factor(
      income_group,
      levels = c(
        "Low Income", "Lower Middle Income", "Upper Middle Income",
        "High Income", "Not Assigned"
      )
    )
  )

plot_sum <- plot_df2 |>
  group_by(income_group, cost_scenario) |>
  summarise(
    n = sum(n, na.rm = TRUE),
    incremental_cost_val   = sum(incremental_cost_val, na.rm = TRUE),
    incremental_cost_lower = sum(incremental_cost_lower, na.rm = TRUE),
    incremental_cost_upper = sum(incremental_cost_upper, na.rm = TRUE),
    daly_averted_val       = sum(daly_averted_val, na.rm = TRUE),
    daly_averted_lower     = sum(daly_averted_lower, na.rm = TRUE),
    daly_averted_upper     = sum(daly_averted_upper, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    icer_val   = incremental_cost_val / daly_averted_val,
    icer_lower = incremental_cost_lower / daly_averted_lower,
    icer_upper = incremental_cost_upper / daly_averted_upper
  )

p <- ggplot(plot_sum, aes(x = income_group, y = icer_val / 1e3, 
       color = cost_scenario)) +
  coord_cartesian(
    ylim = c(min(plot_sum$icer_val / 1e3, na.rm = TRUE), 60)
  ) +
  scale_y_continuous( breaks = seq(0, 52, by = 10)) +
  geom_point(size =4, alpha = 0.8, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = c(1, 5, 50, 100), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.8) +
  annotate("label", x = 4, y = 1,  label = "1k Int$/DALY",  hjust = 0, vjust = -0.3,
           size =3, color = "gray40",
           fill = alpha("white", 0.7), label.size = 0) +
  annotate("label", x = 4, y = 5,  label = "5k Int$/DALY",  hjust = 0, vjust = -0.3,
           size = 3, color = "gray40",
           fill = alpha("white", 0.7), label.size = 0) +
  annotate("label", x = 4, y = 50, label = "50k Int$/DALY", hjust = 0, vjust = -0.6,
           size = 3, color = "gray40",
           fill = alpha("white", 0.7), label.size = 0) +
  geom_vline(
    xintercept = seq(1.5, length(levels(plot_sum$income_group)) - 0.5, by = 1),
    color = "gray90",
    linewidth = 0.4
  ) +
  labs(
    x = NULL,
    y = "ICER\n(1,000 Int$ per DALY averted)",
    color = "Semaglutide Cost Scenario"
  ) +
  theme_classic(base_size = 11) +
  scale_color_lancet() +
  theme(
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.5),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

p

ggsave(
  filename = file.path(out_folder, "FIG_incomelvl_icer.png"),
  plot = p,
  width = 10, height = 8,
  dpi = 600, bg = "white"
)