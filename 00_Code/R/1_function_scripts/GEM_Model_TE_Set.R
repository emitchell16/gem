# -------------------'
# Function to set treatment effect input parameters
# Author: Liz Mitchell
# Date updated: March 2026
# -------------------'

#########' set before calling in 00_Top file:
#########'   te_source = config_sims$glp1$treatment_effect_source


set_treatment_effect <- function(te_source= c("rct", "rwe"), 
                                 value = c("val", "lower", "upper")) {
  te_source <- match.arg(te_source)
  value     <- match.arg(value)
  group_vars <- c("treated_pop", "age_group", "outcome_unit", "outcome")
  
  glp1_impact_data <- read_excel(input_param_fp, sheet = "semaglutide_treatment_effects") |> 
    select(treated_pop,	age_group, outcome_unit,	outcome, 
           outcome_val,	outcome_lower, outcome_upper, source_type) |>
    group_by(across(all_of(group_vars)), source_type) |> 
    summarise(
      outcome_val   = mean(outcome_val,   na.rm = TRUE),
      outcome_lower = mean(outcome_lower, na.rm = TRUE),
      outcome_upper = mean(outcome_upper, na.rm = TRUE),
      n_rows        = n(),
      .groups = "drop"
    )
  if (te_source == "rct") {
    out <- glp1_impact_data |>
      filter(source_type == "rct") |>
      mutate(
        outcome_val = case_when(
          value == "lower" ~ outcome_lower,
          value == "upper" ~ outcome_upper,
          TRUE ~ outcome_val
        ),
        source_used = "rct"
      ) |>
      select(all_of(group_vars), outcome_val, outcome_lower, outcome_upper, n_rows, source_used)
    
    if (value == "lower") {
      message("Applying lower estimated treatment effects from RCT source.")
    } else if (value == "upper") {
      message("Applying upper estimated treatment effects from RCT source.")
    } else {
      message("Applying point estimate treatment effects from RCT source.")
    }
  }
  
  if (te_source == "rwe") {
    rwe_te <- glp1_impact_data |>
      filter(source_type == "rwe") |>
      select(-source_type) |>
      rename(
        rwe_outcome_val   = outcome_val,
        rwe_outcome_lower = outcome_lower,
        rwe_outcome_upper = outcome_upper,
        rwe_n_rows        = n_rows
      )
    
    rct_te <- glp1_impact_data |>
      filter(source_type == "rct") |>
      select(-source_type) |>
      rename(
        rct_outcome_val   = outcome_val,
        rct_outcome_lower = outcome_lower,
        rct_outcome_upper = outcome_upper,
        rct_n_rows        = n_rows
      )
    
    out <- rct_te |>
      left_join(rwe_te, by = group_vars) |>
      mutate(
        outcome_val   = coalesce(rwe_outcome_val, rct_outcome_val),
        outcome_lower = coalesce(rwe_outcome_lower, rct_outcome_lower),
        outcome_upper = coalesce(rwe_outcome_upper, rct_outcome_upper),
        n_rows        = coalesce(rwe_n_rows, rct_n_rows),
        source_used = case_when(
          !is.na(rwe_outcome_val) ~ "rwe",
          !is.na(rct_outcome_val) ~ "rct_imputed",
          TRUE ~ NA_character_
        )
      ) |>
      select(all_of(group_vars), outcome_val, outcome_lower, outcome_upper, n_rows, source_used)
    
    message("Applying RWE treatment effects where available; imputing missing cells with RCT values.")
    if (value %in% c("lower", "upper")) {
      message("For te_source = 'rwe', lower/upper option ignored and point estimates are used.")
    }
  }
  
  out
}
  