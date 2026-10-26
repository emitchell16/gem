# projects/semaglutide_global_access/R/glp1_inputs.R

# description ------
# GLP-1 input tables, loaded per simulation config by hook_load_inputs() (called at the end of set_in_tables.R):
#   treatment dynamics parameters (sheet treatment_dynamics), treatment effects, double-counting calibration
#   multipliers, baseline GLP-1 use, and drug price scenarios (sheets rwe_sema_costs, alt_sema_costs).
# Globals set: model_param, model_param_ci (extended), glp1_impact_data, dc_multiplier_ds, glp1_baseline_use,
#   glp1_cost_inputs, cost_treat_cols, cost_total_cols.

# double-counting calibration multipliers (glp1_double_counting_calibration.R writes the file) ----
load_dc_multipliers <- function(te_type) {
  fp <- file.path(paths$working, paste0("double_counting_calibration_results_", te_type, ".csv"))
  if (!file.exists(fp)) return(data.table(outcome = character(), subgroup = character(), multiplier = numeric(), year = integer()))
  fread(fp)[, .(outcome, subgroup, multiplier, year)]
}
# trial HRs without a calibration row would count the BMI pathway twice
dc_check <- function() {
  te <- as.data.table(glp1_impact_data)[outcome %in% names(dc_outcome_def), .(outcome, subgroup = treated_pop)]
  miss <- te[!dc_multiplier_ds, on = .(outcome, subgroup)]
  if (nrow(miss)) warning("no double-count calibration for: ", paste(miss$outcome, miss$subgroup, sep = "/", collapse = "; "))
  invisible(nrow(miss) == 0)
}

# entry point (overrides the core no-op hook) ----
hook_load_inputs <- function(cs) {
  ge <- .GlobalEnv
  # treatment dynamics join the scalar model parameters and their PSA intervals
  treat_dyn <- read_excel(input_param_fp, sheet = "treatment_dynamics")
  need <- c("wl_resp_cv", "regain_frac_per_year", "disc_prob_not2d", "disc_prob_t2d",
            "restart_prob_not2d", "restart_prob_t2d", "disc_year_frac_on_tx",
            "disc_prob_t2d_later", "disc_prob_not2d_later", "restart_prob_t2d_later", "restart_prob_not2d_later")
  if (!all(need %in% treat_dyn$parameter)) stop("missing treatment_dynamics parameters: ", paste(setdiff(need, treat_dyn$parameter), collapse = ", "))
  assign("model_param", c(get("model_param", ge), setNames(as.numeric(treat_dyn$value), treat_dyn$parameter)), envir = ge)
  ci <- treat_dyn |> select(any_of(c("parameter", "lower", "upper", "psa_dist", "psa_family"))) |>
    filter(!is.na(lower), !is.na(upper), psa_dist %in% c("beta", "lognormal", "normal")) |>
    mutate(across(c(lower, upper), as.numeric))
  assign("model_param_ci", bind_rows(get("model_param_ci", ge), ci), envir = ge)

  assign("glp1_impact_data", set_treatment_effect(cs$glp1$treatment_effect_source), envir = ge)
  assign("dc_multiplier_ds", load_dc_multipliers(cs$glp1$treatment_effect_source), envir = ge)
  assign("glp1_baseline_use", baseline_use_inputs(cs$baseline_use), envir = ge)
  if (isTRUE(cs$costing)) glp1_load_cost_inputs() else {
    assign("glp1_cost_inputs", NULL, envir = ge)
    assign("cost_treat_cols", character(0), envir = ge); assign("cost_total_cols", character(0), envir = ge)
  }
  invisible(NULL)
}

# drug price scenarios (monthly price by country) + price-scenario cost columns ----
glp1_load_cost_inputs <- function() {
if (config_sims$costing ) { 
  if (!exists("glp1_cost_inputs", envir = .GlobalEnv) || is.null(get("glp1_cost_inputs", envir = .GlobalEnv))) {   # NULL left from a costing = FALSE run
    country_groups = read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv"))
    ### real world costs: current prices only (data_year >= price_year_min); other countries imputed below ----
    glp1_rwe_cost_data <- read_excel(input_param_fp, sheet = "rwe_sema_costs") |>
      filter(!is.na(iso3), as.numeric(data_year) >= price_year_min) |>
      select(iso3, unit_price_lcu, unit_price_int_dollar) 
  
    glp1_rwe_cost_data_imp <- heirarchical_imputation(glp1_rwe_cost_data, var1 = "unit_price_int_dollar", var1_imp = "unit_price_int_dollar_imp",
                                                  out_fn = "glp1_rwprices_impute_summary",
                                                  stat = median, levels = c("income", "region", "global"), min_n = 3L)
    rm(glp1_rwe_cost_data )
    ### literature-based costs ----
    glp1_lit_cost_scenarios <- read_excel(input_param_fp, sheet = "alt_sema_costs") |>
      select(est_name, unit_price_low, unit_price_high) |> 
      mutate( "low income" = unit_price_low,
              "high income" = unit_price_high,
              "lower middle income" = round(unit_price_low + (unit_price_high - unit_price_low) / 3,2),
              "upper middle income" = round(  unit_price_low + 2 * (unit_price_high - unit_price_low) / 3,2)) |> 
      select(-unit_price_low, -unit_price_high) |> 
      pivot_longer(cols = -est_name,
        names_to = "income_group",
        values_to = "value"
      ) |>
      pivot_wider(names_from = est_name, values_from = value)
    
    glp1_rwe_cost_data <- glp1_rwe_cost_data_imp |> 
      left_join(glp1_lit_cost_scenarios, by="income_group") 
    
    rm(glp1_rwe_cost_data_imp, glp1_lit_cost_scenarios)
      #### convert usd -> lcu -> ppp international dollar ----
    lcu_per_id = read_xlsx(file.path(paths$data, "for_costing", "ppp_LCU.xlsx")) |>    # LCU per int. dollar (ID) aka implied PPP conversion rate
      select(iso3, c_2021:c_2026) |> 
      rowwise() |>
      mutate(lcu_per_id = {
        vals <- c_across(c_2021:c_2026)
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) NA_real_ else tail(vals, 1)
      }) |>  ungroup() |>  select(-starts_with("c_"))
    lcu_per_usd = read_xlsx(file.path(paths$data, "for_costing", "WB_LCU_per_USD.xlsx")) |> 
      select(iso3="Country Code", lcu_per_usd =  lcu_per_usd_clean)
  
    lcu_usd_imputed <- heirarchical_imputation(lcu_per_usd, var1 = "lcu_per_usd", var1_imp = "lcu_per_usd_imp",
                                                      out_fn = "glp1_lcu-usd_conversionfactor_impute_summary")
    lcu_id_imputed <- heirarchical_imputation(lcu_per_id, var1 = "lcu_per_id", var1_imp = "lcu_per_id_imp",
                                               out_fn = "glp1_lcu-id_conversionfactor_impute_summary")
    
    lcu_id_and_usd_imputed <- lcu_usd_imputed |> select(iso3, lcu_per_usd_imp) |> 
      left_join(lcu_id_imputed |>  select(iso3, lcu_per_id_imp), by = "iso3") |> 
      mutate(usd_per_id_plr = lcu_per_id_imp / lcu_per_usd_imp )  # price level ratio (PLR)
  
      # clean and apply conversion 
          ##' int$ = USD * PLR 
    glp1_cost_inputs <- glp1_rwe_cost_data |> 
      left_join(lcu_id_and_usd_imputed %>% select(iso3, usd_per_id_plr), by = "iso3") |> 
      mutate(
        price_id_sust_costx1 = sust_cost_price_usd * usd_per_id_plr,  # price floor from Barber et al. paper
        price_id_sust_costx3 = sust_cost_pricex3_usd * usd_per_id_plr, # implementation-adjusted sustainable price 
        
        price_id_lowest_mp = lowest_market_price_usd * usd_per_id_plr,
        price_id_generic_entry = generic_entry_price_usd * usd_per_id_plr
      ) |>
      select(iso3, region, price_intl_rwe =  unit_price_int_dollar_imp, starts_with("price")) |> 
      filter(region != "no who region" ) |>  # removes one XOD which isn't iso3 in UN data
      mutate(
        across(starts_with("price"), ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x)
        ))
    if (!file.exists(file.path(paths$out, "costing_region_marketprice_vs_realworld_summary.xlsx"))){ 
      summary_region_mp_vs_rw <- glp1_cost_inputs |> 
        group_by(region) |> 
        summarise(
          mean_price_id_lowest_mp = mean(price_id_lowest_mp, na.rm = TRUE),
          mean_price_intl_rwe     = mean(price_intl_rwe, na.rm = TRUE),
          n_countries = n(),
          .groups = "drop"
        )
      write_xlsx(summary_region_mp_vs_rw, file.path(paths$out, "costing_region_marketprice_vs_realworld_summary.xlsx"))
      rm(summary_region_mp_vs_rw)
    }
    rm(lcu_id_imputed, lcu_usd_imputed, glp1_rwe_cost_data, lcu_per_id, lcu_per_usd,lcu_id_and_usd_imputed)
    
    glp1_cost_inputs <- glp1_cost_inputs |>  select(-region) 
    } else { glp1_cost_inputs <- get("glp1_cost_inputs", envir = .GlobalEnv) }
  # treatment price scenarios carried in person-level cost columns
  price_suffixes  <- sub("^price_", "", grep("^price_", names(glp1_cost_inputs), value = TRUE))
  cost_treat_cols <- paste0("cost_treat_total_", price_suffixes)
  cost_total_cols <- paste0("cost_total_", price_suffixes)
  assign("glp1_cost_inputs", glp1_cost_inputs, envir = .GlobalEnv)
  assign("cost_total_cols", cost_total_cols, envir = .GlobalEnv)
  assign("cost_treat_cols", cost_treat_cols, envir = .GlobalEnv)
  } else assign("glp1_cost_inputs", NULL, envir = .GlobalEnv)
}
