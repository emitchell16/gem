#' GEM simulation input data tables
#' run 1x per simulation
#' -------------------------

##############'
# Helpers ----
input_param_fp <- file.path(paths$data, "GEM_input_parameters.xlsx")
input_params <- read_excel(input_param_fp, sheet = "baseline_dataset") |> select(parameter, value)

# clean input data strings for matching
norm_key <- function(x) {
  x |> stri_replace_all_fixed("\u00A0", " ") |> str_squish() |>  str_to_lower()                           
}

# for costing data, provide estimated values for missing countries/territories
heirarchical_imputation <- function(ds, var1, var1_imp, out_fn) {
  country_groups = read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")) |> 
    filter(who_region != "no who region")
  
  source_col <- paste0(var1_imp, "_source")
  
  ds2 <- ds |>
    right_join( country_groups, by=c("iso3")) |> select(-X) |>
    rename(income_group = world_bank_income_group, region = who_region) |>
    group_by(income_group, region) |>
    mutate(
      income_reg_mean = mean(.data[[var1]], na.rm = TRUE)
    ) |>
    ungroup() |>
    group_by(income_group) |>
    mutate(
      income_mean = mean(.data[[var1]], na.rm = TRUE)
    ) |>
    ungroup() |>
    group_by(region) |>
    mutate(
      region_mean = mean(.data[[var1]], na.rm = TRUE)
    ) |>
    ungroup() |>
    mutate(
      global_mean = mean(.data[[var1]], na.rm = TRUE),
      !!var1_imp := coalesce( 
        .data[[var1]],
        income_reg_mean,
        income_mean,
        region_mean,
        global_mean
      ),
      !!source_col := case_when(
        !is.na(.data[[var1]]) ~ "observed",
        !is.na(income_reg_mean) & is.finite(income_reg_mean) ~ "income_region_mean",
        !is.na(income_mean) & is.finite(income_mean) ~ "income_mean",
        !is.na(region_mean) & is.finite(region_mean) ~ "region_mean",
        TRUE ~ "global_mean")
    ) |> select(-ends_with("_mean"))
  # summary table: 
  summary_table <- ds2 |>
    group_by(region, income_group) |>
    summarise(
      n_countries = n(),
      mean_est = mean(.data[[var1_imp]], na.rm = TRUE),
      median_est = median(.data[[var1_imp]], na.rm = TRUE),
      sd_est = sd(.data[[var1_imp]], na.rm = TRUE),
      min_est = min(.data[[var1_imp]], na.rm = TRUE),
      max_est = max(.data[[var1_imp]], na.rm = TRUE),
      n_imputed = sum(source_col != "observed", na.rm = TRUE),
      prop_imputed = n_imputed / n_countries,
      .groups = "drop"
    ) |> 
    mutate(
      mean_sd = case_when(!is.na(sd_est) ~ sprintf("%.0f (%.0f)", mean_est, sd_est), TRUE ~ sprintf("%.0f", mean_est)),
      median_iqr = sprintf("%.0f [%.0f–%.0f]", median_est, min_est, max_est),
      prop_imputed_pct = sprintf("%.1f%%", 100 * prop_imputed)
    ) |>
    select(region, income_group, n_countries, mean_sd, median_iqr, prop_imputed_pct)
  
  write_xlsx(summary_table, file.path(out_folder, paste0("costing_", out_fn,".xlsx")))
  
  ds2
}

##############'
# Double-counting calibration multipliers -----

te_type = config_sims$glp1$treatment_effect_source
dc_multiplier_ds <- read.csv(
  file.path(paths$working, paste0("double_counting_calibration_results_", te_type, ".csv"))) |>
  select(outcome, subgroup, multiplier, year)

##############'
# Incidence and death parameters ----
all_cohort_param_rates <- cohort_input_data |> 
  select(-starts_with("prevalence_"),  -starts_with("prob_death"), - pop_size) |>
  mutate(incidence_ac_death_upper = pmin(incidence_ac_death_upper, 1))

# Risk modifiers ----
.prop_stroke_ischaemic <- input_params$value[input_params$parameter == "prop_stroke_ischaemic"]

obesity_rr = read_excel(input_param_fp, sheet = "obesity_rr_nejm") |>  select(-source)
obesity_rr_by_age <- obesity_rr |> 
  mutate(across(c(upper, lower, value), as.numeric)) |> 
  rename(rr_condition_with_obesity = RR_morbidity_and_mortality_obesity) |> 
  mutate(rr_key = norm_key(rr_condition_with_obesity)) |>
  group_by(rr_key) |>
  # expand age categories and impute with nearest
  complete(age_low = c(10, 15, 20, 25, 85, 90, 95)) |> 
  arrange(age_low) |> 
  fill(value, upper, lower, .direction = "downup") |> 
  ungroup() |> 
  select(- age_range) |> 
  group_by(age_low) |> 
  mutate(
    # Replace stroke subtypes with weighted average 
    across(c(value, upper, lower), ~ if_else(
      rr_key == norm_key("Ischaemic stroke"),
      .prop_stroke_ischaemic * . + (1 - .prop_stroke_ischaemic) * value[rr_key == norm_key("Hemorrhagic stroke")],
      .)),
    rr_condition_with_obesity =  if_else(rr_key == norm_key("Ischaemic stroke"), norm_key("stroke"), rr_key)
  ) |>
  ungroup() |> 
  filter(rr_key != norm_key("Hemorrhagic stroke")) |> 
  arrange(rr_condition_with_obesity) |> 
  mutate(age_high = age_low + 4) |>
  select(-rr_key)

## * obesity and mortality ---- 
obesity_rr_acdeath_byclass = read_excel(input_param_fp, sheet = "obesity_rr_acdeath_byclass") |>  select(-source)
obesity_rr_acdeath_bysex   = read_excel(input_param_fp, sheet = "obesity_rr_acdeath_bysex") |>  select(-source)
obesity_rr_acdeath_byage   = read_excel(input_param_fp, sheet = "obesity_rr_acdeath_byage") |>  select(-source)
obesity_rr_cvd             = read_excel(input_param_fp, sheet = "obesity_rr_cvd") |>  select(-source)
obesity_rr_ckd             = read_excel(input_param_fp, sheet = "obesity_rr_ckd") |>  select(-source)
## * t2d comorbidities ---- 
t2d_rr                     = read_excel(input_param_fp, sheet = "t2d_rr") |>  select(-source)
## * stroke ----
cvd_rr_stroke_bysex        = read_excel(input_param_fp, sheet = "cvd_rr_stroke_bysex") |>  select(-source)
## * eskd mortality ---- 
min_death_risk_eskd        = input_params$value[input_params$parameter == "min_1yrdeath_risk_eskd"]

##############'
# Disability weights ----
dw_tbl <- fread(file.path(paths$working, "condition_disability_wts.csv"))

#################'
# Sample parameters if PSA --------

rr_tables <- list(
  obesity_rr_by_age = obesity_rr_by_age,
  obesity_rr_acdeath_byclass = obesity_rr_acdeath_byclass,
  obesity_rr_acdeath_bysex = obesity_rr_acdeath_bysex,
  obesity_rr_acdeath_byage = obesity_rr_acdeath_byage,
  obesity_rr_cvd = obesity_rr_cvd,
  obesity_rr_ckd = obesity_rr_ckd,
  t2d_rr = t2d_rr,
  cvd_rr_stroke_bysex = cvd_rr_stroke_bysex
)

if (config_sims$psa) {
  # sample parameter from distribution for PSA runs
  sample_rate <- function(distribution = c("normal", "logit-normal", "log-normal"),
                          value, lower, upper, n = 1) {
    distribution <- match.arg(distribution)
    z <- 1.96 ; eps =1e-8
    if (upper == lower || is.na(upper) || is.na(lower)) return(rep(value, n))
    
    if (distribution == "normal") {
      lo <- min(lower, upper)
      hi <- max(lower, upper)
      sd_val <- (hi - lo) / (2 * z)
      if (is.na(sd_val) || abs(sd_val) < eps) {
        return(rep(value, n))
      }
      
      return(rnorm(n, mean = value, sd = sd_val))
    }
    if (distribution == "log-normal") {
      if (any((lower <= 0 | upper <= 0 | value <= 0))) { 
        warning("sample_rate: non-positive inputs for log-normal.")
        return(rep(value, n))
      }
      mu <- log(value)
      sd_val <- (log(upper) - log(lower)) / (2 * z)
      if (abs(sd_val) < eps) {  
        return(rep(value, n))  
      }
      return(exp(rnorm(n, mean = mu, sd = sd_val)))
    }
    if (distribution == "logit-normal") {
      logit <- function(p) log(p / (1 - p))
      inv_logit <- function(x) 1 / (1 + exp(-x))
      value <- pmin(1-eps, pmax(eps, value))
      upper = pmin(1-eps, upper )
      lower = pmax(eps, lower) 
      mu <- logit(value)
      sd_val <- (logit(upper) - logit(lower)) / (2 * z) 
      if (abs(sd_val) < eps) {  
        return(rep(value, n))  
      }
      return(inv_logit(rnorm(n, mean = mu, sd = sd_val)))
    }
  }

  glp1_impact_data <- read_excel(input_param_fp, sheet = "semaglutide_treatment_effects") |>
    filter(source_type == "rct") |>
    select(treated_pop,	age_group, outcome_unit,	outcome, outcome_val,	outcome_lower, outcome_upper) |>
    group_by(treated_pop, age_group, outcome_unit, outcome)  |>
    summarise(
      outcome_val   = mean(outcome_val,   na.rm = TRUE),
      outcome_lower = mean(outcome_lower, na.rm = TRUE),
      outcome_upper = mean(outcome_upper, na.rm = TRUE),
      n_rows        = dplyr::n(),
      .groups = "drop"
    )
  
  psa_params <- draw_psa_params(
    all_cohort_param_rates,
    dw_tbl,
    glp1_impact_data,
    rr_tables)
  
  all_cohort_param_rates<- psa_params$cohort_params
  dw_tbl                <- psa_params$dw_tbl
  glp1_impact_data      <- psa_params$glp1_impact_data
  rr_tables             <- psa_params$rr_tables
  
  obesity_rr_by_age          <- rr_tables$obesity_rr_by_age
  obesity_rr_acdeath_byclass <- rr_tables$obesity_rr_acdeath_byclass
  obesity_rr_acdeath_bysex   <- rr_tables$obesity_rr_acdeath_bysex
  obesity_rr_acdeath_byage   <- rr_tables$obesity_rr_acdeath_byage
  obesity_rr_cvd             <- rr_tables$obesity_rr_cvd
  obesity_rr_ckd             <- rr_tables$obesity_rr_ckd
  t2d_rr                     <- rr_tables$t2d_rr
  cvd_rr_stroke_bysex        <- rr_tables$cvd_rr_stroke_bysex
}

##############'
# Health event costs -----
#'      sets 'health_event_cost_inputs'

# conventional 3% discount factor
disc_rate <- 0.03

if (config_sims$costing) { 
  if (!exists("health_event_cost_inputs")) {
  # CHE per capita scaling:
  che_ppp_percap = read_xlsx(file.path(paths$data, "for_costing", "WHO_GHED_CHE_percap_PPP_2024.XLSX")) |> 
    select(location, iso3 = code, year, che_ppp_pc) |> 
    filter(year > 2020) |> arrange(iso3, desc(year)) |>
    group_by(iso3) |>
    summarise(che_ppp_pc = che_ppp_pc[which(!is.na(che_ppp_pc))[1]],
      che_ppp_pc_year = year[which(!is.na(che_ppp_pc))[1]],
      .groups = "drop")
  cat("range of years from WHO GHED CHE per capita across countries/territories:\n")
  print(summary(che_ppp_percap$che_ppp_pc_year))
  
  che_ppp_percap_imp <- heirarchical_imputation(che_ppp_percap, var1 = "che_ppp_pc", var1_imp = "che_ppp_pc_imp",
                                                out_fn = "CHD_per capita_ppp_impute_summary")
  # process health cost data:
 
  # country-level cost estimates: (ckd/stroke/eskd) 
  intl_countrylevel_health_costs <-  read_excel(input_param_fp, sheet = "intl_health_costs") |>  
     select(-source, -"source notes", - region, -income_group)
  intl_countrylevel_health_costs_obs <- intl_countrylevel_health_costs |>  
      pivot_wider(
        names_from = event_category,
        values_from = cost_intl,
        names_prefix = "cost_id_")
  
  ref_conds <- intl_countrylevel_health_costs |> filter(iso3 == "USA") |> 
    rename(cost_usd = cost_intl)
  #US cost estimates: (death/t2d/cvd/stroke (acute and chronic))
  usd_health_costs <- read_excel(input_param_fp, sheet = "us_health_costs") |> 
    bind_rows(ref_conds) |> filter(! ( event_category == "stroke" & is.na(source))) # prefer ahrq stroke cost 
    
    # Scale costs for each iso3 using CHE proxy
    ## cost_intl = cost_usd * (ratio of local to us CHE per capita)
    ref_che <- che_ppp_percap_imp |> filter(iso3 == "USA") |> pull(che_ppp_pc_imp)
    
    costs_int <- usd_health_costs |> select(-iso3) |> 
      crossing(
        che_ppp_percap_imp |> select(iso3, che_ppp_pc_imp)
      ) |>
      mutate(
        cost_int_dollar = cost_usd * (che_ppp_pc_imp / ref_che)
      ) |>
      select(iso3, event_category, cost_int_dollar) |>
      pivot_wider(
        names_from = event_category,
        values_from = cost_int_dollar,
        names_prefix = "cost_id_"
      )
    
    # combine sources of health state costs:
    cost_cols_both <- intersect(
      grep("^cost_id", names(costs_int), value = TRUE),
      grep("^cost_id", names(intl_countrylevel_health_costs_obs), value = TRUE)
    )
    
    health_event_cost_inputs <- costs_int |>
      left_join(intl_countrylevel_health_costs_obs, by = "iso3", suffix = c("", "_temp")) |>
      mutate(across(starts_with("cost_id"), as.double))
    for (nm in cost_cols_both) {
      temp_nm <- paste0(nm, "_temp")
      health_event_cost_inputs[[nm]] <- coalesce(
        health_event_cost_inputs[[temp_nm]],
        health_event_cost_inputs[[nm]])
    }
    health_event_cost_inputs <- health_event_cost_inputs |> select(-ends_with("temp"), -cost_id_stroke) 
    
    write_xlsx(health_event_cost_inputs, file.path(paths$working, "health_event_cost_inputs.xlsx"))
    
    rm(cost_cols_both,intl_countrylevel_health_costs_obs, intl_countrylevel_health_costs,
       che_ppp_percap_imp, che_ppp_percap,costs_int, usd_health_costs, ref_che)
  } else {health_event_cost_inputs }
  } else { health_event_cost_inputs = NULL }

##############'
# *** Semaglutide Treatment Effect & Costing Inputs: *** --------
#'        sets 'glp1_impact_data' & 'glp1_cost_inputs'
##############'
## 2. Treatment cost data -----
if (config_sims$costing ) { 
  if (!exists("glp1_cost_inputs")) { 
    country_groups = read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv"))
    ### real world costs (collected Jan-March 2026) ----
    glp1_rwe_cost_data <- read_excel(input_param_fp, sheet = "rwe_sema_costs") |>
      select(iso3, unit_price_lcu, unit_price_int_dollar) 
  
    glp1_rwe_cost_data_imp <- heirarchical_imputation(glp1_rwe_cost_data, var1 = "unit_price_int_dollar", var1_imp = "unit_price_int_dollar_imp",
                                                  out_fn = "semaglutide_rwprices_impute_summary")
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
                                                      out_fn = "sema_lcu-usd_conversionfactor_impute_summary")
    lcu_id_imputed <- heirarchical_imputation(lcu_per_id, var1 = "lcu_per_id", var1_imp = "lcu_per_id_imp",
                                               out_fn = "sema_lcu-id_conversionfactor_impute_summary")
    
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
        across(starts_with("price"), ~ if_else(is.na(.x), mean(.x, na.rm = TRUE), .x)
        ))
    if (!file.exists(file.path(out_folder, "costing_region_marketprice_vs_realworld_summary.xlsx"))){ 
      summary_region_mp_vs_rw <- glp1_cost_inputs |> 
        group_by(region) |> 
        summarise(
          mean_price_id_lowest_mp = mean(price_id_lowest_mp, na.rm = TRUE),
          mean_price_intl_rwe     = mean(price_intl_rwe, na.rm = TRUE),
          n_countries = n(),
          .groups = "drop"
        )
      write_xlsx(summary_region_mp_vs_rw, file.path(out_folder, "costing_region_marketprice_vs_realworld_summary.xlsx"))
      rm(summary_region_mp_vs_rw)
    }
    rm(lcu_id_imputed, lcu_usd_imputed, glp1_rwe_cost_data, lcu_per_id, lcu_per_usd,lcu_id_and_usd_imputed)
    
    glp1_cost_inputs <- glp1_cost_inputs |>  select(-region) 
    } else { glp1_cost_inputs = glp1_cost_inputs} 
  } else { glp1_cost_inputs = NULL}
