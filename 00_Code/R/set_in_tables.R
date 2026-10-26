#' GEM simulation input data tables
#' sourced once per simulation config (run_gem); requires paths, cohort_input_data, config_sims
#' -------------------------

##############'
# Helpers ----
input_param_fp <- file.path(paths$data, "GEM_input_parameters.xlsx")
input_params <- read_excel(input_param_fp, sheet = "baseline_dataset") |> select(parameter, value)
model_param  <- setNames(as.numeric(input_params$value), input_params$parameter)    # named scalar parameters
# PSA intervals / distributions for scalar parameters (rows with lower & upper)
#  an intervention module can append its own parameters in hook_load_inputs()
model_param_ci <- read_excel(input_param_fp, sheet = "baseline_dataset") |>
  select(any_of(c("parameter", "lower", "upper", "psa_dist", "psa_family"))) |>
  filter(!is.na(lower), !is.na(upper), psa_dist %in% c("beta", "lognormal", "normal")) |>
  mutate(across(c(lower, upper), as.numeric))

# clean input data strings for matching
norm_key <- function(x) {
  x |> stri_replace_all_fixed("\u00A0", " ") |> str_squish() |>  str_to_lower()                           
}

# for costing data, provide estimated values for missing countries/territories ----
#  fill order = levels (income_region, income, region, global); group value = stat over observed countries,
#  used only if the group has >= min_n observed values (default: mean, no minimum, as before)
#  drug prices: median of income group, then region, then global; min_n = 3 (robust to outliers)
heirarchical_imputation <- function(ds, var1, var1_imp, out_fn, stat = mean,
                                    levels = c("income_region", "income", "region", "global"), min_n = 1L) {
  country_groups = read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")) |> 
    filter(who_region != "no who region")
  source_col <- paste0(var1_imp, "_source")
  by_level <- list(income_region = c("income_group", "region"), income = "income_group", region = "region", global = character(0))
  grp_val <- function(x) { x <- x[!is.na(x)]; if (length(x) >= min_n) stat(x) else NA_real_ }

  ds2 <- ds |>
    right_join( country_groups, by=c("iso3")) |> select(-X) |>
    rename(income_group = world_bank_income_group, region = who_region)
  ds2[[var1_imp]] <- ds2[[var1]]
  ds2[[source_col]] <- ifelse(is.na(ds2[[var1]]), NA_character_, "observed")
  for (lv in levels) {
    v <- ds2 |> group_by(across(all_of(by_level[[lv]]))) |> mutate(v__ = grp_val(.data[[var1]])) |> ungroup() |> pull(v__)
    fill <- is.na(ds2[[var1_imp]]) & !is.na(v)
    ds2[[var1_imp]][fill]  <- v[fill]
    ds2[[source_col]][fill] <- paste0(lv, "_", deparse(substitute(stat)))
  }
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
      n_imputed = sum(.data[[source_col]] != "observed", na.rm = TRUE),
      prop_imputed = n_imputed / n_countries,
      .groups = "drop"
    ) |> 
    mutate(
      mean_sd = case_when(!is.na(sd_est) ~ sprintf("%.0f (%.0f)", mean_est, sd_est), TRUE ~ sprintf("%.0f", mean_est)),
      median_iqr = sprintf("%.0f [%.0f–%.0f]", median_est, min_est, max_est),
      prop_imputed_pct = sprintf("%.1f%%", 100 * prop_imputed)
    ) |>
    select(region, income_group, n_countries, mean_sd, median_iqr, prop_imputed_pct)
  
  write_xlsx(summary_table, file.path(paths$out, paste0("costing_", out_fn,".xlsx")))
  
  ds2
}

##############'
# Incidence and death parameters ----
#  prevalence kept: GBD condition prevalence (excess mortality = CSMR / prevalence) and
#  NCD-RisC BMI categories (continuous BMI distribution; add_bmi_dist_params() in GEM_Model_BMI_RR.R)
source(file.path(paths$fn_scripts, "GEM_Model_BMI_RR.R"))
all_cohort_param_rates <- cohort_input_data |> 
  select(-starts_with("prob_death"), - pop_size) |>
  mutate(incidence_ac_death_upper = pmin(incidence_ac_death_upper, 1)) |>
  as.data.table() |>
  add_bmi_dist_params()

# Risk modifiers ----
.prop_stroke_ischaemic <- input_params$value[input_params$parameter == "prop_stroke_ischaemic"]

## * BMI (continuous, per-5 kg/m2; joint RRs) ----
#   age-specific: obesity_rr_nejm; other-cause death: bmi_rr_continuous
nejm_age_rr <- load_nejm_age_rr(input_param_fp, p_isch = .prop_stroke_ischaemic)
bmi_rr_tbl  <- load_bmi_rr(input_param_fp)
## * t2d comorbidities ---- 
t2d_rr                     = read_excel(input_param_fp, sheet = "t2d_rr") |>  select(-source)
## * stroke ----
cvd_rr_stroke_bysex        = read_excel(input_param_fp, sheet = "cvd_rr_stroke_bysex") |>  select(-source)
## * scalar model parameters used in risk equations (baseline_dataset sheet) ----
need_params <- c("eskd_inc_prob_total_pop", "eskd_excess_mort",
                 "stroke_recur_prob_yr1", "stroke_recur_prob_after", "cvd_recur_prob_yr1", "cvd_recur_prob_after")
if (any(is.na(model_param[need_params]))) stop("missing baseline_dataset parameters: ",
                                               paste(need_params[is.na(model_param[need_params])], collapse = ", "))

##############'
# Disability weights ----
dw_tbl <- fread(file.path(paths$working, "condition_disability_wts.csv"))

#################'
# PSA: parameters are drawn per run inside run_psa() (GEM_Model_PSA_Functions.R),
#   from the point-estimate tables loaded here (captured by psa_capture_base()); nothing drawn at load
##############'
# Health event costs -----
#'      sets 'health_event_cost_inputs'

# conventional 3% discount factor
disc_rate <- 0.03
# cost price year: health-state costs inflated to cost_year with US CPI-U medical care (sheet us_cpi_medical);
#  drug prices are current list/retail prices (2025-2026; price_year_min), not adjusted
cost_year      <- 2025L
price_year_min <- 2025L
inflate_cost <- function(x, from_year, to_year = cost_year) {
  cpi <- read_excel(input_param_fp, sheet = "us_cpi_medical")
  f <- cpi$cpi[match(to_year, cpi$year)] / cpi$cpi[match(from_year, cpi$year)]
  if (anyNA(f)) stop("us_cpi_medical: missing year(s) ", paste(unique(c(from_year, to_year)[is.na(match(c(from_year, to_year), cpi$year))]), collapse = ", "))
  x * f
}

if (config_sims$costing) { 
  if (!exists("health_event_cost_inputs") || is.null(health_event_cost_inputs)) {   # NULL left from a costing = FALSE run
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
     mutate(cost_intl = inflate_cost(cost_intl, price_year)) |>          # 2022 USD (PPP) -> cost_year
     select(-source, -"source notes", - region, -income_group, -price_year)
  intl_countrylevel_health_costs_obs <- intl_countrylevel_health_costs |>  
      pivot_wider(
        names_from = event_category,
        values_from = cost_intl,
        names_prefix = "cost_id_")
  
  ref_conds <- intl_countrylevel_health_costs |> filter(iso3 == "USA") |> 
    rename(cost_usd = cost_intl)
  #US cost estimates: (death/t2d/cvd/stroke (acute and chronic))
  usd_health_costs <- read_excel(input_param_fp, sheet = "us_health_costs") |> 
    mutate(cost_usd = inflate_cost(cost_usd, price_year)) |> select(-price_year) |>   # 2024 USD -> cost_year
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

##############'
# Intervention inputs (project module) -----
#  cost_total_cols: person-level cost columns that receive health-state costs in addition to cost_conditions_total
#  (none in the core model); hook_load_inputs() may set it and load further input tables
cost_total_cols <- character(0)
hook_load_inputs(config_sims)
