# R/configs.R

# description ------
# This file defines default settings for:
#   1. Input data selection (config_inputs)
#   2. Simulation settings (config_simulation)
# Users can override any parameter by passing values when creating a configuration object.
# Unknown field names stop with an error (a misspelled field would otherwise be added silently and ignored).
# Key settings:
# - Data years (UN, GBD, NCD-RisC)
# - Simulation horizon and random seed
# - Cohort sampling proportion
# - DALY and costing sensitivity analyses
# - Hook for an intervention module (see "Intervention modules" below)
#
# Example:
# input_data_config <- config_inputs(
#   sample_prop = 0.20)
# simulation_config <- config_simulation(
#   sim_years = 10)
#
# Use print_config() to display the active configuration prior to model runs.

# config_inputs() fields ------
#   UN_data_yr, GBD_data_yr, NCDRisC_data_yr   input data years
#   country_inclusion   "max" | "has_rate_data" (drop the 31 territories with UN pop data not in GBD or NCD-RisC)
#   min_age_model       minimum modeled age (18; adults only)
#   sample_prop         sampled fraction of each country-sex-5yr cohort (0-1)
#   min_cohort_n        minimum sampled persons per cohort (weights w = pop / n)
#   iso3_subset         NULL = all countries; character iso3 codes for quick test runs
#   seed, chunk_cap, short_date, in_data_run_date

# config_simulation() fields ------
#   sim_years           time horizon (years)
#   psa                 FALSE = primary (deterministic-input) run; TRUE = probabilistic sensitivity analysis
#   psa_*               PSA settings (GEM_Model_PSA_Functions.R)
#   costing, daly_sensitivity, open_pop_cohort
#   mortality_calibration  TRUE = apply all-cause mortality calibration from validation (calibration_estimates.csv)
#   intervention        NULL = no intervention module; a project wrapper sets a name (e.g. "semaglutide")
#   intervention_arm    TRUE = simulate an intervention arm next to the comparison arm (needs a project module)

# Intervention modules ------
# The core model simulates a population without an intervention. A project adds an intervention (treatment, policy,
# program) by defining the hook_*() functions in GEM_Model_Intervention_Hooks.R in its own files, and by wrapping
# config_simulation() to add its settings. Example: projects/semaglutide_global_access/ (config_semaglutide()).

# helpers -----
# modifyList with name check (recursive for nested lists)
modify_checked <- function(defaults, args, where = "config") {
  bad <- setdiff(names(args), names(defaults))
  if (length(bad)) stop(where, ": unknown field(s) ", paste(bad, collapse = ", "),
                        ". Valid: ", paste(names(defaults), collapse = ", "))
  for (nm in names(args)) {
    if (is.list(defaults[[nm]]) && is.list(args[[nm]])) {
      args[[nm]] <- modify_checked(defaults[[nm]], args[[nm]], paste0(where, "$", nm))
    }
  }
  modifyList(defaults, args)
}
.short_date <- function() paste0(as.numeric(format(Sys.Date(), "%m")), ".", as.numeric(format(Sys.Date(), "%d")))

# define -----
config_inputs <- function(...) {
  defaults <- list(
    UN_data_yr            = 2024L,
    GBD_data_yr           = 2023L,  # life table, disability weights, incidence, prevalence, mortality
    NCDRisC_data_yr       = 2022L,
    country_inclusion     = "max",
    min_age_model = 18L,
    sample_prop   = 0.1,
    min_cohort_n  = 400L,     # min persons per country-sex-5yr cohort (weighted); see initialize_baseline()
    seed = 1L,
    chunk_cap = 5e5,
    short_date = .short_date(),
    in_data_run_date = NULL,    # last run date as needed
    iso3_subset = NULL          # test runs only, e.g. c("USA", "GBR", "BRA", "IND", "NGA")
  )
  cfg <- modify_checked(defaults, list(...), "config_inputs")
  if (cfg$min_age_model < 18) stop("min_age_model must be >= 18 (adult model)")
  cfg
}
# folder tag for a country subset (keeps subset baseline/sim files apart from full runs)
subset_tag <- function(cfg) if (is.null(cfg$iso3_subset)) "" else paste0("_iso-", paste(sort(cfg$iso3_subset), collapse = "-"))

config_simulation <- function(...) {
  defaults <- list(
    seed = 1L,
    short_date = .short_date(),
    sim_years = 5L,
    psa       = FALSE,
    # PSA settings (see GEM_Model_PSA_Functions.R)
    psa_runs          = 100L,    # parameter draws; check psa_convergence.csv
    psa_sample_prop   = 0,       # per-run sample fraction (small samples; noise removed by ANOVA)
    psa_min_cohort_n  = 50L,     # per-run min persons per country-sex-5yr cohort
    psa_batches       = 10L,     # person batches per run for within-run (Monte Carlo) variance
    psa_cost_cv       = 0.25,    # assumed CV for health-state costs (gamma; no CIs available)
    daly_sensitivity = FALSE,    # TRUE to use conservative lower bound
    open_pop_cohort  = FALSE,    # population replacement (e.g. for longer time horizon)
    costing = TRUE,
    mortality_calibration = TRUE,
    intervention     = NULL,     # set by a project wrapper (marker for hook_*() functions)
    intervention_arm = FALSE     # TRUE: intervention arm vs comparison arm
  )
  cfg <- modify_checked(defaults, list(...), "config_simulation")
  cfg$sim_years <- as.integer(cfg$sim_years)
  if (!is.numeric(cfg$sim_years) || cfg$sim_years < 1) stop("sim_years must be a positive integer")
  cfg
}

# print the active configuration before a run
print_config <- function(input_data_config, config_simulation) {
  cat("GEM configuration\n")
  cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("-----------------\n")
  cat("Input data years: UN", input_data_config$UN_data_yr, "| GBD", input_data_config$GBD_data_yr,
      "| NCD-RisC", input_data_config$NCDRisC_data_yr, "\n")
  cat("Countries:", if (is.null(input_data_config$iso3_subset)) "all" else paste(input_data_config$iso3_subset, collapse = ", "), "\n")
  cat("Seed:", config_simulation$seed, "\n")
  cat("Stratified random sample proportion:", input_data_config$sample_prop,
      "| min persons per cohort:", input_data_config$min_cohort_n, "\n")
  cat("Years:", config_simulation$sim_years, "| PSA:", config_simulation$psa, "\n")
  cat("Open cohort/age replacement:", config_simulation$open_pop_cohort, "| costing:", config_simulation$costing, "\n")
  cat("Mortality calibration:", config_simulation$mortality_calibration, "\n")
  if (!is.null(config_simulation$intervention)) {
    cat("Intervention module:", config_simulation$intervention, "| intervention arm:", config_simulation$intervention_arm, "\n")
    if (exists("hook_print_config")) hook_print_config(config_simulation)
  }
}

# scenario labels (output folder names) ----
#  comparison arm: <years>yr_basecase[_tags]; shared by configs with the same horizon and tags
#  intervention arm: hook_scenario_label(cs); a project module sets extra tags in hook_label_tags()
scenario_tags <- function(cs) {
  c(if (!isTRUE(cs$mortality_calibration)) "nocal", hook_label_tags(cs),
    if (isTRUE(cs$costing)) "c")                                       # costed arms never paired with uncosted arms
}
base_label_of <- function(cs) paste(c(paste0(cs$sim_years, "yr_basecase"), scenario_tags(cs)), collapse = "_")
scenario_label <- function(cs) if (isTRUE(cs$intervention_arm)) hook_scenario_label(cs) else base_label_of(cs)

# outcome lists for post-simulation summaries ----
#  prevalence outcomes: weighted mean among living at time point
#  cumulative outcomes: weighted mean over starting cohort (incl. those who died)
#  count outcomes: weighted sums (events, DALYs, costs; discounted)
conditions_prev_outcomes <- c("obesity", "class1_obesity", "class2_obesity", "class3_obesity",
                              "t2d", "cvd", "ckd", "eskd")
conditions_cum_outcomes  <- c("t2d_new", "cvd_new", "ckd_new", "eskd_new", "stroke_ever", "ac_death")
conditions_list_outcomes <- c(conditions_prev_outcomes, conditions_cum_outcomes)
event_count_outcomes     <- c("stroke_n", "cvd_event_n")
daly_list_outcomes       <- c("yll_adj", "yld_adj_total", "yld_t2d_adj_total", "yld_cvd_adj_total", "yld_ckd_adj_total")
cost_list_outcomes       <- c("cost_t2d_tot", "cost_death_tot", "cost_ckd_tot", "cost_eskd_tot", "cost_stroke_tot",
                              "cost_cvd_tot", "cost_conditions_total")   # + intervention cost columns added by a project module when costing
all_outcomes_list        <- c(conditions_list_outcomes, event_count_outcomes, daly_list_outcomes, cost_list_outcomes)
