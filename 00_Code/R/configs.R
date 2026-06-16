# 00_Code/R/configs.R


# description ------
#' Purpose: provide defaults for simulation model configurations
#' Users can override any of the following fields:
#'
#' Core run settings:
#'   seed                Integer RNG seed
#'   sim_years           Simulation time horizon (years)
#'   treatment_enabled   TRUE/FALSE
#'   sample_prop         Proportion of each cohort sampled (0-1)
#'   psa                 TRUE/FALSE
#'
#' Input data years:
#'   UN_data_yr
#'   GBD_data_yr
#'   NCDRisC_data_yr

#'   
#' Model parameters:
#'   min_age_model       Minimum modeled age (>= 12)
#'   country_inclusion   "max" | "has_rate_data" (drop the 31 territories with UN pop data not in GBD or NCD-RisC datasets)
#'
#' GLP-1 settings (nested in $glp1):
#'   coverage_rate       Proportion treated (0-1)
#'   drug_type          "semaglutide"
#'   indication          "all" | "t2d_only" | "obese_only"
#'   treat_strategy      Selection strategy for targeted treatment (currently placeholder)
#'   adherence_vs_rct    Effect scaling factor (0-1)
#'   treatment_effect_source  "rct" | "rwe" for sensitivity 
#'   open_pop_cohort     Whether new people added to cohort each year as they age-in
#'   daly_sensitivity           Whether to use primary or sensitivity daly inputs/equations
#'
#' Example:
#' config <- config_base_glp1(
#'   seed = 5,
#'   treatment_enabled = TRUE,
#'   glp1 = list(coverage_rate = 0.25)
#' )

# define -----
config_inputs <- function(...) {
  defaults <- list(
    UN_data_yr            = 2024L,
    GBD_data_yr           = 2023L,  # life table, disability weights, incidence, prevalence, mortality
    NCDRisC_data_yr       = 2022L,
    country_inclusion     = "max", 
    min_age_model = 12L,
    sample_prop   = 0.1,
    
    seed = 1L,
    chunk_cap = 5e5,
    
    short_date = paste0(as.numeric(format(Sys.Date(), "%m")), ".", as.numeric(format(Sys.Date(), "%d"))) ,
    in_data_run_date = NULL     # last run date as needed
  ) 
  modifyList(defaults, list(...))
  }

config_simulation <- function(...) {
  defaults <- list(
    seed = 1L,
    short_date = paste0(as.numeric(format(Sys.Date(), "%m")), ".", as.numeric(format(Sys.Date(), "%d"))),
    sim_years = 5L,
    psa       = FALSE,
    daly_sensitivity = FALSE,                 # TRUE to use conservative lower bound
    open_pop_cohort = FALSE,                  # population replacement (e.g. for longer time horizon)
    treatment_enabled = TRUE,
    costing = FALSE,
    
    # GLP-1 treatment settings (only used if treatment_enabled=TRUE)
    glp1 = list(
      coverage_rate      = 1,          # 0-1 (e.g., 0.25 = 25% eligible treated)
      indication         = "all",      # "all" | "t2d only" | "obese_only"  
      treat_strategy     = "random",   # "random" | "age_risk" | "t2d_only" | "cvd_risk" etc target populations *(#TODO not currently coded) 
      adherence_vs_rct   = 1,          # 0-1; simple multiplier on effect or uptake                             *(#TODO not currently coded) 
      treatment_effect_source = "rct",  # rct vs rwe   
      drug_type    = "semaglutide",
      treatment_cost_scenario  = "price_intl_rwe"   # current options: "price_id_sust_costx1" | "price_id_sust_costx3" | "price_intl_rwe" (var names in glp1_cost_inputs table)
      )
  )
  
  modifyList(defaults, list(...))
}

print_config <- function(input_data_config, config_simulation) {
  cat("GEM configuration\n")
  cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("-----------------\n")
  cat("Seed:", config_simulation$seed, "\n")
  cat("Stratified random sample proportion:", input_data_config$sample_prop, "\n")
  cat("Years:", config_simulation$sim_years, "\n")
  cat("Open cohort/age replacement:", config_simulation$open_pop_cohort, "\n")
  cat("Treatment enabled:", config_simulation$treatment_enabled, "\n")
  if (config_simulation$treatment_enabled) {
    cat("GLP-1 coverage:", config_simulation$glp1$coverage_rate, "\n")
    cat("Indication:", config_simulation$glp1$indication, "\n")
    cat("Treatment strategy: ", config_simulation$glp1$treat_strategy,  "\n")
    cat("Treatment effect source: ", config_simulation$glp1$treatment_effect_source, "\n")
  }
}