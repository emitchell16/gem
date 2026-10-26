# projects/semaglutide_global_access/R/glp1_config.R

# description ------
# config_semaglutide(): config_simulation() plus GLP-1 settings for the semaglutide global access project.
# Core fields keep their meaning (sim_years, psa, costing, ...); intervention = "semaglutide" switches on the
# GLP-1 hooks. Unknown field names stop with an error.
#
# Example:
# simulation_config <- config_semaglutide(
#   sim_years = 10,
#   glp1 = list(coverage_rate = 0.25))
#
# Fields added here:
#   treatment_enabled   FALSE = comparison arm only (no policy-driven initiation; baseline users continue)
#   baseline_use        prevalent GLP-1 RA use at year 0, both arms: "base" | "low" | "high" (obesity 10/5/15% in HIC) | "none"
#   baseline_use_mode   "maintain" (baseline users stay on drug) | "dynamics" (same discontinuation/restart as new users)
#   glp1 (nested; glp1_access_policy.R):
#     strategy          "random" | "tiered" | "capacity"
#     coverage_rate     0-1: share of eligible with access (random: every tier; capacity: filled top tier first, per country)
#     tiers_open        tiered: priority tiers with access (1 = T2D BMI>=30 ... 6 = BMI 30-34.9)
#     tier3_bmi27_cvd   TRUE: tier 3 also includes no-T2D BMI 27-29.9 with CVD (CDA-AMC, CMS Bridge)
#     age_max_start     no start or restart at this age or older (85: Wegovy SmPC, limited experience >= 85; Inf = none)
#     uptake            "base" | "optimistic" (sheet glp1_uptake) | "full" (all with access start in year 1) | numeric 0-1
#     treatment_effect_source  "rct" | "rwe"
#     drug_type         "semaglutide" (only option currently)
#     treatment_cost_scenario  price column in glp1_cost_inputs

glp1_defaults <- list(
  strategy        = "random",
  coverage_rate   = 1,
  tiers_open      = 1:6,
  tier3_bmi27_cvd = FALSE,
  age_max_start   = 85,
  uptake          = "base",
  treatment_effect_source = "rct",
  drug_type       = "semaglutide",
  treatment_cost_scenario = "price_intl_rwe"
)

config_semaglutide <- function(..., glp1 = list(), treatment_enabled = TRUE,
                               baseline_use = "base", baseline_use_mode = "maintain") {
  cfg <- config_simulation(..., intervention = "semaglutide", intervention_arm = isTRUE(treatment_enabled))
  cfg$treatment_enabled <- isTRUE(treatment_enabled)
  cfg$baseline_use      <- baseline_use
  cfg$baseline_use_mode <- baseline_use_mode
  cfg$glp1              <- modify_checked(glp1_defaults, glp1, "config_semaglutide$glp1")
  if (cfg$glp1$coverage_rate > 1) stop("glp1$coverage_rate is a proportion (0-1), e.g. 0.25 for 25%")
  if (cfg$glp1$drug_type != "semaglutide") stop("only semaglutide is available in the current model")
  if (!is.numeric(cfg$glp1$age_max_start) || cfg$glp1$age_max_start <= 18) stop("glp1$age_max_start must be a number > 18 (Inf = no limit)")
  if (!cfg$baseline_use %in% c("base", "low", "high", "none")) stop("baseline_use must be one of: base, low, high, none")
  if (!cfg$baseline_use_mode %in% c("maintain", "dynamics")) stop("baseline_use_mode must be 'maintain' or 'dynamics'")
  cfg
}
