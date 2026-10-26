# projects/semaglutide_global_access/run_scenarios.R
# description ------
# Scenario set for the global semaglutide access analysis (13 primary configs + 1 PSA config).
# Edit scenarios here. run_gem(): inputs -> baseline sample -> simulations (primary or PSA, per config) -> post-sim tables
# Defaults and options: R/configs.R (core), R/glp1_config.R (GLP-1). Access tiers and uptake: projects/semaglutide_global_access/R/glp1_access_policy.R
# GLP-1 code (hooks into the core model) lives in this project's R/ folder and is passed to run_gem(project_modules = ).
# Run from the repo root (open the project there, or setwd()). Input data: see README.

# load packages + run_gem() --------
.code <- if (dir.exists(here::here("00_Code"))) here::here("00_Code") else here::here()   # repo layout: 00_Code/ or flat
source(file.path(.code, "R", "setup_environment.R"))
load_gem_packages()
source(file.path(.code, "R", "configs.R"))
source(file.path(.code, "R", "run_gem.R"))
proj_modules <- file.path(.code, "projects", "semaglutide_global_access", "R")   # GLP-1 module (hook functions)
source(file.path(proj_modules, "glp1_config.R"))                                 # config_semaglutide()
source_gem_modules(get_paths(), proj_modules)                                    # definitions only; run_gem() sources them again

.TESTRUN     <- FALSE   # TRUE: small test (5 countries, 2 configs); FALSE: all scenarios
run_primary  <- TRUE    # primary scenarios (status quo, tier steps, capacity, sensitivities)
run_psa_runs <- FALSE   # probabilistic sensitivity analysis (long)

# input data + baseline sample (adults 18+) ------
input_data_config <- config_inputs(
  in_data_run_date = "10.20",
  sample_prop   = 0.05,                                  # production runs: 0.03-0.1
  min_cohort_n  = if (.TESTRUN) 50L else 400L,            # test: smaller cohorts
  iso3_subset   = if (.TESTRUN) c("USA", "GBR", "BRA", "TON", "NGA") else NULL)   # test: 5 countries

# ----------------------------------'
#  ********** scenarios: **********
# ----------------------------------'
# status quo (no new initiation; baseline GLP-1 users continue, baseline_use = "base") -----
sim0_5y_config <- config_semaglutide(sim_years = 5L, treatment_enabled = FALSE)

#' primary uptake: sheet glp1_uptake "optimistic" (base x 2.5, Mailhac 2024) -- data-anchored;
#' intervention scenarios assume expanded access raises uptake above real-world levels
upt_primary <- "optimistic"

# A. coverage expansion by priority tier (primary): tiers 1..k open, primary uptake -----
#    1 T2D BMI>=30 | 2 T2D BMI 25-29.9 + CVD/CKD | 3 no T2D BMI>=30 + CVD/CKD | 4 BMI>=40 | 5 BMI 35-39.9 | 6 BMI 30-34.9
tier_configs <- setNames(
  lapply(1:6, function(k) config_semaglutide(glp1 = list(strategy = "tiered", tiers_open = 1:k, uptake = upt_primary))),
  paste0("tier1to", 1:6))
# full coverage (tiers 1-6) at full uptake: hypothetical maximum (all with access start in year 1)
tier_max <- list(tier1to6_full = config_semaglutide(glp1 = list(strategy = "tiered", tiers_open = 1:6, uptake = "full")))

# B. capacity-constrained access (fixed share of eligible, top tier first) -----
access_shares <- c(0.1, 0.25)
capacity_configs <- setNames(lapply(access_shares, function(cv)
  config_semaglutide(glp1 = list(strategy = "capacity", coverage_rate = cv, uptake = upt_primary))),
  paste0("capacity", access_shares * 100))

# C. sensitivity (full coverage, tiers 1-6) -----
sens_configs <- list(
  tier1to6_rwe     = config_semaglutide(glp1 = list(strategy = "tiered", tiers_open = 1:6, uptake = upt_primary,
                                                   treatment_effect_source = "rwe")),
  tier1to6_3y      = config_semaglutide(sim_years = 3L, glp1 = list(strategy = "tiered", tiers_open = 1:6, uptake = upt_primary)),
  # mortality calibration off: status_quo_summary.csv compares 5yr_basecase vs 5yr_basecase_nocal (obesity, deaths, T2D incidence)
  tier1to6_nocal   = config_semaglutide(mortality_calibration = FALSE, glp1 = list(strategy = "tiered", tiers_open = 1:6, uptake = upt_primary))
)

# D. PSA (run separately; long) -----
psa_configs <- list(
  psa_tier1to6 = config_semaglutide(psa = TRUE, glp1 = list(strategy = "tiered", tiers_open = 1:6, uptake = upt_primary))
)

# collect (13 configs; status quo arms added automatically per horizon/calibration) -----
primary_configs <- c(list(sim0_5y = sim0_5y_config), tier_configs, tier_max, capacity_configs, sens_configs)

# run the model -----
configs_to_run <- if (.TESTRUN) primary_configs[c("sim0_5y", "tier1to6")] else primary_configs
print_config(input_data_config, primary_configs[[1]])
if (run_primary)  manifest     <- run_gem(input_data_config, configs_to_run, project_modules = proj_modules)
if (run_psa_runs) manifest_psa <- run_gem(input_data_config, psa_configs, project_modules = proj_modules)
