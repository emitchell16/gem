# 00_Code/scripts/scenario_sem_2026.R

# load packages + get paths --------
source("00_Code/R/setup_environment.R")
load_gem_packages()

# load run_gem() and helpers
source("00_Code/R/run_gem.R")
source("00_Code/R/configs.R")
source("00_Code/R/prep_main_inputs.R")

source(file.path("00_Code/R/", "1_function_scripts", "GEM_Model_1_Baseline_Summary_Stats.R"))

# define configurations ------

        #TODO: edit here to alter default configs: (see:"00_Code/R/configs.R" for defaults + options)

.TESTRUN = TRUE   #set to FALSE to run all scenarios 

# settings for input demographic data 
input_data_config <- config_inputs(
  in_data_run_date = "3.1",
  min_age_model = 12L, 
  sample_prop   = 0.03)

# ----------------------------------'
#  ********** scenarios: ********** 
# ----------------------------------'
# no semaglutide, status quo
sim0_5y_config  <- config_simulation(
  sim_years = 5L,
  treatment_enabled = FALSE)

sim0_3y_config  <- config_simulation(
  sim_years = 3L,
  treatment_enabled = FALSE)
# -----------------------------'
# incremental scale up to 100% access scenarios
access_scenarios = c(0.1, 0.25, 0.5, 1)

scaleup_configs <- lapply(access_scenarios, function(cov) {
  config_simulation(
    treatment_enabled = TRUE,
    glp1 = list(coverage_rate = cov))
})
names(scaleup_configs) <- paste0("sim", access_scenarios * 100, "_5y_config")
sim100_5y_config <- scaleup_configs[["sim100_5y_config"]]

# -----------------------------'
# rwe treatment effect estimates in 100% scenario
sim100rwe_config <- config_simulation(
  treatment_enabled = TRUE,
  glp1 = list(treatment_effect_source = "rwe"))
# -----------------------------'
# 3-year sensitivity in 100% scenario
sim100_3y_config <- config_simulation(
  sim_year = 3L,
  treatment_enabled = TRUE)
# -----------------------------'
# collect
collect_sim_configs <- function(env = parent.frame()) {
  sim_names <- grep("^sim.*config$", ls(envir = env), value = TRUE)
  sim_configs <- mget(sim_names, envir = env)
}
all_configs <- collect_sim_configs()
# -----------------------------'


# run the model -----
configs_to_run <- if (.TESTRUN) {
  all_configs[c("sim0_5y_config", "sim100_5y_config")]
} else {
  all_configs
}

results <- run_gem(input_data_config, configs_to_run)
