# examples/run_example.R

# description ------
# Minimal GEM run of the core model: small sample, five countries, 5 years, no intervention.
# Use it to check that the setup, input data and packages work before a full run.
# An intervention is added through a project module (see projects/semaglutide_global_access/).
# Run from the repo root. Input data (01_Data/) are described in the README.

# setup ------
.code <- if (dir.exists(here::here("00_Code"))) here::here("00_Code") else here::here()   # repo layout: 00_Code/ or flat
source(file.path(.code, "R", "setup_environment.R"))
load_gem_packages()
source(file.path(.code, "R", "configs.R"))
source(file.path(.code, "R", "run_gem.R"))
source_gem_modules(get_paths())                       # core only: hook_*() defaults (no intervention)

# configuration ------
input_data_config <- config_inputs(
  sample_prop  = 0.01,
  min_cohort_n = 50L,
  iso3_subset  = c("USA", "GBR", "BRA", "IND", "NGA"))

sim_configs <- list(core_5y = config_simulation(sim_years = 5L))

print_config(input_data_config, sim_configs$core_5y)

# run ------
# writes the simulation to the local temp folder (out_root); no result tables without an intervention arm
manifest <- run_gem(input_data_config, sim_configs)
print(manifest[, .(config, label, mode)])

# weighted outcomes at year 5, global ------
d <- as.data.table(arrow::open_dataset(file.path(manifest$out_root[1], paste0("sim=", manifest$label[1]), "year=5")) |>
                     dplyr::collect())
print(summarize_weighted(d, c("obesity", "t2d", "cvd", "ckd", "eskd", "ac_death", "yll_adj", "yld_adj_total")))
