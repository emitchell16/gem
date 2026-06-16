### 
#' Title: Validation study replication and comparison codes
#' Date: August 2025
#' Updated: March 2026
#' Author: Liz Mitchell
###

#---------------------- ***1. set up*** -----------------
## to rerun from pre-calibration:
# unlink(file.path(working_folder, "calibration_estimates.csv"))

rm(list = ls())

suppressMessages({
  library(tidyverse)
  library(stringi)
  library(pbapply)
  library(writexl)
  library(readxl)
  library(openxlsx)
  library(R.utils)
  library(rlang)
  library(truncnorm) 
  library(here)     
  # efficiency libraries
  library(data.table)
  library(arrow)
  # for parallel processing:
  library(future.apply)
  plan(multisession, workers = parallel::detectCores() - 2)
  options(future.globals.maxSize = 5 * 1024^3)  
  print(nbrOfFreeWorkers())
  library(progressr)
  # for results plotting
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  library(patchwork)
  library(ggtext)
})

folder_base <- tryCatch({
  normalizePath('C:/Users/Liz/OneDrive - Emory University', mustWork = T)
}, error = function(e1) { 
  tryCatch({
    normalizePath('C:/Users/Liz/OneDrive - Emory', mustWork = T)
  }, error = function(e2) {stop("Both paths failed.")})
})
od_fp <- file.path(folder_base, "/1_Research proj/Hui Shao Group/Global impact of GLP-1RAs") 
shared_folder_root <-file.path(folder_base, "1_Research proj/GEM Model_shared")

code           <- file.path(od_fp, "00_Code_326")
working_folder <- file.path(od_fp, "02_Working_326")
out_folder     <- file.path(od_fp, "03_Out_326")

local_working_folder <- "C:\\Users\\Liz\\Research_Local" 
if (!dir.exists(local_working_folder)) {stop("Directory does not exist: ", local_working_folder)} 

# load initial parameters and helper functions --------------
paths <- list(
  code            = file.path(od_fp, "00_Code_326"),
  r_code          = file.path(shared_folder_root, "00_Code/R"),
  data            = file.path(shared_folder_root, "01_Data"),
  working         = file.path(shared_folder_root, "02_Working"),
  fn_scripts      = file.path(shared_folder_root, "00_Code/R/1_function_scripts"),
  out             = file.path(od_fp, "03_Out_326"), 
  v_code          = file.path(shared_folder_root, "00_Code/R/2_validation_simulation_codes")
)
validation_out_dir <- file.path(paths$out, "validation1"); dir.create(validation_out_dir, recursive = TRUE, showWarnings = FALSE)

prep_inputs_path         <- file.path(paths$r_code, "prep_main_inputs.R")
run_date                 <- paste0(as.numeric(format(Sys.Date(), "%m")), ".", as.numeric(format(Sys.Date(), "%d")))
source(file.path(paths$r_code, "configs.R"))

##############'
# set configs ----
config <- config_inputs(
  sample_prop   = 0.1,
  in_data_run_date = "3.2" )

config_sims <- config_simulation() 
##############'


## validation inputs & outcome data with light cleaning -------
# load data
cohort_input_data <- read.csv(filePath(paths$working,  "pop_cc_cohort.csv")) 

# nts: used for input prevalence rates when missing from validation cohort data
region_income_country_groups <- read.csv(filePath(paths$working, "who_world_bank_country_groups_iso_cleaned.csv"))
conditions_list        <- c("t2d",  "obesity", "cvd", "ckd", "stroke", "eskd",  "ac_death")

## calibration parameter  
if (file.exists(file.path(working_folder, "calibration_estimates.csv"))) {
  acm_calibration <- read.csv(file.path(working_folder, "calibration_estimates.csv"), stringsAsFactors = FALSE)
}else{acm_calibration <- NULL}

validation_data_input   <- read_excel(filePath(paths$data, "for_validation1", "GEM_data_validation_workbook.xlsx"),
                                      sheet = "DataInput")
names(validation_data_input) <- str_replace_all(tolower(names(validation_data_input)), " ", "_")
validation_data_input <- validation_data_input |> 
  mutate(futime = ceiling(as.numeric(`follow-up_years`))) |> 
  rename( 
    europe = european,
    "eastern mediterranean" = eastern_mediterranean,
    "western pacific" = western_pacific,
    "south-east asia" = "south-east_asia"
  )

comparison_outcome_data <- read_excel(filePath(paths$data, "for_validation1", "GEM_data_validation_workbook.xlsx"),
                                      sheet = "ComparisonOutcomes")
names(comparison_outcome_data) <- str_replace_all(tolower(names(comparison_outcome_data)), " ", "_")
comparison_outcome_data <- comparison_outcome_data |> 
  rename(ac_death = death, ckd_new = new_ckd, eskd_new = new_eskd, cvd_new = new_cvd, t2d_new = new_t2d)

# define other input tables
source(file.path(paths$r_code, "set_in_tables.R")) 

# define functions: ----

  # GEM functions and risk equation inputs used in _2.simulation_function.R
  # updated from primary GEM model codes, for the following:
  #' 1. add cvd_death updating in f1
  #' 2. change cohort incidence rates at iso3-age-sex level to be by region-sex-age
  #' 3. add if file exists, read ac_death_calibration.csv and apply      
            
source(file.path(paths$v_code, "GEM_Validate1_replicate_study.R"))        # def gem_simulate() which calls v_initialize_baseline; generate_vresults() which calls v_initialize_baseline() + gem_simulate()
source(file.path(paths$v_code, "GEM_Validate1_baseline_data.R"))          # def v_initialize_baseline()   

source(file.path(paths$v_code, "GEM_Validate1_vUpdate_Cohort_States.R")) 
source(file.path(paths$fn_scripts, "GEM_Model_Risk_Equations.R")) 

source(file.path(paths$v_code, "GEM_Validate1_process_results.R"))        
source(file.path(paths$v_code, "GEM_Validate1_manusript_fig.R"))


#---------------------- ***2. Run function to generate validation results*** -----------------
generate_vresults(vres_outfile = "validation_summaries_1.csv")

# check output:
validation_output_comp_data <- read_csv(file.path(validation_out_dir, "validation_summaries_1.csv"), show_col_types = FALSE)
print(validation_output_comp_data)

#---------------------- Summarize validation results, fig & table -----------------

reg_voutput_wide <- process_results(
  vres_outfile = "validation_summaries_1.csv",
  out_subdir = (file.path(validation_out_dir, "run1--notcalib")))
process_results_death(vres_outfile = "validation_summaries_1.csv", out_subdir = (file.path(validation_out_dir, "run1--notcalibdeath")))

process_results_noMACE(vres_outfile = "validation_summaries_1.csv", out_subdir = (file.path(validation_out_dir, "run1--notcalib-nomace")))

#---------------------- Calibration estimates ------------
log_safe <- function(x) log(pmax(x, .Machine$double.eps))
ac <- reg_voutput_wide %>%
  filter(outcome_label == "All-Cause Mortality") %>%
  filter(is.finite(ref_1y), is.finite(sim_1y), ref_1y > 0, sim_1y > 0)

# power calibration of mortality, weighted by sample size
ac_w_n <- ac |>
  left_join(validation_data_input |>  select(study_id, cohort_id, sample_size), by = c("study_id", "cohort_id")) |>
  mutate(
    w_raw = sample_size,
    w_cap = pmin(w_raw, quantile(w_raw, 0.95, na.rm = TRUE)),  # winsorize
    w     = sqrt(w_cap),                                              # √n weighting
    w     = w / median(w, na.rm = TRUE)
  )

fit_w <- lm(log_safe(ref_1y) ~ log_safe(sim_1y), data = ac_w_n, weights = ac_w_n$w)
ac_alpha <- unname(coef(fit_w)[1])
ac_beta  <- unname(coef(fit_w)[2])
ac_r2    <- summary(fit_w)$r.squared

calib_out <- bind_rows(outcome="ac_death", scale="annual_hazard",
  transform = "log h_ref = alpha + beta log h_sim (weighted by sample size)",
  alpha=ac_alpha,  beta=ac_beta,  r2=NA_real_, n=nrow(ac))
write_csv(calib_out, file.path(working_folder, "calibration_estimates.csv"))

#---------------------- Rerun validation results after applying calibration ---------
acm_calibration <- read.csv(file.path(working_folder, "calibration_estimates.csv"), stringsAsFactors = FALSE)

generate_vresults(vres_outfile = "validation_summaries_2.csv")
process_results(vres_outfile = "validation_summaries_2.csv",out_subdir = (file.path(validation_out_dir, "2powercalib-nwgt")))
process_results_death(vres_outfile = "validation_summaries_2.csv",out_subdir = (file.path(validation_out_dir, "2pc-nwghtdeath")))

process_results_noMACE(vres_outfile = "validation_summaries_2.csv", out_subdir = (file.path(validation_out_dir, "2pc-nwghtdeath-nomace")))

# Panel A — pooled (pre-calibration)
new_outdir <- file.path(out_folder, "Results viz-validation"); dir.create(new_outdir, recursive = TRUE, showWarnings = FALSE)

make_panel_figs(
  vres_pre    = "validation_summaries_1.csv",
  out_path_png = file.path(new_outdir, paste0("panelA_pooled_pre-v3_", config$in_data_run_date ,".png"))
)

# Panel B — ac_death pre vs post
make_panel_figs(
  vres_pre    = "validation_summaries_1.csv",
  vres_post   = "validation_summaries_2.csv",
  out_path_png = file.path(new_outdir, paste0("panelB_ac_death_pre_post-v2_", config$in_data_run_date ,".png")),
  panel = "B"
  )

# Panel C - without MACE; pooled pre vs post
make_panel_figs(
  vres_pre    = "validation_summaries_1.csv",
  vres_post   = "validation_summaries_2.csv",
  out_path_png = file.path(new_outdir, paste0("panelC_MACE_pre_post-v2_", config$in_data_run_date ,".png")),
  panel = "C"
)