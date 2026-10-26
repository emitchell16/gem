### 
#' Title: Validation study replication and comparison codes
#' Simulation uses the main-model update_cohort()
#' Author: Elizabeth Mitchell
###

# 1. set up -----------------
## to rerun from pre-calibration:
# unlink(file.path(paths$working, "calibration_estimates.csv"))

rm(list = ls())

.code <- if (dir.exists(here::here("00_Code"))) here::here("00_Code") else here::here()   # repo layout: 00_Code/ or flat
source(file.path(.code, "R", "setup_environment.R"))
load_gem_packages()

paths <- get_paths()
local_working_folder <- paths$local_temp 
validation_out_dir   <- file.path(paths$out, "validation1"); dir.create(validation_out_dir, recursive = TRUE, showWarnings = FALSE)
prep_inputs_path     <- file.path(paths$r_code, "prep_main_inputs.R")
run_date             <- paste0(as.numeric(format(Sys.Date(), "%m")), ".", as.numeric(format(Sys.Date(), "%d")))
source(file.path(paths$r_code, "configs.R"))

##############'
# set configs ----
config <- config_inputs(sample_prop   = 0.1)
      # can add param in_data_run_date = e.g. "3.2" to use baseline data generated on earlier date for efficiency
config_sims <- config_simulation()
min_age_model <- config$min_age_model
crn_stream_offset <- 0L
##############'


## validation inputs & outcome data with light cleaning -------
# load data
cohort_input_data <- read.csv(filePath(paths$working,  "pop_cc_cohort.csv")) 

# nts: used for input prevalence rates when missing from validation cohort data
region_income_country_groups <- read.csv(filePath(paths$working, "who_world_bank_country_groups_iso_cleaned.csv"))
conditions_list        <- c("t2d",  "obesity", "cvd", "ckd", "stroke", "eskd",  "ac_death")

## calibration parameter  
if (file.exists(file.path(paths$working, "calibration_estimates.csv"))) {
  acm_calibration <- read.csv(file.path(paths$working, "calibration_estimates.csv"), stringsAsFactors = FALSE)
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

# main model modules + input tables ----
source_gem_modules(paths)
source(file.path(paths$r_code, "set_in_tables.R")) 

# validation functions: ----
#  same state updates as the main model (update_cohort); validation-only: region rates, trial-arm HRs, CVD death
source(file.path(paths$v_code, "GEM_Validate1_baseline_data.R"))          # v_initialize_baseline()
source(file.path(paths$v_code, "GEM_Validate1_vUpdate_Cohort_States.R"))  # v_reference_consts(), v_prepare_baseline(), v_update_cohort()
source(file.path(paths$v_code, "GEM_Validate1_replicate_study.R"))        # gem_simulate(), generate_vresults()
source(file.path(paths$v_code, "GEM_Validate1_process_results.R"))        
source(file.path(paths$v_code, "GEM_Validate1_manuscript_fig.R"))


# 2. run function to generate validation results -----------------
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
write_csv(calib_out, file.path(paths$working, "calibration_estimates.csv"))

#---------------------- Rerun validation results after applying calibration ---------
acm_calibration <- read.csv(file.path(paths$working, "calibration_estimates.csv"), stringsAsFactors = FALSE)

generate_vresults(vres_outfile = "validation_summaries_2.csv")
process_results(vres_outfile = "validation_summaries_2.csv",out_subdir = (file.path(validation_out_dir, "2powercalib-nwgt")))
process_results_death(vres_outfile = "validation_summaries_2.csv",out_subdir = (file.path(validation_out_dir, "2pc-nwghtdeath")))

process_results_noMACE(vres_outfile = "validation_summaries_2.csv", out_subdir = (file.path(validation_out_dir, "2pc-nwghtdeath-nomace")))

# Panel A — pooled (pre-calibration)
new_outdir <- file.path(paths$out, "Results viz-validation"); dir.create(new_outdir, recursive = TRUE, showWarnings = FALSE)

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

# Panel C version without MACE and ESKD (ESKD not a main outcome; few, low-incidence points pull the fit line down)
make_panel_figs(
  vres_pre    = "validation_summaries_1.csv",
  vres_post   = "validation_summaries_2.csv",
  out_path_png = file.path(new_outdir, paste0("panelC_noMACE_noESKD_pre_post_", config$in_data_run_date ,".png")),
  panel = "C", exclude = c("MACE", "ESKD")
)

# Manuscript figure (2 panels, double-column size): A mortality before/after; B all outcomes except MACE before/after
make_calib_figure(
  vres_pre  = "validation_summaries_1.csv",
  vres_post = "validation_summaries_2.csv",
  out_path  = file.path(new_outdir, paste0("fig_validation_calibration_", config$in_data_run_date))   # .png + .pdf
)
