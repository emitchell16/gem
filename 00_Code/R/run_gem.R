# 00_code/R/run_gem.R
# 
# ### debug
config_init = input_data_config
# config_sims = sim0_5y_config

run_gem <- function(config_init, config_sims) {
  #***************************
  # -------1. Inputs ----------
  #*************************** 
  paths      <- get_paths()   
  cohort_inputs <- prep_main_inputs(paths, config_init)                  # done 3/1  
  cohort_inputs <- cohort_inputs |> filter(pop_size > 10)     # removes a couple old age cohorts with pop size 0-4
  
  #***************************
  # # -------2. Initialize baseline ----------
  #***************************
  baseline_state <- initialize_baseline(cohort_inputs, paths, config_init) 

        ### summary stats:  ----
  if (!is.null(config_init$in_data_run_date)) { 
    .date = config_init$in_data_run_date} else {.date = config_init$short_date}
  y0_summary_filename  <- paste0("summary_stats_y0fmt_", .date, ".csv")
  
        if (file.exists(file.path(paths$out, y0_summary_filename))) {message("already generated; not re-running")
          } else {
            global_summary_conds <- c("female", "obesity", "class1_obesity", "class2_obesity", "class3_obesity", "cvd", "ckd", "eskd", "t2d")
            global_summary_stats(config =config_init, date=.date, paths)      # from "1_function_scripts\Global_GLP1_Model_1_Baseline_Global_Summary_Stats.R"
            rm(global_summary_conds, y0_summary_filename)
            gc(verbose=F)
          }
    ## apply sampling:  ----
    sample_baseline_fp <- sample_baseline_data() #returns file path location of sampled cohort hive files

  #***************************
  # # -------3. Prep for simulation ---------- 
  #***************************
  #TODO rerun calibration *post* feb 2026 updates
  # # calibration parameters   
  # if (file.exists(file.path(paths$working, "calibration_estimates.csv"))) {
  #   acm_calibration <- read.csv(file.path(paths$working, "calibration_estimates.csv"), stringsAsFactors = FALSE)
  #   } else{ acm_calibration <- NULL}
  
  ## helpers ------
  all_cohort_incidence_rates <- cohort_inputs |> select(-starts_with("prevalence_"), - pop_size) 
  conditions_list <- c("t2d",  "obesity", "cvd", "ckd", "stroke", "eskd",  "ac_death")
  wf_withdate <- file.path(paths$local_temp, paste0("file_sim_", config_init$sample_prop, "subsample_", config_init$short_date)) 
  if (!dir.exists(wf_withdate)) dir.create(wf_withdate, recursive = TRUE)
  
  
  if (config_sims$treatment_enabled == TRUE) { 
    ### Generate treatment impact double counting adjustment calibration
    
    }
  # # unlink(dc_test_fp, recursive = T)
  # dc_test_fp <- file.path(working_folder, paste0("file_sim_dctest0.1", "subsample_", "9.2")) 
  # script_path <- filePath(code, "Global_GLP1_Model_2.2_Calculate_Double_Count_Adjustment.R")
  # 
  # ## Check if result file already exists, don't run:
  # if (dir.exists(dc_test_fp)) {
  #   message("already generated; not re-running")
  # } else {
  #   dir.create(dc_test_fp, recursive = TRUE)
  #   
  #   source(filePath(code, "Global_GLP1_Model_2.1_Risk_Equation_Data_Processing.R"))
  #   source(filePath(code, "Global_GLP1_Model_2.3_Update_Obesity_Status_With_Treatment.R"))
  #   source(filePath(code, "Global_GLP1_Model_3.1_Sim_Yearly_Update.R"))
  #   source(filePath(code, "Global_GLP1_Model_3.2_Sim_State_Transitions.R"))
  #   source(filePath(code, "Global_GLP1_Model_3.2_Risk_Equations.R"))
  #   source(script_path)
  # }
  
  #***************************
  # # -------4. Run simulation ---------- 
  #***************************
  
  # for (c in config_sims) {
  #   sim_results <- run_simulation(baseline_state, inputs, paths, c)
  # }
  # 
  # 
  # # -------5. Summarize output ----------
  # summary <- summarize_outputs(sim_results, paths, config)
  # 
  # return(summary)
  return(NULL)
}





#******NOTE TO SELF:
#*       + flow = edit scenario, run scenario, scenario pulls source() over all modules and then calls run_gem() which puts functions together