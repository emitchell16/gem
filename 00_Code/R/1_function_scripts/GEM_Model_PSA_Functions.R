#' Probabilistic Sensitivity Analysis (PSA) parameter draws and main function
#' March 2026
#' -------------------------

#' draw once at start for each sim run
draw_psa_params <- function(cohort_params, dw_tbl, glp1_impact_data, rr_tables) {
  get_draw <- function(value, lower, upper, distribution) {
    sample_rate(
      distribution = distribution,
      value = value,
      lower = lower,
      upper = upper,
      n = 1L)
  }
  classify_dist <- function(stem) {
    if (grepl("^le_yrs", stem) || grepl("^incidence_obesity", stem) || grepl("^incidence_class", stem) ) {
      "normal"
    } else {
      "logit-normal"
    }
  }
  
  # input data parameters for incidence, death, life expectancy ----
    # draw 1 parameter per group and save into new *_val column
  cohort_params_psa <- as.data.table(copy(cohort_params))
  val_cols <- grep("_val$", names(cohort_params_psa), value = TRUE)
  stems <- sub("_val$", "", val_cols)
  
  stems <- stems[paste0(stems, "_lower") %in% names(cohort_params_psa) &
                   paste0(stems, "_upper") %in% names(cohort_params_psa)]
  
  for (stem in stems) {
    vcol <- paste0(stem, "_val")
    lcol <- paste0(stem, "_lower")
    hcol <- paste0(stem, "_upper")
    dist <- classify_dist(stem)    
    cohort_params_psa[, g__ := .GRP, by = c(vcol, lcol, hcol)]
    
    draws__ <- cohort_params_psa[,
      .(draw__ = get_draw(
        value = get(vcol)[1L],
        lower = get(lcol)[1L],
        upper = get(hcol)[1L],
        distribution = dist
      )), by = g__]
    
    cohort_params_psa[draws__, (vcol) := i.draw__, on = "g__"]
    cohort_params_psa[, g__ := NULL]
  }

  # disability weights
  dw_tbl_psa <- as.data.table(copy(dw_tbl))
  dw_tbl_psa[, mean := mapply(
    function(v, l, u) get_draw(v, l, u, distribution = "logit-normal"),
    mean, lower, upper
  )]
  # treatment effects
  glp1_impact_data_psa <- as.data.table(copy(glp1_impact_data))
  
  glp1_impact_data_psa[, outcome_val := mapply(
    function(v, l, u, outcome) {
      if (is.na(l) || is.na(u)) return(v)
      dist <- if (grepl("bmi", outcome, ignore.case = TRUE)) {
        "normal"
      } else { "log-normal"}
      get_draw(v, l, u, distribution = dist)
    },
    outcome_val, outcome_lower, outcome_upper, outcome
  )]
  #  RR tables
  sample_rr_table <- function(dt) {
    dt_psa <- as.data.table(copy(dt))
    dt_psa[, value := mapply(
      function(v, l, u) get_draw(v, l, u, distribution = "log-normal"),
      value, lower, upper
    )]
    dt_psa
  }
  
  rr_tables_psa <- lapply(rr_tables, sample_rr_table)
  
  list(
    cohort_params = cohort_params_psa,
    dw_tbl = dw_tbl_psa,
    glp1_impact_data = glp1_impact_data_psa,
    rr_tables = rr_tables_psa)
}

run_psa <- function(R, psa_temp_folder, psa_out_folder) {  
  checkpoint_regional <- file.path(psa_out_folder, "regional_runs_partial.rds")
  if (!dir.exists(psa_out_folder)) {dir.create(psa_out_folder, recursive = TRUE, showWarnings = FALSE)}
  
  if (file.exists(checkpoint_regional)) {
    regional_psa <- readRDS(checkpoint_regional)
  } else {regional_psa <- tibble(run = integer())}
  
  done_runs <- unique(regional_psa$run)
  next_run   <- if (length(done_runs)) max(done_runs) + 1 else 1
  
  if (next_run > R) {
    message("All ", R, " PSA runs already completed; skipping.")
  } else{
    for (r in next_run:R) {
      cat("=== PSA run", r, "of", R, "===\n")
      # sample input data
      source(file.path(paths$r_code, "set_in_tables.R")) 
      
      
      # DEBUG: 
      cat("\n  ---- AC mortality incidence parameters: \n ")
      print(summary(all_cohort_param_rates$incidence_ac_death_val))
      
      cat("\n  ---- AC mortality RR parameters -- obesity by class: \n ")
      print(obesity_rr_acdeath_byclass)
      
      # run 2 scenario simulations
      main_sim_file_management(files_yr0 = sample_baseline_files,
                              out_folder_path = wf_withdate,
                              glp1_coverage = 0,
                              sim_years = 5)
      cat("  starting comparison...")
      main_sim_file_management(files_yr0 = sample_baseline_files,
                               out_folder_path = wf_withdate,
                               glp1_coverage = 100,
                               sim_years = 5)
      # pull summary results
      start_time <- Sys.time()
      cat("summarizing post-sim...")
      regional1 <- make_region_impact_table(
                                  5, 100, 
                                  basecase_root   = wf_withdate, 
                                  sim_hive_folder = wf_withdate,
                                  level = "all", psa=psa) |>
        mutate(run = r)
    
      regional_psa <- bind_rows(regional_psa, regional1)
      saveRDS(regional_psa, checkpoint_regional)
    
      # clean up temp files from this iteration
      unlink(psa_temp_folder, recursive = TRUE) 
      dir.create(psa_temp_folder, recursive = TRUE)
      gc(verbose = F)
      
      end_time <- Sys.time()
      elapsed <- difftime(end_time, start_time, units = "mins")
      cat("-- Post-sim summary elapsed time:", round(elapsed,2), "mins\n")
    }
  }
  
  psa_results <- regional_psa |>
    arrange(desc(who_region == "global"), who_region)
  
  # Compute 2.5–97.5 percentile intervals for each region×condition estimate
  psa_intervals <- psa_results |> 
    group_by(who_region, condition) |> 
    summarise(
      across(
        starts_with("outcome_"),
        list(
          low  = ~ quantile(.x, 0.025, na.rm = TRUE),
          high = ~ quantile(.x, 0.975, na.rm = TRUE))
      ),
      across(
        c(abs_diff_pct, rel_diff_pct),
        list(
          low  = ~ quantile(.x, 0.025, na.rm = TRUE),
          high = ~ quantile(.x, 0.975, na.rm = TRUE))
      ),
      .groups = "drop"
    ) |>
    mutate(across(where(is.numeric), ~round(.x, 2)),
           condition = case_when( condition == "class1_obesity"  ~ "obesity_class1",
                                  condition == "class2_obesity"  ~ "obesity_class2",
                                  condition == "class3_obesity"  ~ "obesity_class3",
                                  TRUE ~ condition
           )
    )|> 
    filter(condition != "t2d")
  
  # write out
  out_dir_subfolder <- file.path(psa_out_folder, "final")
  if (!dir.exists(out_dir_subfolder)) { dir.create(out_dir_subfolder, recursive = TRUE)}
  
  write_csv(psa_intervals, file.path(out_dir_subfolder, "psa_intervals_by_region_condition.csv"))
  write_csv(psa_results, file.path(out_dir_subfolder, "psa_results_all.csv"))
  
}
