# -------------------'
# Run 2-year test simulation for calibrating double counting adjustment
# Date updated: Feb. 2026
# -------------------'

# 1. File management wrapper function ----------
dct_main_sim_wrapper <- function(files_yr0 , glp1_access_percentage, out_folder_path) {   
  cat("Start time:", format(Sys.time(), "%c"), "\n")
  start_time <- Sys.time()
  
  on.exit({
    unlink(temp_dir, recursive = TRUE)
    rm(temp_dir)
    gc(verbose = FALSE)
  })
  
  # file directories --------------'
  temp_dir <- file.path(local_working_folder, "temp_simulation")
  unlink(temp_dir, recursive = TRUE) 
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  
  # create final output folder hive system
  group_umbrella_folder <- if (glp1_access_percentage > 0) {
    file.path(out_folder_path, "sim=glp1_access")
  } else {
    file.path(out_folder_path, "sim=basecase")
  }
  if (!dir.exists(group_umbrella_folder)) dir.create(group_umbrella_folder, recursive = TRUE)
  files <- list.files(path = files_yr0, pattern = "\\.parquet$", full.names = TRUE)
  
  # Process each file in files_yr0
  sim_results_by_group <- future_lapply(seq_along(files), function(i) {    
    # check first if cohort file already exists
    sim_endyear_folder <- file.path(group_umbrella_folder, paste0("year=", 2)) 
    if (!dir.exists(sim_endyear_folder)) {dir.create(sim_endyear_folder, recursive = TRUE)  }
    cohort_output_file <- file.path(sim_endyear_folder, paste0("cohort_", i, ".parquet"))
    if (file.exists(cohort_output_file)) {return(NULL)}
    
    temp_dir_i <- file.path(temp_dir, paste0("group=", i))
      if (!dir.exists(temp_dir_i)) dir.create(temp_dir_i, recursive = TRUE)
    group_temp_ds <- read_parquet(files[i]) 
    group_data <- as.data.table(group_temp_ds)
    rm(group_temp_ds)
   
    if (is.na(group_data[1, iso3])) {return(NULL)}
    
    if (nrow(group_data) > config$chunk_cap) {
      nchunks <- ceiling(nrow(group_data) / config$chunk_cap)
      group_chunks <- split(group_data, ceiling(seq_len(nrow(group_data)) / config$chunk_cap))
      cat(group_data[1, iso3], "processing chunks:",  nchunks, "\n" )
      # Process each chunk 
      chunk_results <- lapply(seq_along(group_chunks), function(j) {
        temp_dir_i_chunk <- file.path(temp_dir_i, paste0("tempchunk_", i, "_j", j))  
        if (!dir.exists(temp_dir_i_chunk)) dir.create(temp_dir_i_chunk, recursive = TRUE)
        
        res <- simulate_disease_progression(
          conditions = conditions_list,
          input_data = group_chunks[[j]],
          incidence_rates = all_cohort_param_rates,
          sim_temp_dir = temp_dir_i_chunk,
          glp1_access_percentage = glp1_access_percentage
        )
      })
      
      # Combine chunk results for the group by simulation year:
      group_result <- vector("list", 3L)
      for (year in seq_len(3)) {
        real_year = year - 1
        files_chunk <- lapply(chunk_results, function(res) res[[year]])
        files_vec <- as.character(unlist(files_chunk))
        
        combined_data <- rbindlist(lapply(files_vec, read_parquet), use.names = T, fill = T)
        out_fname <- paste0("sim_year", real_year, "_group", i, ".parquet")
        out_file <- file.path(temp_dir_i, out_fname)
        write_parquet(combined_data, out_file)
        group_result[[year]] <- out_file
        unlink(files_vec)
        rm( combined_data); gc()
      }
      res <- group_result
    } else {  
      res <- simulate_disease_progression(
        conditions = conditions_list,
        input_data = group_data,
        incidence_rates = all_cohort_param_rates,
        sim_temp_dir = temp_dir_i,
        glp1_access_percentage = glp1_access_percentage
      )
      rm(group_data); gc(verbose = FALSE)
    }
    
    # Recombine group-level results to segment level
    final_segment_files <- list()
    
    # write data to folder
    for (year in seq_along(res)) {
      real_year = year - 1
      files_vec <- as.character(res[[year]])
      group_data <- rbindlist(lapply(files_vec, arrow::read_parquet), use.names = T, fill = T)
      year_folder <- file.path(group_umbrella_folder, paste0("year=", real_year))
      if (!dir.exists(year_folder)) dir.create(year_folder, recursive = TRUE)
      group_file <- file.path(year_folder, paste0("cohort_", i,  ".parquet"))
      write_parquet(group_data, group_file)
      final_segment_files[[year]] <- group_file
      
      unlink(files_vec)
      rm(group_data); gc()
    }
    unlink(temp_dir_i, recursive = TRUE)
  }, future.seed = T)
  
  rm(sim_results_by_group); gc()
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  cat("Elapsed time:", elapsed, "mins\n")
}

# 2. Primary Model yearly update function ----------
simulate_disease_progression <- function(conditions = conditions_list, 
                                         input_data,  
                                         incidence_rates, 
                                         sim_temp_dir, # (not OneDrive due to locking issues)
                                         glp1_access_percentage,
                                         indication = "all"
) {
  
  # create y0 t2d-comorbidity correlations for joint probability adjustment
  corr_val_t2ob  <- get_corr(as.integer(input_data$t2d),  as.integer(input_data$obesity))
  corr_val_t2cvd <- get_corr(as.integer(input_data$t2d),  as.integer(input_data$cvd))
  corr_val_t2ckd <- get_corr(as.integer(input_data$t2d),  as.integer(input_data$ckd))
        
    # create correlation variables and assign to every row
    input_data[, `:=`(
      obese_t2d_corr = corr_val_t2ob,
      cvd_t2d_corr   = corr_val_t2cvd,
      ckd_t2d_corr   = corr_val_t2ckd
    )]
  # set placeholders for incident case indicators for non-reversible conditions
  input_data[, c("t2d_new", "cvd_new", "ckd_new", "eskd_new") := 0L] # stroke and ac_death start at 0 so they're effectively already _new
  
  if (glp1_access_percentage > 0) {
      glp1_impact_data <- read_excel(file.path(paths$data, "GEM_input_parameters.xlsx"), sheet = "semaglutide_treatment_effects") |>
        filter(source_type == config_sims$glp1$treatment_effect_source) |> 
        select(treated_pop,	age_group, outcome_unit,	outcome, outcome_val,	outcome_lower, outcome_upper) |> 
        group_by(treated_pop, age_group, outcome_unit, outcome)  |> 
        summarise(
          outcome_val   = mean(outcome_val,   na.rm = TRUE),
          outcome_lower = mean(outcome_lower, na.rm = TRUE),
          outcome_upper = mean(outcome_upper, na.rm = TRUE),
          n_rows        = dplyr::n(),
          .groups = "drop"
        )
      
      assign_glp1_eligible <- function( indication = config_sims$glp1$indication, input_data) {
        if (indication == "all") {
          input_data[(age >= 12) & (obesity == 1), treatment_eligible := 1L]
          input_data[(age >= 18) & (t2d == 1) & (bmi >= 27), treatment_eligible := 1L] # for semaglutide for weight management  
          
        } else if (indication == "t2d only") {
          input_data[(age >= 18) & (t2d == 1) & (bmi >= 27), treatment_eligible := 1L]
        } else if (indication == "obese only") { 
          input_data[(age >= 12) & (obesity == 1), treatment_eligible := 1L]
        }
        return(input_data)
      }
      
      # Set GLP-1 treatment indicators for baseline input data ------------
      input_data[, treatment_eligible := 0L]   
      input_data <- assign_glp1_eligible(indication=indication, input_data = input_data)  # function from model parameters
      input_data[, treated := 0L]
      input_data[treatment_eligible == 1, 
                 treated := as.integer(runif(.N) < glp1_access_percentage/100)]  # make number [0,100] a proportion
    
    input_data[, years_treated := 0L]
  }   
 
  yearly_data <- vector("list", length = 3)
  # Baseline data from function input:
  baseline_file <- file.path(sim_temp_dir, "grp_baseline_cohort.parquet")  
  write_parquet(input_data, baseline_file)
  yearly_data[[1]] <- list(baseline_file)
  
# Run simulation for 2 years -------------
  for (year in seq_len(2)) {  
    prev_year_file <- yearly_data[[year]][[1]] 
    pd <-  as.data.table(read_parquet(prev_year_file, use_memory_map = FALSE))  

    # Get cohort-level incidence rates ------------
    current_iso3 <- pd[1, iso3]
    cohort_params <- incidence_rates[incidence_rates$iso3 == current_iso3, ]

  # call function to update person-level disease statuses and age -------------
    updated_person_data <- update_cohort_dc_test(pd,      
                                       cohort_params,    
                                       conditions,
                                       year,             
                                       glp1_access = glp1_access_percentage,
                                       glp1_impact_data
                                      )  
    # --- Write Updated Data to Disk ---
    year_file <- file.path(sim_temp_dir, paste0("simulated_year_", year, ".parquet"))
    write_parquet(updated_person_data, year_file)
  
    # Save the file path for use in subsequent simulation years.
    yearly_data[[year + 1]] <- list(year_file)
    rm(year_file, pd, updated_person_data); gc(verbose = F)
  }
  return(yearly_data)
}

# 3. Update cohort function (dc_test version) ----------
# For double-counting calibration: GLP-1 access > 0 comparison just reduces obesity in this test

update_cohort_dc_test <- function(pd, cohort_params, conditions, year,
                          glp1_access, glp1_impact_data) {
  # Set up --------------'
  params_df   <- read_excel(file.path(paths$data, "GEM_input_parameters.xlsx"), sheet = "baseline_dataset") |> select(parameter, value) 
  new_state_vals <- copy(pd)                    # Initialize updated year data
  living_idx <- which(new_state_vals$ac_death == 0)
  
  if (glp1_access > 0) {
    treated_idx  <- which(pd[["treated"]][living_idx] == 1)
    new_state_vals[living_idx[treated_idx], years_treated := years_treated + 1] 
  }
  
  # match individuals with previous year (pd) age & sex in cohort_params
  setDT(cohort_params)
  setkey(cohort_params, female, age_low)
  suppressWarnings({   #  suppress "new name" warnings
    matched_params <- cohort_params[pd, 
                                    on        = .(female, age_low = age),
                                    roll      = TRUE,
                                    rollends  = c(TRUE, TRUE)  
    ]
  })  
  
  ####
  ## 2. Update condition incidence rates -------
  ####
  for (cond in c("ckd", "t2d", "stroke")) {  
    if (length(living_idx) == 0) break
    
    # pull person level cohort incidence rates
    inc_val   <- matched_params[[paste0("incidence_", cond, "_val")]]
    
    inc_probs <- inc_val[living_idx]
    inc_probs <- pmin(1, pmax(0, inc_probs))
    inc_probs <- calculate_adjusted_risk(cond, inc_probs, matched_params,
                                         pd, living_idx, inc_val)
    ### draw events -------
    pos_no_cond <- which(new_state_vals[[cond]][living_idx] == 0)
    idx         <- living_idx[pos_no_cond]
    
    if (length(idx) > 0) {
      p <- inc_probs[pos_no_cond]
      new_state_vals[idx, (cond) := rbinom(.N, 1, pmin(1, pmax(0, p)))]
    }
  }
  
  # update cvd:
  # note: stroke modeled as subtype of cvd. draw non-stroke cvd among those without prior cvd
  if (("cvd" %in% conditions) && length(living_idx) > 0) {
    inc_prob_cvd_total_base  <- matched_params[["incidence_cvd_val"]][living_idx]
    inc_probs <- calculate_adjusted_risk("cvd", inc_prob_cvd_total_base, matched_params,
                            pd, living_idx, inc_val)
    
    ### draw events -------
    pos_no_cvd_no_stroke <- which(new_state_vals[["cvd"]][living_idx] == 0 & new_state_vals[["stroke"]][living_idx] == 0)
    idx         <- living_idx[pos_no_cvd_no_stroke]
    if (length(idx) > 0) {
      p  <- inc_probs[pos_no_cvd_no_stroke]
      new_state_vals[idx, cvd := rbinom(.N, 1, pmin(1, pmax(0, p)))]
    }
    # if stroke this year --> cvd this year = 1 
    pos_no_cvd_have_stroke <- which(new_state_vals[["cvd"]][living_idx] == 0 & new_state_vals[["stroke"]][living_idx] == 1)
    idx         <- living_idx[pos_no_cvd_have_stroke]
    if (length(idx) > 0) {
      new_state_vals[idx, cvd := 1L]
    }
  }
  
  ######'  
  ## 3. Update obesity conditions statuses -----------
  #' assumes adjacent only transitions
  ######'  
  # annual transition probabilities for this cohort
  #' class incidence rates can have 3 competing risk steps (down, stay, or up)
  inc_obesity <- pmin(1, pmax(-1, matched_params[living_idx, incidence_obesity_val]))
  inc_c1 <- pmin(1, pmax(-1, matched_params[living_idx, incidence_class1_obesity_val]))
  inc_c2 <- pmin(1, pmax(-1, matched_params[living_idx, incidence_class2_obesity_val]))
  inc_c3 <- pmin(1, pmax(-1, matched_params[living_idx, incidence_class3_obesity_val]))
  
  # Initialize prior and current class (0=non,1,2,3)
  prior_class <- integer(length(living_idx))
  prior_class[pd$class1_obesity[living_idx] == 1] <- 1
  prior_class[pd$class2_obesity[living_idx] == 1] <- 2
  prior_class[pd$class3_obesity[living_idx] == 1] <- 3
  
  current_class <- prior_class
  
  # class 0 -> 1 or stay 
  i0   <- which(prior_class == 0)
  p_up0   <- pmax(0, inc_obesity[i0])
  u0      <- runif(length(i0))
  current_class[i0] <- ifelse(u0 < p_up0, 1, 0)
  
  # — class 1: can go down (to 0), stay, or up (to 2) —
  i1      <- which(prior_class == 1)
  p_down1 <- pmax(0, -inc_c1[i1])   # remission
  p_up1   <- pmax(0,  inc_c2[i1])   # progression
  sum1    <- p_down1 + p_up1
  over1   <- sum1 > 1
  p_down1[over1] <- p_down1[over1] / sum1[over1]
  p_up1[over1]   <- p_up1[over1]   / sum1[over1]
  p_stay1 <- 1 - p_down1 - p_up1
  u1      <- runif(length(i1))
  current_class[i1] <- ifelse(
    u1 < p_down1,            0,
    ifelse(u1 < p_down1+p_stay1, 1, 2)
  )
  
  # — class 2: can go down (to 1), stay, or up (to 3) —
  i2      <- which(prior_class == 2)
  p_down2 <- pmax(0, -inc_c2[i2])
  p_up2   <- pmax(0,  inc_c3[i2])
  sum2    <- p_down2 + p_up2
  over2   <- sum2 > 1
  p_down2[over2] <- p_down2[over2] / sum2[over2]
  p_up2[over2]   <- p_up2[over2]   / sum2[over2]
  p_stay2 <- 1 - p_down2 - p_up2
  u2      <- runif(length(i2))
  current_class[i2] <- ifelse(
    u2 < p_down2,              1,
    ifelse(u2 < p_down2+p_stay2, 2, 3)
  )
  
  # — class 3: can go down (to 2) or stay —
  i3      <- which(prior_class == 3)
  p_down3 <- pmax(0, -inc_c3[i3])
  u3      <- runif(length(i3))
  current_class[i3] <- ifelse(u3 < p_down3, 2, 3)
  
  # Identify who changed status and update data table
  changed_idx <- which(prior_class != current_class)
  new_state_vals[living_idx, `:=`(
    class1_obesity = as.integer(current_class == 1),
    class2_obesity = as.integer(current_class == 2),
    class3_obesity = as.integer(current_class == 3),
    obesity        = as.integer(current_class != 0)
  )]
  
  if (length(changed_idx) > 0) {
    global_idx <- living_idx[changed_idx]
    
    # Update BMI for those who changed class
    new_state_vals[global_idx[current_class[changed_idx] == 1], bmi := runif(.N, 30, 34.99)]
    new_state_vals[global_idx[current_class[changed_idx] == 2], bmi := runif(.N, 35.1, 39.99)]
    new_state_vals[global_idx[current_class[changed_idx] == 3], bmi := runif(.N, 40.1, 65)]
    
    # Update BMI for those who became non-obese
    remit_idx <- changed_idx[prior_class[changed_idx] == 1 & current_class[changed_idx] == 0]
    if (length(remit_idx)) {
      global_remit <- living_idx[remit_idx]
      new_state_vals[global_remit, bmi := runif(.N, 18.5, 29.99)]  # normal - overweight range
    }}
  
  ## 3.1 Update obesity status with GLP-1 treatment ---------
  if (glp1_access > 0 && length(treated_idx) > 0) { 
    treated_abs <- living_idx[treated_idx]
    years_treated_vec <- new_state_vals$years_treated[treated_abs]

      # For individuals with years_treated > 1, override with past year's values:
      restore  <- treated_abs[years_treated_vec > 1]
      if (length(restore) > 0) {
        new_state_vals[restore, `:=`( 
          obesity        = pd$obesity[restore],
          class1_obesity = pd$class1_obesity[restore],
          class2_obesity = pd$class2_obesity[restore],
          class3_obesity = pd$class3_obesity[restore],
          bmi            = pd$bmi[restore]
        )]
      }
      # implement lower BMI in year 1 of treatment only
      first_year_idx <- treated_abs[years_treated_vec == 1 & pd$bmi[treated_abs] >= 27]
      if ( length(first_year_idx) > 0 ) {
        new_state_vals <- update_obesity_transition_glp1(
          new_state_vals, pd, living_idx, first_year_idx, glp1_impact_data)
      } 
    }
  
  ## 3.2 Final alignment of obesity with class indicators ------
  new_state_vals[(class1_obesity == 1) | (class2_obesity == 1) | (class3_obesity == 1), obesity := 1]
  
  ######'
  ## 4. Update ESKD among people with CKD in previous year ---------
  ######'
  eskd_inc_prob_total_pop <- params_df$value[match("eskd_inc_prob_total_pop", params_df$parameter)] 
  ckd_prevalence <-  mean(pd[["ckd"]][living_idx], na.rm = TRUE)       # calculated from previous year data by cohort
  eskd_inc_prob_ckd <- eskd_inc_prob_total_pop/ ckd_prevalence
  
  eskd_update_idx <- which(pd[["ckd"]][living_idx] == 1 & pd[["eskd"]][living_idx] == 0)
  if (length(eskd_update_idx) > 0) {   
    global_idx <- living_idx[eskd_update_idx]
    new_state_vals[global_idx, eskd := rbinom(.N, 1, eskd_inc_prob_ckd)]
  }
  ######'
  ## 5. Incident case indicators updated for non-reversible conditions -------
  ######'
  nr_conds <- c("t2d", "cvd", "ckd", "eskd") 
  for (cond in nr_conds) {
    new_col <- paste0(cond, "_new")
    # if prior value is 0 and new value is 1, flag as 1
    new_state_vals[, (new_col) := ifelse(pd[[cond]] == 0 & get(cond) == 1, 1L, get(new_col))]
  }
  
  ######'
  ## 6. Update death -------
  ######'
  inc_val_death   <- matched_params[["incidence_ac_death_val"]]
  inc_probs_death <- pmin(1, pmax(0, inc_val_death[living_idx]))
  
  ## A. Apply risk adjustment for comorbidities, accounting for those accrued in past year -----
  inc_probs_death <- calculate_adjusted_risk("ac_death", inc_probs_death, matched_params,
                                             new_state_vals, living_idx, inc_val_death)
  
  ## B. Draw deaths ----
  new_state_vals$ac_death[living_idx] <- rbinom(length(living_idx), 1, pmin(1, pmax(0, inc_probs_death)))
  
  ######'
  # 7. Age increased one year among living -----------
  ######'
  new_state_vals[ac_death == 0, age := pmin(age + 1, 99)] # top coded at 99

  return(new_state_vals)
}

 