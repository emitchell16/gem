# -------------------
# Call simulation model and update yearly data
# Author: Liz Mitchell
# Date updated: March 2026
# -------------------'

# helper function to return correlations
get_corr <- function(x, y) {
  if (any(is.na(x)) || any(is.na(y))) {
    cat("Warning: NA values found in input vectors.\n")
  }
  sdx <- sd(x, na.rm = TRUE)
  sdy <- sd(y, na.rm = TRUE)
  # If either sd is NA or 0, return 1
  if (is.na(sdx) || is.na(sdy) || sdx == 0 || sdy == 0) {return(1)}
  r <- cor(x, y, use = "pairwise.complete.obs")  # Pearson correlation = phi coefficient for bin
  if (is.na(r)) {return(0)
  } else {return( pmin(pmax(r, 0), 1))}
}

# GLP1 (semaglutide) treatment indications
assign_glp1_eligible <- function( indication, input_data) {
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

#######'
# Primary Model wrapper function --------------------------------------------------------
#######'
main_sim_file_management <- function(files_yr0,
                                     out_folder_path,      
                                     glp1_coverage,
                                     sim_years,
                                     cs = config_sims) { 
  chunk_size=config$chunk_cap      # if break, make config a function parameter
  
  cat("Start time:", format(Sys.time(), "%c"), "\n")
  start_time <- Sys.time()
  on.exit({gc(verbose = FALSE) })
  
  # file directories --------------
  temp_dir <- file.path(local_working_folder, "temp_simulation")
   unlink(temp_dir, recursive = TRUE) # remove any prior
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  
  # create output folder hive system
  group_umbrella_folder <- if (glp1_coverage > 0) {
    file.path(out_folder_path, paste0("sim=", sim_years, "yr_", glp1_coverage, "glp1"))
  } else {file.path(out_folder_path, paste0("sim=", sim_years, "yr_basecase"))}
  if (!dir.exists(group_umbrella_folder)) dir.create(group_umbrella_folder, recursive = TRUE)
  
  # -----------'
  files <- list.files(path = files_yr0, pattern = "\\.parquet$", full.names = TRUE)
  file_count <- length(files)
  # cat("Total number of files:", file_count, "\n") 
                
                # # small set for testing: ###############
                # files <- head(files, 20)   # debugging runs only
                # cat("Testing with 20 cohort files only\n")
                # # ###########'
                
  # Process each file in files_yr0
  # sim_results_by_group <- pblapply(seq_along(files), function(i) {      # change files to test files and vice versa
  sim_results_by_group <- future_lapply(seq_along(files), function(i) {  # parallelized version, faster 
    # cat("select cohorts only....\n")
    # sim_results_by_group <- pblapply(c(1781, 1782), function(i) {          # catch to rerun a cohort if code breaks in the middle
    # sim_results_by_group <- future_lapply(c(20:40), function(i) { 
   
     # check first if cohort sim file already exists (useful if crashed) -----
    for (yr in c(0, sim_years)) {
      dir.create(file.path(group_umbrella_folder, paste0("year=", yr)), recursive = TRUE, showWarnings = FALSE)
      }
    cohort_output_file <- file.path(group_umbrella_folder, paste0("year=", sim_years), paste0("cohort_", i, ".parquet"))
    if (file.exists(cohort_output_file)) {return(NULL)}
    
    temp_dir_i <- file.path(temp_dir, paste0("group=", i))
    if (!dir.exists(temp_dir_i)) dir.create(temp_dir_i, recursive = TRUE)
    
    group_temp_ds <- read_parquet(files[i])  # for non test run 
    # group_temp_ds <- read_parquet(files[i]) |> head(1e4) ## debug test run only; do first n rows of each file
    group_data <- as.data.table(group_temp_ds)
    if (cs$open_pop_cohort){group_data[, pid := as.character(pid)]} else {group_data[, pid := NULL]}
    rm(group_temp_ds)

        # ## debug test run only #####
        # if (!"iso3" %in% names(group_data) || !any(group_data$iso3 %in% c("ASM", "USA", "MEX"))) {
        #   return(NULL)  # Return early if iso3 is not ASM or USA
        # }
        # #########'
  
    if (is.na(group_data[1, iso3])) {return(NULL) }
    
    # define years of results to keep -------
      idx_in_list   <- c(1, 2)
      real_years    <- c(0, sim_years)
      
    if (nrow(group_data) > chunk_size) {
      N       <- nrow(group_data)
      nchunks <- ceiling(N / chunk_size)
      group_chunks <- split(seq_len(N), ceiling(seq_len(N) / chunk_size))
        # cat(group_data[1, iso3],i, "cohort, processing chunks:",  nchunks, "\n" ) # toggle noisy output
  
      # simulate_disease_progression() over each chunk
      chunk_results <- vector("list", nchunks)
      for (j in seq_len(nchunks)) {
        start <- (j - 1) * chunk_size + 1L
        end   <- min(j * chunk_size, N)
        chunk_data <- group_data[start:end] # slice
        
        temp_dir_i_chunk <- file.path(temp_dir_i, paste0("tempchunk_", i, "_j", j))  # needed so chunk pieces don't overwrite before combined
           if (!dir.exists(temp_dir_i_chunk)) dir.create(temp_dir_i_chunk, recursive = TRUE)
        
        # check if already run
        chunk_file <- file.path(temp_dir_i_chunk, paste0("simulated_year_", sim_years, ".parquet"))
        if (file.exists(chunk_file)) {chunk_results[[j]] <- chunk_file; next}
        
        chunk_results[[j]] <- simulate_disease_progression(input_data = chunk_data, sim_temp_dir = temp_dir_i_chunk,  
                                                           sim_years, glp1_coverage)
      }
      # Combine results for the group by year:
      group_result <- vector("list", length(idx_in_list))  
      
      for (k in seq_along(idx_in_list)) {
        files_chunk <- unlist(lapply(chunk_results, `[[`, idx_in_list[k]))
        files_vec <- as.character(unlist(files_chunk))
        combined_data <- open_dataset(files_vec, format = "parquet") |>  collect()
                         # if IO error: rbindlist(lapply(files_vec, read_parquet),   use.names = T, fill = T)
        out_file <- file.path(temp_dir_i, paste0("sim_year", real_years[k], "_group", i, ".parquet"))
        write_parquet(combined_data, out_file)
        group_result[[k]] <- out_file
        unlink(files_vec)
        rm(combined_data); gc()
      }
      res <- group_result
      
    } else {
      res <- simulate_disease_progression(input_data = group_data, sim_temp_dir = temp_dir_i,
                                          sim_years, glp1_coverage)
      res <- res[idx_in_list]
      rm(group_data); gc(verbose = FALSE)
    }
    
    # Recombine group-level results to segment level
    final_segment_files <- vector("list", length(res))
    for (k in seq_along(res)) {
      files_vec <- res[[k]] 
      group_data <- open_dataset(files_vec, format = "parquet") |>  collect()
      year_folder <- file.path(group_umbrella_folder, paste0("year=", real_years[k]))
        if (!dir.exists(year_folder)) dir.create(year_folder, recursive = TRUE)
      group_file <- file.path(year_folder, paste0("cohort_", i,  ".parquet"))
      write_parquet(group_data, group_file)
      
      final_segment_files[[k]] <- group_file
      unlink(files_vec)
      rm(group_data)
    }
    unlink(temp_dir_i, recursive = TRUE)
   }, future.seed=config$seed)
 # }) # for pblapply instead of future
  
  rm(sim_results_by_group); gc()
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  cat("Elapsed time:", round(elapsed,2), "mins\n")
}

#######'
# Yearly update function --------------------------------------------------------
#######'
simulate_disease_progression <- function(input_data,     # for the sex-country cohort
                                         sim_temp_dir,   # temp directory
                                         sim_years,
                                         glp1_coverage, # (0-100) 
                                         conditions = conditions_list, 
                                         cs = config_sims) {
  
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
  # set placeholders
  input_data[, c(
    "t2d_new", "cvd_new", "ckd_new", "eskd_new",            # note stroke/ac_death/ESKD start at 0
    "stroke_this_year", "death_year")    
    := 0L]
  input_data[, c("yll", "yll_adj", "yld_total", "yld_ckd_total", "yld_cvd_total", "yld_t2d_total",
                 "yld_adj_total", "yld_ckd_adj_total", "yld_cvd_adj_total", "yld_t2d_adj_total")
             := 0.0]
  if (config_sims$costing) {
    input_data[, c("cost_t2d_tot", "cost_death_tot", "cost_ckd_tot", "cost_eskd_tot", "cost_stroke_tot", "cost_cvd_tot", 
                   "cost_conditions_total", cost_treat_cols, cost_total_cols ) 
               := 0.0]
  }
  if (glp1_coverage > 0) {
    # Read in GLP-1 drug treatment impacts data table
    if (cs$glp1$drug_type == "semaglutide") {
      # Set GLP-1 treatment indicators for baseline input data ------------
      input_data[, `:=`(
        treatment_eligible        = 0L,
        treatment_access_checked  = 0L,  # has this person ever been evaluated for access
        treated                   = 0L,
        years_treated             = 0L
      )]
      
      input_data <- assign_glp1_eligible(indication=cs$glp1$indication, input_data = input_data)  # function from model parameters
      input_data[treatment_eligible == 1,`:=`(
                   treated       = as.integer(runif(.N) < glp1_coverage / 100),
                   treatment_access_checked         = 1L
                 )]
    } else {warning("Only semaglutide available in current model")}
  }
  #save accessible temp data file of youngest cohort at baseline to clone in open cohort
  if (cs$open_pop_cohort && min(input_data$age) <= min_age_model) {
    youngest_y0_person_data <- input_data[age == min_age_model] 
    fwrite(youngest_y0_person_data, file.path(sim_temp_dir, "youngest_y0_person_data.csv"))
    rm(youngest_y0_person_data); gc(verbose = F)
  }
  
  ############'
  yearly_data <- vector("list", length = 2) # only save baseline and final years
  
  # Baseline data from function input:
  baseline_file <- file.path(sim_temp_dir, "grp_baseline_cohort.parquet")  
  write_parquet(input_data, baseline_file)
  yearly_data[[1]] <- baseline_file
  
# Run simulation for given number of years -------------
  for (year in seq_len(sim_years)) {
                          # cat("updating year ", year, "\n")  # DEBUG
    # use previous year data as new input
    prev_year_file <- if (year == 1) yearly_data[[1]] else yearly_data[[2]]
    pd <-  as.data.table(read_parquet(prev_year_file))  # only holds one year data in memory at a time

    # Get cohort-level incidence rates ------------
    current_iso3 <- pd[1, iso3]
    cohort_params <- all_cohort_param_rates[all_cohort_param_rates$iso3 == current_iso3, ]   # pass filtered parameters by iso3 (matched to age and sex in update_cohort())
 
  ###########'
  # call update_cohort() -------------
    updated_person_data <- update_cohort(pd,             # input data
                                       cohort_params,    # iso3 matched incidence rates
                                       conditions,
                                       year,             
                                       glp1_access      = glp1_coverage,
                                       glp1_drug_type   = cs$glp1$drug_type,
                                       glp1_impact_data = glp1_impact_data,
                                       doublecount_multiplier_table = dc_multiplier_ds, 
                                       open_cohort      = cs$open_pop_cohort,
                                       sim_temp_dir     = sim_temp_dir,
                                       cal              = acm_calibration , 
                                       dw_sensitivity   = cs$daly_sensitivity
                                       ) 
    # --- Write Updated Data to Disk ---
    year_file <- file.path(sim_temp_dir, paste0("simulated_year_", year, ".parquet"))
    write_parquet(updated_person_data, year_file)
    
    if (year > 1) {
      prev_file <- file.path(sim_temp_dir, paste0("simulated_year_", year - 1, ".parquet"))
      if (file.exists(prev_file)) unlink(prev_file)
    }
      yearly_data[[2]] <- year_file   # overwrite previous years, so yearly_data contains baseline_file data and final year 
    
    rm(year_file, pd, updated_person_data)
    gc(verbose = F)
    Sys.sleep(0.00001) # time to close file 
  }
  return(yearly_data)
}