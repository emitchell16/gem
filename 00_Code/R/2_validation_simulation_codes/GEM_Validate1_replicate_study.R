#' Define functions to replicate input studies with GEM model

# 1. Function to run simulation ---------
gem_simulate <- function(input_data = y0_person_data, sim_years, 
                         apply_treatment = FALSE, te_input = NULL,
                         inc_data, storage_dir) {
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
  
  cat("\n starting sim \n")
  
  # sex-region level correlation variables and assign to every relevant row
  if (!is.data.table(input_data)) setDT(input_data)
  input_data <- copy(input_data) 
  
  input_data[, obese_t2d_corr := get_corr(t2d, obesity), by = .(who_region, female)]
  input_data[, cvd_t2d_corr   := get_corr(t2d, cvd),     by = .(who_region, female)]
  input_data[, ckd_t2d_corr   := get_corr(t2d, ckd),     by = .(who_region, female)]
  
  # set placeholders for incident case indicators for non-reversible conditions
  input_data[, c("t2d_new", "cvd_new", "ckd_new", "eskd_new", 
                 "stroke_this_year", "death_year") := 0L]
  
  yearly_data <- vector("list", length = 1) 
  baseline_file <- file.path(storage_dir, "grp_baseline_cohort.parquet")  
  write_parquet(input_data, baseline_file)
  yearly_data[[1]] <- baseline_file
  
  for (year in seq_len(sim_years)) {
    # use previous year data as new input
    prev_year_file <- yearly_data[[1]] 
    pd <-  as.data.table(read_parquet(prev_year_file))
    
    ###########'
    # call update_cohort() to update person-level disease statuses and age -------------
    updated_person_data <- v_update_cohort(pd,           
                                           cohort_params = inc_data,
                                           year, conditions = conditions_list,
                                           apply_treatment,
                                           treatment_impact_data = te_input
    ) 
    year_file <- file.path(storage_dir, paste0("simulated_year_", year, ".parquet"))
    write_parquet(updated_person_data, year_file)
    if (year == 1 && file.exists(baseline_file)) unlink(baseline_file)
    if (year > 1) {
      prev_file <- file.path(storage_dir, paste0("simulated_year_", year - 1, ".parquet"))
      if (file.exists(prev_file)) unlink(prev_file)
    }
    yearly_data[[1]] <- year_file   # overwrite previous years, results should just be final year
    rm(pd, updated_person_data)
    gc(verbose = F)
  }
  return(yearly_data)  
}

# 2. Function to make baseline data, simulate, and store results ------
generate_vresults <- function(vres_outfile) {
  #---------------------- generate cohort person data, simulate, & summarize outcomes -----------------
  validation_temp_dir <- file.path(working_folder, "validation_temp") 
  dir.create(validation_temp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # delete previous debugging runs:
    unlink(file.path(validation_out_dir, vres_outfile), recursive = T)
    
  # treatment effect cols
  te_cols <- names(validation_data_input)[grepl("^(treatment|hr)", names(validation_data_input))]
  
  # outcome condition list (excludes obesity)
  inc_outcome_conds <- c("ac_death", "ckd_new","cvd_new","eskd_new","mace", "stroke","t2d_new")
  
  res <- pblapply(seq_len(nrow(validation_data_input)), function(i){ 
  # res <- pblapply(seq_len(4), function(i){                                   # testing/debugging only
    # generate baseline person-level data if not already made for this study/cohort ----
    v_input_row_i    <- validation_data_input[i, , drop = FALSE]
    study_id <- v_input_row_i$study_id[[1]]
    cohort_id <- v_input_row_i$cohort_id[[1]]
    unique_cohort_baseline <- isTRUE(v_input_row_i$unique_cohort_baseline[i])
    
    if (unique_cohort_baseline) {
      fname <- sprintf("y0_input_%s_c%s.parquet", study_id, cohort_id)
    } else {  # if cohorts not unique at baseline
      fname <- sprintf("y0_input_%s.parquet", study_id)
    }
    fpath <- file.path(validation_temp_dir, fname)
    if (!file.exists(fpath)) {
      y0 <- v_initialize_baseline(v_input_data = v_input_row_i, paths = paths)
      write_parquet(y0, fpath)
    }
    # set simulation inputs -----
    y0 <- read_parquet(fpath)
    # treatment effect data for this study/cohort 
    te_row <- v_input_row_i |> select(any_of(te_cols))
    has_te <- !all(is.na(unlist(te_row[ , te_cols, drop = FALSE])))
    if (has_te) {
      te_long <- te_row |> 
        rename_with(~ sub("^treatment_hr(\\d+)$", "hr\\1_effect", .x)) |> 
        rename_with(~ sub("^hr(\\d+)_outcome$", "hr\\1_outcome", .x)) |> 
        pivot_longer(
          everything(),
          names_to = c("n", ".value"),
          names_pattern = "^hr(\\d+)_(effect|outcome)$"
        ) |> 
        transmute(
          outcome = tolower(as.character(outcome)),
          effect  = suppressWarnings(as.numeric(effect)),
        ) |> 
        mutate(
          outcome = case_when(
            outcome == "all cause mortality" ~ "ac_death",
            # outcome == "mean change in bmi"  ~ "bmi",
            TRUE ~ outcome
          )    ) |> 
        filter(!is.na(effect), !is.na(outcome), outcome != "")
    }
    
    sim_years <- v_input_row_i$futime[[1]]
    # get region-sex-age level incidence and mortality rates:
    cohort_incidence_rate_data <-  cohort_input_data |> 
      left_join(region_income_country_groups |> select(iso3, who_region), by = "iso3") |> 
      filter(age_high >= v_input_row_i$min_age[[1]], age_low  <= v_input_row_i$max_age[[1]])  %>%
      { 
        iso3_val <- v_input_row_i$iso3[[1]]
        if (!is.na(iso3_val) && nchar(iso3_val) == 3) filter(., iso3 == iso3_val) else .
      } %>%
      group_by(who_region, female, age_low, age_high) |> 
      summarise(across(
        c(where(is.numeric), -pop_size),
        ~ weighted.mean(.x, w = pop_size, na.rm = TRUE)
      ),
      .groups = "drop"
      ) |> 
      select(-starts_with("prevalence_"))
    
    temp_dir_i <- file.path(validation_temp_dir, paste0("v_group=", v_input_row_i$study_id[1], "-", v_input_row_i$cohort_id[1]))
    if (!dir.exists(temp_dir_i)) dir.create(temp_dir_i, recursive = TRUE) 
    # *** run sim *** ----
    sim_results_end_yr <- gem_simulate(
      input_data = y0,
      sim_years = sim_years, 
      apply_treatment = has_te,
      te_input = if (has_te) te_long else NULL,
      inc_data = cohort_incidence_rate_data,
      storage_dir = temp_dir_i
    )
    
    # comparison outcomes for this study-cohort:
    comp_row <- comparison_outcome_data |> 
      filter(.data$study_id == !!study_id, .data$cohort_id == !!cohort_id) |> 
      slice(1)
    
    incidence_time_unit_i <- comp_row$incidence_time_unit[1]
    
    ## calculate outcomes ----
    sim_df <- read_parquet(sim_results_end_yr[[1]]) |> 
      mutate(mace = 
               case_when( cvd_new ==1 | cvd_death == 1L | stroke ==1 ~ 1, TRUE ~ 0L))
    
    # incidence outcomes that exist in comparison data
    inc_keys <- inc_outcome_conds[!vapply(comp_row[inc_outcome_conds], function(x) all(is.na(x)), logical(1))]
    n_person <- nrow(sim_df)
    
    inc_outcome <- bind_rows(
      lapply(c(inc_keys), function(k){
        events <- sum(sim_df[[k]] %in% 1)
        estimates <- if (!grepl("cumulative", incidence_time_unit_i, ignore.case = TRUE)) {
          events / (n_person * sim_years) 
        } else {
          events / n_person # cumulative incidence per person over entire time horizon 
        }
        
        tibble(
          study_id  = study_id,
          cohort_id = cohort_id,
          outcome = k,
          measure = if (!grepl("cumulative", incidence_time_unit_i, TRUE)) "incidence_rate_1py" else "incidence_rate_cumulative", 
          estimate = estimates,
          test_indicator = "simulated"
        )
      })) 
    if (is.null(inc_outcome)) inc_outcome <- tibble()
    
    # add from reference study
    inc_ref <- bind_rows(lapply(inc_keys, function(k){
      ref_val <- suppressWarnings(as.numeric(comp_row[[k]][1]))
      tibble(
        study_id  = study_id,
        cohort_id = cohort_id,
        outcome   = k,
        measure   = if (!grepl("cumulative", incidence_time_unit_i, TRUE))
          "incidence_rate_1py" else "incidence_rate_cumulative",
        estimate  = ref_val,
        test_indicator = "reference_study"
      )
    }))
    if (is.null(inc_ref)) inc_ref <- tibble()
    
    ##combine all outcome tibbles for this study-cohort: ----
    sim_results_ir_ds <- bind_rows(inc_outcome, inc_ref) |> 
      mutate(estimate = round(estimate, 4))
    
    # save output as csv
    
    out_fp <- file.path(validation_out_dir, vres_outfile)
    if (file.exists(out_fp)) {
      old <- read_csv(out_fp, show_col_types = FALSE) 
      sim_results_ir_ds <- bind_rows(old, sim_results_ir_ds) |>  arrange(outcome)
    }
    write_csv(sim_results_ir_ds, out_fp)
    
    sim_results_ir_ds
    
    unlink(temp_dir_i, recursive = T)
    gc(verbose = F)
  })

}