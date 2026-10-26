#' Define functions to replicate input studies with GEM model
#' Simulation uses the main model state updates (v_update_cohort -> update_cohort)

# 1. Function to run simulation (in memory; validation cohorts are small) ---------
gem_simulate <- function(input_data, sim_years, inc_data, ext_te = NULL) {
  pd <- input_data
  for (year in seq_len(sim_years)) pd <- v_update_cohort(pd, inc_data, year, ext_te)
  pd
}

# 2. Function to make baseline data, simulate, and store results ------
generate_vresults <- function(vres_outfile) {
  #---------------------- generate cohort person data, simulate, & summarize outcomes -----------------
  validation_temp_dir <- file.path(paths$working, "validation_temp") 
  dir.create(validation_temp_dir, recursive = TRUE, showWarnings = FALSE)
  consts <- v_reference_consts(validation_temp_dir)     # population rr scaling constants (once)
  
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
    unique_cohort_baseline <- isTRUE(v_input_row_i$unique_cohort_baseline[[1]])
    
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
    iso3_val  <- v_input_row_i$iso3[[1]]
    iso3_study <- if (!is.na(iso3_val) && nchar(iso3_val) == 3) iso3_val else NA_character_
    # region-sex-age rates and prevalences (prevalences: excess mortality and BMI distribution)
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
      as.data.table() |>
      add_bmi_dist_params()
    setkey(cohort_incidence_rate_data, who_region, female, age_low)
    
    # run sim ----
    sim_df <- gem_simulate(
      input_data = v_prepare_baseline(y0, cohort_incidence_rate_data, consts, iso3_study, crn_base = i * 1000L),
      sim_years  = sim_years,
      inc_data   = cohort_incidence_rate_data,
      ext_te     = if (has_te) setNames(te_long$effect, te_long$outcome) else NULL
    )
    
    # comparison outcomes for this study-cohort:
    comp_row <- comparison_outcome_data |> 
      filter(.data$study_id == !!study_id, .data$cohort_id == !!cohort_id) |> 
      slice(1)
    
    incidence_time_unit_i <- comp_row$incidence_time_unit[1]
    
    ## calculate outcomes ----
    #  stroke: any stroke during follow-up; MACE: new CVD (incl. stroke) or CVD death
    sim_df[, `:=`(stroke = stroke_ever, mace = as.integer(cvd_new == 1 | cvd_death == 1L | stroke_ever == 1))]
    
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
    
    gc(verbose = F)
  })

}