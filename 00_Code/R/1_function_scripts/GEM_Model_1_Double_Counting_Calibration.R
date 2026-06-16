# -------------------'
# 00_code/R/1_function_scripts/GEM_Model_1_Double_Counting_Calibration.R
# Run versions of model to estimate indirect TE
# Generate treatment effect (TE) double counting adjustment calibrators for modeling direct and indirect TEs
# -------------------'

## Set up & run the two sims ---------------------
dc_calibration_fn <- function() {

  #  2-year simulation-------
  dct_main_sim_wrapper(files_yr0 = sample_baseline_files,
                       glp1_access_percentage = 0,
                       out_folder_path = dc_test_fp)
  cat("comparison case: ")
  dct_main_sim_wrapper(files_yr0 = sample_baseline_files,
                       glp1_access_percentage = 100,
                       out_folder_path = dc_test_fp)

  # Define outcomes & subgroups -------
  double_count_yvars <- c("ac_death", "stroke", "cvd", "ckd", "eskd", "t2d")
  
    #ac_death and stroke cohorts to match observed HRs:  obese without t2d; t2d with ckd; t2d
    #ckd and cvd cohorts to match observed HRs: t2d; obese without t2d
    #eskd cohorts to match: obese and ckd without t2d; t2d and ckd
    # t2d cohorts to match: obese without t2d
  subgroup_defs <- list(
    ac_death = exprs(
      `obese without t2d`    = (obesity==1 & t2d==0),
      `t2d and ckd`          = (t2d==1 & ckd==1),
      t2d                   = (t2d==1 & ckd==0) # does this need obesity == 0 or is it implied by order?
    ),
    cvd    = exprs(
      t2d                   = (t2d==1),
      `obese without t2d`   = (obesity==1 & t2d==0)
    ),
    eskd   = exprs(
      `obese without t2d` = (obesity==1 & t2d==0),
      `t2d and ckd`       = (t2d==1 & ckd==1)
    ),
    t2d   = exprs(
      `obese without t2d` = (bmi>=25 & t2d==0)
    )
  )
  subgroup_defs$stroke <- subgroup_defs$ac_death
  subgroup_defs$ckd    <- subgroup_defs$cvd
  
  # Read in sims ---------
  ds_treatment    <- open_dataset(filePath(dc_test_fp, "sim=glp1_access"), format="parquet")
  ds_base         <- open_dataset(filePath(dc_test_fp,"sim=basecase"), format="parquet")
  calc_dc_results <- function(ds_base, ds_treatment, double_count_yvars, subgroup_defs, glp1_impact_data) {
    
    map_dfr(double_count_yvars, function(y) {
      cat("processing", y, "\n")
      
      defs <- subgroup_defs[[y]]          # keep as expr list (from exprs())
      y_sym <- sym(y)
      is_death <- (y == "ac_death")
      
      ref_yr  <- if (is_death) 0L else 1L
      eval_yr <- if (is_death) 1L else 2L
      
      # baseline /reference year tables
      base_at_risk <- ds_base %>%
        filter(year == ref_yr, ac_death == 0) %>%
        transmute(
          y0      = !!y_sym, pid,
          obesity, t2d, ckd, bmi)

      trt_at_risk <- ds_treatment %>%
        filter(year == ref_yr, ac_death == 0) %>%
        transmute(y0      = !!y_sym, pid,
          obesity, t2d, ckd, bmi)
      
      # evaluation-year tables 
      base_eval <- ds_base %>%
        filter(year == eval_yr) %>%
        transmute( y1      = !!y_sym, pid)
      
      trt_eval <- ds_treatment %>%
        filter(year == eval_yr) %>%
        transmute(y1      = !!y_sym, pid)
      
      map_dfr(names(defs), function(subname) {
        cat("  subgroup", subname, "\n")
        colexpr <- defs[[subname]]
        
        # from trial data
        hr_obs <- glp1_impact_data  |> 
          filter(outcome == y, treated_pop == subname) %>%
          pull(outcome_val)
        
        # pid sets at risk 
        base_pidset <- base_at_risk %>%
          filter(!!colexpr, y0 == 0) %>%
          select(pid)
        
        trt_pidset <- trt_at_risk %>%
          filter(!!colexpr, y0 == 0) %>%
          select(pid)
      
        
        hazard_base <- base_eval %>%
          semi_join(base_pidset, by = "pid") %>%
          summarise(h = mean(y1, na.rm = TRUE)) %>%
          collect() %>%
          dplyr::pull(h)
        
        hazard_treat <- trt_eval %>%
          semi_join(trt_pidset, by = "pid") %>%
          summarise(h = mean(y1, na.rm = TRUE)) %>%
          collect() %>%
          dplyr::pull(h)
        
        hr_impl      <- hazard_treat / hazard_base
        
        tibble(
          outcome      = y,
          subgroup     = subname,
          year         = eval_yr,
          hazard_base  = hazard_base,
          hazard_treat = hazard_treat,
          HR_implied   = hr_impl,
          HR_observed  = hr_obs,
          multiplier   = HR_observed / HR_implied
        )
      })
    })
  }
  results <- calc_dc_results(ds_base, ds_treatment, double_count_yvars, subgroup_defs, glp1_impact_data)
  
  # write out the double‐count multipliers:
  write.csv(results,
            file.path(paths$working, paste0("double_counting_calibration_results_",te_source , ".csv")),
            row.names = FALSE)
  
  print(results)
}

# scratch process notes ---------
# for each endpoint c("ac_death", "stroke", "cvd", "ckd", "t2d"): 
# ac_death HR_implied evaluated in year 1 compared to year 0 ; stroke, cvd,ckd, eskd, t2d in year 2 compared to year 1
#' run treatment with just lower obesity prevalence 
#' calculate HR_implied = hazard_treated/hazard_control
#' compare to clinical trial HR (HRobserved)
#' compute multiplier = HRobs/HRimplied 
#' apply to yearly HR risk of event in full simulation for every treated year; starting in year 1 for death; and year 2 for other non-fatal conditions 

# calculate multipliers over test sim, then feed into 3.2_ codes and multiply treatment effect for applicable outcomes
# note double-counting adjustment not needed for obesity bc there's not a risk equation pathway 
# external validation and calibration applied only after double counting calibration applied. 
##############################'