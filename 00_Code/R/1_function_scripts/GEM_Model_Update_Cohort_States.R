#' -------------------
# GEM Model - Simulate disease status, treatment, and age updates
# Author: Liz Mitchell
# Date started:  Feb. 2025
# wip updating march 2026-- daly etc
#' -------------------


# helper function ----- 
# get risk and treatment adjusted incidence rates, apply
modify_risk <- function(cond, inc_probs, matched_params, pd, living_idx, inc_val, new_state_vals,
                        glp1_access, doublecount_multiplier_table) {
  ####
  # A. Modify incidence rates with condition-specific risk factors: --------
  ####
  inc_probs <- calculate_adjusted_risk(cond, inc_probs, matched_params,
                                       pd, living_idx, inc_val, new_state_vals)
  ####
  # B. Modify risks based on GLP-1RA treatment ------------
  ####
    if (glp1_access > 0) {
      treated_idx  <- which(pd[["treated"]][living_idx] == 1)
      if (length(living_idx[treated_idx]) > 0){
         if (cond %in% c("stroke", "ac_death")) {
          grp_labels <- case_when(
            pd[["t2d"]][living_idx[treated_idx]] == 1 & pd[["ckd"]][living_idx[treated_idx]] == 1 ~ "t2d and ckd",
            pd[["t2d"]][living_idx[treated_idx]] == 1 ~ "t2d",
            TRUE ~ "obese without t2d"  # non-t2d glp1 treated persons currently or previously obese
          )} else if(cond %in% c("ckd", "cvd")) {
          grp_labels <- case_when(
            pd[["t2d"]][living_idx[treated_idx]] == 1 ~ "t2d",
            pd[["obesity"]][living_idx[treated_idx]] == 1 ~ "obese without t2d",
            TRUE ~ "obese without t2d"   # non-t2d glp1 treated persons currently or previously obese
          )} else if(cond %in% c("t2d")) {
          grp_labels <- case_when(
            pd[["t2d"]][living_idx[treated_idx]] == 0 ~ "obese without t2d",  # include all non-t2d treated persons currently or previously obese or overweight
            TRUE ~ NA_character_
          )}
      for (grp in unique(grp_labels[!is.na(grp_labels)])) {
        # Get indices for treated individuals in this subgroup:
        pos <- which(grp_labels == grp)
        grp_idx_abs <- living_idx[treated_idx][pos]
        grp_idx_rel <- treated_idx[pos]
        ### i. Observed TE multiplier ----
        outcome_data <- glp1_impact_data[glp1_impact_data$outcome == cond &
                                           glp1_impact_data$treated_pop == grp, ]
        if (nrow(outcome_data) == 0) {
          multiplier <- rep(1, length(grp_idx_rel))
          if (grp != "other") {
            warning("outcome data not matched for group", grp)}
        } else {
          multiplier <- rep(outcome_data$outcome_val[1], length(grp_idx_rel))
         }
        ### ii. Calibrated double-counting adjustment multiplier ----
        dc_info <- doublecount_multiplier_table |>
          filter(outcome == cond, subgroup == grp)
        if (nrow(dc_info) == 0) {
          warning("No double‐count entry for ", cond, "/", grp, "; using dc_mult=1")
          dc_mult   <- 1; threshold <- if (cond == "ac_death") 1 else 2
        } else {
          dc_mult   <- dc_info$multiplier[1]; threshold <- dc_info$year[1]
        }
        years_treated_vec <- new_state_vals$years_treated[grp_idx_abs]
        adj <- ifelse(years_treated_vec >= threshold,          # 1 for ac_death, 2 for non-fatal outcomes
                      1 - (1 - multiplier) * dc_mult, multiplier )
        inc_probs[grp_idx_rel] <- inc_probs[grp_idx_rel] * adj
      }
    }}
              ## debugging code (active)----
              if (any(is.na(inc_probs))) {
                cat("NA values found in inc_probs: ", paste(cond, head(inc_probs, 10), collapse = ", "))
                stop("break due to NAs")
              }

              if (any(inc_probs < 0)) {
                neg_count <- sum(inc_probs < 0)
                total_count <- length(inc_probs)
                # Print negative values (up to the first 10)
                neg_values <- head(inc_probs[inc_probs < 0], 10)
                cat("Negative values found in inc_probs for", cond, ":", paste(neg_values, collapse = ", "), "\n")
                cat("Count of negative values:", neg_count, "out of", total_count, "\n")
                stop("break due to negative incidences")
              }
              ###----
  return(inc_probs)
}

# MAIN FUNCTION ------
update_cohort <- function(pd, cohort_params, conditions, year,
                          glp1_access, glp1_drug_type, glp1_impact_data, doublecount_multiplier_table, open_cohort = F,
                          sim_temp_dir, cal = acm_calibration , dw_sensitivity = F) {
  ####
  # 1. Set up --------------
  ####
  new_state_vals <- copy(pd)
  living_idx <- which(new_state_vals$ac_death == 0)
  
  if (glp1_access > 0) {
    treated_idx  <- which(pd[["treated"]][living_idx] == 1)
    new_state_vals[living_idx[treated_idx], years_treated := years_treated + 1]    # years_treated starts at 0, updated at the start of each sim year
    } 

  # match individual age & sex in cohort_params (age increases over modeled years so starting cohort can split)
  setDT(cohort_params); setkey(cohort_params, female, age_low)
  suppressWarnings({  
    matched_params <- cohort_params[pd, 
                                    on        = .(female, age_low = age),
                                    roll      = TRUE,
                                    rollends  = c(TRUE, TRUE)  ]
  })
  
  ####
  # 2. Update condition incidence rates -------
  ####
    ## 2.1. Update ckd, t2d, stroke -----
    for (cond in c("ckd", "t2d", "stroke")) { 
      if (length(living_idx) == 0) break
      vcol <- paste0("incidence_", cond, "_val")
      inc_val   <- matched_params[[vcol]]
      inc_probs <- inc_val[living_idx]
      inc_probs <- pmin(1, pmax(0, inc_probs))

                #'### debug: ----
                if (any(is.na(inc_val))) {
                  stop(sprintf(
                    "cohort data returning NA in inc_val for cond '%s' (iso3=%s, female=%s)",
                    cond, pd$iso3[1], pd$female[1]))}
                #####'

      ### function call ----
      inc_probs <- modify_risk(cond,inc_probs, matched_params,
                               pd, living_idx, inc_val, new_state_vals, glp1_access, doublecount_multiplier_table) 
      
      ### draw events -------
      pos_no_cond <- which(new_state_vals[[cond]][living_idx] == 0)
      idx         <- living_idx[pos_no_cond]
      if (length(idx) > 0) {
        p  <- inc_probs[pos_no_cond]
        new_state_vals[idx, (cond) := rbinom(.N, 1, pmin(1, pmax(0, p)))]
      }
    }
  
  ## 2.2. update CVD after stroke -----
      # note: stroke modeled as subtype of cvd. draw non-stroke cvd among those without prior cvd
   if (("cvd" %in% conditions) && length(living_idx) > 0) {
    inc_prob_cvd_total_base   <- matched_params[["incidence_cvd_val"]][living_idx]
   
    inc_probs <- modify_risk(
      cond = "cvd",
      inc_probs = inc_prob_cvd_total_base,
      matched_params, pd, living_idx,
      inc_val = matched_params[["incidence_cvd_val"]],
      new_state_vals = new_state_vals,
      glp1_access = glp1_access,
      doublecount_multiplier_table = doublecount_multiplier_table
    )
    ### draw events -------
      # cvd incidence among ppl without past cvd or stroke this year
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
      new_state_vals[idx,  cvd := 1]
    }
   }
  
    ######'  
    # 3. Update obesity -----------
    #' assumes adjacent only transitions
    ######'  
      # annual transition probabilities for this cohort
      #' class incidence rates can have 3 competing risk steps (down, stay, or up)
      inc_obesity <- pmin(1, pmax(-1, matched_params[living_idx, incidence_obesity_val]))
      inc_c1 <-pmin(1, pmax(-1, matched_params[living_idx, incidence_class1_obesity_val]))
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
        }
        # ##### debugging --------
        # old_obesity <- pd[living_idx, obesity]
        # new_obesity <- new_state_vals[living_idx, obesity]
        # diff_count <- sum(new_obesity != old_obesity, na.rm = TRUE)
        # if (diff_count > 0) {
        #   cat("DEBUG: Obesity indicator changed for", diff_count, "individuals.\n")
        # }
        # #####'
      }
    
    ######'
      ## 3.1 Update obesity status with GLP-1 treatment ---------
    ######'
      #' for treated, maintain past year weight vars before applying glp-1 update logic
      #' applies when bmi >= 27
    
      if (glp1_access > 0 && length(treated_idx) > 0) { 
        treated_abs <- living_idx[treated_idx]
        years_treated_vec <- new_state_vals$years_treated[treated_abs]
        
        # For individuals with years_treated > 1, restore to past year's values
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
  # 4. Update ESKD among people with CKD --------- 
  ######'
  eskd_inc_prob_total_pop        = input_params$value[input_params$parameter == "eskd_inc_prob_total_pop"]
  ckd_prevalence    <- mean(new_state_vals[["ckd"]][living_idx], na.rm = TRUE)
  eskd_inc_prob_ckd <- eskd_inc_prob_total_pop/ ckd_prevalence
  eskd_inc_prob_ckd <- pmin(1, pmax(0, eskd_inc_prob_ckd))
  eskd_update_idx   <- which(new_state_vals[["ckd"]][living_idx] == 1 & new_state_vals[["eskd"]][living_idx] == 0)
  # base case updating:
  if (glp1_access == 0 && length(eskd_update_idx) > 0) {   
    global_idx <- living_idx[eskd_update_idx]
    new_state_vals[global_idx, eskd := rbinom(.N, 1, eskd_inc_prob_ckd)]
  }
  # intervention simulation updating:
  else if (glp1_access > 0 && length(eskd_update_idx) > 0) {
    treated_in_update <- intersect(eskd_update_idx, treated_idx)
    if (length(treated_in_update) > 0) {
      treated_ckd <- treated_idx[
        pd[["ckd"]][living_idx][treated_idx] == 1 & 
        pd[["eskd"]][living_idx][treated_idx] == 0]
      # Divide treated_ckd into cohorts
          # (1) obese without t2d & (2) t2d
        obese_no_t2d_idx <- treated_ckd[pd[["obesity"]][living_idx][treated_ckd] == 1 &
                                           pd[["t2d"]][living_idx][treated_ckd] == 0]
        t2d_idx <- treated_ckd[pd[["t2d"]][living_idx][treated_ckd] == 1]
        treated_eskd_idx <- union(obese_no_t2d_idx, t2d_idx)
        eskd_update_no_glp1_idx <- setdiff(eskd_update_idx, treated_eskd_idx)
        
        # Update outcomes for each subgroup using apply_eskd function from code file *_2_Risk_Equation_Data_Processing.R
        if (length(obese_no_t2d_idx) > 0) {
          global_idx <- living_idx[obese_no_t2d_idx]
          years_treated_vec <- new_state_vals$years_treated[global_idx]
          
          dc_info <- doublecount_multiplier_table |> 
            filter(outcome == "eskd", subgroup == "obese without t2d") |>   
            slice(1)
        
          adj_events <- apply_eskd(obese_no_t2d_idx, 
                                     glp1_impact_data[glp1_impact_data$outcome == "eskd" & 
                                                      glp1_impact_data$treated_pop == "obese without t2d", ], 
                                     base_prob = rep(eskd_inc_prob_ckd, length(global_idx)),
                                     dc_factor = dc_info$multiplier,   
                                     years_treated_vec,
                                     dc_threshold = dc_info$year
                                   ) 
          new_state_vals[global_idx, eskd := adj_events]
        }
        
        if (length(t2d_idx) > 0) {
          global_idx <- living_idx[t2d_idx]
          years_treated_vec <- new_state_vals$years_treated[global_idx]
          dc_info <- doublecount_multiplier_table |> 
            filter(outcome == "eskd", subgroup == "t2d and ckd") |>   
            slice(1)
          
          adj_events <-apply_eskd(t2d_idx, 
                                 glp1_impact_data[glp1_impact_data$outcome == "eskd" & 
                                                    glp1_impact_data$treated_pop == "t2d and ckd", ], 
                                 base_prob = rep(eskd_inc_prob_ckd, length(global_idx)),
                                 dc_factor = dc_info$multiplier,   
                                 years_treated_vec,
                                 dc_threshold = dc_info$year)
          
          new_state_vals[global_idx, eskd := adj_events]
        }
        if (length(eskd_update_no_glp1_idx) > 0) {  # untreated within update set
          global_idx <- living_idx[eskd_update_no_glp1_idx]
          new_state_vals[global_idx, eskd := rbinom(.N, 1, eskd_inc_prob_ckd)]
        }
      } else if ((length(treated_in_update) == 0)) { # fall back if no one in update is treated
        global_idx <- living_idx[eskd_update_idx]
        new_state_vals[global_idx, eskd := rbinom(.N, 1, eskd_inc_prob_ckd)]
      }
    }
  ####
  # 5. Incident case indicators updated for non-reversible conditions -------
  ####
  nr_conds <- c("t2d", "cvd", "ckd", "eskd") 
  for (cond in nr_conds) {
    new_col <- paste0(cond, "_new")
    # if prior value is 0 and new value is 1, flag as 1
    new_state_vals[, (new_col) := ifelse(pd[[cond]] == 0 & get(cond) == 1, 1L, get(new_col))]
  }
  new_state_vals[, stroke_this_year   := fifelse(pd[["stroke"]] == 0 & stroke == 1, 1L, 0L)]
  
  ####
  # 6. Update death -------
  ####
  vcol <- paste0("incidence_ac_death", "_val")

  inc_val_death   <- matched_params[[vcol]]
  inc_probs_death <- pmin(1, pmax(0, inc_val_death[living_idx]))

  ## mortality calibration from validation step if available: -----
  if (!is.null(cal)) {
    if ("alpha" %in% names(cal) && "beta" %in% names(cal)) {
      idx <- which(tolower(trimws(cal$outcome)) == "ac_death")[1]
      a <- as.numeric(cal$alpha[idx]); b <- as.numeric(cal$beta[idx])
      p <- pmin(pmax(inc_probs_death, 1e-12), 1 - 1e-12)
      h <- -log1p(-p)
      hc <- exp(a) * (h ^ b)
      inc_probs_death <- pmin(pmax(1 - exp(-hc), 1e-12), 1 - 1e-12)
    } 
  } 
  ## 6.1. Apply risk adjustment for comorbidities, accounting for those accrued in past year -----
  inc_probs_death <- modify_risk(cond="ac_death", inc_probs_death, matched_params,
                                 pd=new_state_vals, living_idx, inc_val_death, new_state_vals, glp1_access, doublecount_multiplier_table) 
 
  ## 6.2. Draw deaths ----
  new_state_vals$ac_death[living_idx] <- rbinom(length(living_idx), 1, pmin(1, pmax(0, inc_probs_death)))

  ##########' 
  # 7. Update DALY inputs (YLL + YLD) ---------
  ##########'
  #' METHODS NOTE:
  #'   - multiplicative complement method for combining, aligned with GBD for YLD
  #'          -- purpose: avoid impossible additive disability burden across comorbidities
  #'  - compute YLD at end of update_cohort() based on current-year updated states. this is internally consistent with death based on current-year states
  #'  - for those who die during year, use half cycle factor: survived = 1.0 // died = 0.5  
  
  # conventional 3% discount factor, first year undiscounted
  disc <- 1 / ((1 + disc_rate)^(year-1)) 
  
  ## 7.1 YLL -----
  newly_dead <- which(pd$ac_death == 0 & new_state_vals$ac_death == 1)
  if (length(newly_dead) > 0) {
    yll_vals <- matched_params[["le_yrs_val"]][newly_dead]
    if (any(is.na(yll_vals))) {  # active debug
      stop("Missing remaining life expectancy values for newly deceased")
    }
    new_state_vals[newly_dead, `:=`(
                                    death_year = year,
                                    yll = yll_vals,
                                    yll_adj = yll_vals * disc
    )]
  }
  ## 7.2 YLD -----
  
  # Time alive in this cycle (half-cycle correction for deaths)
    # NOTE: living_idx set before drawing this cycle's deaths, allowing YLD for newly dead with time_alive 0.5
  time_alive <- numeric(nrow(new_state_vals))
  time_alive[living_idx] <- ifelse(new_state_vals$ac_death[living_idx] == 1L, 0.5, 1.0)
  
  dw_t2d_vec  <- assign_dw_t2d(new_state_vals, living_idx, dw_sensitivity)
  dw_cvd_vec  <- assign_dw_cvd(new_state_vals, living_idx, dw_sensitivity)
  dw_ckd_vec  <- assign_dw_ckd(new_state_vals, living_idx, dw_sensitivity)
  dw_stroke_vec <- assign_dw_stroke(new_state_vals, living_idx, dw_sensitivity)
  
  ### Combine & store disability weights -----
    # multiplicative complement of (1 - dw_k)
  dw_combined_vec <- 1 - (
    (1 - dw_t2d_vec) *
      (1 - dw_cvd_vec) *
      (1 - dw_ckd_vec) *
      (1 - dw_stroke_vec)
  )
  dw_combined_vec <- pmin(1, pmax(0, dw_combined_vec))
  
  # current year yld
  yld_vec     <- dw_combined_vec * time_alive[living_idx]
  yld_t2d_vec <- dw_t2d_vec * time_alive[living_idx]
  yld_cvd_vec <- dw_cvd_vec * time_alive[living_idx]
  yld_ckd_vec <- dw_ckd_vec * time_alive[living_idx]
  
  
  # cumulative totals
  new_state_vals[living_idx, yld_total := yld_total + yld_vec]  # unadjusted comparison
  new_state_vals[living_idx, yld_adj_total := yld_adj_total + (yld_vec * disc)]
  
  new_state_vals[living_idx, yld_t2d_adj_total := yld_t2d_adj_total + (yld_t2d_vec * disc)]
  new_state_vals[living_idx, yld_cvd_adj_total := yld_cvd_adj_total + (yld_cvd_vec * disc)]
  new_state_vals[living_idx, yld_ckd_adj_total := yld_ckd_adj_total + (yld_ckd_vec * disc)]
  #-----------------------------------
  ##########'
  #8. Costing ---------
  ##########'
  # - Costs discounted at 3% annually
  # - if deceased in current year, prorate chronic maintenance costs by time_alive (treatment, t2d, cvd, chronic stroke, ckd, eskd)
  
  ## 8.1 Treatment costs ------
  #' annual updates to the following vars:
  if (config_sims$costing) {    
    curr_iso3 <- unique(new_state_vals$iso3)
    if (glp1_access > 0 && length(treated_idx) > 0) {
      treated_abs <- living_idx[treated_idx]
      glp1_row <- glp1_cost_inputs |>
        filter(iso3 == curr_iso3)
      
      for (sfx in price_suffixes) {
        price_col <- paste0("price_", sfx)
        treat_col <- paste0("cost_treat_total_", sfx)
        total_col <- paste0("cost_total_", sfx)
        
        glp1_cost_val <- glp1_row |> pull(.data[[price_col]])
        
        glp1_cost_val <- glp1_cost_val * 12  # annualize monthly cost
        
        c_treated_disc <- (glp1_cost_val * time_alive[treated_abs]) * disc
        
        new_state_vals[[treat_col]][treated_abs] <- new_state_vals[[treat_col]][treated_abs] + c_treated_disc
        
        new_state_vals[[total_col]][treated_abs] <-
          new_state_vals[[total_col]][treated_abs] + c_treated_disc
      }
    }
    
    # helper: add health-state cost to every treatment cost scenario total
    add_to_all_cost_totals <- function(idx, add_cost, state_cost_col = NULL) {
      if (length(idx) == 0) return(NULL)
      if (!is.null(state_cost_col)) {
        new_state_vals[[state_cost_col]][idx] <<-
          new_state_vals[[state_cost_col]][idx] + add_cost
      }
      new_state_vals$cost_conditions_total[idx] <<-
        new_state_vals$cost_conditions_total[idx] + add_cost
      for (tc in cost_total_cols) {
        new_state_vals[[tc]][idx] <<-
          new_state_vals[[tc]][idx] + add_cost
      }
    }
    
    ## 8.2 Health state costs ------
    if (length(living_idx) > 0) { 
      health_costs <- health_event_cost_inputs  |>  filter(iso3 == curr_iso3)
      
      ### t2d ----
      c_t2d_lt10  <- health_costs |>  pull(cost_id_t2d_lt10)
      c_t2d_geq10 <- health_costs |>  pull(cost_id_t2d_geq10)
      
      idx_new  <- living_idx[new_state_vals$t2d_new[living_idx] == 1]
      idx_prev <- living_idx[new_state_vals$t2d_new[living_idx] == 0 & new_state_vals$t2d[living_idx] == 1]
      
      if (length(idx_new) > 0) {
        add_cost <- c_t2d_lt10 * time_alive[idx_new] * disc 
        add_to_all_cost_totals(idx_new, add_cost, "cost_t2d_tot")
      }
      if (length(idx_prev) > 0) {
        add_cost <- ((c_t2d_lt10 + c_t2d_geq10) / 2) * time_alive[idx_prev] * disc
        add_to_all_cost_totals(idx_prev, add_cost, "cost_t2d_tot")
      }
         
      ### cvd, stroke  ----
      # acute / chronic stroke
      c_stroke_acute  <- health_costs |>  pull(cost_id_stroke_acute)
      c_stroke_post   <- health_costs |>  pull(cost_id_stroke_chronic)
     
      newly_stroke <- living_idx[new_state_vals$stroke_this_year[living_idx] == 1L]
      
      if (length(newly_stroke) > 0) { 
        add_cost <- c_stroke_acute * disc
        add_to_all_cost_totals(newly_stroke, add_cost, "cost_stroke_tot")
      }
      post_stroke <- which(pd$stroke == 1 & new_state_vals$stroke == 1)
      if (length(post_stroke) > 0) { 
        add_cost <- c_stroke_post *time_alive[post_stroke] * disc
        add_to_all_cost_totals(post_stroke, add_cost, "cost_stroke_tot")
      }
      # cvd when stroke == 0
      c_cvd  <- health_costs |>  pull(cost_id_cvd)
      idx_cvd_nostroke <- living_idx[new_state_vals$cvd[living_idx] == 1 & new_state_vals$stroke[living_idx] == 0]
      if (length(idx_cvd_nostroke) > 0) {
        add_cost <- c_cvd * time_alive[idx_cvd_nostroke] * disc
        add_to_all_cost_totals(idx_cvd_nostroke, add_cost, "cost_cvd_tot")
      }
      ### ckd, eskd (mutually exclusive costs)   ----
      c_ckd  <- health_costs |>  pull(cost_id_ckd)
      c_eskd <- health_costs |>  pull(cost_id_eskd)
      idx_eskd <- living_idx[new_state_vals$eskd[living_idx] == 1]
      idx_ckd  <- living_idx[new_state_vals$ckd[living_idx] == 1 & new_state_vals$eskd[living_idx] == 0]
      if (length(idx_ckd) > 0) {
        add_cost <- c_ckd * time_alive[idx_ckd] * disc
        add_to_all_cost_totals(idx_ckd, add_cost, "cost_ckd_tot")
      }
      if (length(idx_eskd) > 0) {
        add_cost <- c_eskd  * time_alive[idx_eskd] * disc
        add_to_all_cost_totals(idx_eskd, add_cost, "cost_eskd_tot")
      }
      ### ac_death ---- 
      c_death    <- health_costs |>  pull(cost_id_death)
      newly_dead <- which(pd$ac_death == 0 & new_state_vals$ac_death == 1)
      if (length(newly_dead) > 0) {
        add_cost <- c_death * disc  
        add_to_all_cost_totals(newly_dead, add_cost, "cost_death_tot")
      }
    }
  }  
  
  ##########'
  # 9. Age increased one year among living -----------
  new_state_vals[ac_death == 0, age := pmin(age + 1, 99)] # top coded at 99
  
  # if open cohort, introduce new youngest group
  if (open_cohort && min(matched_params$age_low) <= min_age_model) {
    new_entrants <- fread(file.path(sim_temp_dir, "youngest_y0_person_data.csv")) |>
      mutate(pid = paste0(pid, "_entry", year)) 
    #append to new_state_vals
    new_state_vals[, pid := as.character(pid)]
    new_state_vals <- bind_rows(new_state_vals, new_entrants)
    rm(new_entrants)
  }
    
  # 10. New GLP-1 treatment eligibility and treatment status ------
    if (glp1_access > 0) {
      new_state_vals <- assign_glp1_eligible(indication="all", input_data = new_state_vals)
      # assign treatment among eligible, not treated in past year, don't re-evaluate if they already didn't have access:
      new_state_vals[treatment_eligible == 1 & 
                       treated == 0 & 
                       treatment_access_checked == 0,
                     `:=`(
                       treated                  = as.integer(runif(.N) < glp1_access/100),
                       treatment_access_checked = 1)
                     ]
    }
  return(new_state_vals)
  }
