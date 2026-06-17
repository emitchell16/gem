#' -------------------------
#' Apply risk equations with calculate_adjusted_risk()
#' Date updated: March 2026
#' -------------------------

##*****************
# main risk equation function ------

# new March 2026: stroke drawn before cvd, then cvd equation sets cvd = 1, 
# non-stroke cvd updated based on cvd incidence + risk modifiers

calculate_adjusted_risk <- function(cond, inc_probs, matched_params, pd, living_idx,
                                    inc_val, new_state_vals) {
  # T2D ----
  if (cond == "t2d") {
    obesity_multiplier <- rep(1, length(inc_probs))
    rr_rows <- obesity_rr_by_age[obesity_rr_by_age$rr_condition_with_obesity == "t2d", ]
    for (i in seq_len(nrow(rr_rows))) {
      row <- rr_rows[i, ]
      idx <- which(pd[["obesity"]][living_idx] == 1 &
                     pd[["age"]][living_idx] >= row$age_low &
                     pd[["age"]][living_idx] <= row$age_high)
      if (length(idx) > 0) {obesity_multiplier[idx] <- row$value[1]}
    }
    inc_probs <- pmin(1, inc_probs * obesity_multiplier)
    
  # CVD -----
  #' risk factors: obesity, t2d
  } else if (cond == "cvd") {
    ## 1) obesity by degree -----
      # multiplier is RR of 5-unit increase in BMI>30
    obese_idx <- which(pd[["obesity"]][living_idx] == 1)
    obesity_multiplier <- rep(1, length(inc_probs))
    
    if (length(obese_idx) > 0) {
      bmi_increment = (pd[["bmi"]][living_idx][obese_idx]- 30)/5
      obesity_multiplier[obese_idx] <- obesity_rr_cvd$value[1]
     
      # Raise the risk multiplier to the power of the number of 5-unit increments above 30
      obesity_multiplier[obese_idx] = obesity_multiplier[obese_idx]^ bmi_increment 
  
      inc_probs <- pmin(1, inc_probs * obesity_multiplier)
    }
     
    ## 2) t2d  -------
    t2d_idx     <- which(pd[["t2d"]][living_idx] == 1)
    t2d_multiplier <- rep(1, length(inc_probs))
    if (length(t2d_idx) > 0) {
      t2d_rr_row <- t2d_rr[t2d_rr$condition == "cvd", ]
      t2d_multiplier[t2d_idx] <-t2d_rr_row$value[1]
      inc_probs <-  pmin(1, inc_probs * t2d_multiplier)
      }
    ### joint adjustment ------
    joint_idx <- which(pd[["obesity"]][living_idx] == 1 & pd[["t2d"]][living_idx] == 1)
    
    if (length(joint_idx) > 0 && all(pd$obese_t2d_corr > 0)) {
      obese_t2d_corr <- pd[["obese_t2d_corr"]][living_idx][joint_idx]
      # combine multipliers: weighted blend of product vs max
      joint_mult <- (1 - obese_t2d_corr) * (obesity_multiplier[joint_idx] * t2d_multiplier[joint_idx]) +
                    obese_t2d_corr * pmax(obesity_multiplier[joint_idx], t2d_multiplier[joint_idx])
      # difference joint_mult and current combination is the adjustment factor
      d = (obesity_multiplier[joint_idx] * t2d_multiplier[joint_idx])
      d[d<=  1e-9] <- 1e-9 
      adj_factor <- joint_mult / d
      
      inc_probs[joint_idx] <- pmin(1, inc_probs[joint_idx] * adj_factor)
      }
  
   # CKD -----
  #' risk factors: obesity, t2d
  } else if (cond == "ckd") {  
  
    ## 1) obesity ----
    obesity_multiplier <- rep(1, length(inc_probs))
    obese_idx <- which(pd[["obesity"]][living_idx] == 1)
    if (length(obese_idx) > 0) {
      obesity_multiplier[obese_idx] <- obesity_rr_ckd$value[1]
      inc_probs <-  pmin(1, inc_probs * obesity_multiplier)
      }
    ## 2) t2d ----
    t2d_multiplier <- rep(1, length(inc_probs))
    t2d_idx <- which(pd[["t2d"]][living_idx] == 1)
    
    if (length(t2d_idx) > 0) {
      t2d_rr_row <- t2d_rr[t2d_rr$condition == "ckd", ]
      t2d_multiplier[t2d_idx]  <- t2d_rr_row$value[1]
      inc_probs <-  pmin(1, inc_probs * t2d_multiplier)  
    }
    
    ### joint adjustment ------
    joint_idx <- which(pd[["obesity"]][living_idx] == 1 & pd[["t2d"]][living_idx] == 1)
    
    if (length(joint_idx)>0 && all(pd$obese_t2d_corr > 0)) {
      obese_t2d_corr <- pd[["obese_t2d_corr"]][living_idx][joint_idx]
      # combine multipliers: weighted blend of product vs max
      joint_mult <- (1 - obese_t2d_corr) * (obesity_multiplier[joint_idx] * t2d_multiplier[joint_idx]) +
        obese_t2d_corr * pmax(obesity_multiplier[joint_idx], t2d_multiplier[joint_idx])
      # difference joint_mult and current combination is the adjustment factor
      d = (obesity_multiplier[joint_idx] * t2d_multiplier[joint_idx])
      d[d<=  1e-9] <- 1e-9 
      adj_factor <- joint_mult / d
      inc_probs[joint_idx] <- pmin(1, inc_probs[joint_idx] * adj_factor)
    }
  
  # STROKE ----
  #' risk factors: cvd, obesity, t2d
  } else if (cond == "stroke") {
    ## 1) CVD ----
  
    cvd_multiplier <- rep(1, length(inc_probs))
    cvd_idx <- which(pd[["cvd"]][living_idx] == 1)
    
    if (length(cvd_idx) > 0) {
      cohort_female <- pd[["female"]][living_idx][1]
      rr_val <- if (cohort_female == 1) {
        cvd_rr_stroke_bysex$value[cvd_rr_stroke_bysex$female == 1][1]
      } else {
        cvd_rr_stroke_bysex$value[cvd_rr_stroke_bysex$female == 0][1]
      }
      cvd_multiplier[cvd_idx] <- rr_val
    }
    
    inc_probs <- pmin(1, inc_probs * cvd_multiplier) 
    
    ## 2) obesity ----
    obesity_multiplier <- rep(1, length(inc_probs))
    rr_rows <- obesity_rr_by_age[obesity_rr_by_age$rr_condition_with_obesity == "stroke", ]  
    
    for(i in seq_len(nrow(rr_rows))) {   # for each age group, match risk adjustment
      row <- rr_rows[i, ]
      group_idx <- which(
        pd[["obesity"]][living_idx] == 1 & 
          pd[["age"]][living_idx] >= row$age_low & 
          pd[["age"]][living_idx] <= row$age_high
      )
      if (length(group_idx) > 0) {
        obesity_multiplier[group_idx] <- row$value[1]}
    }
    inc_probs <-  pmin(1, inc_probs * obesity_multiplier)  
  
    ## 3) t2d ----
    t2d_multiplier <- rep(1, length(inc_probs))
    t2d_idx  <- which(pd[["t2d"]][living_idx] == 1)
    
    if (length(t2d_idx) > 0) {
      t2d_rr_row <- t2d_rr[t2d_rr$condition == "stroke", ]
      t2d_multiplier[t2d_idx]  <- t2d_rr_row$value[1]
      inc_probs <-  pmin(1, inc_probs * t2d_multiplier)
    }
    ### joint adjustment ------ 
      ### 1. obesity-t2d 
    joint_idx <- which(pd[["obesity"]][living_idx] == 1 & pd[["t2d"]][living_idx] == 1)
    
    if (length(joint_idx) > 0 && all(pd$obese_t2d_corr > 0)) {
      obese_t2d_corr <- pd[["obese_t2d_corr"]][living_idx][joint_idx]
      # combine multipliers: weighted blend of product vs max
      joint_mult <- (1 - obese_t2d_corr) * (obesity_multiplier[joint_idx] * t2d_multiplier[joint_idx]) +
        obese_t2d_corr * pmax(obesity_multiplier[joint_idx], t2d_multiplier[joint_idx])
      d = (obesity_multiplier[joint_idx] * t2d_multiplier[joint_idx])
      d[d<=  1e-9] <- 1e-9 
      adj_factor <- joint_mult / d
      
      inc_probs[joint_idx] <- pmin(1, inc_probs[joint_idx] * adj_factor)
    }
      ### 2. t2d-cvd
    joint_idx <- which(pd[["cvd"]][living_idx] == 1 & pd[["t2d"]][living_idx] == 1)
    
    if (length(joint_idx) > 0 && all(pd$cvd_t2d_corr > 0)) {
      cvd_t2d_corr <- pd[["cvd_t2d_corr"]][living_idx][joint_idx]
      joint_mult <- (1 - cvd_t2d_corr) * (cvd_multiplier[joint_idx] * t2d_multiplier[joint_idx]) +
        cvd_t2d_corr * pmax(cvd_multiplier[joint_idx], t2d_multiplier[joint_idx])
      d = (cvd_multiplier[joint_idx] * t2d_multiplier[joint_idx])
      d[d<=  1e-9] <- 1e-9 
      adj_factor <- joint_mult / d
      
      inc_probs[joint_idx] <- pmin(1, inc_probs[joint_idx] * adj_factor)
    }
    
  # AC DEATH ---- 
  #' risk factors: obesity, t2d, ckd, cvd 
  } else if (cond == "ac_death") {
    
    ## 1) obesity by age, sex, and class  ----
      # set helper variables and parameters
      obesity_vals <- pd[["obesity"]][living_idx]
      class1_vals <- pd[["class_1_obesity"]][living_idx]
      class2_vals <- pd[["class_2_obesity"]][living_idx]
      class3_vals <- pd[["class_3_obesity"]][living_idx]
      female_vals <- pd[["female"]][living_idx]
      age_vals    <- pd[["age"]][living_idx]
      
      class1_params <- obesity_rr_acdeath_byclass[obesity_rr_acdeath_byclass$condition == "class1_obesity", ]
      class2_params <- obesity_rr_acdeath_byclass[obesity_rr_acdeath_byclass$condition == "class2_obesity", ]
      class3_params <- obesity_rr_acdeath_byclass[obesity_rr_acdeath_byclass$condition == "class3_obesity", ]
    
    # Identify obese and non-obese individuals
    obesity_multiplier <- rep(1, length(inc_probs))
    obese_idx <- which(obesity_vals == 1)
    if (length(obese_idx) > 0) {
      # Class-based multiplier:
      mult_class <- rep(1, length(obese_idx))
      idx_class1 <- obese_idx[ which(class1_vals[obese_idx] == 1) ]
      if (length(idx_class1) > 0) {
        mult_class[match(idx_class1, obese_idx)] <- class1_params$value[1]
        }
      idx_class2 <- obese_idx[ which(class2_vals[obese_idx] == 1) ]
      if (length(idx_class2) > 0) {
          mult_class[match(idx_class2, obese_idx)] <- class2_params$value[1]
      }
      idx_class3 <- obese_idx[ which(class3_vals[obese_idx] == 1) ]
      if (length(idx_class3) > 0) {
          mult_class[match(idx_class3, obese_idx)] <- class3_params$value[1]
      }
      
      # Sex-based multiplier:
      mult_sex <- rep(1, length(living_idx))
      cohort_female <- female_vals[1]
      if (cohort_female == 1) {
          mult_sex[] <- obesity_rr_acdeath_bysex$value[obesity_rr_acdeath_bysex$female == 1][1]
      } else if (cohort_female == 0) {
          mult_sex[] <- obesity_rr_acdeath_bysex$value[obesity_rr_acdeath_bysex$female == 0][1]
      }
      
      # Age group multiplier:
      age_group <- ifelse(age_vals >= 35 & age_vals < 50, "adult",
                          ifelse(age_vals >= 70 & age_vals <= 89, "older adult", NA))
      mult_age <- rep(1, length(living_idx))
      adult_idx <- which(age_group == "adult")
      older_idx <- which(age_group == "older adult")
      
      if (length(adult_idx) > 0) {
        adult_params <- obesity_rr_acdeath_byage[obesity_rr_acdeath_byage$age_cat == "adult", ]
        mult_age[adult_idx] <- adult_params$value[1]
      }
      if (length(older_idx) > 0) {
        older_params <- obesity_rr_acdeath_byage[obesity_rr_acdeath_byage$age_cat == "older_adult", ]        
        mult_age[older_idx] <- older_params$value[1]
      }
      # Combined multiplier:
      obesity_multiplier[obese_idx] <- mult_class  * mult_sex[obese_idx] * mult_age[obese_idx]
    }
    
    ## 2) t2d ---- 
      t2d_multiplier <- rep(1, length(inc_probs))
      t2d_idx <- which(pd[["t2d"]][living_idx] == 1)
      
      if (length(t2d_idx) > 0) { 
        t2d_rr_row <- t2d_rr[t2d_rr$condition == "all_cause_death", ]
        t2d_multiplier[t2d_idx]<- t2d_rr_row$value[1]
      }
    
    ## 3) cvd ----
      # NOTE: stroke death included in GBD cvd-death
      cvd_multiplier <- rep(1, length(inc_probs))
      cvd_idx     <- which(pd[["cvd"]][living_idx] == 1)
      if (length(cvd_idx) > 0)  {
        inc_ac_death <- matched_params[["incidence_ac_death_val"]]
        death_cvd    <- matched_params[["death_cvd_val"]]
        
        # non-cvd death incidence:
        num <- inc_ac_death[living_idx][cvd_idx]
        den <- (inc_ac_death - death_cvd)[living_idx][cvd_idx]
        
        ok <- is.finite(num) & is.finite(den) & den > 1e-12
        cvd_multiplier[cvd_idx][ok] <- num[ok] / den[ok]
        cvd_multiplier[!is.finite(cvd_multiplier) | cvd_multiplier < 0] <- 1
      }
      
    ## 4) ckd ---- 
      ckd_multiplier <- rep(1, length(inc_probs))
      ckd_idx     <- which(pd[["ckd"]][living_idx] == 1)
      if (length(ckd_idx) > 0){
        death_ckd <- matched_params[["death_ckd_val"]]
        inc_ac_death <- matched_params[["incidence_ac_death_val"]]
        # non-ckd death incidence:
        num <- inc_ac_death[living_idx][ckd_idx]
        den <- (inc_ac_death - death_ckd)[living_idx][ckd_idx]
        
        ok <- is.finite(num) & is.finite(den) & den > 1e-12
        ckd_multiplier[ckd_idx][ok] <- num[ok] / den[ok]
        ckd_multiplier[!is.finite(ckd_multiplier) | ckd_multiplier < 0] <- 1
      }
    ## apply risk modifiers ----
    risk_params <- list(obesity_multiplier, t2d_multiplier, cvd_multiplier, ckd_multiplier)
    combined_multiplier <- Reduce(`*`, risk_params)
    inc_probs <- pmax(0, pmin(1, inc_probs * combined_multiplier))
    
    #DEBUG --------
    #debug helper
    check_vec <- function(x, name) {
      bad <- which(is.na(x) | !is.finite(x))
      if (length(bad)) {
        message("[", name, "] bad n=", length(bad),
                " (NA=", sum(is.na(x)), ", Inf=", sum(is.infinite(x)), ")",
                " example idx=", paste(head(bad, 10), collapse=","))
        return(bad)
      }
      integer(0)
    }
    bad_ob <- check_vec(obesity_multiplier, "obesity_multiplier")
    bad_t2d <- check_vec(t2d_multiplier, "t2d_multiplier")
    bad_cvd <- check_vec(cvd_multiplier, "cvd_multiplier")
    bad_ckd <- check_vec(ckd_multiplier, "ckd_multiplier")
    bad_comb <- check_vec(combined_multiplier, "combined_multiplier")
    #############'
    
    ### joint adjustment ------
      ### 1. t2d-obesity ----
      joint_idx <- which(pd$obesity[living_idx] == 1 & pd$t2d[living_idx] == 1)
      if (length(joint_idx) && all(pd$obese_t2d_corr > 0)) {
        obese_t2d_corr <- pd[["obese_t2d_corr"]][living_idx][joint_idx]
        
        # Calculate weighted joint multiplier using the baseline correlation
        joint_mult <- (1 - obese_t2d_corr) * (obesity_multiplier[joint_idx] * t2d_multiplier[joint_idx]) +
          obese_t2d_corr * pmax(obesity_multiplier[joint_idx], t2d_multiplier[joint_idx])
        # Compute adjustment factor relative to the product of multipliers
        d = (obesity_multiplier[joint_idx] * t2d_multiplier[joint_idx])
        d[d<=  1e-9] <- 1e-9 
        adj_factor <- joint_mult / d
        
        inc_probs[joint_idx] <- pmin(1, inc_probs[joint_idx] * adj_factor)
      }
    
      ### 2. t2d-cvd ----
      joint_idx <- which(pd$cvd[living_idx] == 1 & pd$t2d[living_idx] == 1)
    
      if (length(joint_idx) > 0 && all(pd$cvd_t2d_corr > 0)) {
        cvd_t2d_corr <- pd[["cvd_t2d_corr"]][living_idx][joint_idx]
        
        # Calculate weighted joint multiplier using the baseline correlation
        joint_mult <- (1 - cvd_t2d_corr) * (cvd_multiplier[joint_idx] * t2d_multiplier[joint_idx]) +
          cvd_t2d_corr * pmax(cvd_multiplier[joint_idx], t2d_multiplier[joint_idx])
        # Compute adjustment factor relative to the product of multipliers
        d = (cvd_multiplier[joint_idx] * t2d_multiplier[joint_idx])
        d[d<=  1e-9] <- 1e-9 
        adj_factor <- joint_mult / d
        inc_probs[joint_idx] <- pmin(1, inc_probs[joint_idx] * adj_factor)
      }
    
      ### 3. t2d-ckd ----
      joint_idx <- which(pd$ckd[living_idx] == 1 & pd$t2d[living_idx] == 1)
      if (length(joint_idx)>0 && all(pd$ckd_t2d_corr > 0)) {
        ckd_t2d_corr <- pd[["ckd_t2d_corr"]][living_idx][joint_idx]
        joint_mult <- (1 - ckd_t2d_corr) * 
          (t2d_multiplier[joint_idx] * ckd_multiplier[joint_idx]) +
          ckd_t2d_corr * pmax(t2d_multiplier[joint_idx], ckd_multiplier[joint_idx])
        
        d = (t2d_multiplier[joint_idx] * ckd_multiplier[joint_idx])
        d[d<=  1e-9] <- 1e-9 
        adj_factor <- joint_mult / d
        
        inc_probs[joint_idx] <- pmin(1, inc_probs[joint_idx] * adj_factor)
      }
      # Enforce minimum risk of ac_death for ESKD ---- 
      eskd_idx <- which(pd[["ckd"]][living_idx] == 1 & pd[["eskd"]][living_idx] == 1)  # among ckd individuals
      non_eskd_idx <- which(pd[["ckd"]][living_idx] == 1 & pd[["eskd"]][living_idx] == 0)
      original_avg <- mean(inc_probs[ckd_idx], na.rm = TRUE)
      N_total <- length(ckd_idx)
      inc_probs[eskd_idx] <- pmax(inc_probs[eskd_idx], (1-sqrt(min_death_risk_eskd))) 
      
      # Balance CKD non-ESKD to maintain average incidence
      S_eskd <- sum(inc_probs[eskd_idx], na.rm = TRUE)
      S_non  <- sum(inc_probs[non_eskd_idx], na.rm = TRUE)
      f <- (original_avg * N_total - S_eskd) / S_non
      f <- max(f, 0) # no negative incidence 
      inc_probs[non_eskd_idx] <- inc_probs[non_eskd_idx] * f
  }
  return(inc_probs)
}

  
