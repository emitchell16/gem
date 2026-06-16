# -------------------'
# GEM Model - GLP-1 obesity & ESKD adjustment functions
# Last major update: Feb. 2026
# -------------------'

##*****************
# Function to update ESKD risk with treatment
# ESKD incidence rate is not cohort-specific (current data limitation)

apply_eskd <- function(indices, glp1_outcome_impact_row, 
                       base_prob, 
                       dc_factor,
                       years_treated_vec,
                       dc_threshold) {
  multiplier <- rep(glp1_outcome_impact_row$outcome_val[1], length(indices))
  res <- numeric(length(indices))
  
  idx_low  <- which(years_treated_vec < dc_threshold)
  idx_high <- which(years_treated_vec >= dc_threshold)
  
  # For individuals with years_treated < threshold, use full multiplier
  if (length(idx_low) > 0) {
    prob_low <-  pmin(1, base_prob[idx_low] * multiplier[idx_low])
    res[idx_low] <- rbinom(length(idx_low), size = 1, prob = prob_low)
  }
  # For individuals with years_treated >= threshold, apply double counting adjustment
  if (length(idx_high) > 0) {
    adj_mult <- 1 + (multiplier[idx_high] - 1) * dc_factor
    prob_high <- pmin(1, base_prob[idx_high] * adj_mult)
    res[idx_high] <- rbinom(length(idx_high), 1, prob_high)
  }
  return(res)
}

##*****************
# Function to update obesity status based on prop BMI reduction 
# Uses glp1_impact_data$outcome_unit == "bmi (proportion change)" 

update_obesity_transition_glp1 <- function(new_state_vals, 
                                           pd ,
                                           living_idx ,
                                           treated_in_obese,  # (overweight and obese: bmi 27+)
                                           glp1_impact_data = glp1_impact_data) {
  
  setDT(glp1_impact_data)
  # internal helper functions ----- 
  # function to stratify obesity class-group by age group  -----
  age_group_indices <- function(group_index) {
    list(
      adolescent = group_index[pd$age[group_index] < 18],
      adult      = group_index[pd$age[group_index] >= 18]
    )
  }
  
  # function to process each obesity class ------
  update_obesity_vars <- function(class_index, t2d_yn) {
    if (t2d_yn == "no") {
      class_index_no_t2d <- class_index
      
      if (length(class_index_no_t2d) > 0) {
        # age stratify
        age_groups <- age_group_indices(class_index_no_t2d)
        
        # initialize new bmi using current
        new_bmi <- pd[["bmi"]][class_index_no_t2d]
        
        # adolescents: ----
        if (length(age_groups$adolescent) > 0) {
          
          impact_adolescent <- glp1_impact_data[outcome_unit == "bmi (proportion change)" &
              treated_pop == "obese without t2d" & age_group == "adolescent" ]
          delta_adolescent <- rep(impact_adolescent$outcome_val[1], length(age_groups$adolescent))
          
          # find relative positions and update new_bmi
          adolescent_position <- match(age_groups$adolescent, class_index_no_t2d)
          new_bmi[adolescent_position] <- new_bmi[adolescent_position] * (1 + delta_adolescent)
        }
        # adults: ---------
        if (length(age_groups$adult) > 0) {
          impact_adult <- glp1_impact_data[outcome_unit == "bmi (proportion change)"&
              treated_pop == "obese without t2d" & age_group == "adult" ]
          delta_adult <- rep(impact_adult$outcome_val[1], length(age_groups$adult))
          adult_position <- match(age_groups$adult, class_index_no_t2d)
          new_bmi[adult_position] <- new_bmi[adult_position] * (1 + delta_adult)
        }
        
        # Update obesity class, bmi, and overall obesity assignments  --------
        # new index masks
        mask0 <- new_bmi < 30
        mask1 <- new_bmi >= 30 & new_bmi < 35
        mask2 <- new_bmi >= 35 & new_bmi < 40
        mask3 <- new_bmi >= 40
        
        
        # Non-obese transition: 
        rows0 <- class_index_no_t2d[mask0]
        if (length(rows0)) {
          new_state_vals[rows0, `:=`(
            obesity        = 0L,
            class1_obesity = 0L,
            class2_obesity = 0L,
            class3_obesity = 0L,
            bmi            = new_bmi[mask0]
          )]
        }
        
        # For BMI between 30 and 35: Class 1 obesity
        rows1 <- class_index_no_t2d[mask1]
        if (length(rows1)) {
          new_state_vals[rows1, `:=`(
            obesity        = 1L,
            class1_obesity = 1L,
            class2_obesity = 0L,
            class3_obesity = 0L,
            bmi            = new_bmi[mask1]
          )]
        }
        
        # For BMI between 35 and 40: Class 2 obesity
        rows2 <- class_index_no_t2d[mask2]
        if (length(rows2)) {
          new_state_vals[rows2, `:=`(
            obesity        = 1L,
            class1_obesity = 0L,
            class2_obesity = 1L,
            class3_obesity = 0L,
            bmi            = new_bmi[mask2]
          )]
        }
        
        
        # For BMI 40 and above: Class 3 obesity
        rows3 <- class_index_no_t2d[mask3]
        if (length(rows3)) {
          new_state_vals[rows3, `:=`(
            obesity        = 1L,
            class1_obesity = 0L,
            class2_obesity = 0L,
            class3_obesity = 1L,
            bmi            = new_bmi[mask3]
          )]
        }
        
        # return updated values
        return(new_state_vals)
      } 
    } else if (t2d_yn == "yes") {
      class_index_t2d <- class_index
      
      if (length(class_index_t2d) > 0) {
        impact_t2d <- glp1_impact_data[ glp1_impact_data$outcome_unit == "bmi (proportion change)" &
                                        glp1_impact_data$treated_pop == "t2d and obese", ]
        delta_t2d <- rep(impact_t2d$outcome_val[1], length(class_index_t2d))
        
        new_bmi <- pd[["bmi"]][class_index_t2d]
        new_bmi <- new_bmi * ( 1 + delta_t2d)
        
        # Update obesity class, bmi, and overall obesity assignments  --------
        mask0 <- new_bmi < 30
        mask1 <- new_bmi >= 30 & new_bmi < 35
        mask2 <- new_bmi >= 35 & new_bmi < 40
        mask3 <- new_bmi >= 40
        # Non‑obese
        rows0 <- class_index_t2d[mask0]
        if (length(rows0)) {
          new_state_vals[rows0, `:=`(
            obesity        = 0L,
            class1_obesity = 0L,
            class2_obesity = 0L,
            class3_obesity = 0L,
            bmi            = new_bmi[mask0]
          )]
        }
        
        # Class 1
        rows1 <- class_index_t2d[mask1]
        if (length(rows1)) {
          new_state_vals[rows1, `:=`(
            obesity        = 1L,
            class1_obesity = 1L,
            class2_obesity = 0L,
            class3_obesity = 0L,
            bmi            = new_bmi[mask1]
          )]
        }
        
        # Class 2
        rows2 <- class_index_t2d[mask2]
        if (length(rows2)) {
          new_state_vals[rows2, `:=`(
            obesity        = 1L,
            class1_obesity = 0L,
            class2_obesity = 1L,
            class3_obesity = 0L,
            bmi            = new_bmi[mask2]
          )]
        }
        
        # Class 3
        rows3 <- class_index_t2d[mask3]
        if (length(rows3)) {
          new_state_vals[rows3, `:=`(
            obesity        = 1L,
            class1_obesity = 0L,
            class2_obesity = 0L,
            class3_obesity = 1L,
            bmi            = new_bmi[mask3]
          )]
        }
        
        # return updated values
        return(new_state_vals)
      } 
    }
  }
  
  # sample BMI reduction for each class ----------
  
  ## class 1 ------------
  class1_idx <- treated_in_obese[pd[["class1_obesity"]][treated_in_obese] == 1]
  # 2 groups based on t2d status:
  c1_no_t2d_idx <- treated_in_obese[
    pd[["class1_obesity"]][treated_in_obese] == 1 & pd[["t2d"]][treated_in_obese] == 0]
  c1_t2d_idx <- treated_in_obese[
    pd[["class1_obesity"]][treated_in_obese] == 1 & pd[["t2d"]][treated_in_obese] == 1]
  
  # Process group without and with t2d ----
  if (length(c1_no_t2d_idx) > 0) {
    new_state_vals <- update_obesity_vars(c1_no_t2d_idx, t2d_yn = "no")
  }
  if (length(c1_t2d_idx) > 0) {
    new_state_vals <- update_obesity_vars(c1_t2d_idx, t2d_yn = "yes")
  }
  
  ## class 2  ------------
  class2_idx <- treated_in_obese[pd[["class2_obesity"]][treated_in_obese] == 1]
  c2_no_t2d_idx <- treated_in_obese[pd[["class2_obesity"]][treated_in_obese] == 1 & 
                                      pd[["t2d"]][treated_in_obese] == 0]
  c2_t2d_idx <- treated_in_obese[pd[["class2_obesity"]][treated_in_obese] == 1 & 
                                   pd[["t2d"]][treated_in_obese] == 1]
  # Process group without and with t2d ----
  if (length(c2_no_t2d_idx) > 0) {
    new_state_vals <- update_obesity_vars(c2_no_t2d_idx, t2d_yn = "no")
  }
  if (length(c2_t2d_idx) > 0) {
    new_state_vals <- update_obesity_vars(c2_t2d_idx, t2d_yn = "yes")
  }
  
  
  ## class 3:  ------------
  class3_idx <- treated_in_obese[pd[["class3_obesity"]][treated_in_obese] == 1]
  c3_no_t2d_idx <- treated_in_obese[pd[["class3_obesity"]][treated_in_obese] == 1 & 
                                      pd[["t2d"]][treated_in_obese] == 0]
  c3_t2d_idx <- treated_in_obese[pd[["class3_obesity"]][treated_in_obese] == 1 & 
                                   pd[["t2d"]][treated_in_obese] == 1]
  
  # Process group without and with t2d ----
  if (length(c3_no_t2d_idx) > 0) {
    new_state_vals <- update_obesity_vars(c3_no_t2d_idx, t2d_yn = "no")
  }
  if (length(c3_t2d_idx) > 0) {
    new_state_vals <- update_obesity_vars(c3_t2d_idx, t2d_yn = "yes")
  }
  
  #####################'
  ## t2d & overweight: BMI 27-29.9  --------
  t2d_overweight_idx <-  treated_in_obese[pd[["bmi"]][treated_in_obese] >= 27 & pd[["bmi"]][treated_in_obese]< 29.99]
  if (length(t2d_overweight_idx) > 0) {
    new_state_vals <- update_obesity_vars(t2d_overweight_idx, t2d_yn = "yes")
  } 
  
  return(new_state_vals) 
}