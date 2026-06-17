# -------------------'
# Functions to get disability weights by disease state
# Author: Liz Mitchell
# Date updated: March 2026
# -------------------'

##***************** 
## T2D -----
assign_dw_t2d <- function(updated_data, living_idx, sensitivity) { 
  dw_vec <- rep(0, length(living_idx))
  t2d_idx <- which(updated_data[["t2d"]][living_idx] == 1)
  if (length(t2d_idx) == 0) return(dw_vec)
  
  # identify input data rows
  row_uncomp <- dw_tbl[         # uncomplicated baseline dw
    grepl("t2d", condition_logic, ignore.case = TRUE) &
      grepl("generic uncomplicated", healthstate_name, ignore.case = TRUE)]
  row_comp1_retino <- dw_tbl[
    grepl("t2d", condition_logic, ignore.case = TRUE) &
      grepl("retinopathy complications", healthstate_name, ignore.case = TRUE) ]
  row_comp2_neuro <- dw_tbl[
    grepl("t2d", condition_logic, ignore.case = TRUE) &
      grepl("neuropathy and foot", healthstate_name, ignore.case = TRUE)]

  dw_uncomp <- row_uncomp$mean[1]
  dw_retin  <- row_comp1_retino$mean[1]
  dw_neuro  <- row_comp2_neuro$mean[1]
  if (sensitivity) { 
    dw_uncomp <- row_uncomp$lower[1]
    dw_retin  <- row_comp1_retino$lower[1]
    dw_neuro  <- row_comp2_neuro$lower[1]
  }
  
  # complication probabilities
  p_retin <- row_comp1_retino$probability_weight[1] 
  p_neuro <- row_comp2_neuro$probability_weight[1] 
  
  # assign probability-weighted average DW to all with T2D
  dw_t2d_avg <- dw_uncomp + (p_retin * dw_retin) + (p_neuro * dw_neuro)
  dw_vec[t2d_idx] <- dw_t2d_avg
  dw_vec
}
##*****************
## CVD ----
assign_dw_cvd <- function(updated_data, living_idx, sensitivity) { 
 
  dw_vec  <- rep(0, length(living_idx))
  cvd_idx <- which(updated_data[["cvd"]][living_idx] == 1)
  if (length(cvd_idx) == 0) return(dw_vec)

  row_avg <- dw_tbl[         # uncomplicated baseline dw
    grepl("cvd", condition_logic, ignore.case = TRUE) &
      grepl("computed average, cvd", healthstate_name, ignore.case = TRUE)]
  dw_avg  <-row_avg$mean[1]
  dw_vec[cvd_idx] <- dw_avg
  
  if (sensitivity) { 
    # overwrite with lower bound in sensitivity analysis
    dw_lower <- row_avg$lower[1]
    dw_vec[cvd_idx] <- dw_lower
  }
  dw_vec
}

## Stroke ----
assign_dw_stroke <- function(updated_data, living_idx, sensitivity) { 
  dw_vec <- rep(0, length(living_idx))
  stroke_idx <- which(updated_data[["stroke"]][living_idx] == 1)
  if (length(stroke_idx) == 0) return(dw_vec)
  
  chronic_rows <- dw_tbl[
    grepl("^stroke", condition_logic, ignore.case = TRUE) & !grepl("acute stroke", condition_logic, ignore.case = TRUE)]
  acute_rows <- dw_tbl[
    grepl("acute stroke", condition_logic, ignore.case = TRUE)]
  
  w_chronic <- chronic_rows$probability_weight / sum(chronic_rows$probability_weight)
  w_acute   <- acute_rows$probability_weight / sum(acute_rows$probability_weight)
 
  if (!isTRUE(sensitivity)) {
    dw_chronic_avg <- sum(w_chronic * 
                            vapply(seq_len(nrow(chronic_rows)), function(i) {
                              chronic_rows$mean[i]}, numeric(1)), na.rm = TRUE)
    dw_acute_avg <- sum(w_acute * 
                          vapply(seq_len(nrow(acute_rows)), function(i) {
                            acute_rows$mean[i]}, numeric(1)),na.rm = TRUE)
  } else {
    dw_chronic_avg <- sum(w_chronic * chronic_rows$lower, na.rm = TRUE)
    dw_acute_avg   <- sum(w_acute   * acute_rows$lower, na.rm = TRUE)
  }
  
  dw_vec[stroke_idx] <- dw_chronic_avg
  # incident stroke this year: acute
  stroke_this_year_idx <- which(
    updated_data[["stroke"]][living_idx] == 1 & updated_data[["stroke_this_year"]][living_idx] == 1)
  if (length(stroke_this_year_idx) > 0) {
    dw_vec[stroke_this_year_idx] <- dw_acute_avg
  }
  dw_vec
}

##*****************
## CKD & ESKD ---- 
# note: everyone with ESKD also CKD ==1
assign_dw_ckd <- function(updated_data, living_idx, sensitivity) { 
  income_group_map <- fread(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv"))
  dw_vec <- rep(0, length(living_idx))
  ckd_idx <- which(updated_data[["ckd"]][living_idx] == 1)
  if (length(ckd_idx) == 0) return(dw_vec)
  
  ### CKD non-ESKD rows ----
  row_asymp <- dw_tbl[
    grepl("ckd", condition_logic, ignore.case = TRUE) &
      grepl("asymptomatic", healthstate_name, ignore.case = TRUE)]
  row_stage3 <- dw_tbl[
    grepl("ckd", condition_logic, ignore.case = TRUE) &
      grepl("stage3", healthstate_name, ignore.case = TRUE)]
  row_stage4 <- dw_tbl[
    grepl("ckd", condition_logic, ignore.case = TRUE) &
      grepl("stage4", healthstate_name, ignore.case = TRUE)]
  
  ### ESKD rows ----
  row_eskd_stage5 <- dw_tbl[
    grepl("eskd", condition_logic, ignore.case = TRUE) &
      grepl("stage5", healthstate_name, ignore.case = TRUE)]
  row_eskd_transpl <- dw_tbl[
    grepl("eskd", condition_logic, ignore.case = TRUE) &
      grepl("transplant", healthstate_name, ignore.case = TRUE)]
  row_eskd_dialys <- dw_tbl[
    grepl("eskd", condition_logic, ignore.case = TRUE) &
      grepl("dialysis", healthstate_name, ignore.case = TRUE)]
  
  if (!isTRUE(sensitivity)) {
    dw_asymp        <- row_asymp$mean[1]
    dw_stage3       <- row_stage3$mean[1]
    dw_stage4       <- row_stage4$mean[1]
    dw_eskd_stage5  <- row_eskd_stage5$mean[1]
    dw_eskd_transpl <- row_eskd_transpl$mean[1]
    dw_eskd_dialys  <- row_eskd_dialys$mean[1]
  } else {
    dw_asymp        <- row_asymp$lower[1]
    dw_stage3       <- row_stage3$lower[1]
    dw_stage4       <- row_stage4$lower[1]
    dw_eskd_stage5  <- row_eskd_stage5$lower[1]
    dw_eskd_transpl <- row_eskd_transpl$lower[1]
    dw_eskd_dialys  <- row_eskd_dialys$lower[1]
  }
  w_ckd <- c(
    row_asymp$probability_weight[1],
    row_stage3$probability_weight[1],
    row_stage4$probability_weight[1]
  )
  w_ckd <- w_ckd / sum(w_ckd)
  
  dw_ckd_non_eskd_avg <- sum(w_ckd * c(dw_asymp, dw_stage3, dw_stage4))
  dw_vec[ckd_idx] <- dw_ckd_non_eskd_avg
  
  # overwrite for eskd cases
  eskd_idx <- which(updated_data[["eskd"]][living_idx] == 1)
  
  if (length(eskd_idx) > 0) {
    
    # map iso3 -> income group
    iso_vals <- updated_data[["iso3"]][living_idx][eskd_idx]
    income_grp <- income_group_map$world_bank_income_group[match(iso_vals, income_group_map$iso3)]
    
    # treatment probability by income group
    p_treated <- rep(0, length(eskd_idx))
    
    p_treated[grepl("low income", income_grp, ignore.case = TRUE)] <-
      input_params$value[input_params$parameter == "prob_eskd_treated_lic"]
    
    p_treated[grepl("lower middle", income_grp, ignore.case = TRUE)] <-
      input_params$value[input_params$parameter == "prob_eskd_treated_lmic"]
    
    p_treated[grepl("upper middle", income_grp, ignore.case = TRUE)] <-
      input_params$value[input_params$parameter == "prob_eskd_treated_umic"]
    
    p_treated[grepl("high", income_grp, ignore.case = TRUE)] <-
      input_params$value[input_params$parameter == "prob_eskd_treated_hic"]
    
    p_treated[is.na(p_treated)] <- 0
    p_treated <- pmin(1, pmax(0, p_treated))
    
    # treated vs untreated
    treated <- rbinom(length(eskd_idx), 1, p_treated)
    
    # among treated, modality assignment
    p_transpl <- row_eskd_transpl$probability_weight[1]
    p_dialys  <- row_eskd_dialys$probability_weight[1]
    
    mod_w <- c(p_transpl, p_dialys)
    mod_w <- mod_w / sum(mod_w)
    p_transpl <- mod_w[1]
    p_dialys  <- mod_w[2]
    
    # weighted-average ESKD disability weight
    dw_eskd <- (1 - p_treated) * dw_eskd_stage5 +
      p_treated * (p_transpl * dw_eskd_transpl + p_dialys * dw_eskd_dialys)
    
    dw_vec[eskd_idx] <- dw_eskd
  }
  dw_vec
}
