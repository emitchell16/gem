# -------------------'
#' Generate baseline year 0 data for model validation studies
# -------------------'

v_initialize_baseline <- function(v_input_data, paths, config = config_inputs){
  tol <- 1e-12
  cat("\n making baseline for study/cohort", v_input_data$study_id, "/", v_input_data$cohort_id, "\n")
  n <- as.integer(v_input_data$sample_size[[1]]) 
  
  inc_prev_input_params <- read_excel(file.path(paths$data, "GEM_input_parameters.xlsx"), sheet = "baseline_dataset") |>
    select(parameter, value)
  param <- setNames(as.numeric(inc_prev_input_params$value), inc_prev_input_params$parameter)
  
  ## region and income group assignment block----
  regions <- c("africa","americas","eastern mediterrean","europe","south-east asia","western pacific")
  active_regions <- regions[vapply(regions, function(col) isTRUE(as.numeric(v_input_data[[col]][[1]]) == 1), logical(1))]
  
  who_region_vec <- if (identical(v_input_data$iso3[[1]], "MULTI") && length(active_regions) > 0) {
    # random split among the active regions
    sample(active_regions, size = n, replace = TRUE)
  } else if (identical(v_input_data$iso3[[1]], "USA")) {
    rep("americas", n)
  } else if (length(active_regions) >= 1) {
    # pick the single flagged region (if multiple, take the first)
    rep(active_regions[1], n)
  } 
  
  wb_income_group_vec <- if (isTRUE(v_input_data$hic[[1]] == 1) && 
                             (isTRUE(v_input_data$lmic[[1]] == 1))) {
    half1 <- floor(n/2); half2 <- n - half1
    sample(c(rep("hic", half1), rep("lmic", half2)))
  } else if (isTRUE(v_input_data$hic[[1]] == 1)) {
    rep("hic", n)
  } else if (isTRUE(v_input_data$lmic[[1]] == 1)) {
    rep("lmic", n)
  } 
  
  ## initialize person level data table -----
  working_dt <- data.table(
    pid       = seq_len(n),
    study_id  = v_input_data$study_id[[1]],
    cohort_id = v_input_data$cohort_id[[1]],
    
    who_region = who_region_vec,
    wb_income_group = wb_income_group_vec,
    
    age       = round(rtruncnorm(
      n,
      a = v_input_data$min_age[[1]],
      b = v_input_data$max_age[[1]],
      mean = v_input_data$avg_age[[1]],
      sd = v_input_data$sd_age[[1]]
    ), 0),
    female    = rbinom(n, 1, prob = {
      p <- as.numeric(v_input_data$female[[1]])
      if (is.na(p)) p <- 0.5
      p <- min(max(p, 0), 1)
      p
    }),
    ac_death  = 0L,
    cvd_death = 0L,
    mace      = 0L,
    stroke    = 0L
  )
  
  
  # obesity assignment----
      ## note: obesity/bmi assigned first because it is the most reliably reported in validation input sources 
  r_ob <- v_input_data$obesity_rate[[1]]
  r1   <- v_input_data$class_1_obesity_rate[[1]]
  r2   <- v_input_data$class_2_obesity_rate[[1]]
  r3   <- v_input_data$class_3_obesity_rate[[1]]
  bmin <- v_input_data$min_bmi[[1]]
  bmax <- v_input_data$max_bmi[[1]]
  
  if (any(is.na(c(r1, r2,r3)))){
    # if category rates missing, apply truncated normal distribution over bmi data
    working_dt[, bmi  := rtruncnorm(
      n,
      a = v_input_data$min_bmi[[1]],
      b = v_input_data$max_bmi[[1]],
      mean = v_input_data$avg_bmi[[1]],
      sd = v_input_data$sd_bmi[[1]]
    )]  
    
    working_dt[, `:=` (
      class3_obesity = as.integer(bmi >= 40),
      class2_obesity = as.integer(bmi >= 35 & bmi < 40),
      class1_obesity = as.integer(bmi >= 30 & bmi < 35),
      obesity = as.integer(bmi >= 30)
    )]
  } else {
    sum_ob <- r1 + r2 + r3
    if (sum_ob <= 0 && r_ob > 0) stop("Class rates sum to 0 but obesity_rate > 0")
    if (abs(sum_ob - r_ob) > 1e-6) {
      # normalize classes to match total obesity rate
      r1 <- r1 * (r_ob / sum_ob); r2 <- r2 * (r_ob / sum_ob); r3 <- r3 * (r_ob / sum_ob)
    }
    
    n1 <- round(n *  r1); n2 <- round(n * r2);    n3 <- round(n * r3)
    n0 <-  n - (n1 + n2 + n3)
    tot <- n0 + n1 + n2 + n3
    # balance to sum to n exactly 
    if (tot != n) {
      diff <- n - tot
      buckets <- c("n0","n1","n2","n3"); j <- 1
      while (diff != 0) {
        b <- buckets[(j - 1) %% 4 + 1]; val <- get(b)
        if (diff > 0) { assign(b, val + 1); diff <- diff - 1 }
        else if (val > 0) { assign(b, val - 1); diff <- diff + 1 }
        j <- j + 1
      }}
    
    cls <- integer(n)
    idx <- sample.int(n); pos <- 1
    if (n3 > 0) { cls[idx[pos:(pos + n3 - 1)]] <- 3; pos <- pos + n3 }
    if (n2 > 0) { cls[idx[pos:(pos + n2 - 1)]] <- 2; pos <- pos + n2 }
    if (n1 > 0) { cls[idx[pos:(pos + n1 - 1)]] <- 1; pos <- pos + n1 }
    
    # uniform BMI within class ranges, intersected with [bmin, bmax]
    draw_unif <- function(k, a, b) {
      if (k <= 0) return(numeric(0))
      a2 <- max(a, bmin); b2 <- min(b, bmax)
      if (b2 <= a2) return(rep(a2, k))
      runif(k, a2, b2)
    }
    bmi <- numeric(n)
    i0 <- which(cls == 0); if (length(i0)) bmi[i0] <- draw_unif(length(i0), bmin, 29.99)
    i1 <- which(cls == 1); if (length(i1)) bmi[i1] <- draw_unif(length(i1), 30,   34.99)
    i2 <- which(cls == 2); if (length(i2)) bmi[i2] <- draw_unif(length(i2), 35,   39.99)
    i3 <- which(cls == 3); if (length(i3)) bmi[i3] <- draw_unif(length(i3), 40,   bmax)
    
    working_dt[, `:=`(
      bmi             = bmi,
      class3_obesity = as.integer(cls == 3),
      class2_obesity = as.integer(cls == 2),
      class1_obesity = as.integer(cls == 1),
      obesity         = as.integer(cls > 0)
    )]
  }
  
  ## add conditions based on baseline input data or imputed from GBD/NCD-risc studies  -----
  temp_external_cohort_data <- cohort_input_data |> 
    left_join(region_income_country_groups |> select(iso3, who_region), by = "iso3") |> 
    filter(age_high >= v_input_data$min_age[[1]], age_low<= v_input_data$max_age[[1]]) %>%
    { 
      iso3_val <- v_input_data$iso3[[1]]
      if (!is.na(iso3_val) && nchar(iso3_val) == 3) filter(., iso3 == iso3_val) else .
    } %>%
    group_by(who_region, female, age_low, age_high) |> 
    summarise(
      prevalence_t2d_val = weighted.mean(prevalence_t2d_val, w = pop_size, na.rm = TRUE),
      prevalence_ckd_val = weighted.mean(prevalence_ckd_val, w = pop_size, na.rm = TRUE),
      prevalence_cvd_val = weighted.mean(prevalence_cvd_val, w = pop_size, na.rm = TRUE),
      pop_size           = sum(pop_size, na.rm = TRUE),
      .groups = "drop")
  
  ### T2d assignment ----
  rate_val <- v_input_data$t2d_rate[[1]]
  proportion_t2d_obese <- param[["proportion_t2d_obese"]]       #P(obesity|t2d)
  
  if (!is.na(rate_val)) {
    p_total <- pmin(pmax(rate_val, 0), 1)
    if (p_total >= 1 - tol) {working_dt[, t2d := 1L]
    } else if (p_total <= tol) {working_dt[, t2d := 0L]
    } else {
      # increase P(T2D) for obese, decrease for non-obese to maintain overall p_total
      obesity_rate <- mean(working_dt$obesity == 1)
      
      if (obesity_rate > 0 && obesity_rate < 1) {
        # obtain group targets
        p1 <- p_total * proportion_t2d_obese / obesity_rate              # P(T2D|obese)
        p1 <- pmin(pmax(p1, 0), 1)
        p0 <- (p_total - obesity_rate * p1) / (1 - obesity_rate) # P(T2D|non-obese)
        p0 <- pmin(pmax(p0, 0), 1)
        
        working_dt[obesity == 1, t2d := rbinom(.N, 1, p1)]
        working_dt[obesity == 0, t2d := rbinom(.N, 1, p0)]
        
      } else { # all or no obesity:
        working_dt[, t2d := rbinom(.N, 1, p_total)]
      }}
    # missing input data for cohort rate
  } else {
    # non-equi join to impute rate: match age, region, sex
    temp_external_cohort_data <- as.data.table(temp_external_cohort_data)
    wd <- working_dt[, .(pid, female, age, age_low = age, age_high = age)]
    setkey(temp_external_cohort_data, female, age_low, age_high)
    setkey(wd,  female, age_low, age_high)
    join_dt <- foverlaps(wd, temp_external_cohort_data, type = "within", nomatch = 0L)
    prev_col <- "prevalence_t2d_val"
    prob_by_pid <- join_dt[, .(prob = mean(get(prev_col), na.rm = TRUE)), by = pid]
    prob_vec <- prob_by_pid$prob[match(working_dt$pid, prob_by_pid$pid)]
    prob_vec <- pmin(pmax(prob_vec, 0), 1)
    
    # increase P(T2D) for obese, decrease for non-obese to maintain overall p_total
    obesity_rate <- mean(working_dt$obesity == 1)
    p_adj        <- prob_vec
    p_total      <- mean(prob_vec)
    
    if (obesity_rate > 0 && obesity_rate < 1 && any(!is.na(prob_vec))) {
      # obtain group targets
      p1 <- pmin(pmax(p_total * proportion_t2d_obese / obesity_rate, 0), 1)                  # P(T2D | obese)
      p0 <- pmin(pmax((p_total - obesity_rate * p1) / (1 - obesity_rate), 0), 1) # P(T2D | non-obese)
      
      mean_ob  <- mean(prob_vec[working_dt$obesity == 1], na.rm = TRUE)
      mean_non <- mean(prob_vec[working_dt$obesity == 0], na.rm = TRUE)
      delta_plus  <- p1 - mean_ob
      delta_minus <- mean_non - p0
      
      p_adj[working_dt$obesity == 1] <- pmin(pmax(prob_vec[working_dt$obesity == 1] + delta_plus, 0), 1)
      p_adj[working_dt$obesity == 0] <- pmin(pmax(prob_vec[working_dt$obesity == 0] - delta_minus, 0), 1)
      
    } else if (isTRUE(obesity_rate == 1)) {
      p_adj <- pmin(pmax(p_adj, proportion_t2d_obese), 1) # if everyone is obese in cohort, take minimum proportion_t2d_obese from matched imputation rates
    }
    # draw one per row based on that rows probability
    working_dt[, t2d := rbinom(.N, 1, p_adj)] 
  } 
  
  ### Cvd and ckd assignment ----
  p_t2d_input <- v_input_data$t2d_rate[[1]]                 # may be NA
  r_t2d       <- mean(working_dt$t2d == 1, na.rm = TRUE)  #  T2D prevalence
  temp_external_cohort_data <- as.data.table(temp_external_cohort_data)
  proportion_t2d_cvd  <- param[["proportion_t2d_cvd"]]
  proportion_t2d_ckd <- param[["proportion_t2d_ckd"]]
  
  for (cond in c("cvd", "ckd")) {
    rate_col <- paste0(cond, "_rate")
    rate_val <- v_input_data[[rate_col]][[1]]   
    target_p <- if (cond == "cvd") proportion_t2d_cvd else proportion_t2d_ckd
    if (!is.na(rate_val)) {
      p_total <- pmin(pmax(rate_val, 0), 1)
      
      if ( p_total >= 1 - tol) {
        working_dt[, (cond) := 1L]
        next
      }
      # split by T2D when the input T2D rate exists and is < 1
      if (!is.na(p_t2d_input) && p_t2d_input < 1 - tol && r_t2d < 1 - tol && r_t2d > tol) {
        # P(cond | T2D) = target_p; solve P(cond | non-T2D) to keep overall at p_total
        p1 <- target_p
        p0 <- (p_total - r_t2d * p1) / (1 - r_t2d)
        p0 <- min(max(p0, 0), 1)
        
        # vectorized per-row probs
        p_i <- ifelse(working_dt$t2d == 1, p1, p0)
        working_dt[, (cond) := rbinom(.N, 1, p_i)]
      } else {
        # t2d_rate == 1 uniform at p_total
        working_dt[, (cond) := rbinom(.N, 1, p_total)]
      }
    } else {  # no cohort rate; impute from external sources and modify based on t2d status
      wd <- working_dt[, .(pid, female, age, age_low = age, age_high = age)]
      setkey(temp_external_cohort_data, female, age_low, age_high)
      setkey(wd,  female, age_low, age_high)
      join_dt <- foverlaps(wd, temp_external_cohort_data, type = "within", nomatch = 0L)
      prev_col <- paste0("prevalence_", cond, "_val")
      prob_by_pid <- join_dt[, .(prob = mean(get(prev_col), na.rm = TRUE)), by = pid]
      prob_vec <- prob_by_pid$prob[match(working_dt$pid, prob_by_pid$pid)]
      prob_vec <- pmin(pmax(prob_vec, 0), 1)
      
      p_total_ext <- mean(prob_vec, na.rm = TRUE)
      p_adj <- prob_vec
      
      if (r_t2d > tol && r_t2d < 1 - tol && any(!is.na(prob_vec))) {
        p1 <- target_p
        p0_target <- (p_total_ext - r_t2d * p1) / (1 - r_t2d)
        p0_target <- min(max(p0_target, 0), 1)
        
        mean_t2d  <- mean(prob_vec[working_dt$t2d == 1], na.rm = TRUE)
        mean_non  <- mean(prob_vec[working_dt$t2d == 0], na.rm = TRUE)
        delta_plus  <- p1 - mean_t2d
        delta_minus <- mean_non - p0_target
        
        p_adj[working_dt$t2d == 1] <- pmin(pmax(prob_vec[working_dt$t2d == 1] + delta_plus, 0), 1)
        p_adj[working_dt$t2d == 0] <- pmin(pmax(prob_vec[working_dt$t2d == 0] - delta_minus, 0), 1)
      } else if (r_t2d >= 1 - tol) {
        # everyone T2D → ensure at least target_p
        p_adj <- pmax(p_adj, target_p)
      } 
      working_dt[, (cond) := rbinom(.N, 1, p_adj)]
    }
  }
  ### ESKD assignment ----  
  eskd_val <- v_input_data$eskd_rate[[1]]
  if (!is.na(eskd_val)) {
    r_ckd <- mean(working_dt$ckd == 1, na.rm = TRUE)  # P(CKD=1)
    p_ckd <- pmin(pmax(eskd_val / r_ckd, 0), 1)      # P(ESKD=1 | CKD=1)
    working_dt[, eskd := 0L]
    working_dt[ckd == 1, eskd := rbinom(.N, 1, p_ckd)]
  } else { # imputation by income group
    prob_eskd_lmic <- param[["prob_eskd_lmic"]]  
    prob_eskd_hic  <- param[["prob_eskd_hic"]]
    prob_eskd_na   <- param[["prob_eskd_na"]]
    working_dt[, prob_eskd := fifelse(
      wb_income_group == "lmic", prob_eskd_lmic,
      fifelse(wb_income_group == "hic", prob_eskd_hic, prob_eskd_na)
    )]
    working_dt[, eskd := 0L]
    working_dt[ckd == 1, eskd := rbinom(.N, 1, prob_eskd)]
    working_dt[, prob_eskd := NULL]
  }
  return(working_dt)
}