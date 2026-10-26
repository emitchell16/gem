# R/1_function_scripts/GEM_Model_Update_Cohort_States.R
# GEM model: yearly update of disease status, BMI, death, DALYs and costs for one cohort file
# Author: Liz Mitchell
# common random numbers; continuous BMI tracking; stroke as yearly event w/ recurrence; recurrent CVD events
# adults 18+; discounted-only DALYs (YLL annuity)
# intervention logic is not in this file: it enters through hook_*() functions (GEM_Model_Intervention_Hooks.R)

# helpers -----
# common random numbers: uniforms per person x event, seeded by cohort (crn_id), year, stream
#  identical across arms for the same cohort file -> between-arm differences reflect intervention only
#  crn_stream_offset (global, default 0) shifts all streams; run_psa() sets it to the run index
crn_u <- function(crn_id, year, n, k, stream = 1L) {
  off <- if (exists("crn_stream_offset", envir = .GlobalEnv)) get("crn_stream_offset", envir = .GlobalEnv) else 0L
  set.seed(as.integer(((crn_id + 7919L * off) %% 2e5) * 100 + year * 10 + stream))
  matrix(runif(n * k), nrow = n, ncol = k)
}

# stream 1 columns (disease events)
U_CKD <- 1; U_T2D <- 2; U_STROKE <- 3; U_CVD <- 4; U_ESKD <- 5; U_DEATH <- 6
# stream 2 and higher: reserved for intervention modules (crn_u(..., stream = 2L) etc.)

# external multipliers on event probabilities; validation replication of trial arms only ----
#  ext_te: named numeric (outcome -> HR); "mace" applies to cvd and stroke
apply_ext_te <- function(cond, p, ext_te) {
  if (is.null(ext_te)) return(p)
  k <- intersect(c(cond, if (cond %in% c("cvd", "stroke")) "mace"), names(ext_te))
  if (length(k)) pmin(1, p * ext_te[[k[1]]]) else p
}

# risk factors + intervention hook + external multipliers
modify_risk <- function(cond, inc_probs, matched_params, pd, living_idx, inc_val, ictx = NULL, ystate = NULL, ext_te = NULL) {
  inc_probs <- calculate_adjusted_risk(cond, inc_probs, matched_params, pd, living_idx, inc_val)
  inc_probs <- hook_modify_risk(cond, inc_probs, ictx, ystate, pd, living_idx)
  inc_probs <- apply_ext_te(cond, inc_probs, ext_te)
  if (anyNA(inc_probs))     stop("NA values in inc_probs for ", cond)
  if (any(inc_probs < 0))   stop("negative inc_probs for ", cond)
  pmin(1, inc_probs)
}

# MAIN FUNCTION ------
#  ictx: intervention context for this cohort file (hook_country_context()); NULL = no intervention
#  ext_te: named HRs applied to event probabilities (validation replication only)
update_cohort <- function(pd, cohort_params, year, ictx = NULL, open_cohort = FALSE, sim_temp_dir = NULL,
                          cal = acm_calibration, dw_sensitivity = FALSE, ext_te = NULL,
                          costing = isTRUE(config_sims$costing)) {
  ####
  # 1. Set up --------------
  ####
  new_state_vals <- copy(pd)
  living_idx <- which(new_state_vals$ac_death == 0)
  n_all <- nrow(pd)
  U  <- crn_u(pd$crn_id[1], year, n_all, 6, stream = 1L)          # disease events
  new_state_vals[, `:=`(stroke = 0L, cvd_event = 0L)]              # yearly event flags

  # match individual age & sex in cohort_params (age increases over modeled years so starting cohort can split)
  setDT(cohort_params); setkey(cohort_params, female, age_low)
  suppressWarnings({
    matched_params <- cohort_params[pd, on = .(female, age_low = age), roll = TRUE, rollends = c(TRUE, TRUE)]
  })

  ####
  # 1b. Intervention: start-of-year transitions (hook; state at end of last year) --------
  ####
  ystate <- hook_start_of_year(ictx, new_state_vals, pd, living_idx, year)

  # first-event risks from start-of-year state (pd); intervention effects through hook_modify_risk()
  mr <- function(cond, p) modify_risk(cond, p, matched_params, pd, living_idx, NULL, ictx, ystate, ext_te)
  te <- function(cond, p) apply_ext_te(cond, hook_modify_risk(cond, p, ictx, ystate, pd, living_idx), ext_te)
  inc_of <- function(cond) {
    v <- matched_params[[paste0("incidence_", cond, "_val")]]
    if (any(is.na(v))) stop(sprintf("NA incidence for '%s' (iso3=%s, female=%s)", cond, pd$iso3[1], pd$female[1]))
    pmin(1, pmax(0, v[living_idx]))
  }
  if (length(living_idx) > 0) {
  ####
  # 2. Incidence: CKD, T2D (first events among those without condition) -------
  ####
    for (cond in c("ckd", "t2d")) {
      p    <- mr(cond, inc_of(cond))
      ucol <- if (cond == "ckd") U_CKD else U_T2D
      pos  <- which(new_state_vals[[cond]][living_idx] == 0)
      if (length(pos)) {
        idx  <- living_idx[pos]
        new_state_vals[idx, (cond) := as.integer(U[idx, ucol] < p[pos])]
      }
    }
  ####
  # 3. Stroke (yearly event; first or recurrent) --------
  #    first: GBD first-ever stroke rate x joint RR; recurrent: annual prob by years since sim stroke (Mohan 2011)
  ####
    inc_stroke <- inc_of("stroke")
    p_first <- calculate_adjusted_risk("stroke", inc_stroke, matched_params, pd, living_idx)
    ys <- pd$years_since_stroke[living_idx]
    p_recur <- ifelse(ys <= 1, model_param[["stroke_recur_prob_yr1"]], model_param[["stroke_recur_prob_after"]])
    p <- te("stroke", ifelse(pd$stroke_ever[living_idx] == 1, p_recur, p_first))
    s_idx <- living_idx[U[living_idx, U_STROKE] < p]
    if (length(s_idx)) new_state_vals[s_idx, `:=`(stroke = 1L, stroke_ever = 1L, stroke_n = stroke_n + 1L,
                                                  years_since_stroke = 0L, cvd = 1L)]
  ####
  # 4. Non-stroke CVD events (first among cvd==0; recurrent after sim CVD event) --------
  #    base = GBD CVD incidence - stroke incidence (CVD aggregate includes stroke)
  ####
    inc_cvd_ns <- pmax(0, inc_of("cvd") - inc_stroke)
    p_first <- calculate_adjusted_risk("cvd", inc_cvd_ns, matched_params, pd, living_idx)
    yc <- pd$years_since_cvd_event[living_idx]
    p_recur <- ifelse(yc <= 1, model_param[["cvd_recur_prob_yr1"]], model_param[["cvd_recur_prob_after"]])
    first_ok <- pd$cvd[living_idx] == 0
    recur_ok <- pd$cvd_event_ever[living_idx] == 1
    p <- te("cvd", ifelse(first_ok, p_first, ifelse(recur_ok, p_recur, 0)))
    no_stroke_now <- new_state_vals$stroke[living_idx] == 0
    c_idx <- living_idx[no_stroke_now & U[living_idx, U_CVD] < p]
    if (length(c_idx)) new_state_vals[c_idx, `:=`(cvd = 1L, cvd_event = 1L, cvd_event_ever = 1L,
                                                  cvd_event_n = cvd_event_n + 1L, years_since_cvd_event = 0L)]
  ######'
  # 5. BMI: rank tracking + intervention effect --------
  #    bmi = Q_age(u) x (1 + bmi_effect); obesity indicators from BMI
  #    bmi_effect (proportion) is set by hook_bmi_effect(); stays 0 without an intervention
  ######'
    hook_bmi_effect(ictx, ystate, new_state_vals, living_idx)
    new_state_vals[living_idx, bmi := bmi_from_u(bmi_u, matched_params[living_idx]) * (1 + bmi_effect)]
    set_obesity_from_bmi(new_state_vals, living_idx)

  ######'
  # 6. ESKD among people with CKD ---------
  #    per-CKD probability uses baseline CKD prevalence (fixed so intervention effects are not re-normalized)
  ######'
    eskd_p <- pmin(1, model_param[["eskd_inc_prob_total_pop"]] / pmax(new_state_vals$prev0_ckd[living_idx], 1e-6))
    pos    <- which(new_state_vals$ckd[living_idx] == 1 & new_state_vals$eskd[living_idx] == 0)
    if (length(pos)) {
      p_e <- hook_modify_risk("eskd", eskd_p, ictx, ystate, pd, living_idx)
      idx <- living_idx[pos]
      new_state_vals[idx, eskd := as.integer(U[idx, U_ESKD] < pmin(1, p_e[pos]))]
    }
  }
  ####
  # 7. Incident case indicators for non-reversible conditions -------
  ####
  for (cond in c("t2d", "cvd", "ckd", "eskd")) {
    new_col <- paste0(cond, "_new")
    new_state_vals[, (new_col) := fifelse(pd[[cond]] == 0 & get(cond) == 1, 1L, get(new_col))]
  }
  new_state_vals[, stroke_this_year := stroke]

  ####
  # 8. Death -------
  ####
  inc_val_death   <- matched_params[["incidence_ac_death_val"]]
  inc_probs_death <- pmin(1, pmax(0, inc_val_death[living_idx]))
  ## mortality calibration from validation step if available: -----
  if (!is.null(cal) && all(c("alpha", "beta") %in% names(cal))) {
    idx <- which(tolower(trimws(cal$outcome)) == "ac_death")[1]
    a <- as.numeric(cal$alpha[idx]); b <- as.numeric(cal$beta[idx])
    p <- pmin(pmax(inc_probs_death, 1e-12), 1 - 1e-12)
    inc_probs_death <- pmin(pmax(1 - exp(-exp(a) * (-log1p(-p))^b), 1e-12), 1 - 1e-12)
  }
  ## risk adjustment for comorbidities incl. those accrued this year -----
  if (length(living_idx) > 0) {
    inc_probs_death <- modify_risk("ac_death", inc_probs_death, matched_params, new_state_vals, living_idx,
                                   inc_val_death, ictx, ystate, ext_te)
    new_state_vals$ac_death[living_idx] <- as.integer(U[living_idx, U_DEATH] < inc_probs_death)
  }

  ##########'
  # 9. DALYs (discounted, 3%/yr; first year undiscounted) ---------
  #  - YLL: discounted remaining life expectancy (discrete annuity) x discount factor of death year
  #  - YLD: multiplicative complement across conditions; half-cycle for deaths (time alive 0.5)
  ##########'
  disc <- 1 / ((1 + disc_rate)^(year - 1))
  newly_dead <- which(pd$ac_death == 0 & new_state_vals$ac_death == 1)
  if (length(newly_dead) > 0) {
    L <- matched_params[["le_yrs_val"]][newly_dead]
    if (any(is.na(L))) stop("Missing remaining life expectancy values for newly deceased")
    yll_disc <- (1 - (1 + disc_rate)^(-L)) / (1 - 1 / (1 + disc_rate))   # sum_{t=0}^{L-1} (1+r)^-t
    new_state_vals[newly_dead, `:=`(death_year = year, yll_adj = yll_disc * disc)]
  }
  time_alive <- numeric(nrow(new_state_vals))
  time_alive[living_idx] <- ifelse(new_state_vals$ac_death[living_idx] == 1L, 0.5, 1.0)

  dw_t2d_vec    <- assign_dw_t2d(new_state_vals, living_idx, dw_sensitivity)
  dw_cvd_vec    <- assign_dw_cvd(new_state_vals, living_idx, dw_sensitivity)
  dw_ckd_vec    <- assign_dw_ckd(new_state_vals, living_idx, dw_sensitivity)
  dw_stroke_vec <- assign_dw_stroke(new_state_vals, living_idx, dw_sensitivity)   # acute, event year only
  dw_combined_vec <- pmin(1, pmax(0, 1 - (1 - dw_t2d_vec) * (1 - dw_cvd_vec) * (1 - dw_ckd_vec) * (1 - dw_stroke_vec)))
  ta <- time_alive[living_idx] * disc
  new_state_vals[living_idx, `:=`(
    yld_adj_total     = yld_adj_total     + dw_combined_vec * ta,
    yld_t2d_adj_total = yld_t2d_adj_total + dw_t2d_vec * ta,
    yld_cvd_adj_total = yld_cvd_adj_total + dw_cvd_vec * ta,
    yld_ckd_adj_total = yld_ckd_adj_total + dw_ckd_vec * ta)]

  ##########'
  # 10. Costing (discounted 3%/yr) ---------
  # - chronic costs prorated by time alive; intervention costs from hook_costs()
  ##########'
  if (costing) {
    curr_iso3 <- unique(new_state_vals$iso3)
    hook_costs(ictx, ystate, new_state_vals, time_alive, disc)         # intervention costs (e.g. drug)
    # helper: add health-state cost to every intervention cost scenario total (cost_total_cols; none in the core)
    add_to_all_cost_totals <- function(idx, add_cost, state_cost_col = NULL) {
      if (length(idx) == 0) return(NULL)
      if (!is.null(state_cost_col)) {
        new_state_vals[[state_cost_col]][idx] <<- new_state_vals[[state_cost_col]][idx] + add_cost
      }
      new_state_vals$cost_conditions_total[idx] <<- new_state_vals$cost_conditions_total[idx] + add_cost
      for (tc in cost_total_cols) new_state_vals[[tc]][idx] <<- new_state_vals[[tc]][idx] + add_cost
    }
    if (length(living_idx) > 0) {
      health_costs <- health_event_cost_inputs |> filter(iso3 == curr_iso3)
      ### t2d ----
      c_t2d_lt10  <- health_costs |> pull(cost_id_t2d_lt10)
      c_t2d_geq10 <- health_costs |> pull(cost_id_t2d_geq10)
      idx_new  <- living_idx[new_state_vals$t2d_new[living_idx] == 1]
      idx_prev <- living_idx[new_state_vals$t2d_new[living_idx] == 0 & new_state_vals$t2d[living_idx] == 1]
      if (length(idx_new) > 0)  add_to_all_cost_totals(idx_new, c_t2d_lt10 * time_alive[idx_new] * disc, "cost_t2d_tot")
      if (length(idx_prev) > 0) add_to_all_cost_totals(idx_prev, ((c_t2d_lt10 + c_t2d_geq10) / 2) * time_alive[idx_prev] * disc, "cost_t2d_tot")
      ### stroke: acute on event; chronic post-stroke for sim stroke survivors ----
      c_stroke_acute <- health_costs |> pull(cost_id_stroke_acute)
      c_stroke_post  <- health_costs |> pull(cost_id_stroke_chronic)
      ev_stroke <- living_idx[new_state_vals$stroke[living_idx] == 1L]
      if (length(ev_stroke) > 0) add_to_all_cost_totals(ev_stroke, c_stroke_acute * disc, "cost_stroke_tot")
      post_stroke <- living_idx[pd$stroke_ever[living_idx] == 1 & new_state_vals$stroke[living_idx] == 0]
      if (length(post_stroke) > 0) add_to_all_cost_totals(post_stroke, c_stroke_post * time_alive[post_stroke] * disc, "cost_stroke_tot")
      ### cvd (no sim stroke) ----
      c_cvd <- health_costs |> pull(cost_id_cvd)
      idx_cvd <- living_idx[new_state_vals$cvd[living_idx] == 1 & new_state_vals$stroke_ever[living_idx] == 0]
      if (length(idx_cvd) > 0) add_to_all_cost_totals(idx_cvd, c_cvd * time_alive[idx_cvd] * disc, "cost_cvd_tot")
      ### ckd, eskd (mutually exclusive costs) ----
      c_ckd  <- health_costs |> pull(cost_id_ckd)
      c_eskd <- health_costs |> pull(cost_id_eskd)
      idx_eskd <- living_idx[new_state_vals$eskd[living_idx] == 1]
      idx_ckd  <- living_idx[new_state_vals$ckd[living_idx] == 1 & new_state_vals$eskd[living_idx] == 0]
      if (length(idx_ckd) > 0)  add_to_all_cost_totals(idx_ckd,  c_ckd  * time_alive[idx_ckd]  * disc, "cost_ckd_tot")
      if (length(idx_eskd) > 0) add_to_all_cost_totals(idx_eskd, c_eskd * time_alive[idx_eskd] * disc, "cost_eskd_tot")
      ### death ----
      c_death <- health_costs |> pull(cost_id_death)
      if (length(newly_dead) > 0) add_to_all_cost_totals(newly_dead, c_death * disc, "cost_death_tot")
    }
  }

  ##########'
  # 11. End of year: intervention state changes; time counters and aging among living -----------
  alive_end <- new_state_vals$ac_death == 0
  hook_end_of_year(ictx, new_state_vals)
  new_state_vals[alive_end & stroke_ever == 1, years_since_stroke := years_since_stroke + 1L]
  new_state_vals[alive_end & cvd_event_ever == 1, years_since_cvd_event := years_since_cvd_event + 1L]
  new_state_vals[alive_end, age := pmin(age + 1, 99)] # top coded at 99

  # if open cohort, introduce new youngest group 
  if (open_cohort && min(matched_params$age_low) <= min_age_model) {
    new_entrants <- fread(file.path(sim_temp_dir, "youngest_y0_person_data.csv")) |>
      mutate(pid = paste0(pid, "_entry", year))
    new_entrants <- hook_new_entrants(ictx, new_entrants, pd$crn_id[1], year)
    new_state_vals[, pid := as.character(pid)]
    new_state_vals <- rbindlist(list(new_state_vals, new_entrants), use.names = TRUE, fill = TRUE)
    rm(new_entrants)
  }
  new_state_vals
}
