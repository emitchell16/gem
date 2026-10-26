#' Apply risk equations with calculate_adjusted_risk()
#' joint RRs (multiplicative, continuous BMI) + baseline rescaling,
#'   stroke as yearly event, cause-specific excess mortality
#' diabetes-duration HRs are not used (baseline duration unobserved)
#' requires GEM_Model_BMI_RR.R and tables nejm_age_rr, bmi_rr_tbl, t2d_rr, cvd_rr_stroke_bysex, model_param
#' -------------------------
#'
#' METHODS NOTE
#'  incidence (t2d, cvd [non-stroke], ckd, first stroke):
#'    p_i = r_pop(age,sex,iso3) * rr_joint_i * s_cond
#'    rr_joint_i = product of individual RRs (BMI per-5 x T2D x prior CVD for stroke)
#'    s_cond     = N_alive / sum(rr_joint over susceptible), fixed at baseline per cohort
#'                 -> converts population-average rate to reference-group rate (r0 = r_pop / sum p_k RR_k)
#'                    and GBD per-total-population rates to rates among those at risk
#'  recurrent stroke / CVD events: absolute annual probabilities by time since sim event (update_cohort)
#'  all-cause death (cause-deleted, additive excess):
#'    h_i = r_other * rr_bmi_other_i * s_death + sum_k I_k * EMR_k,  k in {cvd, ckd, t2d}
#'    r_other = r_ac - sum CSMR_k;  EMR_k = CSMR_k / prevalence_k (GBD)
#'    CKD split: ESKD excess = param eskd_excess_mort; non-ESKD CKD gets remainder of CKD deaths
#'  baseline constants fixed at y0 so treatment effects are not re-normalized away

rr_cols <- c("age", "bmi", "female", "t2d", "cvd", "ckd", "eskd", "stroke_ever")

# individual joint RR ----
#  d: list of person vectors (rr_cols)
rr_joint_vec <- function(cond, d) {
  n <- length(d$age)
  rr_t2d <- function(nm) {
    v <- t2d_rr$value[t2d_rr$condition == nm][1]
    if (is.na(v)) v <- 1
    ifelse(d$t2d == 1, v, 1)
  }
  if (cond == "t2d") {
    bmi_rr_vec("t2d", d$bmi, d$age, nejm_age_rr, bmi_rr_tbl)
  } else if (cond == "cvd") {
    bmi_rr_vec("cvd", d$bmi, d$age, nejm_age_rr, bmi_rr_tbl) * rr_t2d("cvd")
  } else if (cond == "ckd") {
    bmi_rr_vec("ckd", d$bmi, d$age, nejm_age_rr, bmi_rr_tbl) * rr_t2d("ckd")
  } else if (cond == "stroke") {
    rr_cvd <- cvd_rr_stroke_bysex$value[match(as.integer(d$female), as.integer(cvd_rr_stroke_bysex$female))]
    bmi_rr_vec("stroke", d$bmi, d$age, nejm_age_rr, bmi_rr_tbl) * rr_t2d("stroke") *
      ifelse(d$cvd == 1, rr_cvd, 1)
  } else if (cond == "ac_death") {
    bmi_rr_vec("death_other_nonmodeled", d$bmi, d$age, nejm_age_rr, bmi_rr_tbl)
  } else rep(1, n)
}

# susceptible set for each first-event outcome
at_risk <- function(cond, d) {
  switch(cond,
         t2d    = d$t2d == 0,
         ckd    = d$ckd == 0,
         stroke = d$stroke_ever == 0,
         cvd    = d$cvd == 0,
         rep(TRUE, length(d$age)))
}

# baseline constants for one group (living persons); w: survey weights (NULL = equal) ----
#  population-representative input required (validation uses a reference population sample)
rr_baseline_consts <- function(d, w = NULL) {
  if (is.null(w)) w <- rep(1, length(d$age))
  wm <- function(x, ww = w) sum(ww * x) / sum(ww)
  sc <- function(cond) {
    rr <- rr_joint_vec(cond, d); sus <- at_risk(cond, d)
    s  <- sum(w) / sum(w[sus] * rr[sus])
    if (!is.finite(s) || s <= 0) 1 else s
  }
  ck <- d$ckd == 1
  list(rr_scale_t2d      = sc("t2d"),
       rr_scale_cvd      = sc("cvd"),
       rr_scale_ckd      = sc("ckd"),
       rr_scale_stroke   = sc("stroke"),
       rr_scale_ac_death = 1 / wm(rr_joint_vec("ac_death", d)),
       prev0_cvd = wm(d$cvd == 1), prev0_ckd = wm(ck), prev0_t2d = wm(d$t2d == 1),
       eskd_share0 = if (any(ck)) wm(d$eskd[ck] == 1, w[ck]) else 0)
}
rr_const_names <- c("rr_scale_t2d", "rr_scale_cvd", "rr_scale_ckd", "rr_scale_stroke", "rr_scale_ac_death",
                    "prev0_cvd", "prev0_ckd", "prev0_t2d", "eskd_share0")

# add baseline constants as columns; call once at y0 on a population-representative cohort file
add_rr_baseline_cols <- function(dt) {
  if (!"stroke_ever" %in% names(dt)) dt[, stroke_ever := 0L]
  alive <- dt$ac_death == 0
  k <- rr_baseline_consts(as.list(dt[alive, ..rr_cols]))
  dt[, names(k) := k]
  invisible(dt)
}

# person vectors for living_idx (missing optional columns -> 0)
get_d <- function(pd, idx) {
  setNames(lapply(rr_cols, function(v) if (is.null(pd[[v]])) rep(0, length(idx)) else pd[[v]][idx]), rr_cols)
}
# baseline constant for idx (fallback: compute from current data, no-treatment runs only)
get_const <- function(pd, idx, nm, d) {
  v <- pd[[nm]]
  if (is.null(v)) return(rr_baseline_consts(d)[[nm]])
  v[idx]
}

# main risk equation function ------
#  inc_probs: population rate for cond among living_idx (cvd: non-stroke rate)
calculate_adjusted_risk <- function(cond, inc_probs, matched_params, pd, living_idx,
                                    inc_val = NULL, new_state_vals = NULL) {
  d <- get_d(pd, living_idx)

  # first-event incidence ----
  if (cond %in% c("t2d", "cvd", "ckd", "stroke")) {
    rr <- rr_joint_vec(cond, d)
    s  <- get_const(pd, living_idx, paste0("rr_scale_", cond), d)
    inc_probs <- pmin(1, inc_probs * s * rr)

  # all-cause death ----
  } else if (cond == "ac_death") {
    ac_raw <- inc_val[living_idx]
    cal    <- ifelse(ac_raw > 0, inc_probs / ac_raw, 1)    # carries mortality calibration if applied
    other  <- ac_raw
    excess <- numeric(length(living_idx))
    emr    <- list()
    for (k in c("cvd", "ckd", "t2d")) {
      csmr <- matched_params[[paste0("death_", k, "_val")]][living_idx]
      csmr[is.na(csmr)] <- 0
      prev <- matched_params[[paste0("prevalence_", k, "_val")]]
      prev <- if (is.null(prev)) get_const(pd, living_idx, paste0("prev0_", k), d) else prev[living_idx]
      ok   <- !is.na(prev) & prev > 0
      other    <- other - ifelse(ok, csmr, 0)               # no prevalence -> deaths stay in other-cause
      emr[[k]] <- ifelse(ok, csmr / pmax(prev, 1e-12), 0)
    }
    # CKD: ESKD fixed excess; non-ESKD CKD gets remaining CKD excess
    s_eskd <- get_const(pd, living_idx, "eskd_share0", d)
    emr_eskd <- model_param[["eskd_excess_mort"]]
    emr_ckd_non <- pmax(0, (emr$ckd - s_eskd * emr_eskd) / pmax(1 - s_eskd, 1e-6))
    excess <- ifelse(d$cvd == 1, emr$cvd, 0) +
              ifelse(d$t2d == 1, emr$t2d, 0) +
              ifelse(d$ckd == 1 & d$eskd == 1, emr_eskd, ifelse(d$ckd == 1, emr_ckd_non, 0))
    other <- pmax(0, other)
    s_d   <- get_const(pd, living_idx, "rr_scale_ac_death", d)
    inc_probs <- pmin(1, cal * (other * rr_joint_vec("ac_death", d) * s_d + excess))
  }
  inc_probs
}
