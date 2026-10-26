# projects/semaglutide_global_access/R/glp1_hooks.R

# description ------
# GLP-1 (semaglutide) implementation of the core model's intervention hooks (GEM_Model_Intervention_Hooks.R).
# Sourced after the core (run_gem(..., project_modules = <this R folder>)); each function here replaces the
# no-op default of the same name. Contents:
#   treatment dynamics (initiation, discontinuation, restart), direct effects on event risks, BMI effect,
#   drug costing, scenario labels, run set-up, manifest fields.
# Related files: glp1_access_policy.R, glp1_treatment_effects.R, glp1_double_counting_calibration.R, glp1_inputs.R
# State columns (person level): treated, ever_treated, years_treated, spell_years, off_years, disc_this_year, disc_count,
#   restart_count, baseline_user, wl_resp (individual response), wl_effect (current BMI effect), on_frac,
#   treatment_eligible, has_access, ever_eligible, ever_access, glp1_tier, tier_at_start, u_access.

# random-number streams (core uses stream 1 for disease events) ----
UT_DISC <- 1; UT_RESTART <- 2; UT_INIT <- 3; UT_RESP <- 4     # columns of stream 2 (treatment dynamics)
# stream 4 = person-fixed access uniform (u_access)

# individual weight-loss response (proportion BMI change at full effect)
#  mean by group from glp1_impact_data; SD = |mean| x CV (STEP 1 distribution; input: wl_resp_cv)
draw_wl_resp <- function(u, t2d, glp1_impact_data) {
  get_m <- function(pop, ag = NULL) {
    x <- glp1_impact_data[glp1_impact_data$outcome_unit == "bmi (proportion change)" &
                            glp1_impact_data$treated_pop == pop, ]
    if (!is.null(ag)) x <- x[x$age_group == ag, ]
    as.numeric(x$outcome_val[1])
  }
  m <- ifelse(t2d == 1, get_m("t2d and obese"), get_m("obese without t2d", "adult"))
  r <- qnorm(pmin(pmax(u, 1e-6), 1 - 1e-6), mean = m, sd = abs(m) * model_param[["wl_resp_cv"]])
  pmin(pmax(r, -0.5), 0.1)
}

# GLP-1 direct effect multiplier for treated persons in one subgroup ----
#  years_treated < threshold (BMI pathway not yet acting on this outcome): trial HR
#  years_treated >= threshold: calibrated dc = HR_observed / HR_implied(BMI pathway),
#    so BMI pathway x dc reproduces the trial HR (glp1_double_counting_calibration.R)
#  no calibration row: trial HR in every year (warned once at load in glp1_inputs.R)
te_multiplier <- function(cond, grp, years_treated, dc_tbl) {
  hr <- glp1_impact_data$outcome_val[glp1_impact_data$outcome == cond & glp1_impact_data$treated_pop == grp][1]
  if (is.na(hr)) return(rep(1, length(years_treated)))
  dc <- dc_tbl[dc_tbl$outcome == cond & dc_tbl$subgroup == grp, ]
  if (nrow(dc) == 0) return(rep(hr, length(years_treated)))
  ifelse(years_treated >= dc$year[1], dc$multiplier[1], hr)
}

# apply direct effects to event probabilities of persons on treatment this year ----
#  pd: state defining subgroups; on_frac: 1 full year on drug, 0.5 in discontinuation year
#  -> multiplier 1 - on_frac x (1 - adj); years_treated: this year's count (incl. current year)
apply_te <- function(cond, inc_probs, pd, living_idx, on_frac, years_treated, dc_tbl, direct_te = TRUE) {
  if (!direct_te) return(inc_probs)
  f  <- on_frac[living_idx]
  tr <- which(f > 0)
  if (length(tr) == 0) return(inc_probs)
  t2d_t <- pd$t2d[living_idx[tr]]
  grp <- switch(cond,
    stroke = , ac_death = fifelse(t2d_t == 1 & pd$ckd[living_idx[tr]] == 1, "t2d and ckd",
                                  fifelse(t2d_t == 1, "t2d", "obese without t2d")),
    ckd = , cvd = fifelse(t2d_t == 1, "t2d", "obese without t2d"),
    eskd = fifelse(t2d_t == 1, "t2d and ckd", "obese without t2d"),
    t2d = fifelse(t2d_t == 0, "obese without t2d", NA_character_),
    return(inc_probs))
  for (g in unique(grp[!is.na(grp)])) {
    pos <- tr[which(grp == g)]
    adj <- te_multiplier(cond, g, years_treated[living_idx[pos]], dc_tbl)
    inc_probs[pos] <- inc_probs[pos] * (1 - f[pos] * (1 - adj))
  }
  inc_probs
}

# labels and configuration ------
# tags for settings that change the comparison arm too
hook_label_tags <- function(cs) {
  c(if (!identical(cs$baseline_use, "base")) paste0("bl-", cs$baseline_use),
    if (identical(cs$baseline_use_mode, "dynamics") && !identical(cs$baseline_use, "none")) "bldyn",
    if (!isTRUE(all.equal(cs$glp1$age_max_start, 85))) paste0("agemax", cs$glp1$age_max_start))
}

# intervention label: <years>yr_<access>_<uptake>[_t3ext][_rwe][_nocal][_bl-*][_bldyn][_agemax*][_c]
hook_scenario_label <- function(cs) {
  g <- cs$glp1
  core <- switch(g$strategy,
                 random   = paste0(round(100 * g$coverage_rate), "glp1"),
                 tiered   = paste0("tier", paste(sort(g$tiers_open), collapse = "")),
                 capacity = paste0("cap", round(100 * g$coverage_rate)))
  upt <- if (is.numeric(g$uptake)) paste0("upt", g$uptake) else paste0("upt-", g$uptake)
  paste(c(paste0(cs$sim_years, "yr_", core), upt,
          if (isTRUE(g$tier3_bmi27_cvd)) "t3ext", if (g$treatment_effect_source != "rct") g$treatment_effect_source,
          scenario_tags(cs)),
        collapse = "_")
}

hook_manifest_fields <- function(cs) {
  g <- cs$glp1
  list(treatment = cs$treatment_enabled, strategy = g$strategy, coverage_rate = g$coverage_rate,
       tiers_open = paste(g$tiers_open, collapse = ""), tier3_bmi27_cvd = g$tier3_bmi27_cvd,
       uptake = as.character(g$uptake), te_source = g$treatment_effect_source,
       baseline_use = cs$baseline_use, baseline_use_mode = cs$baseline_use_mode)
}

hook_print_config <- function(cs) {
  cat("Treatment enabled:", cs$treatment_enabled, "\n")
  cat("Baseline GLP-1 use:", cs$baseline_use, "(", cs$baseline_use_mode, ")\n")
  if (isTRUE(cs$treatment_enabled)) {
    g <- cs$glp1
    cat("Access strategy:", g$strategy, "| coverage:", g$coverage_rate, "| tiers open:", g$tiers_open,
        "| uptake:", g$uptake, "\n")
    cat("Treatment effect source:", g$treatment_effect_source, "| age limit for starting:", g$age_max_start, "\n")
  }
}

# run set-up ------
# policy object for the intervention arm
hook_make_policy <- function(cs, baseline_fp) if (isTRUE(cs$treatment_enabled)) make_access_policy(cs, baseline_fp) else NULL

# double-counting calibration (file written once per treatment-effect source)
hook_pre_run <- function(cs, baseline_fp, out_root) {
  if (!isTRUE(cs$treatment_enabled) || suppressWarnings(dc_check())) return(invisible(NULL))
  dc_calibration(baseline_fp, file.path(out_root, "dc_calibration"), cs$glp1$treatment_effect_source, cs)
  assign("dc_multiplier_ds", load_dc_multipliers(cs$glp1$treatment_effect_source), envir = .GlobalEnv)
  dc_check()
  invisible(NULL)
}

# once per arm: baseline-use scaling from the whole baseline sample; inputs captured for the workers
hook_prepare <- function(cs, files_yr0) {
  list(bl = if (identical(cs$baseline_use, "none")) NULL else baseline_use_scale(glp1_baseline_use, files_yr0),
       age_max_start = cs$glp1$age_max_start,
       bl_maintain = !identical(cs$baseline_use_mode, "dynamics"),
       dc_tbl = dc_multiplier_ds)
}

# per cohort file: policy for the country + settings
#  opts$direct_te = FALSE: BMI pathway only (double-counting calibration)
hook_country_context <- function(ctx0, policy, iso3, opts) {
  pol <- if (is.null(policy)) NULL else
    list(a = policy_access(policy, iso3), uptake = policy$uptake, tier3_bmi27_cvd = policy$tier3_bmi27_cvd,
         age_max = policy$age_max)
  c(ctx0, list(pol = pol, direct_te = if (is.null(opts$direct_te)) TRUE else isTRUE(opts$direct_te)))
}

# simulation ------
hook_init_columns <- function(dt, crn_id) {
  dt[, c("treatment_eligible", "has_access", "ever_eligible", "ever_access", "treated", "ever_treated",
         "years_treated", "spell_years", "off_years", "disc_this_year", "disc_count", "restart_count", "baseline_user") := 0L]
  dt[, c("glp1_tier", "tier_at_start") := NA_integer_]
  dt[, c("wl_resp", "wl_effect", "on_frac") := 0.0]
  dt[, u_access := crn_u(crn_id, 0L, .N, 1, stream = 4L)[, 1]]
  invisible(dt)
}

# year-0 state: prevalent users (same persons in every arm), tier and access for reporting
hook_init_baseline <- function(ictx, dt, crn_id, mp) {
  assign_baseline_use(dt, ictx$bl, crn_id, mp)
  pol <- ictx$pol
  if (!is.null(pol)) set_tier_access(dt, pol$a, pol$tier3_bmi27_cvd, pol$age_max)
  invisible(dt)
}

# transitions at start of year (state at end of last year):
#   restart: annual probability among those off treatment, age < age_max_start; first year off vs later years
#   initiation (policy arm only): eligible tier + access + annual uptake, among never treated
#   discontinuation: annual probability among those on treatment; first year of a spell vs later years
#     (T2D: Lim et al. Diabetologia 2025; no T2D: Rodriguez 2025 year 1, later years scaled by the Lim ratio);
#     baseline users (prevalent) use later-year rates; none for baseline users if bl_maintain
#   stop year counts as half a year on drug (cost, direct effect, BMI regain)
hook_start_of_year <- function(ictx, new_state_vals, pd, living_idx, year) {
  n_all <- nrow(pd)
  UT <- crn_u(pd$crn_id[1], year, n_all, 4, stream = 2L)
  set(new_state_vals, j = "disc_this_year", value = 0L)
  pol <- ictx$pol
  on_frac <- numeric(n_all)
  mp_ <- model_param
  if (length(living_idx) > 0 && (!is.null(pol) || any(pd$ever_treated[living_idx] == 1L))) {
    al  <- living_idx
    off <- al[pd$ever_treated[al] == 1L & pd$treated[al] == 0L & pd$age[al] < ictx$age_max_start]
    if (length(off)) {
      new_state_vals[off, off_years := off_years + 1L]                  # 1 = first year after the stop year
      f1 <- new_state_vals$off_years[off] == 1L
      p_rs <- ifelse(pd$t2d[off] == 1, ifelse(f1, mp_[["restart_prob_t2d"]], mp_[["restart_prob_t2d_later"]]),
                                       ifelse(f1, mp_[["restart_prob_not2d"]], mp_[["restart_prob_not2d_later"]]))
      rs <- off[UT[off, UT_RESTART] < p_rs]
      if (length(rs)) new_state_vals[rs, `:=`(treated = 1L, restart_count = restart_count + 1L, spell_years = 0L, off_years = 0L)]
    }
    if (!is.null(pol)) {
      set_tier_access(new_state_vals, pol$a, pol$tier3_bmi27_cvd, pol$age_max)
      cand <- al[new_state_vals$has_access[al] == 1L & pd$ever_treated[al] == 0L]
      if (length(cand)) {
        p_up <- uptake_prob(pd$t2d[cand], pd$female[cand], pd$age[cand], pol$uptake)
        st <- cand[UT[cand, UT_INIT] < p_up]
        if (length(st)) new_state_vals[st, `:=`(treated = 1L, ever_treated = 1L, tier_at_start = glp1_tier,
                                                wl_resp = draw_wl_resp(UT[st, UT_RESP], t2d, glp1_impact_data))]
      }
    }
    on <- al[new_state_vals$treated[al] == 1L]
    if (length(on)) {
      new_state_vals[on, spell_years := spell_years + 1L]                # 1 = first year of this treatment spell
      y1 <- new_state_vals$spell_years[on] == 1L & !(if (is.null(pd$baseline_user)) FALSE else pd$baseline_user[on] %in% 1L)
      p_disc <- ifelse(pd$t2d[on] == 1, ifelse(y1, mp_[["disc_prob_t2d"]], mp_[["disc_prob_t2d_later"]]),
                                        ifelse(y1, mp_[["disc_prob_not2d"]], mp_[["disc_prob_not2d_later"]]))
      if (ictx$bl_maintain) p_disc[pd$baseline_user[on] == 1L] <- 0
      stp <- on[UT[on, UT_DISC] < p_disc]
      on_frac[on] <- 1
      on_frac[stp] <- model_param[["disc_year_frac_on_tx"]]
      new_state_vals[on, years_treated := years_treated + 1L]
      if (length(stp)) new_state_vals[stp, `:=`(disc_this_year = 1L, disc_count = disc_count + 1L)]
    }
  }
  set(new_state_vals, j = "on_frac", value = on_frac)
  # direct effects: not for year-0 baseline users (input rates already reflect existing use);
  #   baseline users keep observed on-drug BMI and drug cost
  te_frac <- on_frac
  if (!is.null(pd$baseline_user)) te_frac[pd$baseline_user == 1L] <- 0
  list(on_frac = on_frac, te_frac = te_frac, yt = new_state_vals$years_treated)
}

# direct treatment effect on event probabilities
hook_modify_risk <- function(cond, p, ictx, ystate, pd, idx) {
  if (is.null(ystate)) return(p)
  apply_te(cond, p, pd, idx, ystate$te_frac, ystate$yt, ictx$dc_tbl, ictx$direct_te)
}

# BMI effect: full individual response while on drug; regain fraction per year off drug (STEP 1 extension);
#   stop year: regain over the half year off drug
hook_bmi_effect <- function(ictx, ystate, new_state_vals, living_idx) {
  if (is.null(ystate)) return(invisible(NULL))
  on_frac <- ystate$on_frac
  if (any(on_frac > 0) || any(new_state_vals$wl_effect[living_idx] != 0)) {
    rg  <- model_param[["regain_frac_per_year"]]
    f   <- on_frac[living_idx]
    on  <- living_idx[f == 1]
    stp <- living_idx[f > 0 & f < 1]
    off <- living_idx[f == 0 & new_state_vals$wl_effect[living_idx] != 0]
    if (length(on))  new_state_vals[on,  wl_effect := wl_resp]
    if (length(stp)) new_state_vals[stp, wl_effect := wl_resp * (1 - rg)^(1 - on_frac)]   # column (this year)
    if (length(off)) new_state_vals[off, wl_effect := wl_effect * (1 - rg)]
  }
  set(new_state_vals, j = "bmi_effect", value = new_state_vals$wl_effect)
  invisible(NULL)
}

# drug cost: monthly price x 12 x time alive x fraction of year on drug x discount; one column per price scenario
hook_costs <- function(ictx, ystate, new_state_vals, time_alive, disc) {
  if (is.null(ystate)) return(invisible(NULL))
  on_frac <- ystate$on_frac
  treated_abs <- which(on_frac > 0)
  if (length(treated_abs) == 0) return(invisible(NULL))
  glp1_row <- glp1_cost_inputs |> filter(iso3 == unique(new_state_vals$iso3))
  for (sfx in sub("^cost_total_", "", cost_total_cols)) {
    glp1_cost_val  <- (glp1_row |> pull(.data[[paste0("price_", sfx)]])) * 12   # annualize monthly cost
    c_treated_disc <- glp1_cost_val * time_alive[treated_abs] * on_frac[treated_abs] * disc
    treat_col <- paste0("cost_treat_total_", sfx); total_col <- paste0("cost_total_", sfx)
    set(new_state_vals, i = treated_abs, j = treat_col, value = new_state_vals[[treat_col]][treated_abs] + c_treated_disc)
    set(new_state_vals, i = treated_abs, j = total_col, value = new_state_vals[[total_col]][treated_abs] + c_treated_disc)
  }
  invisible(NULL)
}

# stops take effect at end of year
hook_end_of_year <- function(ictx, new_state_vals) {
  new_state_vals[disc_this_year == 1L, treated := 0L]
  invisible(NULL)
}

# open cohort: new access uniforms (CRN stream 4 of entry year)
hook_new_entrants <- function(ictx, entrants, crn_id, year) {
  entrants[, u_access := crn_u(crn_id, year, .N, 1, stream = 4L)[, 1]]
  entrants
}

# columns kept in saved simulation files
hook_output_cols <- function() {
  c("treated", "ever_treated", "treatment_eligible", "has_access", "ever_eligible", "ever_access",
    "on_frac", "years_treated", "disc_count", "restart_count", "glp1_tier", "tier_at_start", "baseline_user")
}

# cross-scenario tables and figures (after the per-scenario tables)
hook_postsim <- function(out_root, manifest, res_dir) {
  prim <- manifest[mode == "primary" & intervention_arm == TRUE]
  if (nrow(prim)) compile_scenarios(out_root, prim, res_dir)
  sq <- unique(manifest[mode == "primary" & !is.na(base_label), .(base_label, sim_years)])
  if (nrow(sq)) status_quo_summary(out_root, sq, res_dir)    # baseline use; calibration on/off comparison
  # circle plots by country (res_dir/figures); a plot error does not stop the run
  if (nrow(prim)) tryCatch(plot_gem_circles(out_root, res_dir, manifest = manifest),
                           error = function(e) message("plot_gem_circles failed: ", conditionMessage(e)))
  # tier expansion circle plots when scenarios named tier1to1 ... tier1to6 all ran (same comparison arm)
  if (all(paste0("tier1to", 1:6) %in% prim$config))
    tryCatch(expansion_scenarios_circleplot(out_root, res_dir, manifest = manifest),
             error = function(e) message("expansion_scenarios_circleplot failed: ", conditionMessage(e)))
  invisible(NULL)
}
