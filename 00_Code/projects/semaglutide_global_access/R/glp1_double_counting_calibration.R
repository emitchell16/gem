# -------------------'
# R/1_function_scripts/GEM_Model_1_Double_Counting_Calibration.R
# Treatment-effect double-counting calibration
# paired-person (common random numbers) estimate among persons on treatment; ratio applied directly
# pairing is done per batch file and summed (memory)
# -------------------'
#
# METHODS NOTE
#  Trial HRs include the effect of weight loss. The model also lowers risk through BMI (joint RRs),
#  so applying the trial HR on top would count the BMI pathway twice.
#  1. two 2-year runs, same persons and random numbers: status quo vs all eligible treated from year 1,
#     direct effects off (opts = list(direct_te = FALSE)) -> treatment acts only through BMI; no baseline users
#  2. HR_implied = event probability treated / status quo, for the same persons (pid), among persons on
#     treatment the whole evaluation year, at risk at its start; subgroups as in apply_te()
#     evaluation year: 1 for death (BMI updated before death in year 1); 2 for non-fatal outcomes
#     (year-1 events use start-of-year BMI, before weight loss)
#  3. dc = HR_observed / HR_implied; model applies dc from the evaluation year on (trial HR before),
#     so BMI pathway x dc = trial HR
#  outcome definitions match where apply_te() acts: stroke and non-stroke CVD as yearly events (first + recurrent),
#  T2D/CKD as first events, ESKD among CKD, all-cause death

dc_outcome_def <- list(
  ac_death = list(y = "ac_death", at_risk = quote(TRUE)),
  stroke   = list(y = "stroke",   at_risk = quote(TRUE)),
  cvd      = list(y = "cvd_event", at_risk = quote(TRUE)),
  ckd      = list(y = "ckd",      at_risk = quote(ckd == 0)),
  t2d      = list(y = "t2d",      at_risk = quote(t2d == 0)),
  eskd     = list(y = "eskd",     at_risk = quote(ckd == 1 & eskd == 0)))

# subgroup labels (same rules as apply_te(); ESKD as in update_cohort section 6)
dc_subgroup <- function(cond, t2d, ckd) {
  switch(cond,
         ac_death = , stroke = fifelse(t2d == 1 & ckd == 1, "t2d and ckd", fifelse(t2d == 1, "t2d", "obese without t2d")),
         ckd = , cvd = fifelse(t2d == 1, "t2d", "obese without t2d"),
         t2d  = fifelse(t2d == 0, "obese without t2d", NA_character_),
         eskd = fifelse(t2d == 1, "t2d and ckd", "obese without t2d"))
}

dc_calibration <- function(baseline_fp, dc_dir, te_source = config_sims$glp1$treatment_effect_source,
                           cs = config_sims, min_events = 20) {
  # 1. two runs (direct effects off) ----
  cs_dc <- modifyList(cs, list(sim_years = 2L, treatment_enabled = TRUE, baseline_use = "none",   # untreated comparator
                               glp1 = modifyList(cs$glp1, list(strategy = "random", coverage_rate = 1, uptake = "full"))))
  pol <- make_access_policy(cs_dc, baseline_fp)
  main_sim_file_management(baseline_fp, dc_dir, 2L, policy = NULL, label = "dc_base", cs = cs_dc,
                           keep_years = 0:2, opts = list(direct_te = FALSE))
  main_sim_file_management(baseline_fp, dc_dir, 2L, policy = pol,  label = "dc_treat", cs = cs_dc,
                           keep_years = 0:2, opts = list(direct_te = FALSE))

  cols <- c("pid", "w", "t2d", "ckd", "eskd", "ac_death", "stroke", "cvd_event", "on_frac", "years_treated", "disc_count")
  hr_obs <- as.data.table(glp1_impact_data)[!grepl("bmi", outcome_unit, ignore.case = TRUE),
                                            .(outcome, subgroup = treated_pop, HR_observed = outcome_val)]

  # 2. paired hazards by outcome and subgroup, one batch file at a time (memory) ----
  #' both arms use the same batch files (same baseline input) -> pairs by pid within batch
  one_batch <- function(bf) {
    rd <- function(arm, yr) as.data.table(read_parquet(file.path(dc_dir, paste0("sim=", arm), paste0("year=", yr), bf),
                                                       col_select = all_of(cols)))
    B <- lapply(0:2, function(y) rd("dc_base", y));  names(B) <- 0:2
    Tt <- lapply(0:2, function(y) rd("dc_treat", y)); names(Tt) <- 0:2
    rbindlist(lapply(names(dc_outcome_def), function(y) {
      def <- dc_outcome_def[[y]]; e <- if (y == "ac_death") 1L else 2L; r <- as.character(e - 1L)
      # treated arm: on drug continuously from year 1 through the evaluation year, alive and at risk at its start
      tr <- Tt[[r]][ac_death == 0 & eval(def$at_risk), .(pid, w, grp = dc_subgroup(y, t2d, ckd))]
      on <- Tt[[as.character(e)]][on_frac == 1 & years_treated == e & disc_count == 0, .(pid, y_t = get(def$y))]
      tr <- tr[on, on = "pid", nomatch = 0L]
      # same persons in status quo arm, at risk at start of evaluation year
      bs <- B[[r]][ac_death == 0 & eval(def$at_risk), .(pid)][B[[as.character(e)]][, .(pid, y_b = get(def$y))], on = "pid", nomatch = 0L]
      d  <- tr[bs, on = "pid", nomatch = 0L][!is.na(grp)]
      d[, .(year = e, n = .N, events_base = sum(y_b), sw = sum(w), swyb = sum(w * y_b), swyt = sum(w * y_t)),
        by = .(subgroup = grp)][, outcome := y][]
    }))
  }
  bfs <- list.files(file.path(dc_dir, "sim=dc_treat", "year=2"), pattern = "\\.parquet$")
  res <- rbindlist(lapply(bfs, function(bf) { x <- one_batch(bf); gc(verbose = FALSE); x }))
  res <- res[, .(n = sum(n), events_base = sum(events_base), hazard_base = sum(swyb) / sum(sw),
                 hazard_treat = sum(swyt) / sum(sw)), by = .(outcome, subgroup, year)]
  res <- hr_obs[res, on = .(outcome, subgroup)]
  res[, `:=`(HR_implied = hazard_treat / hazard_base)]
  res[, multiplier := HR_observed / HR_implied]
  ok <- res[is.finite(multiplier) & events_base >= min_events]
  if (nrow(ok) < nrow(res)) warning("dc calibration: ", nrow(res) - nrow(ok),
                                    " outcome-subgroup rows dropped (too few events); trial HR used for those")
  if (any(ok$multiplier > 1)) message("dc > 1 for: ", paste(ok[multiplier > 1, paste(outcome, subgroup)], collapse = "; "),
                                      " (BMI pathway alone exceeds trial effect; dc restores trial HR)")
  out <- ok[, .(outcome, subgroup, year, n, hazard_base, hazard_treat, HR_implied, HR_observed, multiplier)]
  fwrite(out, file.path(paths$working, paste0("double_counting_calibration_results_", te_source, ".csv")))
  unlink(file.path(dc_dir, c("sim=dc_base", "sim=dc_treat")), recursive = TRUE)
  out[]
}
