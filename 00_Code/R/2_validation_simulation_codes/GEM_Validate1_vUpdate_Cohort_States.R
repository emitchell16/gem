#' -------------------
# GEM validation - yearly update using the main model (update_cohort(), GEM_Model_Update_Cohort_States.R)
# Author: Liz Mitchell
# Validation runs the same state updates as the main model:
#   continuous BMI (rank tracking), joint RRs, stroke/CVD events with recurrence, excess mortality.
#   Validation-only additions: rates by WHO region, non-GLP-1 trial-arm HRs (ext_te), CVD death for MACE.
#' -------------------

# reference population constants ----
#  rr_scale_* convert population-average rates to reference-group rates; they must come from the
#  general population (not the study cohort, e.g. an all-T2D cohort would cancel the T2D RR).
#  weighted stratified reference sample (initialize_baseline) -> constants by region (or country), sex, 5-yr age
v_reference_consts <- function(cache_dir, n_min = 200L) {
  fp <- file.path(cache_dir, "v_reference_consts.rds")
  if (file.exists(fp)) return(readRDS(fp))
  ref_dir <- file.path(cache_dir, "reference_sample")
  initialize_baseline(cohort_input_data, paths, modifyList(config, list(sample_prop = 0, min_cohort_n = n_min)),
                      output_dir = ref_dir)
  ref <- as.data.table(open_dataset(ref_dir) |> collect())
  ref <- as.data.table(region_income_country_groups)[, .(iso3, who_region)][ref, on = "iso3"]
  ref[, `:=`(stroke_ever = 0L, age5 = 5L * (age %/% 5L))]
  k_reg <- ref[ac_death == 0, rr_baseline_consts(.SD, w), by = .(who_region, female, age5), .SDcols = rr_cols]
  k_iso <- ref[ac_death == 0, rr_baseline_consts(.SD, w), by = .(iso3, female, age5), .SDcols = rr_cols]
  out <- list(region = k_reg, iso3 = k_iso)
  saveRDS(out, fp)
  unlink(ref_dir, recursive = TRUE)
  out
}

# prepare a validation baseline cohort for update_cohort() ----
#  y0: v_initialize_baseline() output; inc_data: rates + prevalences by who_region-sex-age (bmi dist params added)
v_prepare_baseline <- function(y0, inc_data, consts, iso3_study, crn_base) {
  dt <- copy(as.data.table(y0))
  dt[, `:=`(iso3 = if (is.na(iso3_study)) "VAL" else iso3_study, age5 = 5L * (age %/% 5L))]
  # BMI percentile from the study's BMI within the population distribution (rank tracking keeps it)
  mp <- inc_data[dt, on = .(who_region, female, age_low = age), roll = TRUE, rollends = c(TRUE, TRUE)]
  dt[, bmi_u := bmi_to_u(bmi, mp)]
  init_sim_columns(dt, crn_base)
  dt[, crn_id := crn_base + match(who_region, sort(unique(who_region)))]   # separate streams per region subset
  # population constants: country-specific for single-country studies, else region
  k <- if (!is.na(iso3_study) && iso3_study %in% consts$iso3$iso3) consts$iso3[iso3 == iso3_study, !"iso3"] else consts$region
  by <- intersect(c("who_region", "female", "age5"), names(k))
  dt[k, (rr_const_names) := mget(paste0("i.", rr_const_names)), on = by]
  # ages outside reference strata: nearest stratum in same region/sex
  miss <- which(is.na(dt$rr_scale_t2d))
  if (length(miss)) {
    kk <- copy(k); setkeyv(kk, by)
    fill <- kk[dt[miss, ..by], on = by, roll = "nearest"]
    for (v in rr_const_names) set(dt, miss, v, fill[[v]])
  }
  dt[, age5 := NULL]
  dt
}

# yearly update by WHO region (rates differ by region; main update_cohort expects one rate set) ----
v_update_cohort <- function(pd, inc_data, year, ext_te = NULL, cal = acm_calibration) {
  rbindlist(lapply(split(pd, by = "who_region"), function(d) {
    cp <- inc_data[who_region == d$who_region[1]]
    d2 <- update_cohort(d, cp, year, cal = cal, ext_te = ext_te, costing = FALSE)
    # CVD death (MACE component): share of deaths from CVD in the population (GBD CSMR / all-cause), among newly dead
    nd <- which(d$ac_death == 0 & d2$ac_death == 1)
    if (length(nd)) {
      mp <- cp[d2[nd], on = .(female, age_low = age), roll = TRUE, rollends = c(TRUE, TRUE)]
      p  <- ifelse(mp$incidence_ac_death_val > 0, pmin(1, mp$death_cvd_val / mp$incidence_ac_death_val), 0)
      d2[nd, cvd_death := as.integer(runif(.N) < p)]
    }
    d2
  }), use.names = TRUE)
}
