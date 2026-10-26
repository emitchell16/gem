# R/1_function_scripts/GEM_Model_PSA_Functions.R

# description ------
# Probabilistic sensitivity analysis (PSA) of an intervention arm vs the comparison arm.
# Intervention inputs enter through hook_psa_inputs() (capture) and hook_psa_draw() (draw) in GEM_Model_Intervention_Hooks.R.
#' METHODS NOTE
#' METHODS NOTE
#'  - parameter families: one standard-normal z per family per run, shared by every row in the family;
#'    each row's value = its own uncertainty interval evaluated at that z (Gaussian copula, comonotone within family)
#'    -> uncertainty is not averaged away across country-age-sex strata (Briggs et al. 2012, ISPOR-SMDM TF-6)
#'  - distributions: HR/RR log-normal; rates/probabilities logit-normal (GBD) or beta; DWs beta; costs gamma;
#'    baseline prevalences fixed (population structure); intervention inputs drawn in hook_psa_draw()
#'  - first-order (Monte Carlo) noise: many runs x small samples, independent random numbers per run;
#'    within-run variance from B person batches; parameter variance = between-run var - mean within-run var
#'    (ANOVA; O'Hagan et al. 2007); intervals from run estimates shrunk toward the mean by sqrt(var_param/var_total)

# distribution helpers (vectorized; z shared within family) ----
.z95 <- qnorm(0.975)
psa_se <- function(lower, upper, scale = c("identity", "log", "logit")) {
  scale <- match.arg(scale)
  f <- switch(scale, identity = identity, log = log, logit = qlogis)
  (f(upper) - f(lower)) / (2 * .z95)
}
ok_ci <- function(v, l, u) !is.na(v) & !is.na(l) & !is.na(u) & u > l

# log-normal (HR, RR)
draw_lognormal <- function(v, l, u, z) {
  ok <- ok_ci(v, l, u) & v > 0 & l > 0
  out <- v; out[ok] <- exp(log(v[ok]) + z * psa_se(l[ok], u[ok], "log")); out
}
# logit-normal (small rates/probabilities with asymmetric UIs, e.g. GBD)
draw_logitnormal <- function(v, l, u, z, eps = 1e-9) {
  ok <- ok_ci(v, l, u) & v > 0 & v < 1
  cl <- function(x) pmin(pmax(x, eps), 1 - eps)
  out <- v; out[ok] <- plogis(qlogis(cl(v[ok])) + z * psa_se(cl(l[ok]), cl(u[ok]), "logit")); out
}
# beta (probabilities, disability weights); method of moments from mean and CI-derived SE; quantile at pnorm(z)
draw_beta <- function(v, l, u, z) {
  ok <- ok_ci(v, l, u) & v > 0 & v < 1
  out <- v
  se <- psa_se(l[ok], u[ok]); m <- v[ok]
  k  <- m * (1 - m) / se^2 - 1
  good <- k > 0
  out[ok][good]  <- qbeta(pnorm(z), m[good] * k[good], (1 - m[good]) * k[good])
  out[ok][!good] <- draw_logitnormal(m[!good], l[ok][!good], u[ok][!good], z)   # CI too wide for beta
  out
}
# gamma (costs); SE = cv x mean when no CI
draw_gamma <- function(v, z, cv) {
  ok <- !is.na(v) & v > 0
  out <- v; shape <- 1 / cv^2
  out[ok] <- qgamma(pnorm(z), shape = shape, scale = v[ok] / shape); out
}
# normal (BMI change, life expectancy)
draw_normal <- function(v, l, u, z) {
  ok <- ok_ci(v, l, u)
  out <- v; out[ok] <- v[ok] + z * psa_se(l[ok], u[ok]); out
}

# capture base (point-estimate) inputs once, after set_in_tables.R ----
psa_capture_base <- function() {
  g <- function(x) if (exists(x, envir = .GlobalEnv)) get(x, envir = .GlobalEnv) else NULL
  list(all_cohort_param_rates = copy(as.data.table(g("all_cohort_param_rates"))),
       nejm_age_rr = g("nejm_age_rr"), bmi_rr_tbl = g("bmi_rr_tbl"),
       t2d_rr = g("t2d_rr"), cvd_rr_stroke_bysex = g("cvd_rr_stroke_bysex"),
       dw_tbl = copy(as.data.table(g("dw_tbl"))), ext = hook_psa_inputs(),
       model_param = g("model_param"), model_param_ci = g("model_param_ci"),
       health_event_cost_inputs = g("health_event_cost_inputs"))
}

# one PSA draw: returns list of drawn inputs + table of family z values ----
draw_psa_params <- function(base, cost_cv = 0.25) {
  zlog <- list()
  zf <- function(fam) { z <- rnorm(1); zlog[[fam]] <<- z; z }
  out <- base

  ## 1. GBD rates and life expectancy: one z per measure (shared over all countries, ages, sexes) ----
  cp <- copy(base$all_cohort_param_rates)
  stems <- sub("_val$", "", grep("_val$", names(cp), value = TRUE))
  stems <- stems[paste0(stems, "_lower") %in% names(cp) & paste0(stems, "_upper") %in% names(cp)]
  stems <- stems[grepl("^(incidence_(t2d|cvd|ckd|stroke|ac_death)|death_|le_yrs)", stems)]  # prevalence fixed
  for (s in stems) {
    z <- zf(s); v <- cp[[paste0(s, "_val")]]; l <- cp[[paste0(s, "_lower")]]; u <- cp[[paste0(s, "_upper")]]
    set(cp, j = paste0(s, "_val"),
        value = if (s == "le_yrs") draw_normal(v, l, u, z) else draw_logitnormal(v, l, u, z))
  }
  out$all_cohort_param_rates <- cp

  ## 2. BMI relative risks: one z per outcome, shared across age bands (log-normal) ----
  ar <- base$nejm_age_rr
  for (k in unique(ar$key)) {
    i <- ar$key == k; ar$value[i] <- draw_lognormal(ar$value[i], ar$lower[i], ar$upper[i], zf(paste0("bmi_rr_", k)))
  }
  out$nejm_age_rr <- ar
  br <- base$bmi_rr_tbl
  for (k in unique(br$outcome)) {
    i <- br$outcome == k; br$value[i] <- draw_lognormal(br$value[i], br$lower[i], br$upper[i], zf(paste0("bmi_rr_", k)))
  }
  out$bmi_rr_tbl <- br

  ## 3. other RRs: T2D RR per outcome; prior-CVD stroke RR shared across sexes (log-normal) ----
  tr <- as.data.frame(base$t2d_rr)
  for (k in unique(tr$condition)) {
    i <- tr$condition == k
    tr$value[i] <- draw_lognormal(as.numeric(tr$value[i]), as.numeric(tr$lower[i]), as.numeric(tr$upper[i]), zf(paste0("t2d_rr_", k)))
  }
  out$t2d_rr <- tr
  cs <- as.data.frame(base$cvd_rr_stroke_bysex)
  cs$value <- draw_lognormal(as.numeric(cs$value), as.numeric(cs$lower), as.numeric(cs$upper), zf("cvd_rr_stroke"))
  out$cvd_rr_stroke_bysex <- cs

  ## 4. intervention inputs (project module): out$ext holds the drawn objects; zf(family) gives a shared z ----
  out <- hook_psa_draw(base, out, zf)

  ## 5. disability weights: one z per condition (beta) ----
  dw <- copy(base$dw_tbl)
  dw[, fam__ := tolower(trimws(condition_logic))]
  for (k in unique(dw$fam__)) {
    i <- which(dw$fam__ == k); z <- zf(paste0("dw_", k))
    dw[i, `:=`(mean = draw_beta(mean, lower, upper, z))]
  }
  out$dw_tbl <- dw[, fam__ := NULL]

  ## 6. scalar parameters with intervals (model_param_ci) ----
  mp <- base$model_param; ci <- base$model_param_ci
  if (!is.null(ci) && nrow(ci)) {
    for (fam in unique(ci$psa_family)) {
      z <- zf(fam); r <- ci[ci$psa_family == fam, ]
      for (j in seq_len(nrow(r))) {
        p <- r$parameter[j]
        mp[[p]] <- switch(r$psa_dist[j],
          beta      = draw_beta(mp[[p]], r$lower[j], r$upper[j], z),
          lognormal = draw_lognormal(mp[[p]], r$lower[j], r$upper[j], z),
          normal    = draw_normal(mp[[p]], r$lower[j], r$upper[j], z),
          mp[[p]])
      }
    }
  }
  out$model_param <- mp

  ## 7. health-state costs: one z per cost category, shared across countries (gamma; cv assumed) ----
  if (!is.null(base$health_event_cost_inputs)) {
    hc <- as.data.table(copy(base$health_event_cost_inputs))
    for (cc in grep("^cost_id_", names(hc), value = TRUE)) set(hc, j = cc, value = draw_gamma(hc[[cc]], zf(cc), cost_cv))
    out$health_event_cost_inputs <- hc
  }
  out$z <- data.table(family = names(zlog), z = unlist(zlog))
  out
}

# assign drawn (or base) inputs to the global environment used by the model
psa_apply <- function(p) {
  for (nm in c("all_cohort_param_rates", "nejm_age_rr", "bmi_rr_tbl", "t2d_rr", "cvd_rr_stroke_bysex",
               "dw_tbl", "model_param", "health_event_cost_inputs")) {
    if (!is.null(p[[nm]])) assign(nm, p[[nm]], envir = .GlobalEnv)
  }
  for (nm in names(p$ext)) assign(nm, p$ext[[nm]], envir = .GlobalEnv)      # intervention inputs
  invisible(NULL)
}

# run summary: weighted outcomes by region, full sample and by person batch ----
#  batch = pid mod B; batch totals scaled by B so batch and full-sample estimates share units
psa_summarize_run <- function(base_dir, int_dir, sim_years, B = 10, outcome_conditions = all_outcomes_list) {
  cg <- as.data.table(read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")))[, .(iso3, who_region)]
  rd <- function(p) {
    d <- as.data.table(open_dataset(normalizePath(p), format = "parquet") |> collect())
    d[, batch := as.integer(as.numeric(gsub("[^0-9]", "", as.character(pid))) %% B) + 1L]
    cg[d, on = "iso3"]
  }
  summ <- function(d) {
    one <- function(x, scale_w = 1) {
      x <- copy(x)[, w := w * scale_w]
      as.data.table(summarize_weighted(x, outcome_conditions))
    }
    full <- rbind(one(d)[, who_region := "global"],
                  rbindlist(lapply(split(d, d$who_region), function(x) one(x)[, who_region := x$who_region[1]])))
    full[, batch := 0L]
    bt <- rbindlist(lapply(split(d, d$batch), function(x) {
      b <- x$batch[1]
      rbind(one(x, B)[, who_region := "global"],
            rbindlist(lapply(split(x, x$who_region), function(y) one(y, B)[, who_region := y$who_region[1]])))[, batch := b]
    }))
    rbind(full, bt, fill = TRUE)
  }
  long <- function(s, arm) {
    melt(s, id.vars = c("who_region", "batch"), measure.vars = intersect(c(outcome_conditions, "daly_total_adj"), names(s)),
         variable.name = "condition", value.name = arm)
  }
  b <- summ(rd(file.path(base_dir, paste0("year=", sim_years))))
  i <- summ(rd(file.path(int_dir,  paste0("year=", sim_years))))
  if (all(c("yll_adj", "yld_adj_total") %in% names(b))) {
    b[, daly_total_adj := yll_adj + yld_adj_total]; i[, daly_total_adj := yll_adj + yld_adj_total]
  }
  out <- long(b, "base")[long(i, "int"), on = c("who_region", "batch", "condition")]
  out[, diff := int - base]
  out[]
}

# main PSA loop ----
#  R runs; each run: new parameter draw, new small weighted baseline sample, independent random numbers,
#  both arms simulated, summarized by region and batch, outputs deleted (no reuse across runs)
#  intervention policy rebuilt per run from that run's baseline (hook_make_policy())
#' psa_out_folder: checkpoint (local temp); final_dir: final csv files (default <psa_out_folder>/final; run_gem: 03_Out/results_<tag>/psa)
#' intervention prices are not sampled (scenario variable); ICER intervals: see project module
run_psa <- function(cs, label, psa_temp_folder, psa_out_folder, R = cs$psa_runs, B = cs$psa_batches,
                    final_dir = file.path(psa_out_folder, "final")) {
  dir.create(psa_out_folder, recursive = TRUE, showWarnings = FALSE)
  checkpoint <- file.path(psa_out_folder, "psa_runs_partial.rds")
  res <- if (file.exists(checkpoint)) readRDS(checkpoint) else data.table(run = integer())
  zs  <- if (file.exists(sub("partial", "z", checkpoint))) readRDS(sub("partial", "z", checkpoint)) else data.table(run = integer())
  base_inputs <- psa_capture_base()
  on.exit({psa_apply(base_inputs); assign("crn_stream_offset", 0L, envir = .GlobalEnv)}, add = TRUE)
  psa_cfg <- modifyList(config, list(sample_prop = cs$psa_sample_prop, min_cohort_n = cs$psa_min_cohort_n))
  base_label <- base_label_of(cs)

  for (r in setdiff(seq_len(R), unique(res$run))) {
    t0 <- Sys.time()
    run_dir <- file.path(psa_temp_folder, paste0("run_", r))
    unlink(run_dir, recursive = TRUE); dir.create(run_dir, recursive = TRUE)
    # parameter draw (reproducible per run)
    set.seed(cs$seed * 100000 + r)
    p <- draw_psa_params(base_inputs, cost_cv = cs$psa_cost_cv)
    psa_apply(p)
    # independent first-order noise per run: new baseline sample + new event random numbers
    assign("crn_stream_offset", r, envir = .GlobalEnv)
    bfp <- initialize_baseline(cohort_input_data, paths, modifyList(psa_cfg, list(seed = cs$seed * 1000 + r)),
                               output_dir = file.path(run_dir, "baseline"))
    main_sim_file_management(bfp, run_dir, cs$sim_years, NULL, base_label, cs)
    main_sim_file_management(bfp, run_dir, cs$sim_years, hook_make_policy(cs, bfp), label, cs)
    s <- psa_summarize_run(file.path(run_dir, paste0("sim=", base_label)),
                           file.path(run_dir, paste0("sim=", label)), cs$sim_years, B)
    res <- rbind(res, s[, run := r], fill = TRUE)
    zs  <- rbind(zs, p$z[, run := r], fill = TRUE)
    saveRDS(res, checkpoint); saveRDS(zs, sub("partial", "z", checkpoint))
    unlink(run_dir, recursive = TRUE)
    cat("PSA run", r, "of", R, "-", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
  }
  psa_intervals <- psa_summarize(res, B)
  conv <- psa_convergence(res)
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  ff <- function(nm) file.path(final_dir, paste0(nm, "_", label, ".csv"))
  fwrite(psa_intervals, ff("psa_intervals_by_region_condition"))
  fwrite(res[batch == 0], ff("psa_results_all_runs"))
  fwrite(conv, ff("psa_convergence"))
  fwrite(zs, ff("psa_family_z_draws"))
  message("PSA results: ", final_dir)
  invisible(psa_intervals)
}

# ANOVA-corrected summaries (O'Hagan et al. 2007) ----
#  est_r = full-sample estimate per run; within_r = var(batch estimates) / B
#  var_param = max(0, var(est_r) - mean(within_r)); intervals from mean + k (est_r - mean), k = sqrt(var_param / var(est_r))
psa_summarize <- function(res, B) {
  full <- res[batch == 0]
  wv <- res[batch > 0, .(within_diff = var(diff) / B, within_int = var(int) / B, within_base = var(base) / B),
            by = .(run, who_region, condition)]
  d <- wv[full, on = c("run", "who_region", "condition")]
  q <- function(x) as.list(quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
  d[, {
    adj <- function(est, within) {
      vt <- var(est); vp <- max(0, vt - mean(within, na.rm = TRUE)); k <- if (vt > 0) sqrt(vp / vt) else 0
      m <- mean(est); list(mean = m, lo = quantile(m + k * (est - m), 0.025, names = FALSE),
                           hi = quantile(m + k * (est - m), 0.975, names = FALSE),
                           mc_share = if (vt > 0) 1 - vp / vt else NA_real_)
    }
    a_d <- adj(diff, within_diff); a_b <- adj(base, within_base); a_i <- adj(int, within_int)
    list(n_runs = .N,
         base_mean = a_b$mean, base_low = a_b$lo, base_high = a_b$hi,
         int_mean = a_i$mean,  int_low = a_i$lo,  int_high = a_i$hi,
         diff_mean = a_d$mean, diff_low = a_d$lo, diff_high = a_d$hi,
         diff_low_unadj = q(diff)[[1]], diff_high_unadj = q(diff)[[2]],
         mc_var_share_diff = a_d$mc_share)
  }, by = .(who_region, condition)][order(who_region != "global", who_region, condition)]
}

# running mean and percentiles of the between-arm difference by number of runs (convergence check)
psa_convergence <- function(res) {
  full <- res[batch == 0][order(run)]
  full[, {
    n <- seq_along(diff)
    list(n_runs = n, running_mean = cumsum(diff) / n,
         running_p025 = sapply(n, function(k) quantile(diff[1:k], 0.025, names = FALSE)),
         running_p975 = sapply(n, function(k) quantile(diff[1:k], 0.975, names = FALSE)))
  }, by = .(who_region, condition)]
}
