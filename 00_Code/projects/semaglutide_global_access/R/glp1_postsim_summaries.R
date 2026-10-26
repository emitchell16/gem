# projects/semaglutide_global_access/R/glp1_postsim_summaries.R
# description ------
# Post-simulation summaries for the GLP-1 access analysis: access and coverage, distribution of benefit,
# scenario comparison tables (vs each scenario's own status quo), tier increments, status quo summary.
# Uses the manifest written by run_gem() and the core summaries in GEM_Model_Postsim_Functions.R.

# Distributional results (person-level pairing across arms; valid with common random numbers) ----------
#  DALYs averted per person = base - intervention (discounted), joined by pid
#  returns: weighted DALYs averted by group; share of total benefit accruing to top 10% of treated persons
#  group_vars options: "bmi_class0", "wb_income", "who_region", "adherence" (never / persistent / stopped / restarted)
benefit_distribution <- function(base_path, int_path, group_vars = c("bmi_class0", "wb_income", "adherence")) {
  cols <- c("pid", "iso3", "w", "yll_adj", "yld_adj_total", "bmi_class0",
            "ever_treated", "disc_count", "restart_count", "years_treated")
  b <- open_dataset(normalizePath(base_path)) |> select(any_of(cols)) |> collect() |> as.data.table()
  i <- open_dataset(normalizePath(int_path))  |> select(any_of(cols)) |> collect() |> as.data.table()
  d <- b[, .(pid, iso3, w, bmi_class0, daly_b = yll_adj + yld_adj_total)][
         i[, .(pid, ever_treated, disc_count, restart_count, years_treated, daly_i = yll_adj + yld_adj_total)], on = "pid"]
  d[, daly_averted := daly_b - daly_i]
  d[, adherence := fcase(ever_treated == 0, "never treated",
                         disc_count == 0, "persistent",
                         restart_count > 0, "stopped and restarted",
                         default = "stopped")]
  cg <- as.data.table(read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")))[, .(iso3, who_region, wb_income = world_bank_income_group)]
  d <- cg[d, on = "iso3"]
  by_group <- lapply(intersect(group_vars, names(d)), function(g)
    d[, .(pop = sum(w), dalys_averted = sum(w * daly_averted), dalys_averted_per_1000 = 1000 * sum(w * daly_averted) / sum(w)),
      by = .(group = get(g))][, group_var := g][])
  # share of net benefit to top 10% of treated (weighted; can exceed 1 when some treated have net DALY loss)
  tr <- d[ever_treated == 1][order(-daly_averted)]
  tr[, cum_w := cumsum(w) / sum(w)]
  top10_share <- tr[cum_w <= 0.10, sum(w * daly_averted)] / tr[, sum(w * daly_averted)]
  list(by_group = rbindlist(by_group, fill = TRUE), top10_share_treated = top10_share)
}

# Access and treatment coverage of the eligible population (intervention arm) ----
#  year 0 (= access at start of year 1): eligible, with access, by tier
#  final year: ever eligible / ever with access / ever treated over the horizon (incl. those who died),
#              eligible with access at start of final year, treated during final year, person-years on treatment
#  levels: global, who_region, wb_income, iso3; tier-specific rows at year 0
#  baseline users (prevalent use at year 0, both arms) reported separately; *_new = started under the policy
access_summary <- function(sim_hive_folder, scen_label, sim_years) {
  cg <- as.data.table(read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")))[, .(iso3, who_region, wb_income = world_bank_income_group)]
  rd <- function(yr, cols) {
    d <- as.data.table(open_dataset(file.path(sim_hive_folder, paste0("sim=", scen_label), paste0("year=", yr))) |>
                         select(any_of(cols)) |> collect())
    if (!"baseline_user" %in% names(d)) d[, baseline_user := 0L]      # older runs without baseline users
    cg[d, on = "iso3"]
  }
  y0 <- rd(0, c("iso3", "w", "glp1_tier", "treatment_eligible", "has_access", "baseline_user"))
  yT <- rd(sim_years, c("iso3", "w", "treatment_eligible", "has_access", "ever_eligible", "ever_access", "ever_treated",
                        "on_frac", "years_treated", "baseline_user"))
  lv <- list(global = NULL, who_region = "who_region", wb_income = "wb_income", iso3 = "iso3")
  one <- function(by) {
    a <- y0[, .(adults_y0 = sum(w), eligible_y0 = sum(w * treatment_eligible), access_y0 = sum(w * has_access),
                baseline_users_y0 = sum(w * baseline_user), eligible_baseline_users_y0 = sum(w * treatment_eligible * baseline_user)), by = by]
    b <- yT[, .(eligible_end = sum(w * treatment_eligible), access_end = sum(w * has_access),
                ever_eligible = sum(w * ever_eligible), ever_access = sum(w * ever_access),
                ever_treated = sum(w * ever_treated), ever_treated_new = sum(w * ever_treated * (1L - baseline_user)),
                treated_final_year = sum(w * (on_frac > 0)), person_years_treated = sum(w * years_treated)), by = by]
    out <- if (is.null(by)) cbind(a, b) else a[b, on = by]
    out[, `:=`(pct_adults_baseline_use = 100 * baseline_users_y0 / adults_y0,
               pct_eligible_baseline_use_y0 = 100 * eligible_baseline_users_y0 / eligible_y0,
               pct_eligible_access_y0 = 100 * access_y0 / eligible_y0,
               pct_ever_eligible_access = 100 * ever_access / ever_eligible,
               pct_ever_eligible_treated = 100 * ever_treated / ever_eligible,
               pct_ever_eligible_treated_new = 100 * ever_treated_new / ever_eligible,
               pct_eligible_treated_final_year = 100 * treated_final_year / eligible_end)]
    out[, `:=`(level = if (is.null(by)) "global" else by, group = if (is.null(by)) "global" else as.character(get(by)))]
    if (!is.null(by)) out[, (by) := NULL]
    out
  }
  overall <- rbindlist(lapply(lv, one), use.names = TRUE)
  tiers <- y0[!is.na(glp1_tier), .(eligible_y0 = sum(w), access_y0 = sum(w * has_access)), by = .(tier = glp1_tier)][order(tier)]
  tiers[, `:=`(tier_label = glp1_tiers[tier], pct_of_eligible = 100 * eligible_y0 / sum(eligible_y0),
               pct_access = 100 * access_y0 / eligible_y0)]
  list(overall = overall[, scenario := scen_label][], by_tier = tiers[, scenario := scen_label][])
}

# Combine scenarios: access, outcomes and standard effect measures; one row per scenario x geography x outcome ----
#  outcomes: weighted totals in final-year files (DALYs discounted; deaths and new cases cumulative over the horizon)
#  effect measures (scenario vs its own status quo arm; baseline users are in both arms):
#    averted, relative_change_pct            population impact
#    per_1000_treated, nnt                   per person started under the policy (ever_treated_new); NNT = treated / averted
#    per_1000_py                             per person-year on treatment added by the policy (stop year counted as 1 year)
#    per_100k_adults                         per 100,000 adults at year 0
#  increments: tier step k vs tier step k-1 (same uptake, effect source, horizon, status quo) -> marginal gain of each tier
#  files: all_scenarios_access.csv, all_scenarios_outcomes.csv, scenario_increments.csv, results_summary_fmt.xlsx
compile_scenarios <- function(out_root, manifest, res_dir,
                              outcomes = c("daly_total_adj", "yll_adj", "yld_adj_total", "ac_death", "t2d_new",
                                           "cvd_new", "ckd_new", "eskd_new", "stroke_n", "cvd_event_n")) {
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- as.data.table(manifest)
  cg <- as.data.table(read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")))[, .(iso3, who_region, wb_income = world_bank_income_group)]
  vars <- c(outcomes, "years_treated")
  tot <- function(label, yrs) {
    d <- cg[as.data.table(open_dataset(file.path(out_root, paste0("sim=", label), paste0("year=", yrs))) |>
                            select(any_of(c("iso3", "w", "yll_adj", "yld_adj_total", vars))) |> collect()), on = "iso3"]
    d[, daly_total_adj := yll_adj + yld_adj_total]
    rbindlist(lapply(list(global = NULL, who_region = "who_region", wb_income = "wb_income", iso3 = "iso3"), function(by) {
      x <- d[, lapply(.SD, function(v) sum(w * v)), by = by, .SDcols = intersect(vars, names(d))]
      x[, `:=`(level = if (is.null(by)) "global" else by, group = if (is.null(by)) "global" else as.character(get(by)))]
      if (!is.null(by)) x[, (by) := NULL]
      x
    }))
  }
  acc <- list(); res <- list()
  for (k in seq_len(nrow(manifest))) {
    m <- manifest[k]
    a <- access_summary(out_root, m$label, m$sim_years)
    acc[[k]] <- a$overall
    fwrite(a$by_tier, file.path(res_dir, paste0("access_by_tier_", m$label, ".csv")))
    bl <- if ("base_label" %in% names(m) && !is.na(m$base_label)) m$base_label else paste0(m$sim_years, "yr_basecase")
    b <- melt(tot(bl, m$sim_years), id.vars = c("level", "group"), value.name = "base")
    i <- melt(tot(m$label, m$sim_years), id.vars = c("level", "group"), value.name = "int")
    res[[k]] <- b[i, on = c("level", "group", "variable")][, `:=`(averted = base - int, scenario = m$label)]
  }
  acc <- rbindlist(acc); res <- rbindlist(res)
  # person-years on treatment added by the policy
  py <- res[variable == "years_treated", .(scenario, level, group, py_added = int - base)]
  res <- res[variable != "years_treated"]
  res <- acc[, .(scenario, level, group, adults_y0, ever_treated_new, pct_ever_eligible_access, pct_eligible_access_y0)][
    res, on = c("scenario", "level", "group")]
  res <- py[res, on = c("scenario", "level", "group")]
  res[, `:=`(relative_change_pct = fifelse(base > 0, 100 * (int - base) / base, NA_real_),
             averted_per_1000_treated = fifelse(ever_treated_new > 0, 1000 * averted / ever_treated_new, NA_real_),
             averted_per_1000_py = fifelse(py_added > 0, 1000 * averted / py_added, NA_real_),
             averted_per_100k_adults = 1e5 * averted / adults_y0,
             nnt = fifelse(averted > 0, ever_treated_new / averted, NA_real_))]
  # scenario descriptions and increments between tier steps ----
  desc <- scenario_descriptions(manifest)
  res <- desc[res, on = "scenario"]
  inc <- scenario_increments(res, desc)
  fwrite(acc, file.path(res_dir, "all_scenarios_access.csv"))
  fwrite(res, file.path(res_dir, "all_scenarios_outcomes.csv"))
  if (nrow(inc)) fwrite(inc, file.path(res_dir, "scenario_increments.csv")) else
    message("scenario_increments.csv not written: needs >= 2 tier steps (e.g. tier1to5 + tier1to6) in the run")
  write_results_summary(res, inc, file.path(res_dir, "results_summary_fmt.xlsx"))
  invisible(list(access = acc, outcomes = res, increments = inc))
}

# readable scenario description and grouping keys from the manifest ----
scenario_descriptions <- function(manifest) {
  m <- copy(as.data.table(manifest))
  for (v in c("strategy", "tiers_open", "uptake", "te_source", "coverage_rate", "base_label"))
    if (!v %in% names(m)) m[, (v) := NA]
  m[, n_tiers := nchar(as.character(tiers_open))]
  m[, description := as.character(label)]
  m[strategy %in% "tiered",   description := paste0("Tiers 1", fifelse(n_tiers > 1, paste0("-", n_tiers), ""), ", ", uptake, " uptake")]
  m[strategy %in% "capacity", description := paste0("Capacity ", round(100 * as.numeric(coverage_rate)), "% of eligible, ", uptake, " uptake")]
  m[strategy %in% "random",   description := paste0("Random ", round(100 * as.numeric(coverage_rate)), "% access, ", uptake, " uptake")]
  m[!is.na(te_source) & te_source != "rct", description := paste0(description, ", real-world effects")]
  m[sim_years != 5, description := paste0(description, ", ", sim_years, "-year horizon")]
  if ("mortality_calibration" %in% names(m)) m[mortality_calibration %in% FALSE, description := paste0(description, ", no mortality calibration")]
  m[, step_group := paste(strategy, uptake, te_source, sim_years, base_label, sep = "|")]
  m[, .(scenario = label, description, strategy, n_tiers, step_group)]
}

# marginal gain of each tier step: tiers 1..k vs tiers 1..k-1 (same step group) ----
scenario_increments <- function(res, desc) {
  steps <- desc[strategy %in% "tiered"][order(step_group, n_tiers)]
  if (!nrow(steps)) return(data.table())
  steps[, prev := shift(scenario), by = step_group]
  steps <- steps[!is.na(prev)]
  if (!nrow(steps)) return(data.table())
  cols <- c("scenario", "level", "group", "variable", "averted", "ever_treated_new", "py_added")
  x <- res[, ..cols]
  out <- x[steps[, .(scenario, prev, description)], on = "scenario"]
  out <- x[, .(prev = scenario, level, group, variable, averted_prev = averted, treated_prev = ever_treated_new,
               py_prev = py_added)][out, on = c("prev", "level", "group", "variable")]
  out[, `:=`(marginal_averted = averted - averted_prev,
             marginal_treated = ever_treated_new - treated_prev,
             marginal_py = py_added - py_prev)]
  out[, `:=`(marginal_per_1000_treated = fifelse(marginal_treated > 0, 1000 * marginal_averted / marginal_treated, NA_real_),
             marginal_per_1000_py = fifelse(marginal_py > 0, 1000 * marginal_averted / marginal_py, NA_real_))]
  out[, .(scenario, description, compared_with = prev, level, group, variable, marginal_treated, marginal_py,
          marginal_averted, marginal_per_1000_treated, marginal_per_1000_py)]
}

# formatted headline tables (main paper / eSupplement) ----
#  sheet 1 global: each scenario vs no intervention; sheet 2 tier-step increments; sheets 3-4 by WHO region / income group
write_results_summary <- function(res, inc, fp) {
  key <- c(ac_death = "Deaths", daly_total_adj = "DALYs (discounted)", t2d_new = "New T2D", cvd_new = "New CVD",
           ckd_new = "New CKD", eskd_new = "New ESKD", stroke_n = "Strokes")
  fmt <- function(d) {
    d <- d[variable %in% names(key)]
    d[, Outcome := factor(key[as.character(variable)], levels = key)]
    d[order(description, Outcome), .(
      Scenario = description, Region = group, Outcome,
      `Persons started on treatment` = round(ever_treated_new),
      `Person-years on treatment added` = round(py_added),
      `Status quo` = round(base), Intervention = round(int), Averted = round(averted),
      `Relative change (%)` = round(relative_change_pct, 2),
      `Averted per 1,000 treated` = round(averted_per_1000_treated, 2),
      `Averted per 1,000 person-years` = round(averted_per_1000_py, 2),
      `Averted per 100,000 adults` = round(averted_per_100k_adults, 1),
      NNT = round(nnt))]
  }
  sheets <- list(global = fmt(res[level == "global"])[, Region := NULL][],
                 who_region = fmt(res[level == "who_region"]),
                 income_group = fmt(res[level == "wb_income"]))
  if (nrow(inc)) {
    i <- inc[level == "global" & variable %in% names(key)]
    sheets$tier_increments <- i[, .(Scenario = description, `Compared with` = compared_with,
                                     Outcome = key[as.character(variable)],
                                     `Additional persons treated` = round(marginal_treated),
                                     `Additional person-years` = round(marginal_py),
                                     `Additional averted` = round(marginal_averted),
                                     `Per 1,000 additionally treated` = round(marginal_per_1000_treated, 2),
                                     `Per 1,000 additional person-years` = round(marginal_per_1000_py, 2))]
  }
  write_xlsx(sheets, fp)
  invisible(sheets)
}

# Status quo arms: baseline GLP-1 use (year 0) and main outcomes (final year) ----
#  compares arms that differ in baseline use or mortality calibration (e.g. 5yr_basecase vs 5yr_basecase_nocal)
#  obesity prevalence among living at end; cumulative mortality and T2D incidence (T2D-free at year 0) over the horizon
status_quo_summary <- function(out_root, arms, res_dir) {
  cg <- as.data.table(read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")))[, .(iso3, who_region, wb_income = world_bank_income_group)]
  rd <- function(lab, yr, cols) {
    d <- as.data.table(open_dataset(file.path(out_root, paste0("sim=", lab), paste0("year=", yr))) |> select(any_of(cols)) |> collect())
    if (!"baseline_user" %in% names(d)) d[, baseline_user := 0L]
    cg[d, on = "iso3"]
  }
  out <- rbindlist(lapply(seq_len(nrow(arms)), function(k) {
    lab <- arms$base_label[k]; yrs <- arms$sim_years[k]
    y0 <- rd(lab, 0, c("pid", "iso3", "w", "t2d", "bmi", "baseline_user"))
    yT <- rd(lab, yrs, c("pid", "iso3", "w", "ac_death", "obesity", "t2d_new"))
    yT[y0[, .(pid, t2d0 = t2d)], t2d0 := i.t2d0, on = "pid"]
    rbindlist(lapply(list(global = NULL, who_region = "who_region", wb_income = "wb_income", iso3 = "iso3"), function(by) {
      a <- y0[, .(adults_y0 = sum(w),
                  pct_baseline_use = 100 * sum(w * baseline_user) / sum(w),
                  pct_baseline_use_t2d = 100 * sum(w * baseline_user * t2d) / sum(w * t2d),
                  pct_baseline_use_obese_not2d = 100 * sum(w * baseline_user * (t2d == 0 & bmi >= 30)) / sum(w * (t2d == 0 & bmi >= 30))), by = by]
      b <- yT[, .(pct_obesity_end = 100 * sum(w * (ac_death == 0) * obesity) / sum(w * (ac_death == 0)),
                  pct_died = 100 * sum(w * ac_death) / sum(w),
                  pct_t2d_incident = 100 * sum(w * t2d_new * (t2d0 %in% 0)) / sum(w * (t2d0 %in% 0))), by = by]
      x <- if (is.null(by)) cbind(a, b) else a[b, on = by]
      x[, `:=`(level = if (is.null(by)) "global" else by, group = if (is.null(by)) "global" else as.character(get(by)))]
      if (!is.null(by)) x[, (by) := NULL]
      x
    }))[, `:=`(arm = lab, sim_years = yrs)]
  }))
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(out, file.path(res_dir, "status_quo_summary.csv"))
  invisible(out)
}

# baseline eligibility by priority tier (weighted; year 0 of the intervention arm); called per scenario by generate_results_tables() ----
hook_results_extra <- function(baseline_ds_path, out_folder, iso_regions) {
  elig_fp <- file.path(out_folder, "treatment_eligibility_by_tier.xlsx")
  if (file.exists(elig_fp)) return(invisible(NULL))
  d <- as.data.table(open_dataset(normalizePath(baseline_ds_path)) |> select(iso3, w, glp1_tier) |> collect())
  d <- as.data.table(iso_regions)[d, on = "iso3"][is.na(who_region), who_region := "no who region"]
  tb <- rbind(d[, .(Region = "global", pop = sum(w)), ], d[, .(pop = sum(w)), by = .(Region = who_region)])
  el <- rbind(d[!is.na(glp1_tier), .(Region = "global", eligible = sum(w)), by = glp1_tier],
              d[!is.na(glp1_tier), .(eligible = sum(w)), by = .(Region = who_region, glp1_tier)])
  el <- tb[el, on = "Region"][, `:=`(tier = glp1_tiers[glp1_tier], pct_adults = 100 * eligible / pop)][order(Region, glp1_tier)]
  write_xlsx(el[, .(Region = str_to_title(Region), tier, eligible, pct_adults)], elig_fp)
  invisible(NULL)
}
