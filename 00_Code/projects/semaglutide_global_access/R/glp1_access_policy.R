#' -------------------
#' GEM - GLP-1 eligibility tiers, access policy and uptake
#' tiered eligibility, access and uptake for GLP-1 RA initiation
#' -------------------
#'
#' METHODS NOTE
#'  eligibility (adults 18+): union of 6 mutually exclusive priority tiers, evaluated each year on current state
#'    1  T2D, BMI >= 30                                    (WHO EML 2025; ADA; KDIGO)
#'    2  T2D, BMI 25-29.9, CVD or CKD                      (WHO EML 2025; KDIGO; ADA)
#'    3  no T2D, BMI >= 30, CVD or CKD                     (SELECT; PBAC; CDA-AMC; CMS)
#'       option tier3_bmi27_cvd: also BMI 27-29.9 with CVD (CDA-AMC; CMS Bridge)
#'    4  no T2D/CVD/CKD, BMI >= 40                         (proxy for obesity + moderate risk; NHS TA1026 phasing)
#'    5  no T2D/CVD/CKD, BMI 35-39.9                       (CMS Bridge; NICE TA875)
#'    6  no T2D/CVD/CKD, BMI 30-34.9                       (WHO 2025 obesity guideline)
#'    CVD includes stroke history (stroke sets cvd = 1)
#'  age limit: no start or restart at age >= glp1$age_max_start (default 85; Inf = none)
#'    Wegovy SmPC 4.2/4.4: "Therapeutic experience in patients >=85 years of age is limited"
#'    persons already on treatment continue past the limit; year-0 baseline users are not affected
#'  access: person-fixed uniform u_access (CRN stream 4, year 0); access if u_access < a[tier]
#'    -> access is stable over time, re-evaluated when tier changes, and nested across scenarios
#'    strategies: "random"   a = coverage_rate for every tier (comparability with earlier % access runs)
#'                "tiered"   a = 1 for tiers_open, 0 otherwise (coverage expansion steps)
#'                "capacity" per country, coverage_rate x eligible (year 0) filled top tier first; a fixed over time
#'  uptake: annual probability that an eligible person with access and never treated starts treatment
#'    by indication (T2D vs obesity), sex, age (sheet glp1_uptake; scenarios "base", "optimistic"; "full" = 1)
#'  continuation: once started, treatment continues/stops/restarts per treatment_dynamics, regardless of later tier
#'  baseline use: prevalent GLP-1 RA users at year 0 in both arms (status quo and policy), sheet glp1_baseline_use
#'    T2D: country estimate, else adult-population-weighted mean of estimates in the same income pool
#'         (HIC excl. US; UMIC + LMIC); LIC 0. Within country: NCHS 2024 age x BMI pattern rescaled to the target
#'    obesity (BMI >= 30, no T2D): HIC 10% (low 5%, high 15%); other income groups 0
#'    observed baseline BMI is on drug: off-drug BMI percentile = BMI / (1 + weight-loss response)
#'    no direct treatment effect for baseline users (input rates already reflect existing use); drug cost kept
#'    mode "maintain": baseline users stay on treatment (prevalence is a steady state of starts and stops);
#'         "dynamics": same discontinuation/restart as new users

glp1_tiers <- c("1: T2D, BMI>=30", "2: T2D, BMI 25-29.9 + CVD/CKD", "3: no T2D, BMI>=30 + CVD/CKD",
                "4: BMI>=40", "5: BMI 35-39.9", "6: BMI 30-34.9")
glp1_strategies <- c("random", "tiered", "capacity")

# tier from state vectors (NA = not eligible) ----
glp1_tier_vec <- function(age, bmi, t2d, cvd, ckd, alive, tier3_bmi27_cvd = FALSE, age_max = Inf) {
  cmb <- cvd == 1 | ckd == 1
  t3  <- t2d == 0 & cmb & bmi >= 30
  if (tier3_bmi27_cvd) t3 <- t3 | (t2d == 0 & cvd == 1 & bmi >= 27)
  fcase(!(alive & age >= 18 & age < age_max), NA_integer_,
        t2d == 1 & bmi >= 30,       1L,
        t2d == 1 & bmi >= 25 & cmb, 2L,
        t3,                         3L,
        t2d == 0 & bmi >= 40,       4L,
        t2d == 0 & bmi >= 35,       5L,
        t2d == 0 & bmi >= 30,       6L,
        default = NA_integer_)
}

# set glp1_tier, treatment_eligible, has_access (and ever_ flags) by reference ----
#  a: access probability by tier (length 6)
set_tier_access <- function(dt, a, tier3_bmi27_cvd = FALSE, age_max = Inf) {
  tier <- glp1_tier_vec(dt$age, dt$bmi, dt$t2d, dt$cvd, dt$ckd, dt$ac_death == 0, tier3_bmi27_cvd, age_max)
  acc  <- as.integer(!is.na(tier) & dt$u_access < a[pmax(1L, tier)])
  dt[, `:=`(glp1_tier = tier, treatment_eligible = as.integer(!is.na(tier)), has_access = acc)]
  dt[treatment_eligible == 1L, ever_eligible := 1L]
  dt[has_access == 1L, ever_access := 1L]
  invisible(dt)
}

# capacity fill: share 'target' of eligible (weighted N_k) gets access, highest tier first ----
fill_tiers <- function(N_k, target) {
  cap <- target * sum(N_k); a <- numeric(length(N_k))
  for (k in seq_along(N_k)) {
    if (N_k[k] <= 0) { a[k] <- as.numeric(cap > 0); next }   # empty tier: open only if capacity remains
    a[k] <- min(1, max(0, cap) / N_k[k]); cap <- cap - a[k] * N_k[k]
  }
  a
}

# uptake table for a scenario (sheet glp1_uptake), or numeric annual probability ----
load_uptake <- function(uptake, fp = input_param_fp) {
  if (is.numeric(uptake)) {
    if (length(uptake) != 1 || uptake < 0 || uptake > 1) stop("numeric uptake must be one value in 0-1")
    return(uptake)
  }
  if (identical(uptake, "full")) return(1)
  u <- as.data.table(read_excel(fp, sheet = "glp1_uptake"))[in_use == "yes"]
  u <- u[scenario == uptake, .(indication, female = as.integer(female), age_low = as.numeric(age_low),
                              age_high = as.numeric(age_high), value = as.numeric(value))]
  if (!nrow(u)) stop("no glp1_uptake rows for scenario '", uptake, "'")
  u
}

# annual uptake probability per person; later rows override earlier (general -> specific) ----
uptake_prob <- function(t2d, female, age, upt) {
  if (is.numeric(upt)) return(rep(upt, length(age)))
  ind <- ifelse(t2d == 1, "t2d", "obesity")
  p <- rep(NA_real_, length(age))
  for (r in seq_len(nrow(upt))) {
    i <- ind == upt$indication[r] & age >= upt$age_low[r] & age <= upt$age_high[r] &
         (is.na(upt$female[r]) | female == upt$female[r])
    p[i] <- upt$value[r]
  }
  if (anyNA(p)) stop("glp1_uptake rows do not cover all indication-sex-age combinations")
  p
}

# baseline tier counts (weighted) from year-0 sample files ----
baseline_tier_counts <- function(baseline_fp, tier3_bmi27_cvd = FALSE, age_max = Inf) {
  d <- as.data.table(open_dataset(baseline_fp, format = "parquet") |>
                       select(iso3, w, age, bmi, t2d, cvd, ckd, ac_death) |> collect())
  d[, tier := glp1_tier_vec(age, bmi, t2d, cvd, ckd, ac_death == 0, tier3_bmi27_cvd, age_max)]
  d[!is.na(tier), .(N = sum(w)), by = .(iso3, tier)]
}

# access policy for one simulation config (NULL = no treatment) ----
#  a_default: access by tier for all countries; a_iso3: country-specific (capacity strategy)
make_access_policy <- function(cs, baseline_fp = NULL) {
  if (!isTRUE(cs$treatment_enabled)) return(NULL)
  g <- cs$glp1
  if (!g$strategy %in% glp1_strategies) stop("glp1$strategy must be one of: ", paste(glp1_strategies, collapse = ", "))
  if (g$strategy != "tiered" && (g$coverage_rate <= 0 || g$coverage_rate > 1)) stop("glp1$coverage_rate must be in (0, 1]")
  if (g$strategy == "tiered" && !all(g$tiers_open %in% 1:6)) stop("glp1$tiers_open must be within 1:6")
  pol <- list(strategy = g$strategy, tier3_bmi27_cvd = isTRUE(g$tier3_bmi27_cvd), age_max = g$age_max_start,
              uptake = load_uptake(g$uptake), a_iso3 = NULL,
              a_default = switch(g$strategy,
                                 random   = rep(g$coverage_rate, 6),
                                 tiered   = as.numeric(1:6 %in% g$tiers_open),
                                 capacity = rep(0, 6)))
  if (g$strategy == "capacity") {
    if (is.null(baseline_fp)) stop("capacity strategy needs the baseline sample folder")
    cnt  <- baseline_tier_counts(baseline_fp, pol$tier3_bmi27_cvd, pol$age_max)
    full <- CJ(iso3 = unique(cnt$iso3), tier = 1:6)
    full[cnt, N := i.N, on = .(iso3, tier)][is.na(N), N := 0]
    pol$a_iso3 <- full[order(iso3, tier), .(tier, a = fill_tiers(N, g$coverage_rate)), by = iso3]
  }
  pol
}

# access vector for one country
policy_access <- function(pol, cc) {
  if (is.null(pol$a_iso3)) return(pol$a_default)
  x <- pol$a_iso3[pol$a_iso3$iso3 == cc]
  if (nrow(x) == 6) x$a[order(x$tier)] else rep(0, 6)     # country absent from year-0 eligible counts
}

# baseline GLP-1 RA use inputs for a scenario ("base" | "low" | "high"; "none" -> NULL) ----
#  returns targets by iso3 (t2d, obesity shares), NCHS relative pattern by age and BMI
baseline_use_inputs <- function(scen = "base", fp = input_param_fp) {
  if (identical(scen, "none")) return(NULL)
  if (!scen %in% c("base", "low", "high")) stop("baseline_use must be one of: base, low, high, none")
  b <- as.data.table(read_excel(fp, sheet = "glp1_baseline_use"))[in_use == "yes"]
  num <- c("age_low", "age_high", "bmi_low", "bmi_high", "value")
  b[, (num) := lapply(.SD, as.numeric), .SDcols = num]
  cg <- cg <- as.data.table(read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")))[
    who_region != "no who region", .(iso3, inc = world_bank_income_group)]
  pop <- as.data.table(cohort_input_data)[, .(pop = sum(pop_size)), by = iso3]
  cg <- pop[cg, on = "iso3"][is.na(pop), pop := 0]
  cg[, pool := fcase(inc == "high income", "hic", inc %in% c("upper middle income", "lower middle income"), "mic", default = "lic")]
  # T2D: country value, else pooled (adult-population-weighted; US excluded from HIC pool), LIC from sheet (0)
  cg[b[type == "country" & indication == "t2d"], t2d := i.value, on = "iso3"]
  pooled <- cg[!is.na(t2d) & iso3 != "USA" & pool != "lic", .(t2d_pool = weighted.mean(t2d, pmax(pop, 1))), by = pool]
  cg[pooled, t2d_pool := i.t2d_pool, on = "pool"]
  lic <- b[type == "income" & indication == "t2d" & income_group == "low income", value][1]
  cg[, t2d_src := fifelse(!is.na(t2d), "country", fifelse(pool == "lic", "lic_assumption", paste0("pooled_", pool)))]
  cg[is.na(t2d), t2d := fifelse(pool == "lic", lic, t2d_pool)]
  # obesity without T2D: by income group and scenario
  ob <- b[type == "income" & indication == "obesity" & scenario %in% c(scen, "all"), .(inc = income_group, obesity = value)]
  cg[ob, obesity := i.obesity, on = "inc"][is.na(obesity), obesity := 0]
  pat <- b[type == "pattern" & indication == "t2d"]
  list(scenario = scen,
       targets = cg[, .(iso3, inc, pop, t2d, t2d_src, obesity)],
       age = pat[!is.na(age_low)][order(age_low), .(lo = age_low, r = value)],
       bmi = pat[!is.na(bmi_low)][order(bmi_low), .(lo = bmi_low, r = value)])
}

# relative use among adults with T2D by age and BMI (NCHS pattern) ----
bl_pattern <- function(age, bmi, bl) bl$age$r[findInterval(age, bl$age$lo)] * bl$bmi$r[findInterval(bmi, bl$bmi$lo)]

# country scaling so the weighted mean among adults with T2D equals the target (whole baseline sample, all files) ----
baseline_use_scale <- function(bl, baseline_fp) {
  if (is.null(bl)) return(NULL)
  d <- as.data.table(open_dataset(baseline_fp, format = "parquet") |> select(iso3, w, age, bmi, t2d) |> collect()) |> setindex(NULL)
  d <- d[t2d == 1]
  k <- d[, .(t2d_k = sum(w) / sum(w * bl_pattern(age, bmi, bl))), by = iso3]
  bl$targets <- k[bl$targets, on = "iso3"][is.na(t2d_k), t2d_k := 1]
  bl
}

# assign baseline users at year 0 (by reference); same persons in every arm (CRN stream 4, columns 2-3) ----
#  bl: baseline_use_scale() output; mp: matched cohort parameters (rows aligned with dt) for the BMI distribution
assign_baseline_use <- function(dt, bl, crn_id, mp) {
  if (is.null(bl)) return(invisible(dt))
  tg <- bl$targets[iso3 == dt$iso3[1]]
  if (!nrow(tg)) return(invisible(dt))
  p  <- numeric(nrow(dt))
  t2 <- which(dt$t2d == 1 & dt$ac_death == 0)
  if (length(t2) && tg$t2d > 0) p[t2] <- pmin(1, tg$t2d * tg$t2d_k * bl_pattern(dt$age[t2], dt$bmi[t2], bl))
  ob <- which(dt$t2d == 0 & dt$bmi >= 30 & dt$ac_death == 0)
  p[ob] <- tg$obesity
  u  <- crn_u(crn_id, 0L, nrow(dt), 3, stream = 4L)      # col 1 = u_access (same draws), 2 = use, 3 = response
  bu <- which(u[, 2] < p)
  if (length(bu)) {
    wl <- draw_wl_resp(u[bu, 3], dt$t2d[bu], glp1_impact_data)
    wl <- pmax(wl, dt$bmi[bu] / bmi_max - 1)               # off-drug BMI within the BMI support (keeps observed BMI)
    dt[bu, `:=`(baseline_user = 1L, treated = 1L, ever_treated = 1L, years_treated = 1L, wl_resp = wl, wl_effect = wl)]
    dt[bu, bmi_u := bmi_to_u(bmi / (1 + wl_resp), mp[bu])]
  }
  invisible(dt)
}
