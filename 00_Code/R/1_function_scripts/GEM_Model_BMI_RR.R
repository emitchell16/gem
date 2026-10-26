# -------------------'
# GEM - continuous BMI: distributions (baseline + yearly tracking) and relative risks
# inputs (GEM_input_parameters.xlsx): obesity_rr_nejm, bmi_rr_continuous
#         NCD-RisC BMI category prevalences in cohort params (prevalence_bmi_*_val, prevalence_class*_obesity_val)
# adults 18+ only; bmi_to_u()
# -------------------'

# METHODS NOTE
#  BMI distribution: piecewise quantile function per country-age-sex
#    - category probabilities exact (NCD-RisC <18.5, 18.5-20, 20-25, 25-30, 30-35, 35-40, >=40)
#    - within-category shape from lognormal fit (probit(F) ~ ln BMI, least squares on category cut points)
#    - bounded [bmi_min, bmi_max]; RR cap harmonized with bmi_max
#  BMI over time: rank (percentile) tracking; each person keeps percentile u, BMI = Q_age(u)
#    -> reproduces NCD-RisC prevalence by age; individual class transitions occur as person ages
#  RR: RR_i = RR5_age ^ (max(0, min(BMI, bmi_max) - ref) / 5), ref = 22.5 (TMREL 20-25)

bmi_min <- 13
bmi_max <- 60                                   # draw upper bound = RR cap (assumes log-linear RR to 60)
bmi_cuts <- c(18.5, 20, 25, 30, 35, 40)
bmi_cat_cols <- c("prevalence_bmi_lt18_val", "prevalence_bmi_18_20_val", "prevalence_bmi_20_25_val",
                  "prevalence_bmi_25_30_val", "prevalence_class1_obesity_val",
                  "prevalence_class2_obesity_val", "prevalence_class3_obesity_val")

# BMI distribution parameters per cohort row ----
#  adds F1..F6 (cumulative prob at cuts) and bmi_ln_mu, bmi_ln_sd
add_bmi_dist_params <- function(dt) {
  setDT(dt)
  P <- as.matrix(dt[, ..bmi_cat_cols])
  P[is.na(P) | P < 0] <- 0
  rs <- rowSums(P); rs[rs <= 0] <- 1
  P <- P / rs                                             # categories sum to 1
  Fc <- t(apply(P, 1, cumsum))[, 1:6, drop = FALSE]       # F at 18.5..40
  colnames(Fc) <- paste0("bmi_F", 1:6)
  dt[, (colnames(Fc)) := as.data.table(Fc)]
  # lognormal fit: qnorm(F) = (ln c - mu) / sd  ->  ln c = mu + sd * qnorm(F)
  x  <- qnorm(pmin(pmax(Fc, 1e-4), 1 - 1e-4))
  lc <- matrix(log(bmi_cuts), nrow(x), 6, byrow = TRUE)
  xm <- rowMeans(x); lm_ <- rowMeans(lc)
  sd_ <- rowSums((x - xm) * (lc - lm_)) / pmax(rowSums((x - xm)^2), 1e-8)
  sd_ <- pmin(pmax(sd_, 0.05), 0.6)
  dt[, `:=`(bmi_ln_sd = sd_, bmi_ln_mu = lm_ - sd_ * xm)]
  invisible(dt)
}

# BMI from percentile u, given person-aligned param columns (from matched_params) ----
bmi_from_u <- function(u, mp) {
  Fk <- cbind(0, as.matrix(mp[, paste0("bmi_F", 1:6), with = FALSE]), 1)   # n x 8
  lo <- c(bmi_min, bmi_cuts); hi <- c(bmi_cuts, bmi_max)
  k  <- rowSums(u >= Fk[, 2:7, drop = FALSE]) + 1                        # interval 1..7
  idx <- cbind(seq_along(u), k)
  F0 <- Fk[idx]; F1 <- Fk[cbind(seq_along(u), k + 1)]
  v  <- ifelse(F1 > F0, (u - F0) / (F1 - F0), 0.5)
  a  <- plnorm(lo[k], mp$bmi_ln_mu, mp$bmi_ln_sd)
  b  <- plnorm(hi[k], mp$bmi_ln_mu, mp$bmi_ln_sd)
  out <- qlnorm(a + pmin(pmax(v, 0), 1) * (b - a), mp$bmi_ln_mu, mp$bmi_ln_sd)
  pmin(pmax(out, lo[k]), hi[k])
}

# inverse of bmi_from_u(): percentile of a given BMI in the person's distribution (validation baselines)
bmi_to_u <- function(bmi, mp) {
  Fk <- cbind(0, as.matrix(mp[, paste0("bmi_F", 1:6), with = FALSE]), 1)
  lo <- c(bmi_min, bmi_cuts); hi <- c(bmi_cuts, bmi_max)
  b  <- pmin(pmax(bmi, bmi_min), bmi_max)
  k  <- findInterval(b, c(bmi_min, bmi_cuts), rightmost.closed = TRUE)
  k  <- pmin(pmax(k, 1L), 7L)
  a  <- plnorm(lo[k], mp$bmi_ln_mu, mp$bmi_ln_sd)
  bb <- plnorm(hi[k], mp$bmi_ln_mu, mp$bmi_ln_sd)
  v  <- ifelse(bb > a, (plnorm(b, mp$bmi_ln_mu, mp$bmi_ln_sd) - a) / (bb - a), 0.5)
  F0 <- Fk[cbind(seq_along(b), k)]; F1 <- Fk[cbind(seq_along(b), k + 1)]
  pmin(pmax(F0 + pmin(pmax(v, 0), 1) * (F1 - F0), 1e-6), 1 - 1e-6)
}

# obesity indicators from BMI (adults; WHO cut-points)
set_obesity_from_bmi <- function(dt, idx) {
  b <- dt$bmi[idx]
  dt[idx, `:=`(obesity        = as.integer(b >= 30),
               class1_obesity = as.integer(b >= 30 & b < 35),
               class2_obesity = as.integer(b >= 35 & b < 40),
               class3_obesity = as.integer(b >= 40))]
  invisible(dt)
}

# relative risks ----
# age-specific per-5 RRs (Afshin 2017 Table S2), keyed t2d / stroke / cvd / ckd
load_nejm_age_rr <- function(fp, p_isch = 0.65) {             # p_isch = prop_stroke_ischaemic
  d <- read_excel(fp, sheet = "obesity_rr_nejm") |>
    transmute(key = tolower(trimws(RR_morbidity_and_mortality_obesity)),
              age_low = as.numeric(age_low),
              across(c(value, lower, upper), as.numeric)) |>
    filter(!is.na(value))
  # stroke: ischaemic/haemorrhagic weighted on log scale (bounds with bounds)
  stroke <- d |>
    filter(key %in% c("ischaemic stroke", "hemorrhagic stroke")) |>
    mutate(w = if_else(key == "ischaemic stroke", p_isch, 1 - p_isch)) |>
    group_by(age_low) |>
    summarise(across(c(value, lower, upper), ~ exp(sum(w * log(.x)))), .groups = "drop") |>
    mutate(key = "stroke")
  bind_rows(d, stroke) |>
    mutate(key = recode(key, "ischaemic heart disease" = "cvd", "ckd due to other causes" = "ckd")) |>
    filter(key %in% c("t2d", "stroke", "cvd", "ckd")) |>
    arrange(key, age_low)
}

load_bmi_rr <- function(fp) {
  read_excel(fp, sheet = "bmi_rr_continuous") |>
    filter(!is.na(value)) |>
    mutate(across(c(ref_bmi, value, lower, upper), as.numeric))
}

rr_per5 <- function(bmi, rr5, ref = 22.5, cap = bmi_max) {
  exp(log(rr5) * pmax(0, pmin(bmi, cap) - ref) / 5)
}

#  cond: "t2d","stroke","cvd","ckd" (age-specific) or "death_other_nonmodeled"
bmi_rr_vec <- function(cond, bmi, age, age_tbl, rr_tbl, ref = 22.5) {
  bmi[is.na(bmi)] <- ref
  if (cond %in% age_tbl$key) {
    tb  <- age_tbl[age_tbl$key == cond, ]
    rr5 <- tb$value[pmax(1L, findInterval(age, tb$age_low))]   # <first band -> first band
    out <- rr_per5(bmi, rr5, ref)
    if (cond == "ckd") out[age < min(tb$age_low)] <- 1         # no source RR <35
  } else {
    tb  <- rr_tbl[rr_tbl$outcome == cond, ]
    if (nrow(tb) == 0) stop("no RR rows for ", cond)
    out <- rr_per5(bmi, tb$value[1], tb$ref_bmi[1])
  }
  out                                           # modeled ages 18+; youngest Afshin band applies at 18-24
}

# PSA draws for BMI RRs: draw_psa_params() in GEM_Model_PSA_Functions.R (one z per outcome, shared across ages)
