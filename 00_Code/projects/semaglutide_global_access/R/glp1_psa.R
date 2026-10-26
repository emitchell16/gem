# projects/semaglutide_global_access/R/glp1_psa.R

# description ------
# GLP-1 parts of the probabilistic sensitivity analysis (core: GEM_Model_PSA_Functions.R):
#   hook_psa_inputs(): GLP-1 inputs captured with the base inputs
#   hook_psa_draw():   one z per outcome shared across subgroups (HR log-normal; BMI change normal)
#   psa_icer_fixed_drug(): ICER intervals with drug price fixed (post-PSA)

hook_psa_inputs <- function() list(glp1_impact_data = get("glp1_impact_data", envir = .GlobalEnv))

hook_psa_draw <- function(base, out, zf) {
  te <- as.data.frame(base$ext$glp1_impact_data)
  for (k in unique(te$outcome)) {
    i <- te$outcome == k; z <- zf(paste0("te_", k))
    is_bmi <- grepl("bmi", te$outcome_unit[i], ignore.case = TRUE)
    v <- te$outcome_val[i]; l <- te$outcome_lower[i]; u <- te$outcome_upper[i]
    te$outcome_val[i] <- ifelse(is_bmi, draw_normal(v, pmin(l, u), pmax(l, u), z), draw_lognormal(v, l, u, z))
  }
  out$ext$glp1_impact_data <- te
  out
}

# ICER intervals with drug cost fixed (post-PSA; no rerun) ----
#' drug price = scenario variable, not sampled. Per PSA run r and region:
#'   incremental cost_r = condition cost change_r (cost_conditions_total, PSA) + incremental drug cost (fixed, primary run)
#'   ICER_r = incremental cost_r / DALYs averted_r (daly_total_adj, PSA)
#' fixed drug cost by region and price = sum over countries of incremental cost(price) - incremental cost(tc_zero)
#'   (cost_effect_by_country.csv, primary run of the same scenario, ce_config)
#' res: PSA run results (psa_results_all_runs*.csv, or psa_runs_partial.rds with batch rows)
#'   batch rows present: DALYs averted and condition cost change shrunk toward their mean per component (as psa_summarize())
#'     to remove Monte Carlo noise; otherwise run values used as is (wider, conservative)
#' outputs: icer (per region x price: mean and 2.5-97.5% of DALYs averted, condition cost, total incremental cost, ICER;
#'   share of runs with DALYs averted > 0 and cost-saving) and ceac (share of runs with net monetary benefit > 0 per WTP)
#' limitation: drug-cost variation from time on treatment (mortality, discontinuation) not included
psa_icer_fixed_drug <- function(res, ce_csv, ce_config = "tier1to6", B = 10,
                                prices = c("id_sust_costx1", "id_sust_costx3", "intl_rwe"),
                                wtp = c(0, 100, 250, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000), out_dir = NULL) {
  res <- as.data.table(res)[, who_region := tolower(trimws(who_region))]
  if (!"batch" %in% names(res)) res[, batch := 0L]
  f <- res[batch == 0 & condition %in% c("daly_total_adj", "cost_conditions_total")]
  if (uniqueN(f$condition) < 2) stop("PSA results need daly_total_adj and cost_conditions_total")
  # Monte Carlo noise removed per component when batch rows exist
  if (any(res$batch > 0)) {
    wv <- res[batch > 0 & condition %in% unique(f$condition), .(within = var(diff) / B), by = .(run, who_region, condition)]
    f <- wv[f, on = c("run", "who_region", "condition")]
    f[, diff := {
      vt <- var(diff); vp <- max(0, vt - mean(within, na.rm = TRUE)); k <- if (vt > 0) sqrt(vp / vt) else 0
      mean(diff) + k * (diff - mean(diff))
    }, by = .(who_region, condition)]
    adj <- TRUE
  } else adj <- FALSE
  x <- dcast(f, run + who_region ~ condition, value.var = "diff")
  x[, `:=`(dalys_averted = -daly_total_adj, cond_cost = cost_conditions_total)]
  # fixed incremental drug cost by region, primary run
  ce <- fread(ce_csv)[config == ce_config & cost_scenario %in% c("tc_zero", prices)]
  ce[, who_region := tolower(trimws(who_region))]
  w <- dcast(ce, iso3 + who_region ~ cost_scenario, value.var = "incremental_cost_val")
  dr <- melt(w, id.vars = c("iso3", "who_region", "tc_zero"), measure.vars = intersect(prices, names(w)),
             variable.name = "price", value.name = "inc")[, .(drug = sum(inc - tc_zero, na.rm = TRUE)), by = .(who_region, price)]
  dr <- rbind(dr, dr[, .(who_region = "global", drug = sum(drug)), by = price])
  miss <- setdiff(unique(x$who_region), dr$who_region)
  if (length(miss)) message("psa_icer_fixed_drug: no drug cost for region(s): ", paste(miss, collapse = ", "))
  y <- dr[x, on = "who_region", allow.cartesian = TRUE, nomatch = NULL][, inc_cost := cond_cost + drug]
  q <- function(v, p) quantile(v, p, na.rm = TRUE, names = FALSE)
  icer <- y[, {
    pos <- dalys_averted > 0; r <- inc_cost[pos] / dalys_averted[pos]
    .(n_runs = .N, mc_adjusted = adj,
      dalys_averted_mean = mean(dalys_averted), dalys_averted_low = q(dalys_averted, 0.025), dalys_averted_high = q(dalys_averted, 0.975),
      cond_cost_mean = mean(cond_cost), cond_cost_low = q(cond_cost, 0.025), cond_cost_high = q(cond_cost, 0.975),
      drug_cost_fixed = drug[1],
      inc_cost_mean = mean(inc_cost), inc_cost_low = q(inc_cost, 0.025), inc_cost_high = q(inc_cost, 0.975),
      icer = if (mean(dalys_averted) > 0) mean(inc_cost) / mean(dalys_averted) else NA_real_,
      icer_low = if (any(pos)) q(r, 0.025) else NA_real_, icer_high = if (any(pos)) q(r, 0.975) else NA_real_,
      share_dalys_pos = mean(pos), share_cost_saving = mean(pos & inc_cost < 0))
  }, by = .(who_region, price)][order(who_region != "global", who_region, price)]
  ceac <- y[, .(wtp = wtp, p_cost_effective = sapply(wtp, function(l) mean(l * dalys_averted - inc_cost > 0))),
            by = .(who_region, price)]
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(icer, file.path(out_dir, "psa_icer_fixed_drug_by_region.csv"))
    fwrite(ceac, file.path(out_dir, "psa_ceac_fixed_drug_by_region.csv"))
    message("PSA ICER (drug cost fixed): ", out_dir)
  }
  list(icer = icer, ceac = ceac)
}
