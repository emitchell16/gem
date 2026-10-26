# R/1_function_scripts/GEM_Model_Postsim_Plots.R
#
# Circle plots: one bar per country, grouped by WHO region or World Bank income group
#   plot_gem_circles()               one scenario vs its status quo: A obesity reduction (pp), B deaths averted
#   expansion_scenarios_circleplot() tiers 1..6 scenarios: bar = tiers 1-6 effect, stacked by tier step;
#                                    A by WHO region, B by income group; one figure per measure;
#                                    plus obesity + deaths figures, one for WHO region and one for income group
#   scales: shared across scenarios (and across panels in a figure) so bars can be compared
#   group labels around the plot (name in bold, group colour; detail in black): adults (year 0), [status quo vs intervention],
#   effect, change as a share of the status quo burden
#   works for test runs (iso3_subset) and full runs
#
#   cost-effectiveness (costing on): plot_ce_plane(), plot_country_icer(), plot_ce_max_price(), plot_ce_price_curve()
#
# two steps: data from the simulation outputs -> csv in res_dir/figures; plots from the csv
#   gem_plot_data(out_root, res_dir)  needs the sim=<label> folders (local temp); writes the csv files
#   plot functions with out_root = NULL read only the csv files (03_Out) -> plots on any computer
#   see projects/*/postsim_result_formatting.Rmd for an example of building all figures from the csv files
#
# use: plot_gem_circles(out_root, res_dir)                 # called by run_gem(): builds csv, then plots
#      plot_gem_circles(NULL, res_dir)                     # from csv only
#      expansion_scenarios_circleplot(NULL, res_dir)
# outputs (res_dir/figures): circle_iso3_<label>.png/.csv, expansion_<measure>.png, expansion_scenarios_by_level.csv,
#   expansion_obesity_deaths_<who_region|wb_income>.png, cost_effect_by_country.csv, ce_*.png, bar_*.png

# weighted sums per country for one arm (final year) ----
#' additive -> any grouping is a sum; obesity among living at end; deaths, DALYs over the horizon
#' costs (when costing was on): cost_total_<price scenario> and cost_conditions_total, discounted, summed with weights
arm_sums <- function(out_root, label, sim_years) {
  ds <- open_dataset(file.path(out_root, paste0("sim=", label), paste0("year=", sim_years)))
  cc <- grep("^cost_total_|^cost_conditions_total$", names(ds), value = TRUE)
  cols <- c("iso3", "w", "ac_death", "obesity", "yll_adj", "yld_adj_total", "ever_treated", "baseline_user", cc)
  d <- as.data.table(ds |> select(any_of(cols)) |> collect())
  for (v in c("baseline_user", "ever_treated", "yll_adj", "yld_adj_total")) if (!v %in% names(d)) d[, (v) := 0]
  out <- d[, .(adults = sum(w), alive = sum(w * (ac_death == 0)), obese = sum(w * (ac_death == 0) * obesity),
               deaths = sum(w * ac_death), dalys = sum(w * (yll_adj + yld_adj_total)),
               treated_new = sum(w * ever_treated * (1L - baseline_user))), by = iso3]
  if (length(cc)) out <- out[d[, lapply(.SD, function(v) sum(w * v)), by = iso3, .SDcols = cc], on = "iso3"]
  out[]
}

# country groups (fixed labels) ----
gem_country_groups <- function() {
  cg <- as.data.table(read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")))
  canon <- c("high income", "upper middle income", "lower middle income", "low income", "no wb income group")
  inc <- canon[match(gsub("[^a-z]", "", tolower(cg$world_bank_income_group)), gsub(" ", "", canon))]
  cg[, .(iso3, who_region, wb_income = fifelse(is.na(inc), world_bank_income_group, inc))]
}

# scenario vs status quo by level (iso3, who_region, wb_income, global) ----
sums_to_deltas <- function(b, i, cg) {
  x <- b[i, on = "iso3", nomatch = NULL]                         # i.* = intervention arm
  x <- cg[x, on = "iso3"]
  if ("i.adults" %in% names(x)) x[, i.adults := NULL]
  sc <- c("adults", "alive", "obese", "deaths", "dalys", "i.alive", "i.obese", "i.deaths", "i.dalys", "i.treated_new")
  agg <- function(by) {
    y <- x[, lapply(.SD, sum), by = by, .SDcols = sc]
    y[, `:=`(level = if (is.null(by)) "global" else by, group = if (is.null(by)) "global" else as.character(get(by)))]
    if (!is.null(by)) y[, (by) := NULL]
    y
  }
  y <- rbindlist(lapply(list(NULL, "who_region", "wb_income", "iso3"), agg), use.names = TRUE)
  y[level == "iso3", c("who_region", "wb_income") := cg[.SD, on = c(iso3 = "group"), .(who_region, wb_income)]]
  y[, `:=`(obesity_base = 100 * obese / alive, obesity_int = 100 * i.obese / i.alive)]
  y[, `:=`(obesity_reduction_pp = obesity_base - obesity_int,
           obesity_rel_change_pct = -100 * (obesity_base - obesity_int) / obesity_base,
           deaths_averted = deaths - i.deaths, dalys_averted = dalys - i.dalys)]
  y[, `:=`(deaths_averted_per_100k = 1e5 * deaths_averted / adults,
           deaths_averted_per_1000_treated = fifelse(i.treated_new > 0, 1000 * deaths_averted / i.treated_new, NA_real_),
           mortality_rel_change_pct = fifelse(deaths > 0, -100 * deaths_averted / deaths, NA_real_),
           dalys_averted_per_100k = 1e5 * dalys_averted / adults,
           daly_rel_change_pct = fifelse(dalys > 0, -100 * dalys_averted / dalys, NA_real_))]
  setnames(y, c("deaths", "dalys", "i.treated_new"), c("deaths_base", "dalys_base", "treated_new"))
  y[, c("alive", "obese", "i.alive", "i.obese", "i.deaths", "i.dalys") := NULL][]
}

# cost-effectiveness by country: incremental cost and DALYs averted per price scenario ----
#' cost_total_<sfx>: health-state + drug cost at that price; tc_zero: health-state costs only (drug free)
#' NULL when costing was off
ce_from_sums <- function(b, i, cg) {
  cc <- intersect(grep("^cost_total_|^cost_conditions_total$", names(b), value = TRUE), names(i))
  if (!length(cc)) return(NULL)
  x <- cg[b[i, on = "iso3", nomatch = NULL], on = "iso3"]
  rbindlist(lapply(cc, function(cn) x[, .(iso3, who_region, income_group = wb_income, n = adults,
    cost_scenario = if (cn == "cost_conditions_total") "tc_zero" else sub("^cost_total_", "", cn),
    incremental_cost_val = get(paste0("i.", cn)) - get(cn), daly_averted_val = dalys - i.dalys)]))[
    , icer_val := fifelse(daly_averted_val > 0, incremental_cost_val / daly_averted_val, NA_real_)][]
}

# country-level deltas for one scenario (kept for earlier callers) ----
iso3_deltas <- function(out_root, label, base_label, sim_years) {
  sums_to_deltas(arm_sums(out_root, base_label, sim_years), arm_sums(out_root, label, sim_years),
                 gem_country_groups())[, `:=`(scenario = label, base_label = base_label)][]
}

# WHO region colors (fixed per region, same in test and full runs); other names grey ----
region_palette <- function(regs) {
  pal <- c("africa" = "#046C8A", "americas" = "#E1BD6D", "eastern mediterranean" = "#D62728",
           "europe" = "#35278a", "south-east asia" = "#E377C2", "western pacific" = "#ABdded")
      # NOTE other color hex codes from wesanderson pals: #5B1A18 maroon #5BBCD6 sky blue #7FA48E sage; tableau blue #1F77B4
  out <- pal[tolower(regs)]; out[is.na(out)] = "grey50"
  setNames(unname(out), regs)
}

# darker shade for text and arcs ----
#' light colors (yellow, light blue) darkened more so text stays readable on white
darken <- function(cols, f = 0.8) {
  m <- grDevices::col2rgb(cols)
  lum <- colSums(m * c(0.299, 0.587, 0.114))
  m <- sweep(m, 2, ifelse(lum > 150, 0.6, f), `*`)
  setNames(grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255), names(cols))
}

# labels ----
grp_name <- function(x) {
  x <- tools::toTitleCase(x)
  x[x == "No Wb Income Group"] <- "Not classified"
  sub(" Income$", " income", sub(" Middle ", " middle ", x))
}
fmt_pop <- function(x) ifelse(x >= 1e9, sprintf("%.2fB", x / 1e9), ifelse(x >= 1e6, sprintf("%.0fM", x / 1e6), sprintf("%.0fk", x / 1e3)))
fmt_n <- function(x) formatC(round(x), format = "d", big.mark = ",")
# group label lines: adults (year 0) / [status quo vs intervention] / effect / relative change vs status quo
#' relative change = averted / status quo total (deaths, DALYs) or reduction / status quo prevalence (obesity):
#'   a % of the status quo burden, not a % change of the per-100k number
#' sq, int: status quo and intervention levels (obesity %, deaths per 100k); NULL = line left out
rel_txt <- function(rel, kind) {
  what <- switch(kind, pp = "prevalence", dalys = "DALYs", "deaths")
  if (kind == "pp") sprintf("%.1f%% %s %s", abs(rel), ifelse(rel <= 0, "lower", "higher"), what) else
    sprintf("%.1f%% %s %s", abs(rel), ifelse(rel <= 0, "fewer", "more"), what)
}
sq_int_txt <- function(sq, int, kind) {
  if (kind == "pp") sprintf("Status quo: %.1f%%\nIntervention: %.1f%%", sq, int) else
    paste0("Status quo: ", fmt_n(sq), " per 100k\nIntervention: ", fmt_n(int), " per 100k")
}
grp_lines <- function(val, rel, adults, kind, sq = NULL, int = NULL) {
  l1 <- switch(kind,
    pp      = ifelse(val >= 0, sprintf("%.1f pp lower", val), sprintf("%.1f pp higher", -val)),
    deaths  = paste(fmt_n(val), "per 100k averted"),
    dalys   = paste(fmt_n(val), "per 100k averted"),
    per1000 = sprintf("%.1f per 1,000 treated averted", val))
  paste0("Population: ", fmt_pop(adults), " adults\n", if (!is.null(sq)) paste0(sq_int_txt(sq, int, kind), "\n"),
         l1, "\n", rel_txt(rel, kind))
}
measure_info <- function(m, horizon) switch(m,
  deaths  = list(col = "deaths_averted_per_100k", rel = "mortality_rel_change_pct", unit = "", kind = "deaths",
                 title = paste0("Deaths averted per 100,000 adults over ", horizon, " years"), short = "deaths averted per 100,000 adults"),
  dalys   = list(col = "dalys_averted_per_100k",  rel = "daly_rel_change_pct",      unit = "", kind = "dalys",
                 title = paste0("DALYs averted per 100,000 adults over ", horizon, " years"), short = "DALYs averted per 100,000 adults"),
  obesity = list(col = "obesity_reduction_pp",    rel = "obesity_rel_change_pct",   unit = "pp", kind = "pp",
                 title = paste0("Reduction in obesity prevalence at year ", horizon, " (percentage points)"), short = "reduction in obesity prevalence (pp)"))
global_txt <- function(g, mi) {
  sprintf("%s (%s)", if (mi$kind == "pp") sprintf("%.2f pp lower", g[[mi$col]]) else paste(fmt_n(g[[mi$col]]), "per 100,000 averted"),
          rel_txt(g[[mi$rel]], mi$kind))
}

# one circular bar panel ----
#' dl: long data (iso3, grp, fill_key, v); bars stack outward in fill_key order; lim = shared c(lo, hi)
#' gl: group labels (grp, name, lines); grp_cols: arc and text colors by group; fill_vals: fill scale
#' country labels at the bar end; size scales with the number of countries
#' panel_in: approximate panel size (inches) in the saved figure; sets the country label size so labels
#'   are no taller than the space per bar (no overlap)
circle_panel <- function(dl, gl, lim, unit, title, grp_cols, fill_vals, show_fill_legend = FALSE, panel_in = 9) {
  tot <- dl[, .(tot = sum(v)), by = .(iso3, grp)]
  setorder(tot, grp, -tot)
  gap <- 2L
  tot[, id := seq_len(.N) + (rleid(grp) - 1L) * gap]
  n <- max(tot$id) + gap
  tot[, ang := 90 - 360 * (id - 0.5) / n]
  tot[, `:=`(hj = fifelse(ang < -90, 1, 0), ang = fifelse(ang < -90, ang + 180, ang))]
  dl <- tot[, .(iso3, id)][dl, on = "iso3"]
  if (is.factor(tot$grp)) gl <- copy(gl)[, grp := factor(grp, levels = levels(tot$grp))]
  lo <- min(0, lim[1]); hi <- max(0, lim[2]); rng <- max(hi - lo, 1e-9)
  arcs <- tot[, .(x = c(min(id) - 0.45, max(id) + 0.45), y = lo - 0.04 * rng), by = grp]
  # group labels outside, horizontal; side from angle; bold name over the detail lines (same block height)
  gpos <- tot[, .(xm = mean(range(id))), by = grp][gl, on = "grp", nomatch = NULL]
  gpos[, th := 2 * pi * (xm - 0.5) / n]
  gpos[, `:=`(hj = fifelse(sin(th) > 0.2, 0, fifelse(sin(th) < -0.2, 1, 0.5)),
              vj = fifelse(cos(th) > 0.6, 0, fifelse(cos(th) < -0.6, 1, 0.5)),
              name_lab = paste0(name, strrep("\n", lengths(strsplit(lines, "\n")))), line_lab = paste0("\n", lines))]
  brks <- pretty(c(lo, hi), n = 4); brks <- brks[brks >= lo & brks <= hi]
  blab <- if (unit == "pp") paste(brks, "pp") else formatC(brks, format = "d", big.mark = ",")
  # country labels on one ring outside all bars (more space per label than at the bar ends); thin leader lines
  #   coord_polar draws the circle at 80% of the panel (R mm); ring radius and label size solved together:
  #   cap height + gap (~1.05 x size) <= arc per bar at the ring; ring + label length + group labels fit within R
  R <- 0.4 * panel_in * 25.4
  y_lo <- lo - 0.18 * rng; y_lab <- hi + 0.03 * rng; A <- y_lab - y_lo + 0.02 * rng
  frac <- 0.8
  for (k in 1:3) {
    iso_size <- max(1.4, min(3.4, (2 * pi * R * frac / n) / 1.05))
    q <- min(0.4, 2.6 * iso_size / R)                              # label length (3 letters + pad) / R
    Ld <- q * A / (1 - q)
    frac <- (y_lab - y_lo) / (A + Ld)
  }
  y_hi <- y_lab + Ld + 0.02 * rng
  gy <- y_hi + 0.06 * rng                                           # group labels just outside the label ring
  lead <- tot[, .(id, y0 = pmax(tot, 0) + 0.01 * rng, y1 = y_lab - 0.005 * rng)][y1 > y0]

  ggplot() +
    geom_hline(yintercept = brks, colour = "grey88", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
    geom_col(data = dl, aes(x = id, y = v, fill = fill_key), width = 0.9, position = position_stack(reverse = TRUE)) +
    geom_path(data = arcs, aes(x = x, y = y, colour = grp, group = grp), linewidth = 1.1, lineend = "round") +
    geom_segment(data = lead, aes(x = id, xend = id, y = y0, yend = y1), colour = "grey85", linewidth = 0.15) +
    geom_text(data = tot, aes(x = id, y = y_lab, label = iso3, angle = ang, hjust = hj, colour = grp),
              size = iso_size, vjust = 0.5) +
    geom_label(data = data.frame(x = n - gap / 2 + 0.5, y = brks, lab = blab), aes(x = x, y = y, label = lab),
               size = 3, colour = "grey30", fill = "white", label.size = 0, label.padding = unit(0.08, "lines"), alpha = 0.9) +
    geom_text(data = gpos, aes(x = xm, y = gy, label = name_lab, hjust = hj, vjust = vj, colour = grp),
              size = 3.6, lineheight = 0.95, fontface = "bold") +
    geom_text(data = gpos, aes(x = xm, y = gy, label = line_lab, hjust = hj, vjust = vj),
              size = 3.6, lineheight = 0.95, colour = "black") +
    scale_fill_manual(values = fill_vals, name = NULL,
                      guide = if (show_fill_legend) guide_legend(ncol = 3, byrow = TRUE) else "none") +
    scale_colour_manual(values = darken(grp_cols), guide = "none") +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_lo, y_hi), oob = scales::oob_keep) +
    coord_polar(start = 0, clip = "off") +
    labs(title = title) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14, margin = margin(b = 14)),
          plot.margin = margin(30, 85, 30, 85), legend.position = "bottom")
}

# data step: simulation outputs -> csv (needs the sim=<label> folders) ----
#' circle_configs: one csv per scenario (circle_iso3_<label>.csv); tier_configs: expansion_scenarios_by_level.csv;
#' ce_configs: cost_effect_by_country.csv (costing on). NULL skips that part. Status quo arms read once.
gem_plot_data <- function(out_root, res_dir, manifest = fread(file.path(out_root, "scenario_manifest.csv")),
                          circle_configs = c("tier1to6", "tier1to6_full"), tier_configs = paste0("tier1to", 1:6),
                          ce_configs = c("tier1to6", "tier1to6_full")) {
  prim <- as.data.table(manifest)[mode == "primary" & treatment == TRUE]
  fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  cg <- gem_country_groups()
  sums <- list()
  get_sums <- function(lab, yrs) { if (is.null(sums[[lab]])) sums[[lab]] <<- arm_sums(out_root, lab, yrs); sums[[lab]] }
  tag <- function(d, m) d[, `:=`(config = m$config, scenario = m$label, base_label = m$base_label, sim_years = m$sim_years)][]
  for (cf in intersect(circle_configs, prim$config)) {
    m <- prim[config == cf]
    d <- tag(sums_to_deltas(get_sums(m$base_label, m$sim_years), get_sums(m$label, m$sim_years), cg), m)
    fwrite(d, file.path(fig_dir, paste0("circle_iso3_", m$label, ".csv")))
  }
  if (length(tier_configs) && all(tier_configs %in% prim$config)) {
    m <- prim[match(tier_configs, config)]
    if (uniqueN(m$base_label) != 1) stop("tier scenarios must share one status quo arm")
    dd <- rbindlist(lapply(seq_len(nrow(m)), function(k)
      tag(sums_to_deltas(get_sums(m$base_label[k], m$sim_years[k]), get_sums(m$label[k], m$sim_years[k]), cg), m[k])[, step := k]))
    fwrite(dd, file.path(fig_dir, "expansion_scenarios_by_level.csv"))
  } else if (length(tier_configs)) message("gem_plot_data: not all tier scenarios in manifest; expansion data skipped")
  ce <- rbindlist(lapply(intersect(ce_configs, prim$config), function(cf) {
    m <- prim[config == cf]
    x <- ce_from_sums(get_sums(m$base_label, m$sim_years), get_sums(m$label, m$sim_years), cg)
    if (!is.null(x)) tag(x, m)
  }))
  if (nrow(ce)) fwrite(ce, file.path(fig_dir, "cost_effect_by_country.csv")) else
    message("gem_plot_data: no cost columns (costing off); cost-effectiveness data skipped")
  invisible(fig_dir)
}

# one scenario vs status quo (A obesity, B deaths; by WHO region) ----
#' out_root given: data built from the sims first; out_root = NULL: from res_dir/figures/circle_iso3_*.csv only
#' configs: default tier1to6 and tier1to6_full; scales shared across all configs plotted
#' mort_unit: "per_100k" or "per_1000_treated"
plot_gem_circles <- function(out_root = NULL, res_dir, configs = c("tier1to6", "tier1to6_full"),
                             mort_unit = c("per_100k", "per_1000_treated"), manifest = NULL) {
  mort_unit <- match.arg(mort_unit)
  fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(out_root)) {
    if (is.null(manifest)) manifest <- fread(file.path(out_root, "scenario_manifest.csv"))
    prim <- as.data.table(manifest)[mode == "primary" & treatment == TRUE]
    if (!nrow(prim)) { message("plot_gem_circles: no primary intervention scenarios"); return(invisible(NULL)) }
    if (!any(configs %in% prim$config)) configs <- prim$config[1]
    gem_plot_data(out_root, res_dir, manifest, circle_configs = configs, tier_configs = NULL, ce_configs = NULL)
  }
  fs <- list.files(fig_dir, pattern = "^circle_iso3_.*\\.csv$", full.names = TRUE)
  all_d <- rbindlist(lapply(fs, fread), fill = TRUE)
  if (!nrow(all_d) || !"config" %in% names(all_d)) stop("no circle_iso3_*.csv with a config column in ", fig_dir,
                                                         " (run gem_plot_data() where the sims are)")
  configs <- configs[configs %in% all_d$config]
  if (!length(configs)) stop("none of the configs in ", fig_dir)
  dd <- lapply(configs, function(cf) all_d[config == cf]); names(dd) <- configs
  mo <- measure_info("obesity", dd[[1]]$sim_years[1])
  md <- measure_info("deaths", dd[[1]]$sim_years[1])
  if (mort_unit == "per_1000_treated") { md$col <- "deaths_averted_per_1000_treated"; md$kind <- "per1000" }
  mlab <- if (mort_unit == "per_100k") "per 100,000 adults" else "per 1,000 persons started on treatment"
  all_c <- rbindlist(dd)[level == "iso3" & !is.na(who_region) & who_region != "no who region"]
  lim_o <- range(c(0, all_c[[mo$col]]), na.rm = TRUE)             # shared across scenarios
  lim_m <- range(c(0, all_c[[md$col]]), na.rm = TRUE)

  for (cf in configs) {
    d <- dd[[cf]]; m <- d[1, .(label = scenario, base_label, sim_years)]
    cty <- d[level == "iso3" & !is.na(who_region) & who_region != "no who region"]
    reg <- d[level == "who_region" & group != "no who region"]
    regs <- sort(unique(cty$who_region)); rc <- region_palette(regs)
    g <- d[level == "global"]
    reg[, `:=`(sq_d = 1e5 * deaths_base / adults, int_d = 1e5 * (deaths_base - deaths_averted) / adults)]
    mk <- function(mi, unit, lim, title) {
      dl <- cty[!is.na(get(mi$col)), .(iso3 = group, grp = who_region, fill_key = who_region, v = get(mi$col))]
      k <- if (mi$kind == "pp") c("obesity_base", "obesity_int") else c("sq_d", "int_d")
      gl <- reg[, .(grp = group, name = grp_name(group),
                    lines = grp_lines(get(mi$col), get(mi$rel), adults, mi$kind, get(k[1]), get(k[2])))]
      circle_panel(dl, gl, lim, unit, paste0(title, "\nGlobal: ", global_txt(g, mi)), rc, rc)
    }
    p1 <- mk(mo, "pp", lim_o, paste0("A. Reduction in obesity prevalence at year ", m$sim_years))
    p2 <- mk(md, "", lim_m, paste0("B. Deaths averted over ", m$sim_years, " years, ", mlab))
    p <- (p1 | p2) + plot_annotation(
      title = paste0(cf, " vs status quo (", m$label, " vs ", m$base_label, ")"),
      caption = paste("Bars: countries, ordered by WHO region and effect; outward = benefit. Around each region: adults at year 0;",
                      "obesity prevalence at year", m$sim_years, "or deaths per 100,000 adults over", m$sim_years, "years, status quo vs intervention;",
                      "effect; change as a share of the status quo burden. Scales are the same for all scenarios plotted.",
                      if (nrow(cty) < 20) "Test run: few countries and small cohorts; results are noisy." else ""),
      theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5)))
    fp <- file.path(fig_dir, paste0("circle_iso3_", m$label, ".png"))
    ggsave(fp, p, width = 24, height = 12, dpi = 300, bg = "white")
    message("circle plot: ", fp)
  }
  invisible(fig_dir)
}

# tiers 1..k expansion: stacked tier steps by country ----
#' configs: tier scenarios in order (tiers 1, 1-2, ..., 1-6), same status quo arm
#' measures: "deaths", "dalys", "obesity" -> one two-panel figure each (A WHO region, B income group)
#' obesity_deaths: also two figures, A obesity + B deaths, one by WHO region and one by income group
#' step k = scenario k - scenario k-1; countries without a WB income group are not shown in income panels
tier_key <- c("Tier 1: T2D, BMI \u226530",
              "+ Tier 2: T2D, BMI 25\u201329.9 with CVD or CKD",
              "+ Tier 3: no T2D, BMI \u226530 with CVD or CKD",
              "+ Tier 4: no T2D, BMI \u226540",
              "+ Tier 5: no T2D, BMI 35\u201339.9",
              "+ Tier 6: no T2D, BMI 30\u201334.9")
tier_cols <- c("#0B1F3A", "#1F3B66", "#355C8C", "#5B82B0", "#8FAED0", "#C6D7EA")   # dark navy -> light

#' out_root given: data built from the sims first; out_root = NULL: from res_dir/figures/expansion_scenarios_by_level.csv
expansion_scenarios_circleplot <- function(out_root = NULL, res_dir, configs = paste0("tier1to", 1:6),
                                           measures = c("deaths", "dalys", "obesity"), obesity_deaths = TRUE,
                                           manifest = NULL) {
  fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(out_root)) {
    if (is.null(manifest)) manifest <- fread(file.path(out_root, "scenario_manifest.csv"))
    gem_plot_data(out_root, res_dir, manifest, circle_configs = NULL, tier_configs = configs, ce_configs = NULL)
  }
  fp0 <- file.path(fig_dir, "expansion_scenarios_by_level.csv")
  if (!file.exists(fp0)) stop("missing ", fp0, " (run gem_plot_data() where the sims are)")
  dd <- fread(fp0)
  m <- data.table(step = sort(unique(dd$step)))
  m[, sim_years := if ("sim_years" %in% names(dd)) dd$sim_years[1] else 5L]
  last <- dd[step == max(step)]
  keys <- tier_key[seq_len(nrow(m))]; fills <- setNames(tier_cols[seq_len(nrow(m))], keys)
  io <- c("high income", "upper middle income", "lower middle income", "low income")      # not classified left out

  # step increments by country + scale shared by all panels of one measure
  prep <- function(ms) {
    mi <- measure_info(ms, m$sim_years[1])
    inc <- dd[level == "iso3", .(iso3 = group, who_region, wb_income, step, val = get(mi$col))][order(iso3, step)]
    inc[, v := val - shift(val, fill = 0), by = iso3]
    inc[, fill_key := factor(keys[step], levels = keys)]
    pos <- inc[, .(p = sum(pmax(v, 0)), n = sum(pmin(v, 0))), by = iso3]
    fwrite(inc, file.path(fig_dir, paste0("expansion_", ms, "_by_country_step.csv")))   # plotted values
    list(mi = mi, inc = inc, lim = c(min(0, pos$n), max(pos$p)))
  }
  panel <- function(x, by, title, panel_in = 9) {
    mi <- x$mi
    if (by == "who_region") {
      d <- x$inc[!is.na(who_region) & who_region != "no who region", .(iso3, grp = who_region, fill_key, v)]
    } else {
      d <- x$inc[wb_income %in% io, .(iso3, grp = factor(wb_income, levels = intersect(io, unique(wb_income))), fill_key, v)]
    }
    gl <- last[level == by & group %in% unique(as.character(d$grp)),
               .(grp = group, name = grp_name(group), lines = grp_lines(get(mi$col), get(mi$rel), adults, mi$kind))]
    lv <- if (is.factor(d$grp)) levels(d$grp) else sort(unique(d$grp))
    circle_panel(d, gl, x$lim, mi$unit, title, setNames(rep("grey45", length(lv)), lv), fills, show_fill_legend = TRUE,
                 panel_in = panel_in)
  }
  note <- function(nc) paste0("Each bar: one country, tiers 1\u2013", nrow(m), " total; segments: added effect of each tier step. ",
                              "Around each group: adults at year 0; tiers 1\u2013", nrow(m), " effect; change as a share of the status quo burden. ",
                              "Countries without a World Bank income group are not shown in income panels.",
                              if (nc < 20) " Test run: few countries; results are noisy." else "")
  g <- last[level == "global"]
  fin <- function(p, title, subtitle, nc) p + plot_layout(guides = "collect") +
    plot_annotation(title = title, subtitle = subtitle, caption = note(nc),
                    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                                  plot.subtitle = element_text(hjust = 0.5, colour = "grey25", size = 12))) &
    theme(legend.position = "bottom", legend.text = element_text(size = 11))
  ms_all <- unique(c(measures, if (obesity_deaths) c("obesity", "deaths")))
  px <- lapply(setNames(ms_all, ms_all), prep)
  nc <- uniqueN(px[[1]]$inc$iso3)

  for (ms in measures) {
    x <- px[[ms]]
    p <- fin(panel(x, "who_region", "A. By WHO region") | panel(x, "wb_income", "B. By World Bank income group"),
             x$mi$title, paste0("Coverage expanded by priority tier (tiers 1 to ", nrow(m), ") vs status quo. Global, tiers 1\u2013",
                                nrow(m), ": ", global_txt(g, x$mi)), nc)
    fp <- file.path(fig_dir, paste0("expansion_", ms, ".png"))
    ggsave(fp, p, width = 22, height = 12.5, dpi = 300, bg = "white")
    message("expansion circle plot: ", fp)
  }
  if (obesity_deaths) {
    o <- px[["obesity"]]; d <- px[["deaths"]]
    for (by in c("who_region", "wb_income")) {
      gn <- if (by == "who_region") "WHO region" else "World Bank income group"
      p <- fin(panel(o, by, "A. Reduction in obesity prevalence (pp)") | panel(d, by, "B. Deaths averted per 100,000 adults"),
               paste0("Obesity and deaths by priority-tier expansion (tiers 1 to ", nrow(m), ") vs status quo, by ", gn, ", ",
                      m$sim_years[1], " years"),
               paste0("Global, tiers 1\u2013", nrow(m), ": obesity ", global_txt(g, o$mi), "; deaths ", global_txt(g, d$mi)), nc)
      fp <- file.path(fig_dir, paste0("expansion_obesity_deaths_", by, ".png"))
      ggsave(fp, p, width = 22, height = 12.5, dpi = 300, bg = "white")
      message("expansion circle plot: ", fp)
    }
  }
  invisible(fig_dir)
}

# ==== cost-effectiveness plots (R version of GEM_cost_effect_plot.py) ====
#' data: res_dir/figures/cost_effect_by_country.csv (gem_plot_data(), costing on)
#' costs in Int$ (discounted 3%/yr), DALYs discounted; incremental = scenario - status quo
#' countries without a World Bank income group are left out (reported)
ce_scenario_labels <- c(id_sust_costx1 = "1. Low cost-based prices", id_sust_costx3 = "2. Reduced market prices",
                        intl_rwe = "3. Current market prices", id_lowest_mp = "Lowest market price",
                        id_generic_entry = "Generic entry price", tc_zero = "Zero treatment cost")
ce_ig_order  <- c("low income", "lower middle income", "upper middle income", "high income")
ce_ig_cols   <- c("Low income" = "#009E73", "Lower middle income" = "#AA490C",
                  "Upper middle income" = "#E69D00", "High income" = "#0072B2")        # Wong palette
ce_sc_cols   <- c(intl_rwe = "#ff7f0e", id_sust_costx1 = "#1f77b4", id_sust_costx3 = "#2ca02c",
                  tc_zero = "#9467bd", id_lowest_mp = "#8c564b", id_generic_entry = "#e377c2")

# willingness-to-pay thresholds by income group (Int$ per DALY averted) ----
#' fixed values per World Bank income group, in place of one high threshold for all countries;
#'   opportunity-cost based estimates are lower than 1-3x GDP per capita (Woods et al. 2016, Value Health 19:929-35)
#' change here or pass wtp_ig = c(...) to the plot functions
ce_wtp_ig <- c("low income" = 500, "lower middle income" = 1000, "upper middle income" = 5000, "high income" = 5000)
wtp_lab <- function(w) paste0("Int$ ", formatC(w, format = "d", big.mark = ","), "/DALY")
#' extra threshold lines for chosen price scenarios only: named numeric, e.g. c(intl_rwe = 50000)
#'   (high-income HTA convention; shown on the current market price panel); set in the Rmd
ce_wtp_extra <- c(intl_rwe = 50000)
ce_wtp_extra_lab <- "high-income convention"
# caption note: uptake assumption (tiers 1-6, primary uptake), from sheet glp1_uptake
ce_uptake_note <- function(uptake = "optimistic", uptake_fp = NULL, horizon = 5) {
  u <- uptake_txt(uptake, uptake_fp)
  paste0("Uptake assumption: all six priority tiers open; each year ", if (is.na(u)) "a fixed share" else sub(" per year", "", u),
         " of eligible adults with access start semaglutide (by indication, sex and age; 2.5 \u00d7 real-world uptake, Mailhac 2024). ",
         horizon, "-year horizon; costs and DALYs discounted 3% per year.")
}
# wrap caption lines to the figure width (about 13 characters per inch at caption size)
ce_wrap <- function(x, fig_w) paste(unlist(lapply(strsplit(x, "\n")[[1]], strwrap, width = floor(13 * fig_w))), collapse = "\n")
# caption note: countries left out by name
ce_excl_note <- function(exclude) if (length(exclude)) paste0("Left out (outlier): ", paste(exclude, collapse = ", "), ".") else NULL
# horizontal threshold segments, one per income-group column (y in 1,000 Int$)
wtp_segs <- function(wtp_ig, lv) {
  w <- data.table(ig_disp = factor(grp_name(names(wtp_ig)), levels = lv), y = unname(wtp_ig) / 1e3)[!is.na(ig_disp)]
  w[, `:=`(x0 = as.integer(ig_disp) - 0.4, x1 = as.integer(ig_disp) + 0.4, lab = wtp_lab(y * 1e3))][]
}

ce_read <- function(res_dir, config, scenarios, exclude = NULL) {
  fp <- file.path(res_dir, "figures", "cost_effect_by_country.csv")
  if (!file.exists(fp)) stop("missing ", fp, " (costing on + gem_plot_data() where the sims are)")
  cfg <- config
  d <- fread(fp)[config == cfg & cost_scenario %in% scenarios & !iso3 %in% exclude]
  d[, income_group := wb_norm_ce(income_group)]
  drop <- d[!income_group %in% ce_ig_order, unique(iso3)]
  if (length(drop)) message("CE plots: ", length(drop), " countries without a World Bank income group left out: ",
                            paste(sort(drop), collapse = ", "))
  d <- d[income_group %in% ce_ig_order]
  d[, `:=`(ig_disp = factor(grp_name(income_group), levels = grp_name(ce_ig_order)),
           sc_disp = factor(ce_scenario_labels[cost_scenario], levels = ce_scenario_labels[scenarios]))][]
}
wb_norm_ce <- function(x) {
  canon <- c(ce_ig_order, "no wb income group")
  out <- canon[match(gsub("[^a-z]", "", tolower(x)), gsub(" ", "", canon))]
  ifelse(is.na(out), x, out)
}
ce_theme <- function() theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom", legend.box = "vertical")

# 1. CE plane: incremental cost vs DALYs averted by country, one panel per price scenario ----
#' WTP lines: income-group thresholds (ce_wtp_ig); share_y = FALSE: own y-axis per panel (market prices >> cost-based)
plot_ce_plane <- function(res_dir, config = "tier1to6", scenarios = c("id_sust_costx1", "id_sust_costx3", "intl_rwe"),
                          share_y = FALSE, wtp_ig = ce_wtp_ig, wtp_extra = NULL, exclude = NULL, caption_extra = NULL,
                          out_file = paste0("ce_plane_", config, ".png")) {
  d <- ce_read(res_dir, config, scenarios, exclude)
  uw <- sort(unique(wtp_ig), decreasing = TRUE)
  wl <- sapply(uw, function(v) paste0(wtp_lab(v), " (", paste(grp_name(names(wtp_ig)[wtp_ig == v]), collapse = ", "), ")"))
  wtp <- CJ(sc_disp = factor(levels(d$sc_disp), levels = levels(d$sc_disp)), wtp = uw)
  wtp[, wtp_lab := wl[match(wtp, uw)]]
  ex <- wtp_extra[names(wtp_extra) %in% scenarios]                   # extra lines: listed scenario panels only
  if (length(ex)) {
    el <- paste0(wtp_lab(ex), " (", ce_wtp_extra_lab, ")")
    wtp <- rbind(wtp, data.table(sc_disp = factor(ce_scenario_labels[names(ex)], levels = levels(d$sc_disp)), wtp = unname(ex), wtp_lab = el))
    wl <- c(unique(el), wl)
  }
  wtp[, wtp_lab := factor(wtp_lab, levels = unique(wl))]
  p <- ggplot(d, aes(x = daly_averted_val / 1e6, y = incremental_cost_val / 1e9)) +
    geom_hline(yintercept = 0, linewidth = 0.4) + geom_vline(xintercept = 0, linewidth = 0.4) +
    geom_abline(data = wtp, aes(slope = wtp / 1000, intercept = 0, linetype = wtp_lab), colour = "grey40", linewidth = 0.5) +
    geom_point(aes(size = n, colour = ig_disp), alpha = 0.85) +
    ggrepel::geom_text_repel(aes(label = iso3), size = 2.4, colour = "grey25", segment.colour = "grey60",
                             segment.size = 0.3, max.overlaps = 12, min.segment.length = 0.2, seed = 42) +
    facet_wrap(~ sc_disp, nrow = 1, scales = if (share_y) "fixed" else "free_y") +
    scale_colour_manual(values = ce_ig_cols, name = "Income group", drop = TRUE) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "longdash", "dotdash")[seq_along(unique(wl)) + !length(ex)],
                          name = "Willingness to pay") +
    guides(linetype = guide_legend(nrow = 2)) +
    scale_size_area(max_size = 7, guide = "none") +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "DALYs averted (millions)", y = "Incremental cost (billion Int$)",
         title = "Cost-effectiveness plane by drug price scenario",
         caption = paste(c(paste("Dots: countries (size = adults). Lines: cost per DALY averted equal to the willingness-to-pay threshold.",
                                 if (!share_y) "Y-axis differs by panel." else "", ce_excl_note(exclude)), caption_extra), collapse = "\n")) +
    ce_theme()
  fp <- file.path(res_dir, "figures", out_file); fw <- 3.6 * length(scenarios) + 1.5
  p$labels$caption <- ce_wrap(p$labels$caption, fw)
  ggsave(fp, p, width = fw, height = 7, dpi = 300, bg = "white")
  message("CE plane: ", fp); invisible(p)
}

# 2. country ICERs by income group, one panel per price scenario ----
#' ICER in 1,000 Int$ per DALY averted; below icer_k_floor and above outlier_thresh left out (reported)
plot_country_icer <- function(res_dir, config = "tier1to6", scenarios = c("id_sust_costx1", "id_sust_costx3", "intl_rwe"),
                              icer_k_floor = -20, outlier_thresh = 200, label_icer_k = 25, label_n = 5e8,
                              wtp_ig = ce_wtp_ig, wtp_extra = NULL, exclude = NULL, caption_extra = NULL,
                              out_file = paste0("ce_country_icer_", config, ".png")) {
  d <- ce_read(res_dir, config, scenarios, exclude)[!is.na(icer_val)]
  d[, icer_k := icer_val / 1e3]
  out <- d[icer_k > outlier_thresh | icer_k <= icer_k_floor]
  if (nrow(out)) { message("Country ICER plot: left out (ICER <= ", icer_k_floor, "k or > ", outlier_thresh, "k):")
                   print(out[, .(n = .N, max_icer_k = round(max(icer_k))), by = .(sc_disp, ig_disp)]) }
  d <- d[icer_k > icer_k_floor & icer_k <= outlier_thresh]
  set.seed(42); d[, xj := as.integer(ig_disp) + runif(.N, -0.15, 0.15)]
  ws <- wtp_segs(wtp_ig, levels(d$ig_disp))
  ex <- wtp_extra[names(wtp_extra) %in% scenarios]                   # extra lines: listed scenario panels only
  exd <- data.table(sc_disp = factor(ce_scenario_labels[names(ex)], levels = levels(d$sc_disp)), y = unname(ex) / 1e3,
                    lab = paste0(wtp_lab(ex), " (", ce_wtp_extra_lab, ")"))
  p <- ggplot(d, aes(x = xj, y = icer_k)) +
    geom_segment(data = ws, aes(x = x0, xend = x1, y = y, yend = y), inherit.aes = FALSE,
                 linetype = "dashed", colour = "grey35", linewidth = 0.5) +
    geom_hline(data = exd, aes(yintercept = y), linetype = "longdash", colour = "grey20", linewidth = 0.5) +
    geom_text(data = exd, aes(x = Inf, y = y, label = lab), inherit.aes = FALSE, hjust = 1.02, vjust = -0.5,
              size = 2.8, colour = "grey20") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_point(aes(size = n, colour = cost_scenario), alpha = 0.7) +
    ggrepel::geom_text_repel(data = d[icer_k > label_icer_k | n > label_n], aes(label = iso3),
                             size = 2.4, colour = "grey25", segment.size = 0.3, max.overlaps = 15, seed = 42) +
    facet_wrap(~ sc_disp, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq_along(levels(d$ig_disp)), labels = levels(d$ig_disp)) +
    scale_colour_manual(values = ce_sc_cols, guide = "none") +
    scale_size_area(max_size = 6, guide = "none") +
    labs(x = NULL, y = "ICER (1,000 Int$ per DALY averted)", title = "Country ICER by income group",
         caption = paste(c(paste("Dots: countries (size = adults).", if (length(scenarios) > 1) "Y-axis differs by panel." else "",
                                 ce_excl_note(exclude)),
                           paste0("Dashed line in each income group: willingness-to-pay threshold (",
                                  paste(grp_name(names(wtp_ig)), wtp_lab(wtp_ig), sep = " ", collapse = "; "), ")."),
                           caption_extra), collapse = "\n")) +
    ce_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  fp <- file.path(res_dir, "figures", out_file); fw <- max(8, 4.5 * length(scenarios))
  p$labels$caption <- ce_wrap(p$labels$caption, fw)
  ggsave(fp, p, width = fw, height = 6.3, dpi = 300, bg = "white")
  message("Country ICER plot: ", fp); invisible(p)
}

# ==== maximum cost-effective drug price (price-threshold analysis) ====
#' drug cost is linear in price (cost = monthly price x 12 x discounted time on treatment), so per country:
#'   incremental drug cost at price p = (incremental cost at p) - (incremental cost with drug free, tc_zero)
#'   max cost-effective price p* : condition cost change + p* x treated time = WTP x DALYs averted
#'   r* = p* / current market price = (WTP x DALYs averted - incremental condition cost) / incremental drug cost (intl_rwe)
#'   r_x3, r_x1 = reduced-market and cost-based price / current market price (= ratio of their incremental drug costs)
#' no price table needed: all ratios come from cost_effect_by_country.csv (tc_zero, intl_rwe, id_sust_costx3, id_sust_costx1)
#' r* <= 0: not cost-effective at any price (no net DALYs averted, or condition costs rise more than WTP x DALYs)
#' small-outlier rule (same for every figure): a country is left out when it is small (adults below the 25th
#'   percentile of countries, small_q) AND an outlier (r* <= 0, or log10 r* outside the quartiles +/- 1.5 x IQR, iqr_k)
#'   small countries have small cohorts -> noisy DALYs averted; left-out countries are listed (message + csv)
ce_price_lv <- c("Current market price", "Reduced market price", "Low cost-based price", "None of the 3 prices")
ce_price_cols <- setNames(c("#0B1F3A", "#5B82B0", "#E69D00", "#B2182B"), ce_price_lv)

ce_price_data <- function(res_dir, configs = c("tier1to6", "tier1to6_full"), wtp_ig = ce_wtp_ig,
                          small_q = 0.25, iqr_k = 1.5, exclude = NULL) {
  fp <- file.path(res_dir, "figures", "cost_effect_by_country.csv")
  if (!file.exists(fp)) stop("missing ", fp, " (costing on + gem_plot_data() where the sims are)")
  need <- c("tc_zero", "intl_rwe", "id_sust_costx3", "id_sust_costx1")
  d <- fread(fp)[config %in% configs & cost_scenario %in% need & !iso3 %in% exclude]
  if (!all(need %in% d$cost_scenario)) stop("cost_effect_by_country.csv needs price scenarios: ", paste(need, collapse = ", "))
  d[, income_group := wb_norm_ce(income_group)]
  w <- dcast(d[income_group %in% ce_ig_order], config + iso3 + who_region + income_group + n + daly_averted_val ~ cost_scenario,
             value.var = "incremental_cost_val")
  w[, `:=`(wtp = unname(wtp_ig[income_group]), drug = intl_rwe - tc_zero)]
  w <- w[drug > 0]                                                   # no added treatment: no price question
  w[, r_star := fifelse(daly_averted_val > 0, (wtp * daly_averted_val - tc_zero) / drug, 0)]
  w[, `:=`(r_x3 = (id_sust_costx3 - tc_zero) / drug, r_x1 = (id_sust_costx1 - tc_zero) / drug)]
  w[, ce_at := factor(fcase(r_star >= 1, ce_price_lv[1], r_star >= r_x3, ce_price_lv[2], r_star >= r_x1, ce_price_lv[3],
                            default = ce_price_lv[4]), levels = ce_price_lv)]
  # small-outlier rule, per config
  w[, small := n < quantile(n, small_q), by = config]
  w[, lr := fifelse(r_star > 0, log10(r_star), NA_real_)]
  w[, c("q1", "q3") := as.list(quantile(lr, c(0.25, 0.75), na.rm = TRUE)), by = config]
  w[, outlier := r_star <= 0 | lr < q1 - iqr_k * (q3 - q1) | lr > q3 + iqr_k * (q3 - q1)]
  w[, omit := small & outlier]
  w[, `:=`(ig_disp = factor(grp_name(income_group), levels = grp_name(ce_ig_order)),
           cfg_disp = factor(fifelse(config == "tier1to6_full", "Tiers 1\u20136, full uptake", "Tiers 1\u20136, optimistic uptake"),
                             levels = c("Tiers 1\u20136, optimistic uptake", "Tiers 1\u20136, full uptake")))]
  om <- w[omit == TRUE]
  if (nrow(om)) {
    message("Price plots: left out ", uniqueN(om$iso3), " small outlier countries (adults < ", 100 * small_q,
            "th percentile and r* <= 0 or outside quartiles +/- ", iqr_k, " IQR): ", paste(sort(unique(om$iso3)), collapse = ", "))
    fwrite(om[, .(config, iso3, income_group, adults = n, r_star, daly_averted_val)],
           file.path(res_dir, "figures", "ce_price_left_out.csv"))
  }
  w[, c("lr", "q1", "q3") := NULL][]
}

# 1. maximum cost-effective price by country (% of current market price), by income group ----
#' dot = country (size = adults); colour = lowest of the 3 modelled prices needed for cost-effectiveness
#' short lines in each income group: median reduced-market and cost-based price (% of market); 100% = current market
plot_ce_max_price <- function(res_dir, configs = c("tier1to6", "tier1to6_full"), wtp_ig = ce_wtp_ig, label_top = 12,
                              out_file = "ce_max_price_by_country.png", caption_extra = NULL, exclude = NULL, ...) {
  w <- ce_price_data(res_dir, configs, wtp_ig, exclude = exclude, ...)
  one <- uniqueN(w$cfg_disp) == 1                                    # one config: no panel title, no config name
  w <- w[omit == FALSE]
  lo <- 0.3                                                          # floor (%): r* <= 0 shown here as "none"
  w[, y := pmax(100 * r_star, lo)]
  set.seed(42); w[, xj := as.integer(ig_disp) + runif(.N, -0.22, 0.22)]
  ref <- w[, .(y = 100 * c(median(r_x3), median(r_x1)), what = c("Reduced market", "Low cost-based")), by = .(cfg_disp, ig_disp)]
  ref[, `:=`(x0 = as.integer(ig_disp) - 0.38, x1 = as.integer(ig_disp) + 0.38)]
  lab <- w[order(-n)][, head(.SD, label_top), by = cfg_disp]
  n_ce <- w[, .(txt = sprintf("%s%d of %d countries cost-effective at current market price; %d at reduced; %d at cost-based",
                              if (one) "" else paste0(cfg_disp[1], ": "), sum(r_star >= 1), .N, sum(r_star >= r_x3), sum(r_star >= r_x1))),
            by = cfg_disp]
  lv <- levels(w$ig_disp)
  p <- ggplot(w, aes(x = xj, y = y)) +
    geom_hline(yintercept = 100, colour = "grey20", linewidth = 0.5) +
    annotate("text", x = 0.5, y = 100, label = "Current market price", hjust = 0, vjust = -0.5, size = 3, colour = "grey20") +
    geom_segment(data = ref, aes(x = x0, xend = x1, y = y, yend = y, linetype = what), inherit.aes = FALSE,
                 colour = "grey35", linewidth = 0.5) +
    geom_point(aes(size = n, colour = ce_at), alpha = 0.8) +
    ggrepel::geom_text_repel(data = lab, aes(label = iso3), size = 2.6, colour = "grey15", segment.size = 0.3,
                             min.segment.length = 0.1, seed = 42) +
    facet_wrap(~ cfg_disp, nrow = 1) +
    scale_y_log10(breaks = c(lo, 1, 3, 10, 30, 100, 300, 1000), labels = c("\u2264 0", "1%", "3%", "10%", "30%", "100%", "300%", "1,000%")) +
    scale_x_continuous(breaks = seq_along(lv), labels = lv, limits = c(0.45, length(lv) + 0.55)) +
    scale_colour_manual(values = ce_price_cols, drop = FALSE, name = "Cost-effective at") +
    scale_linetype_manual(values = c("Reduced market" = "dashed", "Low cost-based" = "dotted"), name = "Median modelled price") +
    scale_size_area(max_size = 9, guide = "none") +
    guides(colour = guide_legend(override.aes = list(size = 4), order = 1)) +
    labs(x = NULL, y = "Maximum cost-effective price (% of current market price, log scale)",
         title = "Maximum price at which semaglutide is cost-effective, by country",
         subtitle = paste(n_ce$txt, collapse = "\n"),
         caption = paste0("Willingness to pay per DALY averted: ", paste(grp_name(names(wtp_ig)), wtp_lab(wtp_ig), sep = " ", collapse = "; "),
                          ".\nDots: countries (size = adults; largest labelled). \u2264 0: not cost-effective at any price. ",
                          "Small outlier countries left out (adults < 25th percentile and outlying maximum price); see ce_price_left_out.csv.",
                          if (length(exclude)) paste0(" ", ce_excl_note(exclude)) else "",
                          if (length(caption_extra)) paste0("\n", paste(caption_extra, collapse = "\n")) else "")) +
    ce_theme() + theme(axis.text.x = element_text(angle = 20, hjust = 1), plot.subtitle = element_text(hjust = 0.5, size = 10))+
    (if (one) theme(strip.text = element_blank()) else NULL)
  fp <- file.path(res_dir, "figures", out_file); fw <- if (one) 9 else 6 * length(unique(w$cfg_disp)) + 1
  p$labels$caption <- ce_wrap(p$labels$caption, fw)
  ggsave(fp, p, width = fw, height = 7.4, dpi = 300, bg = "white")
  message("Max price plot: ", fp); invisible(p)
}

# 2. share of adults and of countries where cost-effective, by drug price ----
#' curve: share with r* >= price (price as % of current market price); vertical lines: current market and the
#'   adult-weighted median reduced-market and cost-based prices
plot_ce_price_curve <- function(res_dir, configs = c("tier1to6", "tier1to6_full"), wtp_ig = ce_wtp_ig,
                                out_file = "ce_price_curve.png", exclude = NULL, ...) {
  w <- ce_price_data(res_dir, configs, wtp_ig, exclude = exclude, ...)[omit == FALSE]
  x <- 10^seq(-1, 3, length.out = 400)                                # 0.1% to 1,000% of market price
  cv <- w[, .(x = x, `Share of adults (%)` = sapply(x, function(v) 100 * sum(n[100 * r_star >= v]) / sum(n)),
              `Share of countries (%)` = sapply(x, function(v) 100 * mean(100 * r_star >= v))), by = cfg_disp]
  cv <- melt(cv, id.vars = c("cfg_disp", "x"), variable.name = "measure", value.name = "share")
  wmed <- function(v, wt) { o <- order(v); v[o][which(cumsum(wt[o]) >= sum(wt) / 2)[1]] }
  ref <- data.table(x = c(100, 100 * wmed(w$r_x3, w$n), 100 * wmed(w$r_x1, w$n)),
                    what = c("Current market", "Reduced market (median)", "Low cost-based (median)"))
  ref[, what := factor(what, levels = what)]
  p <- ggplot(cv, aes(x = x, y = share, colour = cfg_disp)) +
    geom_vline(data = ref, aes(xintercept = x, linetype = what), colour = "grey35", linewidth = 0.5) +
    geom_step(linewidth = 0.9, direction = "vh") +
    facet_wrap(~ measure, nrow = 1) +
    scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000), labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%", "100%", "300%", "1,000%")) +
    scale_y_continuous(breaks = seq(0, 100, 20)) + coord_cartesian(ylim = c(0, 100)) +
    scale_colour_manual(values = c("#5B82B0", "#0B1F3A"), name = NULL) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Modelled price") +
    labs(x = "Semaglutide price (% of current market price, log scale)", y = "Cost-effective (%)",
         title = "Where semaglutide is cost-effective, by drug price",
         caption = paste0("Cost-effective: maximum cost-effective price \u2265 price.\nWillingness to pay per DALY averted: ",
                          paste(grp_name(names(wtp_ig)), wtp_lab(wtp_ig), sep = " ", collapse = "; "),
                          ".\nSmall outlier countries left out (adults < 25th percentile and outlying maximum price); see ce_price_left_out.csv.")) +
    ce_theme() + theme(legend.box = "horizontal")
  fp <- file.path(res_dir, "figures", out_file)
  ggsave(fp, p, width = 12, height = 5.8, dpi = 300, bg = "white")
  message("Price curve plot: ", fp); invisible(p)
}

# ==== group-level bar charts (R version of the GEM_bar_chart.py layout) ====
#' data (03_Out only): all_scenarios_outcomes.csv (deaths, new T2D/CVD/CKD, DALYs; base and scenario totals by level)
#'   + figures/expansion_scenarios_by_level.csv and circle_iso3_<full>.csv (obesity prevalence at year horizon)
#' rates: 5-year cumulative incidence = events / adults at year 0 (%); obesity = prevalence among living (%);
#'   DALYs per 100,000 adults. Primary scenarios only (optimistic tiers 1..6, full uptake tiers 1-6; RCT effects)
#' v1 plot_bar_overall(): status quo vs tiers 1-6 optimistic vs tiers 1-6 full, per group (global first);
#'   status quo rate printed in the status quo bar
#' scenario names show the uptake used: optimistic = range of annual uptake probabilities in sheet glp1_uptake
#'   (GEM_input_parameters.xlsx; uptake_fp), full = 100%
#' v2 plot_bar_expansion(): reduction vs status quo, stacked by tier step (same navy scale as the expansion circle plots)
#' by = "who_region" or "wb_income" (groups without a WHO region / WB income group left out)
bar_outcome_info <- data.table(
  outcome = c("ac_death", "t2d_new", "cvd_new", "ckd_new", "eskd_new", "obesity", "daly_total_adj"),
  title   = c("All-cause mortality", "New T2D", "New CVD", "New CKD", "New ESKD", "Obesity", "DALYs"),
  axis    = c("5-year cumulative incidence (% of adults)", "5-year cumulative incidence (% of adults)",
              "5-year cumulative incidence (% of adults)", "5-year cumulative incidence (% of adults)",
              "5-year cumulative incidence (% of adults)", "Prevalence at year 5 (% of living adults)",
              "DALYs per 100,000 adults (discounted)"),
  unit    = c("pp", "pp", "pp", "pp", "pp", "pp", "per100k"))
bar_cols <- c("Status quo" = "#9AA3AD", "Tiers 1\u20136, optimistic uptake" = "#355C8C", "Tiers 1\u20136, full uptake" = "#0B1F3A")

# uptake text from the input file: "5\u201320% per year" (range over indication, sex, age); NA if not found
uptake_txt <- function(uptake = "optimistic", fp = NULL) {
  if (is.null(fp)) fp <- if (exists("input_param_fp")) input_param_fp else
    tryCatch(file.path(get_paths()$data, "GEM_input_parameters.xlsx"), error = function(e) NA_character_)
  u <- tryCatch(as.data.table(readxl::read_excel(fp, sheet = "glp1_uptake")), error = function(e) NULL)
  if (is.null(u)) return(NA_character_)
  v <- 100 * as.numeric(u[in_use == "yes" & scenario == uptake, value])
  if (!length(v)) return(NA_character_)
  r <- range(v); f <- function(x) as.character(signif(x, 2))
  if (diff(r) < 1e-9) paste0(f(r[1]), "% per year") else paste0(f(r[1]), "\u2013", f(r[2]), "% per year")
}
# scenario names and short label prefixes with the uptake used
bar_scen_names <- function(uptake_fp = NULL) {
  u <- uptake_txt("optimistic", uptake_fp)
  opt <- if (is.na(u)) "optimistic uptake" else paste("uptake", u)
  list(names = c("Status quo", paste0("Tiers 1\u20136, ", opt), "Tiers 1\u20136, full uptake (100%)"),
       short = c(if (is.na(u)) "Optimistic" else sub(" per year", "/yr", u), "100%"), opt = opt)
}
# x-axis title under each panel (free scales differ by outcome): strip = outcome + unit
bar_strip <- function(title, axis) paste0(title, "\n", axis)

bar_data <- function(res_dir, by = c("who_region", "wb_income")) {
  by <- match.arg(by)
  o <- fread(file.path(res_dir, "all_scenarios_outcomes.csv"))
  o <- o[level %in% c("global", by)]
  if (by == "wb_income") o[level == "wb_income", group := wb_norm_ce(group)]
  key <- function(desc) fcase(desc == "Tiers 1-6, full uptake", "full",
                              grepl("^Tiers 1(-[2-6])?, optimistic uptake$", desc), "opt", default = NA_character_)
  o[, sc := key(description)]
  o <- o[!is.na(sc)][, step := fifelse(sc == "full", 99L, as.integer(n_tiers))]
  o[, scale := fifelse(variable == "daly_total_adj", 1e5, 100)]
  r <- rbind(o[step == 1L, .(level, group, outcome = variable, step = 0L, value = scale * base / adults_y0)],
             o[, .(level, group, outcome = variable, step, value = scale * int / adults_y0)])
  # obesity prevalence from the plot csv files
  fd <- file.path(res_dir, "figures")
  ex <- if (file.exists(file.path(fd, "expansion_scenarios_by_level.csv"))) fread(file.path(fd, "expansion_scenarios_by_level.csv")) else NULL
  fu <- rbindlist(lapply(list.files(fd, "^circle_iso3_.*\\.csv$", full.names = TRUE), fread), fill = TRUE)
  fu <- if (nrow(fu) && "config" %in% names(fu)) fu[config == "tier1to6_full"] else NULL
  ob <- rbind(if (!is.null(ex)) rbind(ex[step == 1L, .(level, group, step = 0L, value = obesity_base)],
                                      ex[, .(level, group, step, value = obesity_int)]),
              if (!is.null(fu) && nrow(fu)) fu[, .(level, group, step = 99L, value = obesity_int)])
  if (!is.null(ob)) r <- rbind(r, ob[level %in% c("global", by)][, outcome := "obesity"], use.names = TRUE)
  if (by == "wb_income") r[level == "wb_income", group := wb_norm_ce(group)]
  ok <- if (by == "who_region") r$group != "no who region" else r$group %in% c("global", ce_ig_order)
  r <- unique(r[ok], by = c("group", "outcome", "step"))
  lv <- if (by == "who_region") c("global", sort(setdiff(unique(r$group), "global"))) else c("global", rev(ce_ig_order))
  r[, group := factor(group, levels = rev(intersect(lv, unique(group))))][]  # rev: global at the top
}
bar_group_cols <- function(groups, by) {
  g <- setdiff(as.character(groups), "global")
  cols <- if (by == "who_region") darken(region_palette(g)) else setNames(rep("grey30", length(g)), g)
  c(cols, global = "black")
}
# colored group labels: element_text with one color per label (ggplot2 notes this is not official; note muffled)
bar_save <- function(fp, p, ...) withCallingHandlers(ggsave(fp, p, ...), warning = function(w)
  if (grepl("Vectorized input", conditionMessage(w))) invokeRestart("muffleWarning"))
lab_text <- function(cols) withCallingHandlers(element_text(face = "bold", colour = cols), warning = function(w)
  if (grepl("Vectorized input", conditionMessage(w))) invokeRestart("muffleWarning"))
bar_fmt <- function(x, unit) ifelse(unit == "pp", ifelse(abs(x) < 1, sprintf("%.2f pp", x), sprintf("%.1f pp", x)),
                                    formatC(round(x), format = "d", big.mark = ","))

# v1: status quo vs tiers 1-6 optimistic vs full ----
plot_bar_overall <- function(res_dir, by = c("who_region", "wb_income"),
                             outcomes = c("ac_death", "t2d_new", "cvd_new", "ckd_new", "obesity", "daly_total_adj"),
                             out_file = NULL, uptake_fp = NULL) {
  by <- match.arg(by)
  sn <- bar_scen_names(uptake_fp); cols <- setNames(unname(bar_cols), sn$names)
  r <- bar_data(res_dir, by)[outcome %in% outcomes & step %in% c(0L, 6L, 99L)]
  r[, scen := factor(fcase(step == 0L, sn$names[1], step == 6L, sn$names[2], step == 99L, sn$names[3]),
                     levels = rev(sn$names))]                         # reversed: status quo bar on top within a group
  oi <- bar_outcome_info[outcome %in% outcomes]
  r <- oi[r, on = "outcome"][, title := factor(bar_strip(title, axis), levels = bar_strip(oi$title, oi$axis))]
  # status quo rate inside its bar (other bars: empty label, same dodge)
  r[, sq_lab := fifelse(step == 0L, fifelse(unit == "pp", ifelse(value < 1, sprintf("%.2f%%", value), sprintf("%.1f%%", value)),
                                            formatC(round(value), format = "d", big.mark = ",")), "")]
  # labels: change vs status quo for each scenario (absolute; relative)
  w <- dcast(r, title + unit + group ~ step, value.var = "value")
  setnames(w, c("0", "6", "99"), c("sq", "opt", "full"), skip_absent = TRUE)
  for (v in c("opt", "full")) if (!v %in% names(w)) w[, (v) := NA_real_]
  w[, lab := paste0(ifelse(is.na(opt), "", paste0(sn$short[1], ": ", bar_fmt(opt - sq, unit), " (", sprintf("%+.1f%%", 100 * (opt - sq) / sq), ")")),
                    ifelse(is.na(full), "", paste0("\n", sn$short[2], ": ", bar_fmt(full - sq, unit), " (", sprintf("%+.1f%%", 100 * (full - sq) / sq), ")")))]
  w[, x := pmax(sq, opt, full, na.rm = TRUE)]
  gc <- bar_group_cols(levels(r$group), by)
  p <- ggplot(r, aes(y = group, x = value, fill = scen)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(label = sq_lab, group = scen), position = position_dodge(width = 0.8), hjust = 1.08, size = 2.4,
              colour = "black") +
    geom_text(data = w, aes(y = group, x = x, label = lab), inherit.aes = FALSE, hjust = -0.05, size = 2.5,
              lineheight = 0.9, colour = "grey20") +
    facet_wrap(~ title, scales = "free_x", ncol = min(3, nrow(oi))) +
    scale_fill_manual(values = cols, breaks = names(cols), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.85)), labels = scales::label_comma()) +
    scale_y_discrete(labels = function(x) grp_name(x)) +
    labs(x = NULL, y = NULL, title = paste("Status quo vs tiers 1\u20136 access, by", if (by == "who_region") "WHO region" else "World Bank income group"),
         caption = paste("Deaths and new T2D/CVD/CKD: 5-year cumulative incidence (% of adults at year 0); obesity: prevalence at year 5;",
                         "DALYs per 100,000 adults.\nIn the status quo bar: status quo rate. Right of bars: change vs status quo, absolute (relative).")) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
          axis.line.x = element_line(colour = "grey30", linewidth = 0.4), axis.ticks.x = element_line(colour = "grey30"),
          axis.text.x = element_text(colour = "grey20"),
          axis.text.y = lab_text(unname(gc[levels(r$group)])),
          plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
  fp <- file.path(res_dir, "figures", if (is.null(out_file)) paste0("bar_overall_", by, ".png") else out_file)
  bar_save(fp, p, width = 15, height = 4 + 0.55 * length(levels(r$group)) * ceiling(nrow(oi) / 3), dpi = 300, bg = "white")
  message("bar chart: ", fp); invisible(p)
}

# v2: reduction vs status quo stacked by tier step (optimistic uptake) ----
plot_bar_expansion <- function(res_dir, by = c("who_region", "wb_income"),
                               outcomes = c("ac_death", "t2d_new", "cvd_new", "ckd_new", "obesity", "daly_total_adj"),
                               out_file = NULL, uptake_fp = NULL) {
  by <- match.arg(by)
  sn <- bar_scen_names(uptake_fp)
  r <- bar_data(res_dir, by)[outcome %in% outcomes & step <= 6L][order(group, outcome, step)]
  r[, v := shift(value) - value, by = .(group, outcome)]                   # reduction added by step k
  sq <- r[step == 0L, .(group, outcome, sq = value)]
  r <- r[step >= 1L]
  keys <- tier_key[seq_len(max(r$step))]
  r[, fill_key := factor(keys[step], levels = keys)]
  oi <- copy(bar_outcome_info[outcome %in% outcomes])
  oi[, axis := fifelse(unit == "pp", paste0("Reduction in ", tolower(substr(axis, 1, 1)), substr(sub(" \\(.*\\)$", "", axis), 2, 999), " (pp)"),
                       "DALYs averted per 100,000 adults")]
  r <- oi[r, on = "outcome"][, title := factor(bar_strip(title, axis), levels = bar_strip(oi$title, oi$axis))]
  tot <- r[, .(tot = sum(v), pos = sum(pmax(v, 0))), by = .(title, unit, group, outcome)][sq, on = c("group", "outcome"), nomatch = NULL]
  sq_f <- function(x, unit) ifelse(unit == "pp", ifelse(x < 1, sprintf("%.2f%%", x), sprintf("%.1f%%", x)), formatC(round(x), format = "d", big.mark = ","))
  tot[, lab := paste0(bar_fmt(tot, unit), " (", sprintf("%.1f%%", 100 * tot / sq), " of SQ ", sq_f(sq, unit), ")")]
  gc <- bar_group_cols(levels(r$group), by)
  p <- ggplot(r, aes(y = group, x = v, fill = fill_key)) +
    geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_text(data = tot, aes(y = group, x = pos, label = lab), inherit.aes = FALSE, hjust = -0.08, size = 2.7, colour = "grey20") +
    facet_wrap(~ title, scales = "free_x", ncol = min(3, nrow(oi))) +
    scale_fill_manual(values = setNames(tier_cols[seq_along(keys)], keys), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.6)), labels = scales::label_comma()) +
    scale_y_discrete(labels = function(x) grp_name(x)) +
    guides(fill = guide_legend(ncol = 3, byrow = TRUE)) +
    labs(x = NULL, y = NULL, title = paste0("Reduction vs status quo by priority-tier expansion (", sn$opt, "), by ",
                                           if (by == "who_region") "WHO region" else "World Bank income group"),
         caption = paste("Reduction in percentage points (deaths, new T2D/CVD/CKD: 5-year cumulative incidence; obesity: prevalence at year 5)",
                         "or DALYs averted per 100,000 adults.\nSegments: added reduction of each tier step. Labels: tiers 1\u20136 total (% of the status quo (SQ) rate; SQ rate).")) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
          axis.line.x = element_line(colour = "grey30", linewidth = 0.4), axis.ticks.x = element_line(colour = "grey30"),
          axis.text.x = element_text(colour = "grey20"),
          axis.text.y = lab_text(unname(gc[levels(r$group)])),
          plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
  fp <- file.path(res_dir, "figures", if (is.null(out_file)) paste0("bar_expansion_", by, ".png") else out_file)
  bar_save(fp, p, width = 15, height = 4.5 + 0.45 * length(levels(r$group)) * ceiling(nrow(oi) / 3), dpi = 300, bg = "white")
  message("bar chart: ", fp); invisible(p)
}
