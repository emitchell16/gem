# R/1_function_scripts/GEM_Model_Simulate.R
# Author: Liz Mitchell
# GEM model: call simulation and update yearly data
# scenario label folders; keep_years; mortality calibration on/off by config
# intervention state and effects come from hook_*() functions (GEM_Model_Intervention_Hooks.R)

#######'
# Primary Model wrapper function --------------------------------------------------------
#  files_yr0: baseline sample folder (initialize_baseline); policy: object from hook_make_policy() or NULL (no intervention)
#  opts: named list passed to hook_country_context() (e.g. switches for an intervention module)
#  output hive: <out_folder_path>/sim=<label>/year=<y>/cohort_<i>.parquet
#######'
# columns saved in sim=<label> outputs (kept small for disk space) ----
#' what post-sim tables, PSA summaries and plots read, plus hook_output_cols(); other (working)
#' columns stay in memory only. Absent columns are skipped. Costs: all cost_* columns when costing is on
#' options(gem.keep_all_cols = TRUE): save every column (debugging, model tests)
sim_output_cols <- function(nm) {
  if (isTRUE(getOption("gem.keep_all_cols"))) return(nm)
  keep <- unique(c("pid", "iso3", "w", "age", "female", "bmi", "bmi_class0",
                   all_outcomes_list, "obesity", "t2d", "cvd", "ckd", "eskd", "ac_death", "stroke_ever", "stroke_n", "cvd_event_n",
                   "stroke", "cvd_event", hook_output_cols()))
  intersect(nm, c(keep, grep("^cost_", nm, value = TRUE)))
}
.pq_codec <- function() if (arrow::codec_is_available("zstd")) "zstd" else "snappy"

main_sim_file_management <- function(files_yr0, out_folder_path, sim_years, policy = NULL, label,
                                     cs = config_sims, keep_years = c(0, sim_years), opts = list(),
                                     chunk_size = config$chunk_cap) {
  cat("Start time:", format(Sys.time(), "%c"), "| sim =", label, "\n")
  start_time <- Sys.time()
  on.exit(gc(verbose = FALSE))

  temp_dir <- file.path(local_working_folder, paste0("temp_simulation_", label))
  unlink(temp_dir, recursive = TRUE)
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  group_umbrella_folder <- file.path(out_folder_path, paste0("sim=", label))
  for (yr in keep_years) dir.create(file.path(group_umbrella_folder, paste0("year=", yr)), recursive = TRUE, showWarnings = FALSE)
  files <- list.files(path = files_yr0, pattern = "\\.parquet$", full.names = TRUE)
  # batches: consecutive cohort files of one country up to chunk_size persons -> fewer tasks, reads and writes;
  #   memory per task stays <= chunk_size (larger single files are chunked as before). index cached next to baseline folder
  idx_fp <- paste0(normalizePath(files_yr0, mustWork = FALSE), "_file_index.csv")   # next to (not in) the parquet folder
  fi <- if (file.exists(idx_fp)) fread(idx_fp) else NULL
  if (is.null(fi) || !setequal(fi$f, basename(files))) {
    fi <- rbindlist(lapply(files, function(f) {
      x <- read_parquet(f, col_select = "iso3")$iso3
      data.table(f = basename(f), iso3 = if (length(x)) as.character(x[1]) else NA_character_, n = length(x)) }))
    fwrite(fi, idx_fp)
  }
  fi <- fi[match(basename(files), f)][, i := seq_len(.N)][!is.na(iso3) & n > 0]
  setorder(fi, iso3, i)
  bid <- integer(nrow(fi)); b <- 0L; cc0 <- ""; acc <- Inf
  for (r in seq_len(nrow(fi))) {
    if (fi$iso3[r] != cc0 || acc + fi$n[r] > chunk_size) { b <- b + 1L; acc <- 0; cc0 <- fi$iso3[r] }
    acc <- acc + fi$n[r]; bid[r] <- b
  }
  groups <- split(fi$i, bid)
  # intervention set-up from the whole baseline sample (identical across arms); NULL without an intervention
  ctx0 <- hook_prepare(cs, files_yr0)

  future_lapply(seq_along(groups), function(i) {
    out_name <- paste0("batch_", i, ".parquet")
    # skip batches already written (restart after crash)
    if (file.exists(file.path(group_umbrella_folder, paste0("year=", max(keep_years)), out_name))) return(NULL)
    temp_dir_i <- file.path(temp_dir, paste0("group=", i))
    dir.create(temp_dir_i, recursive = TRUE, showWarnings = FALSE)
    options(datatable.auto.index = FALSE, datatable.use.index = FALSE)   # workers do not inherit options
    group_data <- rbindlist(lapply(groups[[i]], function(f) {
      d <- as.data.table(read_parquet(files[f])); setindex(d, NULL)     # drop stale index from parquet
      d[, cohort_file := f] }), use.names = TRUE, fill = TRUE)
    if (nrow(group_data) == 0) return(NULL)
    if (cs$open_pop_cohort) group_data[, pid := as.character(pid)]   # pid kept: pairs persons across arms (CRN)

    # chunks: same crn_id per chunk in every arm -> identical random numbers across arms
    N <- nrow(group_data)
    chunks <- split(seq_len(N), ceiling(seq_len(N) / chunk_size))
    res <- lapply(seq_along(chunks), function(j) {
      dir_j <- file.path(temp_dir_i, paste0("chunk_", j)); dir.create(dir_j, showWarnings = FALSE)
      simulate_disease_progression(group_data[chunks[[j]]], dir_j, sim_years, policy, cs,
                                   crn_id = i * 100 + if (length(chunks) > 1) j else 0,
                                   keep_years = keep_years, ctx0 = ctx0, opts = opts)
    })
    rm(group_data)
    for (yr in as.character(keep_years)) {
      fv <- vapply(res, `[[`, character(1), yr)
      # mmap = FALSE: file released after read (Windows lock left temp stragglers)
      d  <- rbindlist(lapply(fv, read_parquet, mmap = FALSE), use.names = TRUE, fill = TRUE)
      d  <- d[, sim_output_cols(names(d)), with = FALSE]                  # saved columns only
      write_parquet(d, file.path(group_umbrella_folder, paste0("year=", yr), out_name), compression = .pq_codec())
      rm(d)
    }
    unlink(temp_dir_i, recursive = TRUE)
    NULL
  }, future.seed = config$seed)

  unlink(temp_dir, recursive = TRUE)
  cat("Elapsed time:", round(difftime(Sys.time(), start_time, units = "mins"), 2), "mins\n")
  invisible(group_umbrella_folder)
}

#######'
# Yearly update function --------------------------------------------------------
#######'
# BMI percentile for baseline files without bmi_u (legacy): u uniform within person's current class interval
add_bmi_u_if_missing <- function(input_data, mp) {
  if ("bmi_u" %in% names(input_data)) return(input_data)
  Fk <- cbind(0, as.matrix(mp[, paste0("bmi_F", 1:6), with = FALSE]), 1)
  k  <- ifelse(input_data$class3_obesity == 1, 7, ifelse(input_data$class2_obesity == 1, 6,
         ifelse(input_data$class1_obesity == 1, 5, NA)))
  lo <- ifelse(is.na(k), 0, Fk[cbind(seq_len(nrow(Fk)), pmax(k, 1))])
  hi <- ifelse(is.na(k), Fk[, 5], Fk[cbind(seq_len(nrow(Fk)), pmin(k + 1, 8))])
  input_data[, bmi_u := lo + runif(.N) * (hi - lo)]
  input_data
}

# simulation state columns (placeholders); intervention columns from hook_init_columns(); shared with validation ----
init_sim_columns <- function(dt, crn_id, costing = FALSE) {
  dt[, c("t2d_new", "cvd_new", "ckd_new", "eskd_new", "stroke", "stroke_this_year", "death_year",
         "stroke_ever", "stroke_n", "years_since_stroke",
         "cvd_event", "cvd_event_ever", "cvd_event_n", "years_since_cvd_event") := 0L]
  dt[, c("yll_adj", "yld_adj_total", "yld_ckd_adj_total", "yld_cvd_adj_total", "yld_t2d_adj_total",
         "bmi_effect") := 0.0]
  dt[, crn_id := crn_id]
  hook_init_columns(dt, crn_id)
  if (costing) {
    dt[, c("cost_t2d_tot", "cost_death_tot", "cost_ckd_tot", "cost_eskd_tot", "cost_stroke_tot", "cost_cvd_tot",
           "cost_conditions_total", cost_total_cols) := 0.0]
  }
  invisible(dt)
}

simulate_disease_progression <- function(input_data,     # for the sex-country cohort
                                         sim_temp_dir,   # temp directory
                                         sim_years,
                                         policy = NULL,  # hook_make_policy() object; NULL = no intervention
                                         cs = config_sims,
                                         crn_id = 1,             # common-random-number id (cohort file x chunk)
                                         keep_years = c(0, sim_years),
                                         ctx0 = NULL,            # from hook_prepare()
                                         opts = list()) {        # passed to hook_country_context()
  setDT(input_data)
  current_iso3  <- input_data[1, iso3]
  cohort_params <- as.data.table(all_cohort_param_rates[all_cohort_param_rates$iso3 == current_iso3, ])
  setkey(cohort_params, female, age_low)
  mp0 <- cohort_params[input_data, on = .(female, age_low = age), roll = TRUE, rollends = c(TRUE, TRUE)]

  # baseline BMI from percentile (rank tracking; GEM_Model_BMI_RR.R) ----
  set.seed(as.integer((crn_id %% 2e5) * 100 + 9))
  input_data <- add_bmi_u_if_missing(input_data, mp0)
  if (!"w" %in% names(input_data)) input_data[, w := 1 / config$sample_prop]   # legacy files (uniform sampling fraction)
  input_data[, bmi := bmi_from_u(bmi_u, mp0)]
  set_obesity_from_bmi(input_data, seq_len(nrow(input_data)))
  input_data[, `:=`(bmi0 = bmi,
                    bmi_class0 = fcase(class3_obesity == 1, "class3", class2_obesity == 1, "class2",
                                       class1_obesity == 1, "class1", bmi >= 25, "overweight", default = "normal"))]
  init_sim_columns(input_data, crn_id, isTRUE(cs$costing))
  # y0 joint-RR scaling constants + baseline prevalences (GEM_Model_Risk_Equations.R); identical across arms
  #   per original cohort file (country-sex-5yr cohort) when files are batched
  if ("cohort_file" %in% names(input_data) && length(unique(input_data$cohort_file)) > 1) {
    if (!"stroke_ever" %in% names(input_data)) input_data[, stroke_ever := 0L]
    for (g in unique(input_data$cohort_file)) {
      r  <- which(input_data$cohort_file == g); al <- r[input_data$ac_death[r] == 0]
      kk <- rr_baseline_consts(as.list(input_data[al, ..rr_cols]))
      for (nm in names(kk)) set(input_data, i = r, j = nm, value = kk[[nm]])
    }
  } else add_rr_baseline_cols(input_data)
  # intervention context for this cohort file; year-0 intervention state (reporting; effects start in year 1)
  ictx <- hook_country_context(ctx0, policy, current_iso3, opts)
  hook_init_baseline(ictx, input_data, crn_id, mp0)

  # youngest single-year age at baseline cloned as entrants in open cohort
  if (cs$open_pop_cohort && min(input_data$age) <= min_age_model) {
    fwrite(input_data[age == min_age_model], file.path(sim_temp_dir, "youngest_y0_person_data.csv"))
  }

  out <- list()
  write_year <- function(d, yr) {
    f <- file.path(sim_temp_dir, paste0("simulated_year_", yr, ".parquet"))
    write_parquet(d, f); out[[as.character(yr)]] <<- f
  }
  if (0 %in% keep_years) write_year(input_data, 0)
  pd <- input_data; rm(input_data)
  for (year in seq_len(sim_years)) {
    pd <- update_cohort(pd, cohort_params, year, ictx = ictx,
                        open_cohort = cs$open_pop_cohort, sim_temp_dir = sim_temp_dir,
                        cal = if (isTRUE(cs$mortality_calibration)) acm_calibration else NULL,
                        dw_sensitivity = cs$daly_sensitivity, costing = isTRUE(cs$costing))
    if (year %in% keep_years) write_year(pd, year)
  }
  out     # named list of parquet files by year
}
