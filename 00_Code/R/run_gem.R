# R/run_gem.R
#
# run_gem(): config-driven pipeline (called from a scenario script, e.g. examples/run_example.R)
#   1. inputs (prep_main_inputs) -> 2. weighted baseline sample (initialize_baseline) + summary stats
#   3. per simulation config: input tables; comparison arm and, with a project module, the intervention arm
#      (shared comparison arm per horizon), or PSA (run_psa)
#   4. scenario manifest + post-simulation tables
#   comparison arm label = base_label_of(cs): shared by configs with the same horizon and tags
#   project_modules: folder or R files with hook_*() functions (GEM_Model_Intervention_Hooks.R); NULL = core only
#
# outputs (out_root = <local temp>/gem_sims_<sample>_<run_tag>):
#   sim=<label>/year=<y>/batch_<i>.parquet      person-level results (year 0 and final year)
#   scenario_manifest.csv                       config name -> label, mode, settings
#   psa_<label>/psa_runs_partial.rds            PSA checkpoint (restart resumes)
#   paths$out/results_<run_tag>/psa/*_<label>.csv  PSA intervals, convergence, run-level results
#   paths$out/results_<run_tag>/                post-simulation tables (by scenario label)

#' run_tag: output folder tag; default = in_data_run_date (fixed) so a restart on a later day resumes in the same folder
run_gem <- function(config_init, sim_configs,
                    run_tag = if (!is.null(config_init$in_data_run_date)) config_init$in_data_run_date else config_init$short_date,
                    postsim = TRUE, project_modules = NULL) {
  # sound alert on exit: 3 = finished (fanfare), 1 = error or interrupt (ping)
  done <- FALSE
  on.exit(if (requireNamespace("beepr", quietly = TRUE)) try(beepr::beep(if (done) 3 else 1), silent = TRUE), add = TRUE)
  if (inherits(sim_configs, "list") && !is.null(sim_configs$sim_years)) sim_configs <- list(sim = sim_configs)  # single config
  if (is.null(names(sim_configs)) || any(names(sim_configs) == "")) stop("sim_configs must be a named list of config_simulation() objects")

  # 0. environment (globals used by model modules) ------
  paths <- get_paths()
  assign("paths", paths, envir = .GlobalEnv)
  assign("config", config_init, envir = .GlobalEnv)
  assign("local_working_folder", paths$local_temp, envir = .GlobalEnv)
  assign("min_age_model", config_init$min_age_model, envir = .GlobalEnv)
  assign("crn_stream_offset", 0L, envir = .GlobalEnv)
  source_gem_modules(paths, project_modules)
  cal_fp <- file.path(paths$working, "calibration_estimates.csv")    # mortality calibration (validation step)
  assign("acm_calibration", if (file.exists(cal_fp)) read.csv(cal_fp, stringsAsFactors = FALSE) else NULL, envir = .GlobalEnv)
  if (!file.exists(cal_fp)) message("No calibration_estimates.csv: all-cause mortality uncalibrated (run validation first)")

  # 1. inputs ------
  cohort_inputs <- prep_main_inputs(paths, config_init) |>
    filter(pop_size > 10, age_high >= config_init$min_age_model)     # adults; drops near-empty oldest cohorts
  if (!is.null(config_init$iso3_subset)) {                            # test runs: selected countries only
    cohort_inputs <- filter(cohort_inputs, iso3 %in% config_init$iso3_subset)
    if (!nrow(cohort_inputs)) stop("iso3_subset: no matching countries")
  }
  assign("cohort_input_data", cohort_inputs, envir = .GlobalEnv)

  # 2. weighted baseline sample + summary stats ------
  baseline_fp <- initialize_baseline(cohort_inputs, paths, config_init)
  .date <- if (!is.null(config_init$in_data_run_date)) config_init$in_data_run_date else config_init$short_date
  if (!file.exists(file.path(paths$out, paste0("summary_stats_y0fmt_", .date, ".csv")))) {
    global_summary_stats(config = config_init, date = .date, paths, ds_path = baseline_fp)
  }
  out_root <- file.path(paths$local_temp, paste0("gem_sims_f", config_init$sample_prop, "_n", config_init$min_cohort_n, subset_tag(config_init), "_", run_tag))
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  # 3. simulations by config ------
  manifest <- list()
  for (nm in names(sim_configs)) {
    cs <- sim_configs[[nm]]
    assign("config_sims", cs, envir = .GlobalEnv)
    source(file.path(paths$r_code, "set_in_tables.R"), local = globalenv())
    label <- scenario_label(cs)
    cat("\n==== ", nm, " -> ", label, if (cs$psa) " (PSA)" else "", " ====\n", sep = "")

    hook_pre_run(cs, baseline_fp, out_root)                           # e.g. calibration needed by the intervention
    if (cs$psa) {
      if (!isTRUE(cs$intervention_arm)) stop(nm, ": PSA compares comparison vs intervention arm; set intervention_arm = TRUE")
      run_psa(cs, label, psa_temp_folder = file.path(paths$local_temp, paste0("psa_temp_", label)),
              psa_out_folder = file.path(out_root, paste0("psa_", label)),
              final_dir = file.path(paths$out, paste0("results_", run_tag), "psa"))
    } else {
      main_sim_file_management(baseline_fp, out_root, cs$sim_years, NULL, base_label_of(cs), cs)
      if (isTRUE(cs$intervention_arm)) {
        main_sim_file_management(baseline_fp, out_root, cs$sim_years, hook_make_policy(cs, baseline_fp), label, cs)
      }
    }
    manifest[[nm]] <- data.table(config = nm, label = label, base_label = base_label_of(cs), mode = if (cs$psa) "psa" else "primary",
                                 sim_years = cs$sim_years, intervention_arm = isTRUE(cs$intervention_arm), costing = cs$costing,
                                 mortality_calibration = cs$mortality_calibration,
                                 run_tag = run_tag, out_root = out_root, finished = format(Sys.time(), "%Y-%m-%d %H:%M"))
    extra <- hook_manifest_fields(cs)
    if (length(extra)) manifest[[nm]] <- cbind(manifest[[nm]], as.data.table(extra))
  }
  manifest <- rbindlist(manifest)
  mf <- file.path(out_root, "scenario_manifest.csv")
  if (file.exists(mf)) manifest <- unique(rbind(fread(mf)[!config %in% manifest$config], manifest, fill = TRUE), by = "config")
  fwrite(manifest, mf)

  # 4. post-simulation tables ------
  if (postsim) {
    res_dir <- file.path(paths$out, paste0("results_", run_tag))
    prim <- manifest[mode == "primary" & intervention_arm == TRUE]
    for (k in seq_len(nrow(prim))) {
      generate_results_tables(out_root, prim$sim_years[k], prim$label[k], file.path(res_dir, prim$label[k]),
                              base_label = prim$base_label[k])
    }
    hook_postsim(out_root, manifest, res_dir)          # cross-scenario tables and figures (project module)
  }
  done <- TRUE
  invisible(manifest)
}
