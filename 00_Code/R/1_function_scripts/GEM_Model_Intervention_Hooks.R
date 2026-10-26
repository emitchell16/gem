# R/1_function_scripts/GEM_Model_Intervention_Hooks.R

# description ------
# Hook functions that let a project add an intervention (treatment, policy or program) to the core model.
# The core calls each hook at a fixed point. The defaults below do nothing, so the core alone runs a
# no-intervention simulation. A project overrides a hook by defining a function with the same name in its own
# module files, sourced after the core (run_gem(..., project_modules = <folder>) or source_gem_modules(paths, <folder>)).
# Worked example: projects/semaglutide_global_access/R/glp1_hooks.R
#
# Design notes:
# - hooks are called by name so parallel workers (future) export them with the core functions
# - ictx = intervention context for one cohort file (country); ystate = state for one simulated year;
#   both are opaque to the core and defined by the project
# - hooks that change person-level data must do it by reference (data.table set() or :=)
# - the config marker cfg$intervention (non-NULL) shows that an intervention module is active;
#   cfg$intervention_arm = TRUE runs an intervention arm next to the comparison arm

# labels and configuration ------
hook_label_tags <- function(cs) character(0)                  # extra tags for the comparison-arm label
hook_scenario_label <- function(cs) paste(c(paste0(cs$sim_years, "yr_intervention"), scenario_tags(cs)), collapse = "_")
hook_manifest_fields <- function(cs) list()                   # scalar fields added to the scenario manifest
hook_print_config <- function(cs) invisible(NULL)             # extra lines for print_config()

# inputs and run set-up ------
hook_load_inputs <- function(cs) invisible(NULL)              # end of set_in_tables.R: load intervention inputs (globals)
hook_pre_run <- function(cs, baseline_fp, out_root) invisible(NULL)   # before the simulations of one config (e.g. calibration)
hook_make_policy <- function(cs, baseline_fp) NULL            # policy object for the intervention arm (NULL = none)

# simulation ------
hook_prepare <- function(cs, files_yr0) NULL                  # once per arm, before parallel work; returns ctx0
hook_country_context <- function(ctx0, policy, iso3, opts) NULL   # per cohort file; returns ictx
hook_init_columns <- function(dt, crn_id) invisible(dt)       # add intervention state columns (by reference)
hook_init_baseline <- function(ictx, dt, crn_id, mp) invisible(dt)    # year-0 intervention state (e.g. prevalent users)
hook_start_of_year <- function(ictx, new_state_vals, pd, living_idx, year) NULL   # start of year; returns ystate
hook_modify_risk <- function(cond, p, ictx, ystate, pd, idx) p        # change event probabilities p (aligned with idx)
hook_bmi_effect <- function(ictx, ystate, new_state_vals, living_idx) invisible(NULL)   # set column bmi_effect (proportion)
hook_costs <- function(ictx, ystate, new_state_vals, time_alive, disc) invisible(NULL)  # intervention costs (costing on)
hook_end_of_year <- function(ictx, new_state_vals) invisible(NULL)    # end-of-year state changes
hook_new_entrants <- function(ictx, entrants, crn_id, year) invisible(entrants)   # open cohort: new entrants
hook_output_cols <- function() character(0)                   # extra columns kept in saved simulation files

# post-simulation ------
hook_results_extra <- function(baseline_ds_path, out_folder, iso_regions) invisible(NULL)   # extra tables per scenario
hook_postsim <- function(out_root, manifest, res_dir) invisible(NULL)  # cross-scenario tables and figures

# probabilistic sensitivity analysis ------
hook_psa_inputs <- function() list()                          # named list of intervention input objects (globals) to capture
hook_psa_draw <- function(base, out, zf) out                  # draw them: read base$ext, write out$ext; zf(family) gives a shared z
