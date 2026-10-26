# R/setup_environment.R
# package loading, project paths and model module sourcing
#
# Core model: R/. Project modules (intervention-specific code, hooks, scripts): projects/<name>/ (R/ for module files).
# Expected project layout (root found with here::here(); add an empty .here file or an .Rproj at the repo root):
#   R/                     model code (this folder)
#   projects/              project-specific scenario and reporting scripts
#   01_Data/               input data files you provide (see README)
#   02_Working/            intermediate files written by the model
#   03_Out/                results written by the model
# If a 00_Code/ folder exists at the root, code is read from 00_Code/R instead.
# Data, working and output folders can sit elsewhere: set environment variables GEM_DATA_DIR, GEM_WORKING_DIR, GEM_OUT_DIR
#   (for example in .Renviron or Sys.setenv() before get_paths()). GEM_ROOT overrides the root found by here().
# Large simulation files go to a local temp folder: GEM_TEMP_DIR environment variable if set,
#   otherwise <LOCALAPPDATA or ~/.local/share>/GEM_Local_Temp.

load_gem_packages <- function() {
  pkgs <- c(
    "tidyverse", "stringi", "pbapply", "writexl", "readxl",
    "openxlsx", "R.utils", "rlang", "truncnorm",
    "here",
    "data.table", "arrow",
    "future.apply", "progressr",
    # for results plotting:
    "ggplot2", "ggthemes", "ggrepel", "patchwork", "ggtext"
  )
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Packages not installed: ", paste(missing_pkgs, collapse = ", "),
         "\nInstall with: install.packages(c(\"", paste(missing_pkgs, collapse = "\", \""), "\"))")
  }
  for (p in pkgs) library(p, character.only = TRUE)
  # optional: sound alert when run_gem() finishes or stops
  if (requireNamespace("beepr", quietly = TRUE)) library(beepr)

  plan(multisession, workers = max(1, parallel::detectCores() - 3))
  options(future.globals.maxSize = 5 * 1024^3)
  options(arrow.parquet.thrift_string_size_limit = 5e9)
  # parquet files keep data.table index attributes, which are invalid after reading; do not create or use indices
  options(datatable.auto.index = FALSE, datatable.use.index = FALSE)
}

get_paths <- function() {
  local_temp_folder <- Sys.getenv("GEM_TEMP_DIR", unset = "")
  if (local_temp_folder == "") {
    base_dir <- Sys.getenv("LOCALAPPDATA", unset = "")
    if (base_dir == "") base_dir <- file.path(path.expand("~"), ".local", "share")   # macOS / Linux
    local_temp_folder <- file.path(base_dir, "GEM_Local_Temp")
  }
  dir.create(local_temp_folder, recursive = TRUE, showWarnings = FALSE)

  root <- Sys.getenv("GEM_ROOT", unset = "")
  root <- if (root == "") here() else normalizePath(root, winslash = "/", mustWork = TRUE)
  dir_or <- function(var, default) { x <- Sys.getenv(var, unset = ""); if (x == "") default else x }
  code <- if (dir.exists(file.path(root, "00_Code"))) file.path(root, "00_Code") else root
  paths <- list(
    root            = root,
    code            = code,
    r_code          = file.path(code, "R"),
    fn_scripts      = file.path(code, "R", "1_function_scripts"),
    scripts         = file.path(code, if (dir.exists(file.path(code, "scripts"))) "scripts" else "projects"),
    data            = dir_or("GEM_DATA_DIR", file.path(root, "01_Data")),
    working         = dir_or("GEM_WORKING_DIR", file.path(root, "02_Working")),
    out             = dir_or("GEM_OUT_DIR", file.path(root, "03_Out")),
    local_temp      = local_temp_folder,
    v_code          = file.path(code, "R", "2_validation_simulation_codes")
  )
  for (d in c(paths$working, paths$out)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(paths$data)) warning("Data folder not found: ", paths$data,
                                       "\nPut input files there or set GEM_DATA_DIR (see README).", call. = FALSE)
  return(invisible(paths))
}

# model function modules (order: helpers before users) ----
#  project_modules: NULL (core model only), a folder, or R files of a project module; sourced after the core so that
#  hook_*() functions defined there replace the no-op defaults (GEM_Model_Intervention_Hooks.R)
source_gem_modules <- function(paths, project_modules = NULL) {
  for (f in c("GEM_Model_Intervention_Hooks.R", "GEM_Model_BMI_RR.R", "GEM_Model_Risk_Equations.R",
              "GEM_Model_Disability_Estimation.R", "GEM_Model_Update_Cohort_States.R", "GEM_Model_Simulate.R",
              "GEM_Model_PSA_Functions.R", "GEM_Model_1_Baseline_Summary_Stats.R",
              "GEM_Model_Postsim_Functions.R", "GEM_Model_Postsim_Results.R")) {
    source(file.path(paths$fn_scripts, f), local = globalenv())
  }
  source(file.path(paths$r_code, "prep_main_inputs.R"), local = globalenv())
  if (!is.null(project_modules)) {
    fs <- if (length(project_modules) == 1 && dir.exists(project_modules))
      sort(list.files(project_modules, pattern = "\\.[Rr]$", full.names = TRUE)) else project_modules
    if (!length(fs) || !all(file.exists(fs))) stop("project_modules: files not found: ", paste(project_modules, collapse = ", "))
    for (f in fs) source(f, local = globalenv())
  }
  invisible(TRUE)
}
