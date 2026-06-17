# -------------------'
# 00_Code/R/setup_environment.R
# -------------------'
load_gem_packages <- function() {
  pkgs <- c(
    "tidyverse", "stringi", "pbapply", "writexl", "readxl",
    "openxlsx", "R.utils", "rlang", "truncnorm",
    "here",
    "data.table", "arrow",
    "future.apply", "progressr",
    # for validation results plotting:
    "ggplot2", "ggthemes", "ggrepel", "patchwork", "ggtext"
  )
  
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package not installed: ", p,
           ". Please install it with install.packages('", p, "')")
    }
    library(p, character.only = TRUE)
  }

  plan(multisession, workers = max(1, parallel::detectCores() - 3))
  options(future.globals.maxSize = 5 * 1024^3)
  options(arrow.parquet.thrift_string_size_limit = 5e9) 
}

get_paths <- function(){
  local_temp_folder <- file.path(Sys.getenv("LOCALAPPDATA"), "GEM_Local_Temp")
  if(!dir.exists(local_temp_folder)) {dir.create(
    local_temp_folder, recursive = TRUE, showWarnings = FALSE)}  

  root <- here()  
  paths <- list(
    root            = root,
    code            = file.path(root, "00_Code"),
    r_code          = file.path(root, "00_Code", "R"),
    fn_scripts      = file.path(root, "00_Code", "R", "1_function_scripts"),
    scripts         = file.path(root, "00_Code", "scripts"),
    hpc             = file.path(root, "00_Code", "HPC"),
    data            = file.path(root, "01_Data"),
    working         = file.path(root, "02_Working"),
    out             = file.path(root, "03_Out"),
    local_temp      = local_temp_folder,
    v_code          = file.path(root, "00_Code", "R", "2_validation_simulation_codes")
    )
  return(invisible(paths))
}

