# projects/semaglutide_global_access/run_circleplots.R
#
# Standalone: circle plots from finished simulation outputs (no simulation rerun)
#   1. plot_gem_circles(): tiers 1-6 optimistic and full uptake (obesity, deaths; by WHO region), shared scales
#   2. expansion_scenarios_circleplot(): tiers 1..6 stacked by tier step; deaths, DALYs, obesity
#      (A by WHO region, B by World Bank income group) + 4-panel (WHO / income x obesity / deaths)
#   outputs: 03_Out/results_<run_tag>/figures/
#   needs: sim=<label> folders and scenario_manifest.csv in the local temp run folder

# setup -----
.code <- if (dir.exists(here::here("00_Code"))) here::here("00_Code") else here::here()   # repo layout: 00_Code/ or flat
source(file.path(.code, "R", "setup_environment.R")); load_gem_packages(); source(file.path(.code, "R", "configs.R"))
paths <- get_paths()
source_gem_modules(paths, file.path(.code, "projects", "semaglutide_global_access", "R"))   # core + GLP-1 module (plots)

run_tag <- "10.20"                                          # results_<run_tag> and gem_sims_..._<run_tag>
cand <- list.dirs(paths$local_temp, recursive = FALSE)
cand <- cand[grepl(paste0("^gem_sims_.*_", gsub(".", "\\.", run_tag, fixed = TRUE), "$"), basename(cand))]
if (!length(cand)) stop("no gem_sims_*_", run_tag, " folder in ", paths$local_temp)
out_root <- cand[which.max(file.info(cand)$mtime)]          # newest if several (e.g. test and full run)
res_dir  <- file.path(paths$out, paste0("results_", run_tag))
cat("Run folder:", out_root, "\nResults:", res_dir, "\n")
man <- fread(file.path(out_root, "scenario_manifest.csv"))
print(man[, .(config, label, base_label, mode)])
