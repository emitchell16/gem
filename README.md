# Global Epidemiology Multimorbidity (GEM) Microsimulation Model

GEM is an open microsimulation model of obesity, type 2 diabetes (T2D), cardiovascular disease (CVD, including stroke), and chronic kidney disease (CKD) in adults aged 18 years and older. It covers all countries with UN population data and reports deaths, incident disease, DALYs and costs.

*This readme file was generated on 2025-11-04 and last updated on 2026-10-26 by Elizabeth Mitchell.*

# GENERAL INFORMATION

**Title of Package:** Global Epidemiology Multimorbidity (GEM) Microsimulation Model

## Author Information

**Name:** Elizabeth Mitchell
**ORCID:** 0009-0005-7176-7465
**Institution:** Emory University
**Address:** Atlanta, GA
**Email:** [estaton@emory.edu](mailto:estaton@emory.edu) / [ec.mitchell3@gmail.com](mailto:ec.mitchell3@gmail.com)

## Co-Author Information

**Name:** Hui Shao
**ORCID:** 0000-0002-4088-546
**Institution:** Emory University
**Address:** Atlanta, GA
**Email:** [hui.shao@emory.edu](mailto:hui.shao@emory.edu)

## Overview

* Date of data collection: January 2025 - May 2026; Parameters last updated: Oct 2026
* Type of data: published estimates, cohort studies, epidemiological and population data
* Date of initial model development: 2025-2026
* Information about EM & HS funding sources that supported the project: T32HL130025 (NHLBI), P30DK111024 (NIDDK), R01DK133465 (NIDDK), and 1R01DK143456 (NIDDK)

---

# SHARING/ACCESS INFORMATION

* Licenses/restrictions placed on the data: Input data from public sources, some of which have individual access requirements
* Links to publications that cite or use the data: https://diabetesjournals.org/diabetes/article/74/Supplement_1/2043-LB/159453
* Links to other publicly accessible locations of the data:

### Key Sources

1. United Nations Population Division. Home Page | Data Portal [Internet]. Population Division Data Portal. Available from:
   https://population.un.org/dataportal/home?df=f48f8e56-4fa5-4a22-bbd4-8b066cb8b1b9

2. Global Burden of Disease Collaborative Network. Global Burden of Disease Study 2023 (GBD 2023). Seattle, United States: Institute for Health Metrics and Evaluation (IHME), 2025. Available from: https://www.healthdata.org/research-analysis/gbd-data

      * Not provided in this repository in accordance with the applicable data use agreement. To reproduce the analyses, obtain the GBD 2023 data directly and save the files with the names given in "Input data" below
      * To download data from GBD:
   https://vizhub.healthdata.org/gbd-results/

   
3. NCD Risk Factor Collaboration (NCD-RisC). National Body-Mass Index Data Download. Available from: https://ncdrisc.org/data-downloads-adiposity.html
      * Not provided in this repository in accordance with the applicable data use agreement. To reproduce the analyses, obtain the NCD-RisC 2024 data directly and save the files with the names given in "Input data" below

        
#### FOR OTHER DATA INPUTS AND CITATIONS SEE FILE "GEM_documentation_workbook.xlsx"


## Repository layout

| Folder | Content |
|---|---|
| `00_Code/R/` | Core GEM microsimulation: inputs, baseline sample, yearly update, DALYs, costing, PSA, post-simulation tables. No intervention-specific code. |
| `00_Code/R/1_function_scripts/GEM_Model_Intervention_Hooks.R` | Hook functions. Defaults do nothing, so the core alone simulates a population without an intervention. |
| `00_Code/R/2_validation_simulation_codes/` | Validation against published studies. |
| `00_Code/examples/` | `run_example.R`: small core-only test run. |
| `00_Code/projects/<name>/` | One folder per application. A project adds its intervention by overriding the hooks. |
| `00_Code/projects/semaglutide_global_access/` | Worked example: semaglutide access policies (`R/glp1_*.R` modules, scenario and plot scripts). |
| `01_Data/` | Input data (see "Input data"). Only the two files marked "this repository" are tracked. |
| `02_Working/`, `03_Out/` | Intermediate files and results written by the model (created on first run, not tracked). |

---
## Requirements

- R (tested with R version 4.5.2). Packages are checked by `load_gem_packages()`, which prints an install command for any missing package.
- Set the working directory to the repository root (the folder with the empty `.here` file). Paths are found with `here::here()`.
- Large simulation files are written to a local temp folder (`GEM_TEMP_DIR`, default `GEM_Local_Temp` under `LOCALAPPDATA` or `~/.local/share`).
- Optional environment variables: `GEM_DATA_DIR`, `GEM_WORKING_DIR`, `GEM_OUT_DIR`, `GEM_ROOT` to place folders elsewhere.
- Full run of 5% random sample, 5-year: ~40 min. Example run: ~5 min.
---

# DATA & FILE OVERVIEW
## Input data

Raw GBD, NCD-RisC and UN files are not included (data use agreements). Download each file and save it in `01_Data/` with the name below. `<year>` is the data year set in `config_inputs()`.

| File | Source |
|---|---|
| `GEM_input_parameters.xlsx` | this repository |
| `who_world_bank_country_groups_iso.csv` | this repository |
| `un_pop_data_<year>.csv` | UN World Population, UN Data Portal |
| `IHME-GBD_prev_inc_male_<year>.csv`, `IHME-GBD_prev_inc_female_<year>.csv` | IHME GBD |
| `IHME-GBD_mortality_by_cause_<year>.csv` | IHME GBD |
| `IHME-GBD_life_expectancy_<year>.csv` | IHME GBD |
| `IHME_GBD_disability_weights_<year>.xlsx` | IHME GBD |
| `NCD_RisC_Lancet_2024_BMI_male_age_specific_country.csv` and the female file | NCD-RisC |
| `for_costing/ppp_LCU.xlsx`, `for_costing/WB_LCU_per_USD.xlsx` | World Bank (costing only) |

#### GEM_input_parameters.xlsx

The `in_use` column marks the rows the code reads.

Core model:

| Sheet | Content |
|---|---|
| `baseline_dataset` | baseline parameters; `lower`, `upper`, `psa_dist`, `psa_family` define PSA distributions |
| `obesity_rr_nejm` | age-specific risk ratios per 5 BMI units |
| `bmi_rr_continuous` | continuous-BMI risk ratios |
| `t2d_rr`, `cvd_rr_stroke_bysex` | T2D and stroke risk ratios |
| `us_health_costs`, `intl_health_costs`, `us_cpi_medical` | disease event costs (inflated to 2025 prices) |
| `PPP conversion rate_LCUperID` | purchasing power parity factors |

Semaglutide project (`00_Code/projects/semaglutide_global_access/`):

| Sheet | Content |
|---|---|
| `semaglutide_treatment_effects` | treatment effects on BMI and disease outcomes (trial and real-world) |
| `treatment_dynamics` | discontinuation and restart rates by treatment year |
| `glp1_uptake` | annual start rates by indication, sex, age (base, optimistic, full) |
| `glp1_baseline_use` | GLP-1 RA use at year 0 by country, indication, age |
| `rwe_sema_costs`, `alt_sema_costs` | semaglutide price scenarios (current prices; missing countries imputed by income group, WHO region, global median) |

#### Country groups file

`who_world_bank_country_groups_iso.csv` columns: `iso3` (ISO 3166-1 alpha-3), `who_region` (WHO region), `world_bank_income_group` (World Bank income group).

#### GEM_data_validation_workbook.xlsx

Contains reference study information and extracted data used in model validation.


*Note: Users can download and use different years of population and epidemiological data to model different years. It is recommended to use the same file naming conventions and specify the `_data_yr` variables in the config function.*
---

---
# QUICK START

### Model overview

- Annual state updates for BMI, T2D, CVD (including stroke), CKD, ESKD, and death, using joint risk ratios and common random numbers across scenarios.
- Weighted stratified baseline sample by country, sex and 5-year age group (weight = population / sample size).
- DALYs (discounted), health-state costs, probabilistic sensitivity analysis (PSA).
- Interventions (treatment, policy, program) are added by project modules through hook functions. The core model contains no intervention-specific code.

Set `GEM_DATA_DIR` if the input data are not in `01_Data/`.

### Run the core model

```r
# from the repository root, with input data in 01_Data/
source("00_Code/examples/run_example.R")
```

Settings are defined in `00_Code/R/configs.R`:

```r
input_data_config <- config_inputs(sample_prop = 0.20)
simulation_config <- config_simulation(sim_years = 10)
print_config(input_data_config, simulation_config)
manifest <- run_gem(input_data_config, list(my_scenario = simulation_config))
```

`run_gem()` reads inputs, draws a weighted baseline sample, and simulates each scenario.
With a project module it also simulates an intervention arm against the shared comparison arm and writes result tables to `03_Out/results_<run tag>/`.

Typical order for a full analysis: validation and mortality calibration (`00_Code/R/2_validation_simulation_codes/GEM_Model_Validation_00_Main.R`), primary runs, PSA, results and figures.

### Add an intervention

The core calls hook functions at fixed points of the yearly update. Defaults do nothing. A project overrides the hooks it needs, in its own files.

| Hook group | Purpose |
|---|---|
| `hook_load_inputs`, `hook_pre_run` | Load intervention inputs; run set-up (for example calibration) |
| `hook_make_policy`, `hook_prepare`, `hook_country_context` | Build the policy and the per-country context |
| `hook_init_columns`, `hook_init_baseline` | Intervention state columns; year-0 state |
| `hook_start_of_year`, `hook_modify_risk`, `hook_bmi_effect`, `hook_costs`, `hook_end_of_year`, `hook_new_entrants` | Yearly update: transitions, effects on event risks and BMI, costs, state changes |
| `hook_label_tags`, `hook_scenario_label`, `hook_manifest_fields`, `hook_print_config` | Scenario labels, manifest, config printout |
| `hook_output_cols`, `hook_results_extra`, `hook_postsim` | Saved columns; extra tables; cross-scenario tables and figures |
| `hook_psa_inputs`, `hook_psa_draw` | Intervention inputs in the PSA |

Steps:
1. Create `00_Code/projects/<name>/R/` with files that define the hooks you need. Signatures and defaults: `00_Code/R/1_function_scripts/GEM_Model_Intervention_Hooks.R`.
2. Wrap `config_simulation()` to add your settings. Set `intervention` (a name) and `intervention_arm = TRUE`.
3. Run `run_gem(input_config, sim_configs, project_modules = "00_Code/projects/<name>/R")`.

Hooks that change person-level data must do so by reference (`data.table::set()` or `:=`).

Worked example: `00_Code/projects/semaglutide_global_access/` (GLP-1 access policies).

```r
source("00_Code/projects/semaglutide_global_access/run_scenarios.R")
```

Settings use `config_semaglutide()`:

```r
# after source("00_Code/projects/semaglutide_global_access/R/glp1_config.R")
simulation_config <- config_semaglutide(sim_years = 10, glp1 = list(coverage_rate = 0.25))
```

### Outputs

Core: person-level simulation files in the local temp folder (`sim=<label>/year=<y>/`), `scenario_manifest.csv`, and per-scenario tables in `03_Out/results_<run tag>/<label>/`.

Semaglutide project (in `03_Out/results_<run tag>/`, added by its module): `all_scenarios_access.csv`, `all_scenarios_outcomes.csv`, `status_quo_summary.csv`, `results_summary_fmt.xlsx`, `figures/`, and `psa/`.
---


---
## Citation

Mitchell E, Shao H. Global Epidemiology Multimorbidity (GEM) Microsimulation Model. Version 1.0. 2026. https://github.com/emitchell16/gem
Publication: <DOI tba>

## License

MIT License (see `LICENSE`). Third-party input data (IHME GBD, NCD-RisC, UN, World Bank) are not included and remain subject to their own terms.

## Disclaimer

GEM is a research tool. Results depend on input data and assumptions and are not clinical advice. The software is provided without warranty.

## Issues and contact

Report bugs or questions through GitHub Issues, or email estaton@emory.edu.