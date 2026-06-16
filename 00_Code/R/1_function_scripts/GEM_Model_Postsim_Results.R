# -------------------
# Apply post simulation analysis 
# Author: Liz Mitchell
# Date updated: March 2026
# -------------------'

generate_results_tables <- function(sim_hive_folder, sim_years, access_percent, .out_folder, 
                                    conditions_list, se_ci, psa, skip_country_lvl = FALSE ) {

# input helpers -------
  base_label <- paste0(sim_years, "yr_basecase")
  comp_label <- paste0(sim_years,  "yr_", access_percent, "glp1")
  
    if (dir.exists( file.path(sim_hive_folder, paste0("sim=", base_label), "year=0"))) { 
      basecase_root <- sim_hive_folder
    } else if (dir.exists(  file.path(working_folder, paste0("sim=", base_label), "year=0"))) { basecase_root <- working_folder 
    } else {stop("Neither path exists\n")}
  
  
  baseline_ds_path     <-  file.path(basecase_root, paste0("sim=", base_label), "year=0")
  end_basecase_ds_path <- file.path(basecase_root, paste0("sim=", base_label),  paste0("year=", sim_years))
  
  baseline_ds_path_treatment <- file.path(sim_hive_folder, paste0("sim=", comp_label), "year=0")
  end_intervention_ds_path <- file.path(sim_hive_folder,   paste0("sim=", comp_label), paste0("year=", sim_years))
  
  scale_factor_forN = 1 / config$sample_prop
  un_pop <- cohort_input_data |> group_by(iso3) |> summarise(pop_size = sum(pop_size, na.rm = TRUE)*scale_factor_forN, .groups = "drop")
      #5/7/26 added
  
  country_groups2 <- country_groups |> rename(wb_income = world_bank_income_group)
  iso_regions <- country_groups2 |>  select(iso3, who_region)
  
  if (!file.exists(file.path(paths$working, "who_region_with_pop.xlsx")))  {
    # aggregate population by who region
    who_region_with_pop <- iso_regions |> 
      left_join(un_pop, by = "iso3") |> 
      mutate(who_region = if_else(is.na(who_region), "no who region", who_region)) |> 
      group_by(who_region)  %>% 
      summarise(pop_size = sum(pop_size, na.rm = TRUE), .groups = "drop") |> 
      arrange(who_region) %>% 
      bind_rows(
        tibble(
          who_region = "global",
          pop_size   = sum(.$pop_size, na.rm = TRUE) )
      ) |> 
      mutate(pop_size = round(pop_size))
  
    # save to working_folder who_region_with_pop.xlsx (used in graphing)
    write_xlsx(who_region_with_pop, file.path(paths$working, "who_region_with_pop.xlsx"))
  } else {
    who_region_with_pop <- read.xlsx(file.path(paths$working, "who_region_with_pop.xlsx"))
  }
  
  # dictionary mapping who_region to pop size
  region_pop <- setNames(who_region_with_pop$pop_size, who_region_with_pop$who_region)
  
###############'
## Results: Baseline treatment eligible groups table -------
###############'

  # cat("* NOTE: baseline treatment eligibility summary part commented out.\n")
 if (sample_p == 0.1 && !file.exists(file.path(.out_folder, "treatment_distributions_table.xlsx"))) {
   ds_temp<- open_dataset(normalizePath(baseline_ds_path_treatment), format = "parquet")

   summarise_treatment <- function(df) {
     df |>
       summarize(
         n                  = n(),
         n_t2d_bmi27_18y    = sum(if_else(age >= 18 & t2d == 1 & bmi >= 27, 1L, 0L, missing = 0L)),
         n_obesity_12y      = sum(if_else(age >= 12 & bmi >= 30,                  1L, 0L, missing = 0L)),
         n_obesity_wt2d_12y = sum(if_else(age >= 12 & t2d == 1 & bmi >= 30,       1L, 0L, missing = 0L)),
         
         n_geq18 = sum(if_else(age >= 18, 1L, 0L)),
         n_geq45 =  sum(if_else(age >= 45, 1L, 0L)),
       ) |>
       mutate(
         pct_t2d_bmi27_18y    = 100 * n_t2d_bmi27_18y    / n,
         pct_obesity_12y      = 100 * n_obesity_12y      / n,
         pct_obesity_wt2d_12y = 100 * n_obesity_wt2d_12y / n,
         
         scale_factor_forN = 1 / config$sample_prop
       )
   }

   eligible_global <- summarise_treatment(ds_temp) |>
     mutate(Region = "global") |>
     select(Region, everything()) |>
     collect()

   eligible_by_region <- ds_temp |>
     left_join(iso_regions, by = "iso3") |>
     mutate(who_region = if_else(is.na(who_region), "no who region", who_region)) |>
     group_by(who_region) |>
     summarise_treatment() |>
     rename(Region = who_region) |>
     select(Region, everything()) |>
     collect()

   eligible_all <- bind_rows(eligible_global, eligible_by_region) |>
     mutate(Region = str_to_title(Region))

   write_xlsx(eligible_all, file.path(.out_folder, "treatment_distributions_table.xlsx"))
   rm(ds_temp, eligible_global, eligible_by_region, eligible_all); gc()
 }
  
 ###############'
 ## Results: global & regional summary -------
 ###############'
  combined_fn <- paste0("regional_result_table_", sim_years,"_", access_percent ,".xlsx")
  if ( !file.exists(file.path(.out_folder, combined_fn))) {
    regional_result <- make_region_impact_table(sim_years, access_percent, basecase_root, level = "all", se_ci=se_ci, psa=psa)
    write_xlsx(regional_result, file.path(.out_folder, combined_fn))
    if (sim_years == 5 & access_percent == 100 & sample_p == 0.1){
      age18_regional_result <- make_region_impact_table(sim_years, access_percent, basecase_root,level = "all", filter_min_age = "age>=18")
      age45_regional_result <- make_region_impact_table(sim_years, access_percent, basecase_root, level = "all",filter_min_age = "age>=45")

      write_xlsx(age18_regional_result, file.path(.out_folder, paste0("18_regional_", sim_years, "_", access_percent, ".xlsx")))
      write_xlsx(age45_regional_result, file.path(.out_folder, paste0("45_regional_", sim_years, "_", access_percent, ".xlsx")))
        }
    gc(verbose = F)
    cat("finished regional and combined. ")
  } else{ # load previous:
    cat("- loading existing region result tables \n")
    regional_result       <- read_xlsx(file.path(.out_folder, combined_fn))
    if (sim_years == 5 & access_percent == 100 & sample_p == 0.1){
    age18_regional_result <- read_xlsx(file.path(.out_folder, paste0("18_regional_", sim_years, "_", access_percent, ".xlsx")))
    age45_regional_result <- read_xlsx(file.path(.out_folder, paste0("45_regional_", sim_years, "_", access_percent, ".xlsx")))
    }}

  ###############'
  ## Format combined region/global results ------
  ###############'
  region_result_fmt <- format_region_result(regional_result, who_region_with_pop, .out_folder, se_ci=se_ci)

  if (exists("age18_regional_result") & exists("age45_regional_result")) {
    format_region_result(age18_regional_result, who_region_with_pop, .out_folder, "_age18", se_ci)
    format_region_result(age45_regional_result, who_region_with_pop, .out_folder, "_age45", se_ci)
  }

##############################
# Results: formatted baseline means with N   -------
##############################
 has_ci <- all(c("Low (prev)", "High (prev)") %in% names(region_result_fmt))
 if (!file.exists(file.path(.out_folder,  paste0("baseline_rates_with_n_fmt", "_ci_", has_ci,".xlsx")))) {
   baseline_fmt <- region_result_fmt  |>
     select(any_of(c("Region", "Condition", "Baseline prevalence", "Low (prev)", "High (prev)")))  |>
     group_by(Region) |>
     filter(!Condition %in% c("Stroke", "All-cause mortality"))  |>
     mutate(
       `(95% CI)` = if (has_ci) {
         case_when(
           !is.na(`Low (prev)`) & !is.na(`High (prev)`) ~
             sprintf("(%.2f to %.2f)", `Low (prev)`, `High (prev)`),
           TRUE ~ "")
       } else {
         ""},
       N = round(`Baseline prevalence` / 100 * region_pop[tolower(as.character(Region))]) %>%
         as.integer() %>%
         format(big.mark = ",")
     ) |>
     ungroup() |>
     select(-any_of(c("Low (prev)", "High (prev)"))) |>
     rename(`Baseline prevalence (%)` = `Baseline prevalence`)

   write_xlsx(baseline_fmt, file.path(.out_folder, paste0("baseline_rates_with_n_fmt", "_ci_", has_ci,".xlsx")))
 }


 ###############'
 ## Results: country level -----
 ###############'
   if (skip_country_lvl == T) {return()}

   sim_result_iso3 <- build_sim_result_iso3(.out_folder, sim_years, access_percent,
                                            basecase_root = basecase_root, sim_hive_folder= sim_hive_folder)

   # country-level results workbook:
    wb <- createWorkbook()
    ## 1) condition-rate summaries used for region + count tabs -----

    summary_long <- sim_result_iso3 %>%
      select(iso3, sim, starts_with("rate__")) %>%
      pivot_longer(-c(iso3, sim), names_to = "outcome", values_to = "value") %>%
      mutate(outcome = sub("^rate__", "", outcome)) %>%
      pivot_wider(names_from = sim, values_from = value) %>%
      mutate(
        base_pct     = round(100 * .data[[base_label]], 3),
        int_pct      = round(100 * .data[[comp_label]], 3),
        abs_diff_pct = round(100 * (.data[[comp_label]] - .data[[base_label]]), 3),
        rel_diff_pct = round(
          if_else(.data[[base_label]] == 0, NA_real_,
                  100 * (.data[[comp_label]] - .data[[base_label]]) / .data[[base_label]]),
          3
        )
      ) %>%
      arrange(iso3, outcome)

    summary_wide <- summary_long %>%
      select(iso3, outcome, base_pct, int_pct, abs_diff_pct, rel_diff_pct) %>%
      pivot_wider(names_from = outcome,
                  values_from = c(base_pct, int_pct, abs_diff_pct, rel_diff_pct), names_sep = "__")

   sum_wide_bypop <- summary_wide %>%
     left_join(un_pop,   by = "iso3") %>%
     left_join(country_groups2, by = "iso3") %>%
     mutate(
       across(starts_with("abs_diff_pct__"),
              ~ (.x / 100) * pop_size, .names = "abs_diff_count__{.col}"
       )
     ) %>%
     rename_with(~ sub("^abs_diff_count__abs_diff_pct__", "abs_diff_count__", .x),
                 starts_with("abs_diff_count__abs_diff_pct__"))
   #add sheet
    addWorksheet(wb, "raw_summary_country_diffs")
     writeData(wb, "raw_summary_country_diffs", sum_wide_bypop)

   ## region totals ordered by top obesity reduction -------
   region_totals_all <- bind_rows(lapply(conditions_list_outcomes, function(o) {
     col <- paste0("abs_diff_count__", o)
     sum_wide_bypop %>%
       group_by(who_region) %>%
       summarise(lives_averted = sum(pmax(0, - .data[[col]]), na.rm = TRUE), .groups = "drop") %>%
       mutate(outcome = o)
   })) %>%
     pivot_wider(names_from = outcome, values_from = lives_averted,
                 names_prefix = "lives_averted__") %>%
     arrange(desc(`lives_averted__obesity`))

    #add sheet
   addWorksheet(wb, "region_cond_count_tab")
   writeData(wb, "region_cond_count_tab", region_totals_all)

   # top 10 countries for conditions -----
   top10_helper <- function(df, outcome, group_var = NULL, n = 10) {
     count_col <- paste0("abs_diff_count__", outcome)
     pct_col   <- paste0("abs_diff_pct__",   outcome)

     if (!is.null(group_var)) {
       df <- df %>% group_by(.data[[group_var]])
     }

     res <- df %>%
       slice_min(order_by = .data[[count_col]], n = n, with_ties = FALSE) %>%
       ungroup()

     if (!is.null(group_var)) {
       res <- res %>% arrange(.data[[group_var]], .data[[count_col]])
     } else {
       res <- res %>% arrange(.data[[count_col]])}

     cols <- unique(c(group_var, "iso3", "who_region", "wb_income", "pop_size", count_col, pct_col))
     res %>% select(all_of(cols))
   }

   outcomes <- c("obesity", "t2d_new", "ac_death")
   groups   <- list(overall = NULL, by_region = "who_region", by_income = "wb_income")

   for (oc in outcomes) {
     for (nm in names(groups)) {
       gv  <- groups[[nm]]
       res <- top10_helper(sum_wide_bypop, oc, group_var = gv)

       sheet_name <- paste0("top10_", nm, "_", oc)
       addWorksheet(wb, sheet_name)
       writeData(wb, sheet_name, res)
     }
   }

   # full country results for select outcomes
   focus_outcomes <- c("obesity", "ac_death", "ckd_new", "cvd_new", "t2d_new")

   summary_rates_long <- summary_long %>%
     filter(outcome %in% focus_outcomes) %>%
     select(iso3, outcome, base_pct, int_pct, abs_diff_pct, rel_diff_pct) %>%
     left_join(un_pop, by = "iso3") %>%
     mutate(counts_averted = round((abs_diff_pct / 100) * pop_size, 0))

   country_res_fmt <- summary_rates_long %>%
     pivot_wider(
       names_from = outcome,
       values_from = c(base_pct, int_pct, abs_diff_pct, rel_diff_pct, counts_averted),
       names_sep = "__"
     ) |>
     select(
       iso3,
       base_pct__obesity,  int_pct__obesity,  abs_diff_pct__obesity,  rel_diff_pct__obesity,  counts_averted__obesity,
       base_pct__ac_death, int_pct__ac_death, abs_diff_pct__ac_death, rel_diff_pct__ac_death, counts_averted__ac_death,
       base_pct__ckd_new,  int_pct__ckd_new,  abs_diff_pct__ckd_new,  rel_diff_pct__ckd_new,  counts_averted__ckd_new,
       base_pct__cvd_new,  int_pct__cvd_new,  abs_diff_pct__cvd_new,  rel_diff_pct__cvd_new,  counts_averted__cvd_new,
       base_pct__t2d_new,  int_pct__t2d_new,  abs_diff_pct__t2d_new,  rel_diff_pct__t2d_new,  counts_averted__t2d_new
     )

   addWorksheet(wb, "country_cond_appx")
   writeData(wb, "country_cond_appx", country_res_fmt)

   # 2) Daly/costing summary tabs -----
   make_country_compare <- function(df, outcomes, digits = 3) {
     df %>%
       filter(sim %in% c(base_label, comp_label)) %>%
       select(iso3, sim, all_of(outcomes)) %>%
       pivot_longer(
         cols = all_of(outcomes),
         names_to = "outcome",
         values_to = "value"
       ) %>%
       pivot_wider(names_from = sim, values_from = value) %>%
       mutate(
         base_val  = round(.data[[base_label]], digits),
         int_val   = round(.data[[comp_label]], digits),
         abs_diff  = round(.data[[comp_label]] - .data[[base_label]], digits),
         pct_diff  = round(
           if_else(.data[[base_label]] == 0, NA_real_,
                   100 * (.data[[comp_label]] - .data[[base_label]]) / .data[[base_label]]),
           digits
         )
       ) %>%
       select(iso3, outcome, base_val, int_val, abs_diff, pct_diff)
   }

   ##  top 10 countries by % DALY difference -----
   top10_daly_pct <- make_country_compare(
     sim_result_iso3,
     outcomes = c("daly_total_adj"),
     digits = 3
   ) %>%
     left_join(un_pop, by = "iso3") %>%
     left_join(country_groups2, by = "iso3") %>%
     arrange(pct_diff) %>%
     slice_head(n = 10) %>%
     transmute(
       iso3, who_region, wb_income, pop_size,
       daly_base = base_val,
       daly_int = int_val,
       daly_abs_diff = abs_diff,
       daly_pct_diff = pct_diff
     )

   addWorksheet(wb, "top10_pct_daly")
   writeData(wb, "top10_pct_daly", top10_daly_pct)

   ##  all countries: DALY / YLL / YLD -----
   country_daly_appx <- make_country_compare(
     sim_result_iso3,
     outcomes = c("daly_total_adj", "yll_adj", "yld_adj_total"),
     digits = 3
   ) %>%
     pivot_wider(
       names_from = outcome,
       values_from = c(base_val, int_val, abs_diff, pct_diff),
       names_sep = "__"
     ) %>%
     left_join(country_groups2 %>% select(iso3, who_region, wb_income), by = "iso3") %>%
     select(
       iso3, who_region, wb_income,
       base_val__daly_total_adj, int_val__daly_total_adj, abs_diff__daly_total_adj, pct_diff__daly_total_adj,
       base_val__yll_adj,        int_val__yll_adj,        abs_diff__yll_adj,        pct_diff__yll_adj,
       base_val__yld_adj_total,  int_val__yld_adj_total,  abs_diff__yld_adj_total,  pct_diff__yld_adj_total
     ) %>%
     arrange(iso3)

   addWorksheet(wb, "country_daly_yll_yld")
   writeData(wb, "country_daly_yll_yld", country_daly_appx)

   treat_cols <- names(sim_result_iso3)[grepl("^cost_treat_total_", names(sim_result_iso3))]
   cost_suffixes <- sub("^cost_treat_total_", "", treat_cols)

   ## treatment cost for each country -----
   country_treatment_cost_appx <- map_dfr(cost_suffixes, function(sfx) {

     outcome_nm <- paste0("cost_treat_total_", sfx)

     make_country_compare(
       sim_result_iso3,
       outcomes = outcome_nm,
       digits = 3
     ) |>
       mutate(cost_suffix = sfx)
   }) |>
     left_join(country_groups2 |> select(iso3, who_region, wb_income), by = "iso3") |>
     left_join(un_pop |> select(iso3, pop_size), by = "iso3") |>
     transmute(
       iso3, who_region, wb_income, pop_size, cost_suffix,
       treatment_cost_base = base_val,
       treatment_cost_int  = int_val,
       treatment_cost_int_per_capita = if_else(pop_size > 0, int_val / pop_size, NA_real_),
       treatment_cost_abs_diff = abs_diff,
       treatment_cost_pct_diff = pct_diff
     ) |>
     tidyr::pivot_wider(
       names_from  = cost_suffix,
       values_from = c(
         treatment_cost_base,
         treatment_cost_int,
         treatment_cost_int_per_capita,
         treatment_cost_abs_diff,
         treatment_cost_pct_diff
       ),
       names_glue = "{.value}_{cost_suffix}"
     ) |>
     arrange(iso3)

   addWorksheet(wb, "country_treatment_cost")
   writeData(wb, "country_treatment_cost", country_treatment_cost_appx)


   ## health cost + total cost for each country -----
   cond_costs <- make_country_compare(
     sim_result_iso3, outcomes = "cost_conditions_total", digits = 3
   )

   total_costs <- map_dfr(cost_suffixes, function(sfx) {

     total_outcome_nm <- paste0("cost_total_", sfx)

     make_country_compare(
       sim_result_iso3,
       outcomes = total_outcome_nm,
       digits = 3
     )
   })

   country_costs_appx <- bind_rows(cond_costs, total_costs) |>
     pivot_wider(
       names_from = outcome,
       values_from = c(base_val, int_val, abs_diff, pct_diff),
       names_sep = "__"
     ) |>
     left_join(country_groups2 |> select(iso3, who_region, wb_income), by = "iso3") |>
     left_join(un_pop |> select(iso3, pop_size), by = "iso3") |>
     arrange(iso3)

   addWorksheet(wb, "country_health_total_cost")
   writeData(wb, "country_health_total_cost", country_costs_appx)

   ## save wb -------
   saveWorkbook(wb,
     file = file.path(.out_folder, "detailed_country_region_results_tabs.xlsx"),
     overwrite = TRUE
   )

}
  

    