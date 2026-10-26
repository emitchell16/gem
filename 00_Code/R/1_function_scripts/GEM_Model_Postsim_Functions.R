# -------------------'
# R/1_function_scripts/GEM_Model_Postsim_Functions.R
# Post-simulation summary functions (core): weighted outcomes by arm for global, WHO region and country
# Author: Liz Mitchell
# survey weights (w = pop / n); prevalence outcomes among living only; DALYs discounted (yll_adj + yld_adj_total)
# arms identified by label (scenario_label(), base_label_of())
# -------------------'

# weighted summaries on an arrow dataset / data frame (grouped or not) ----
#  prevalence outcomes: among living; cumulative outcomes: over starting cohort; counts: weighted sums
summarize_weighted <- function(ds, outcome_conditions) {
  outcome_conditions <- intersect(outcome_conditions, names(ds))            # only columns present
  prev_o <- intersect(conditions_prev_outcomes, outcome_conditions)
  cum_o  <- intersect(conditions_cum_outcomes,  outcome_conditions)
  sum_o  <- intersect(c(event_count_outcomes, daly_list_outcomes, cost_list_outcomes,
                        grep("^cost_(treat_)?total_", outcome_conditions, value = TRUE)), outcome_conditions)
  ds |>
    mutate(alive__ = if_else(ac_death == 0, 1, 0)) |>
    summarize(sample_n = n(), pop_w = sum(w), pop_alive_w = sum(w * alive__),
              across(all_of(prev_o), ~ sum(w * alive__ * .x, na.rm = TRUE), .names = "{.col}__num"),
              across(all_of(cum_o),  ~ sum(w * .x, na.rm = TRUE), .names = "{.col}__num"),
              across(all_of(sum_o),  ~ sum(w * .x, na.rm = TRUE), .names = "{.col}"),
              .groups = "drop") |>
    collect() |>
    mutate(across(all_of(paste0(prev_o, "__num")), ~ .x / pop_alive_w),
           across(all_of(paste0(cum_o,  "__num")), ~ .x / pop_w)) |>
    rename_with(~ sub("__num$", "", .x), ends_with("__num"))
}

# Define global and regional result functions ----------
make_region_impact_table <- function(sim_years, scen_label, basecase_root,
                              sim_hive_folder = basecase_root,
                              outcome_conditions = all_outcomes_list,
                              filter_min_age = expr(age >= 0),
                              level = c("global", "all"),
                              se_ci = FALSE, psa = FALSE,
                              base_label = paste0(sim_years, "yr_basecase")) {
  # inputs -----
  level <- match.arg(level)
  sim_order <- c("0_base", "end_base", "end_int")
  comp_label  <- scen_label
  if (is.character(filter_min_age)) {filter_min_age <- parse_expr(filter_min_age)}
  need_region <- identical(level, "all")
  
  # cache of partition summaries inside the run folder (keyed by arm label and age filter)
  cache_dir <- file.path(sim_hive_folder, "_postsim_cache"); dir.create(cache_dir, showWarnings = FALSE)
  age_tag <- gsub("[^0-9a-z]", "", tolower(as_label(filter_min_age)))

  
  if (need_region) {
    country_groups <- read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")) |>
      as_tibble() |> select(iso3, who_region)
  }
  outcome_conditions_final <- unique(c(outcome_conditions, "daly_total_adj"))
  
  # helpers -----
  summarize_one_partition <- function(ds_path, sim_type, year_val, grouping_vars = character(0)) {
    
    geo_tag <- if (length(grouping_vars) == 0) "global" else "region"
    arm <- if (sim_type == "int") comp_label else base_label
    fn_str <- file.path(cache_dir, paste0("RES_", arm, "_", age_tag, "_", year_val, "_", geo_tag, ".rds"))
    if (psa | !file.exists(fn_str)) { 
      ds_part <- open_dataset(normalizePath(ds_path), format = "parquet") |>
        select(any_of(c("iso3", "age", "w", "ac_death", outcome_conditions))) |>
        filter(!!filter_min_age)
      if (need_region) {ds_part <- ds_part |> left_join(country_groups, by = "iso3")}
      if (length(grouping_vars) > 0) ds_part <- ds_part |> group_by(across(all_of(grouping_vars)))
      out <- summarize_weighted(ds_part, outcome_conditions) |>
        mutate(sim_type = sim_type, year = year_val)
      
      rm(ds_part);gc(verbose = FALSE)
      saveRDS(out, fn_str)
      } 
    out <- readRDS(fn_str) |>  as_tibble()
    }
  
  summarize_from_partition_summaries <- function(summary_tbl, grouping_vars = character(0)) {
    
    if (!("yll_adj" %in% names(summary_tbl))) summary_tbl$yll_adj <- NA_real_
    if (!("yld_adj_total" %in% names(summary_tbl))) summary_tbl$yld_adj_total <- NA_real_
    
    summary_tbl <- summary_tbl |>
      mutate(
        daly_total_adj = yll_adj + yld_adj_total,
        sim_label = case_when(
          year == 0 ~ paste0("0_", sim_type),
          year == sim_years ~ paste0("end_", sim_type),
          TRUE ~ NA_character_
        )
      )
    
    if (length(grouping_vars) == 0) {
      summary_tbl <- summary_tbl |>
        mutate(who_region = "global")
      grouping_vars <- "who_region"
    }
    
    pivot_id_cols <- c(grouping_vars, "condition", "sim_label", "sample_n", "outcome_value")
    
    pivoted_temp <- summary_tbl |>
      pivot_longer(
        cols = any_of(outcome_conditions_final),
        names_to = "condition",
        values_to = "outcome_value"
      ) |>
      select(all_of(pivot_id_cols)) |>
      pivot_wider(
        names_from = sim_label,
        values_from = c(sample_n, outcome_value),
        names_glue = "{.value}_{sim_label}"
      )
    
    if (!se_ci) {
      result_main <- pivoted_temp |>
        mutate(
          abs_diff_pct = outcome_value_end_int - outcome_value_end_base,
          rel_diff_pct = abs_diff_pct / outcome_value_end_base
        ) |>
        mutate(
          across(matches("^(outcome_value_|abs_)"),        # proportions -> %; counts (DALYs, costs, events) unchanged
                 ~ if_else(condition %in% conditions_list_outcomes, round(.x * 100, 2), .x)),
          rel_diff_pct = round(rel_diff_pct * 100, 2)
        )|>
        select(
          condition, who_region,
          any_of(paste0("outcome_value_", sim_order)),
          abs_diff_pct, rel_diff_pct
        )
    } else {
      result_main <- pivoted_temp |>
        mutate(
          se_0_base   = sqrt(outcome_value_0_base * (1 - outcome_value_0_base) / sample_n_0_base),
          se_end_base = sqrt(outcome_value_end_base * (1 - outcome_value_end_base) / sample_n_end_base),
          se_end_int  = sqrt(outcome_value_end_int * (1 - outcome_value_end_int) / sample_n_end_int),
          
          conf_low_0_base   = outcome_value_0_base - 1.96 * se_0_base,
          conf_high_0_base  = outcome_value_0_base + 1.96 * se_0_base,
          conf_low_end_base = outcome_value_end_base - 1.96 * se_end_base,
          conf_high_end_base = outcome_value_end_base + 1.96 * se_end_base,
          conf_low_end_int  = outcome_value_end_int - 1.96 * se_end_int,
          conf_high_end_int = outcome_value_end_int + 1.96 * se_end_int,
          
          abs_diff_pct = outcome_value_end_int - outcome_value_end_base,
          se_abs_diff = sqrt(se_end_base^2 + se_end_int^2),
          abs_conf_low = abs_diff_pct - 1.96 * se_abs_diff,
          abs_conf_high = abs_diff_pct + 1.96 * se_abs_diff,
          
          rel_diff_pct = abs_diff_pct / outcome_value_end_base,
          se_rel_diff = se_abs_diff / outcome_value_end_base,
          rel_conf_low = rel_diff_pct - 1.96 * se_rel_diff,
          rel_conf_high = rel_diff_pct + 1.96 * se_rel_diff
        ) |>
        mutate(
          across(
            matches("^(outcome_value_|se_|conf_low_|conf_high_|rel_|abs_)"),
            ~ round(.x * 100, 2)
          )
        ) |>
        select(
          condition, who_region,
          any_of(unlist(lapply(sim_order, function(s) {
            c(
              paste0("outcome_value_", s),
              paste0("conf_low_", s),
              paste0("conf_high_", s)
            )
          }))),
          abs_diff_pct, abs_conf_low, abs_conf_high,
          rel_diff_pct, rel_conf_low, rel_conf_high
        )
    }
    result_main
  }
  
  # partition paths ----
  part_info <- tribble(
    ~sim_type, ~year,      ~path,
    "base",    0,          file.path(basecase_root, paste0("sim=", base_label), "year=0"),
    "base",    sim_years,  file.path(basecase_root, paste0("sim=", base_label), paste0("year=", sim_years)),
    "int",     0,          file.path(sim_hive_folder, paste0("sim=", comp_label), "year=0"),
    "int",     sim_years,  file.path(sim_hive_folder, paste0("sim=", comp_label), paste0("year=", sim_years))
  )
  
  # chunk summaries -----
  global_parts <- vector("list", nrow(part_info))
  for (i in seq_len(nrow(part_info))) {
    global_parts[[i]] <- summarize_one_partition(
      ds_path = part_info$path[i],
      sim_type = part_info$sim_type[i],
      year_val = part_info$year[i],
      grouping_vars = character(0)
    )
    gc(verbose = FALSE)
  }
  summary_global <- bind_rows(global_parts)
  rm(global_parts)
  gc(verbose = FALSE)
  
  result_global <- summarize_from_partition_summaries(summary_global, character(0))
  rm(summary_global)
  gc(verbose = FALSE)
  cat("global done...")
  
  if (need_region) {
    region_parts <- vector("list", nrow(part_info))
    for (i in seq_len(nrow(part_info))) {
      region_parts[[i]] <- summarize_one_partition(
        ds_path = part_info$path[i],
        sim_type = part_info$sim_type[i],
        year_val = part_info$year[i],
        grouping_vars = "who_region"
      )
      gc(verbose = FALSE)
    }
    summary_region <- bind_rows(region_parts)
    rm(region_parts)
    gc(verbose = FALSE)
    
    result_region <- summarize_from_partition_summaries(summary_region, "who_region")
    rm(summary_region)
    gc(verbose = FALSE)
    
    result_main <- bind_rows(result_global, result_region)
    rm(result_region)
    
    cat("region done...")
  } else {
    result_main <- result_global
  }
  rm(result_global)
  gc(verbose = FALSE)
  
  # baseline demographics and optional 95% CI -----
  if (!psa) {
    summarize_baseline_one <- function(ds_path, grouping_vars = character(0)) {
      ds_part <- open_dataset(normalizePath(ds_path), format = "parquet") |>
        select(any_of(c("iso3", "age", "w", "female", "ckd", "cvd", "t2d", "eskd"))) |>
        filter(!!filter_min_age)
      
      if (need_region) {
        ds_part <- ds_part |> left_join(country_groups, by = "iso3")
      }
      out <- ds_part %>%
        {
          if (length(grouping_vars) > 0) {
            group_by(., across(all_of(grouping_vars)))
          } else {
            .
          }
        } %>%
        summarize(
          sample_n = n(), w__ = sum(w),
          across(
            all_of(c("age", "female", "ckd", "cvd", "t2d", "eskd")),
            ~ sum(w * .x, na.rm = TRUE),
            .names = "{.col}"
          ),
          .groups = "drop"
        ) %>%
        collect() |>
        mutate(across(all_of(c("age", "female", "ckd", "cvd", "t2d", "eskd")), ~ .x / w__)) |>
        select(-w__)
      
      rm(ds_part)
      gc(verbose = FALSE)
      
      if (length(grouping_vars) == 0) {
        out <- out %>%
          mutate(who_region = "global")
      }
      
      out <- out %>%
        mutate(
          age = round(age, 2),
          across(c(female, ckd, cvd, t2d, eskd), ~ round(.x * 100, 2))
        ) %>%
        rename_with(
          ~ paste0(.x, "_baseline"),
          c("ckd", "cvd", "t2d", "eskd")
        ) %>%
        pivot_longer(
          cols = c(age, female, ckd_baseline, cvd_baseline, t2d_baseline, eskd_baseline),
          names_to = "condition",
          values_to = "outcome_value_0_base"
        )
      
      if (se_ci) {
        out <- out %>%
          mutate(
            se = case_when(
              condition %in% c("ckd_baseline", "cvd_baseline", "t2d_baseline", "eskd_baseline") ~
                sqrt((outcome_value_0_base / 100) * (1 - outcome_value_0_base / 100) / sample_n) * 100,
              TRUE ~ NA_real_
            ),
            conf_low_0_base = round(outcome_value_0_base - 1.96 * se, 2),
            conf_high_0_base = round(outcome_value_0_base + 1.96 * se, 2)
          ) %>%
          select(who_region, condition, outcome_value_0_base, conf_low_0_base, conf_high_0_base)
      } else {
        out <- out %>%
          select(who_region, condition, outcome_value_0_base)
      }
      as_tibble(out)
    }
    
    baseline_global <- summarize_baseline_one(
      ds_path = file.path(sim_hive_folder, paste0("sim=", comp_label), "year=0"),
      grouping_vars = character(0)
    )
    
    if (need_region) {
      baseline_region <- summarize_baseline_one(
        ds_path = file.path(sim_hive_folder, paste0("sim=", comp_label), "year=0"),
        grouping_vars = "who_region"
      )
      baseline_demog <- bind_rows(baseline_global, baseline_region)
      rm(baseline_region)
    } else {
      baseline_demog <- baseline_global
    }
    result_main <- bind_rows(result_main, baseline_demog)
    rm(baseline_demog, baseline_global)
    gc(verbose = FALSE)
  }
  
  result_main |>
    arrange(desc(who_region == "global"), condition, who_region)
}

######'
## format: -------------
format_region_result <- function(df, who_region_with_pop, .out_folder, suffix = "", se_ci) {
  if (se_ci) { 
    base_cols <- c("outcome_value_0_base","conf_low_0_base","conf_high_0_base")    
  } else { 
    base_cols <- c("outcome_value_0_base")}

  temp_baseline <- df  |> 
    filter(str_ends(condition, "_baseline"))  |> 
    mutate(condition = str_remove(condition, "_baseline")) 
  
  # merge baseline prevalence for irreversible conds
  region_df <- df  |> 
    filter(condition != "t2d") |> 
    mutate(condition = str_remove(condition, "_new")) |> 
    left_join(temp_baseline %>% select(who_region, condition, all_of(base_cols)),
              by = c("who_region","condition"), suffix = c("", "_cc"))
  repl_cond <- function(x) if_else(is.na(x) | x == 0, TRUE, FALSE, missing = TRUE) 
  for (col in base_cols) {
    cc <- paste0(col, "_cc")
    region_df[[col]] <- if_else(
      repl_cond(region_df[[col]]),
      region_df[[cc]],
      region_df[[col]]
    )
    region_df[[cc]] <- NULL 
  }
  region_df <- region_df |>  select(-ends_with("_cc"))

  region_df <- region_df |> 
    filter(!condition %in% c("female","age"), !str_ends(condition, "_baseline"))
  
  # add population and impacted lives
  region_df <- region_df |>
    mutate(
      abs_diff_frac = if_else(condition %in% conditions_list_outcomes, abs_diff_pct / 100, NA_real_),
      pop = who_region_with_pop$pop_size[match(who_region, who_region_with_pop$who_region)],
      impacted_millions = if_else(!is.na(abs_diff_frac),
                                  round(pop * abs_diff_frac / 1e6, 2),
                                  NA_real_),
      `Total impacted lives (millions)` = if_else(!is.na(impacted_millions),
                                                  format(impacted_millions, nsmall = 2),
                                                  NA_character_)
    ) |>
    select(-abs_diff_frac, -impacted_millions)
  # rename cols
  region_df <- region_df |> 
    mutate(who_region = str_to_title(who_region)) |> 
    rename(
      Region = who_region,
      Condition = condition,
      `End-of-horizon value, status quo` = outcome_value_end_base,
      `End-of-horizon value, intervention` = outcome_value_end_int,
      `Absolute difference (percentage points)` = abs_diff_pct,
      `Relative difference (% change)` = rel_diff_pct,
      `Baseline prevalence` = outcome_value_0_base
    )
  round_cols <- c("End-of-horizon value, status quo", "End-of-horizon value, intervention", "Relative difference (% change)", 
                  "Absolute difference (percentage points)" )
  if (se_ci) {
    region_df <- region_df |> 
      rename(
        `Low (AC)` = abs_conf_low, `High (AC)` = abs_conf_high,
        `Low (RC)` = rel_conf_low, `High (RC)` = rel_conf_high,
        `Low (prev)` = conf_low_0_base, `High (prev)` = conf_high_0_base
      )
    round_cols <- c(round_cols, "Low (RC)", "High (RC)", "Low (AC)", "High (AC)" )
  }
  
  region_df <- region_df %>% mutate(across(all_of(round_cols), ~ round(as.numeric(.x), 2)))
  
  # mapping of condition labels
  mapping <- c(
    ac_death = "All-cause mortality",
    class1_obesity = "Obesity class 1",
    class2_obesity = "Obesity class 2",
    class3_obesity = "Obesity class 3",
    obesity = "Obesity", eskd = "ESKD",
    cvd = "CVD", stroke = "Stroke",
    ckd = "CKD", t2d = "T2D"
    )
  region_df$Condition <- recode(region_df$Condition, !!!mapping)
  
  # reorder
  col_order <- c(
      "Region", "Condition",
      "End-of-horizon value, status quo", "End-of-horizon value, intervention",
      "Absolute difference (percentage points)",
      "Relative difference (% change)",
      "Total impacted lives (millions)",
      "Baseline prevalence"
    )
  
  if (se_ci) {
    col_order <- c(
      "Region", "Condition",
      "End-of-horizon value, status quo", "End-of-horizon value, intervention",
      "Absolute difference (percentage points)", "Low (AC)", "High (AC)",
      "Relative difference (% change)", "Low (RC)", "High (RC)",
      "Total impacted lives (millions)",
      "Baseline prevalence", "Low (prev)", "High (prev)"
    )
  }
  region_df <- region_df |>  select(any_of(col_order)) |> 
    arrange(Condition, Region)

  write_xlsx(region_df, file.path(.out_folder, paste0("region_result_fmt", suffix, ".xlsx")))
  return(region_df)
}

# Define country level result function ----------

build_sim_result_iso3 <- function(.out_folder, sim_years, scen_label,
                                  basecase_root, sim_hive_folder,  overwrite = FALSE,
                                  base_label = paste0(sim_years, "yr_basecase")) {
  out_fp <- file.path(.out_folder, paste0("sim_result_iso3_", scen_label, ".xlsx"))
  
  if (file.exists(out_fp) && !overwrite) {
    message("Reading existing sim_result_iso3.xlsx")
    return(read_xlsx(out_fp))
  }
  
  comp_label <- scen_label
  
  sim_info <- tribble(
    ~sim,        ~sim_path,
    base_label,  file.path(basecase_root, paste0("sim=", base_label), paste0("year=", sim_years)),
    comp_label,  file.path(sim_hive_folder, paste0("sim=", comp_label), paste0("year=", sim_years))
  )
  
  # weighted: rate__ = prevalence among living / cumulative over cohort; count__ = weighted persons
  summarize_one_sim_iso3 <- function(ds_path, sim_label) {
    message(" Summarizing: ", sim_label)
    ds_part <- open_dataset(normalizePath(ds_path), format = "parquet") |>
      select(any_of(c("iso3", "w", "ac_death", all_outcomes_list,
                      grep("^cost_(treat_)?total_", names(open_dataset(normalizePath(ds_path))), value = TRUE))))
    out <- summarize_weighted(ds_part |> group_by(iso3), c(all_outcomes_list,
                              grep("^cost_(treat_)?total_", names(ds_part), value = TRUE))) |>
      rename_with(~ paste0("rate__", .x), any_of(conditions_list_outcomes)) |>
      mutate(across(starts_with("rate__"), ~ .x * (if (sub("^rate__", "", cur_column()) %in% conditions_prev_outcomes)   # scalar test; vectors by country
                                                          pop_alive_w else pop_w), .names = "{sub('rate__', 'count__', .col)}"),
             daly_total_adj = yll_adj + yld_adj_total,
             sim = sim_label, .before = 1) |>
      rename(n = sample_n)
    rm(ds_part); gc(verbose = FALSE)
    out
  }
  
  sim_result_list <- vector("list", nrow(sim_info))
  
  for (i in seq_len(nrow(sim_info))) {
    sim_result_list[[i]] <- summarize_one_sim_iso3(
      ds_path = sim_info$sim_path[i],
      sim_label = sim_info$sim[i]
    )
    gc(verbose = FALSE)
  }
  
  sim_result_iso3 <- bind_rows(sim_result_list)
  rm(sim_result_list)
  gc(verbose = FALSE)
  
  write_xlsx(sim_result_iso3, out_fp)
  
  sim_result_iso3
}
