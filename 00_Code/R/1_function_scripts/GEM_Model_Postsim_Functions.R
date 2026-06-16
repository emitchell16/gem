# -------------------'
# Post simulation summary functions 
# Author: Liz Mitchell
# Date updated: March 2026
# -------------------'
# Define global and regional result functions ----------
make_region_impact_table <- function(sim_years, access_percent, basecase_root,
                              sim_hive_folder = wf_withdate,
                              outcome_conditions = all_outcomes_list,
                              filter_min_age = expr(age >= 0),
                              level = c("global", "all"),
                              se_ci = FALSE, psa = FALSE) {
  # inputs -----
  level <- match.arg(level)
  sim_order <- c("0_base", "end_base", "end_int")
  base_label  <- paste0(sim_years, "yr_basecase")
  comp_label  <- paste0(sim_years, "yr_", access_percent, "glp1")
  scale_factor <- 1 / config$sample_prop
  if (is.character(filter_min_age)) {filter_min_age <- parse_expr(filter_min_age)}
  need_region <- identical(level, "all")
  
  # tag for intervention summaries:
  te_tag <- case_when(
    str_detect(tolower(sim_hive_folder), "lower") ~ "lower",
    str_detect(tolower(sim_hive_folder), "upper") ~ "upper",
    str_detect(tolower(sim_hive_folder), "rwe")   ~ "rwe",
    TRUE ~ NA_character_
  )
  age_tag <- case_when(
    as_label(filter_min_age) == "age >= 18" ~ "age18up",
    as_label(filter_min_age) == "age >= 45" ~ "age45up",
    TRUE ~ NA_character_)
  tag_parts <- c(age_tag, te_tag, paste0(access_percent, "glp1"))
  int_cache_tag <- paste(tag_parts[!is.na(tag_parts)], collapse = "_")

  
  if (need_region) {
    country_groups <- read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv")) |>
      as_tibble() |> select(iso3, who_region)
  }
  outcome_conditions_final <- unique(c(outcome_conditions, "daly_total", "daly_total_adj"))
  
  # helpers -----
  summarize_one_partition <- function(ds_path, sim_type, year_val, grouping_vars = character(0)) {
    
    geo_tag <- if (length(grouping_vars) == 0) "global" else "region"
    fn_str <- if (sim_type == "int") {
      paste0("RES_temp_", sim_type, "_", int_cache_tag, "_", year_val, "_", geo_tag, ".rds")
    } else {
      paste0("RES_temp_", sim_type, "_", age_tag, "_", year_val, "_", geo_tag, ".rds")
    }
    if (psa | !file.exists(file.path(local_working_folder, fn_str))) { 
      ds_part <- open_dataset(normalizePath(ds_path), format = "parquet") |>
        select(any_of(c("iso3", "age", outcome_conditions))) |>
        filter(!!filter_min_age)
      if (need_region) {ds_part <- ds_part |> left_join(country_groups, by = "iso3")}
      
      out <- ds_part %>%
        {if (length(grouping_vars) > 0) {
            group_by(., across(all_of(grouping_vars)))
          } else {.}
        } %>%
        summarize(sample_n = n(),
                  across(
                    all_of(intersect(conditions_list_outcomes, outcome_conditions)),
                    ~ mean(.x, na.rm = TRUE),
                    .names = "{.col}"),
                  across(
                    all_of(intersect(daly_list_outcomes, outcome_conditions)),
                    ~ sum(.x, na.rm = TRUE) * scale_factor,
                    .names = "{.col}"),
                  across(
                    all_of(intersect(cost_list_outcomes, outcome_conditions)),
                    ~ sum(.x, na.rm = TRUE) * scale_factor,
                    .names = "{.col}"),
                  .groups = "drop"
                  ) %>%
        collect() |>
        mutate(
          sim_type = sim_type,
          year = year_val)
      
      rm(ds_part);gc(verbose = FALSE)
      saveRDS(out, file.path(local_working_folder, fn_str))
      } 
    out <- readRDS(file.path(local_working_folder, fn_str)) |>  as_tibble()
    }
  
  summarize_from_partition_summaries <- function(summary_tbl, grouping_vars = character(0)) {
    
    if (!("yll" %in% names(summary_tbl))) summary_tbl$yll <- NA_real_
    if (!("yld_total" %in% names(summary_tbl))) summary_tbl$yld_total <- NA_real_
    if (!("yll_adj" %in% names(summary_tbl))) summary_tbl$yll_adj <- NA_real_
    if (!("yld_adj_total" %in% names(summary_tbl))) summary_tbl$yld_adj_total <- NA_real_
    
    summary_tbl <- summary_tbl |>
      mutate(
        daly_total = yll + yld_total,
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
          across(
            matches("^(outcome_value_|rel_|abs_)") & 
              !matches("daly|yll|yld|cost"),  # exclude counts 
            ~ round(.x * 100, 2)
          )
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
  cat("***global done...")
  
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
    
    cat("***region done...")
  } else {
    result_main <- result_global
  }
  rm(result_global)
  gc(verbose = FALSE)
  
  # baseline demographics and optional 95% CI -----
  if (!psa) {
    summarize_baseline_one <- function(ds_path, grouping_vars = character(0)) {
      ds_part <- open_dataset(normalizePath(ds_path), format = "parquet") |>
        select(any_of(c("iso3", "age", "treated", "female", "ckd", "cvd", "t2d", "eskd"))) |>
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
          sample_n = n(),
          across(
            all_of(c("age", "treated", "female", "ckd", "cvd", "t2d", "eskd")),
            ~ mean(.x, na.rm = TRUE),
            .names = "{.col}"
          ),
          .groups = "drop"
        ) %>%
        collect()
      
      rm(ds_part)
      gc(verbose = FALSE)
      
      if (length(grouping_vars) == 0) {
        out <- out %>%
          mutate(who_region = "global")
      }
      
      out <- out %>%
        mutate(
          age = round(age, 2),
          across(c(treated, female, ckd, cvd, t2d, eskd), ~ round(.x * 100, 2))
        ) %>%
        rename_with(
          ~ paste0(.x, "_baseline"),
          c("treated", "ckd", "cvd", "t2d", "eskd")
        ) %>%
        pivot_longer(
          cols = c(age, treated_baseline, female, ckd_baseline, cvd_baseline, t2d_baseline, eskd_baseline),
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
      `Year 5 Value, Status Quo` = outcome_value_end_base,
      `Year 5 Value, Intervention` = outcome_value_end_int,
      `Absolute difference (percentage points)` = abs_diff_pct,
      `Relative difference (% change)` = rel_diff_pct,
      `Baseline prevalence` = outcome_value_0_base
    )
  round_cols <- c("Year 5 Value, Status Quo", "Year 5 Value, Intervention", "Relative difference (% change)", 
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
    ) #TODO update with costing/daly
  region_df$Condition <- recode(region_df$Condition, !!!mapping)
  
  # reorder
  col_order <- c(
      "Region", "Condition",
      "Year 5 Value, Status Quo", "Year 5 Value, Intervention",
      "Absolute difference (percentage points)",
      "Relative difference (% change)",
      "Total impacted lives (millions)",
      "Baseline prevalence"
    )
  
  if (se_ci) {
    col_order <- c(
      "Region", "Condition",
      "Year 5 Value, Status Quo", "Year 5 Value, Intervention",
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

build_sim_result_iso3 <- function(.out_folder, sim_years, access_percent, 
                                  basecase_root, sim_hive_folder,  overwrite = FALSE) {
  out_fp <- file.path(.out_folder, paste0("sim_result_iso3_",access_percent ,"glp1.xlsx"))
  
  if (file.exists(out_fp) && !overwrite) {
    message("Reading existing sim_result_iso3.xlsx")
    return(read_xlsx(out_fp))
  }
  
  scale_factor <- 1 / config$sample_prop
  base_label <- paste0(sim_years, "yr_basecase")
  comp_label <- paste0(sim_years, "yr_", access_percent, "glp1")
  
  sim_info <- tribble(
    ~sim,        ~sim_path,
    base_label,  file.path(basecase_root, paste0("sim=", base_label), paste0("year=", sim_years)),
    comp_label,  file.path(sim_hive_folder, paste0("sim=", comp_label), paste0("year=", sim_years))
  )
  
  summarize_one_sim_iso3 <- function(ds_path, sim_label) {
    message(" Summarizing: ", sim_label)
    
    ds_part <- open_dataset(normalizePath(ds_path), format = "parquet") |>
      select(any_of(c("iso3",
                      conditions_list_outcomes,
                      daly_list_outcomes,
                      cost_list_outcomes)))
    out <- ds_part |>
      group_by(iso3) |>
      summarize(
        n = n(),
        across(all_of(conditions_list_outcomes),
          ~ sum(if_else(.x == 1, 1L, 0L, missing = 0L), na.rm = TRUE),
          .names = "count__{.col}"),
        across(all_of(conditions_list_outcomes),
          ~ mean(.x, na.rm = TRUE),
          .names = "rate__{.col}" ),
        across(all_of(daly_list_outcomes),
          ~ sum(.x, na.rm = TRUE) * scale_factor,
          .names = "{.col}"),
        across(all_of(cost_list_outcomes),
          ~ sum(.x, na.rm = TRUE) * scale_factor,
          .names = "{.col}" ),
        .groups = "drop") |>
      collect() |>
      mutate(
        sim = sim_label,
        .before = 1)
    
    out <- out |>
      mutate(
        daly_total = yll + yld_total,
        daly_total_adj = yll_adj + yld_adj_total)
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
