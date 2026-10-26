# R/1_function_scripts/GEM_Model_1_Baseline_Summary_Stats.R
##############
# Code to compute global summary statistics for year 0 (weighted stratified sample; w = pop / n; adults 18+)
# Uses arrow to read over individual cohort files
##############

global_summary_stats <- function(config, date, paths, ds_path,
                                 global_summary_conds = c("female", "obesity", "class1_obesity", "class2_obesity",
                                                          "class3_obesity", "cvd", "ckd", "eskd", "t2d")) {

  ds <- open_dataset(ds_path, format = "parquet")
    # summarize without loading the full dataset into memory
  need_cols <- unique(c("w", "age", "female", "t2d", "bmi", "cvd", "ckd", "obesity", global_summary_conds))
   summary_stats <- ds |>
     select(all_of(need_cols)) |>
      summarise(
        sample_n = n(),
        people   = sum(w),
        age_mean = sum(w * age) / sum(w),
        age_m2   = sum(w * age^2) / sum(w),
        bmi_mean = sum(w * bmi) / sum(w),
        across(all_of(global_summary_conds), ~ sum(w * .x) / sum(w), .names = "p_{col}")
      ) |>
      collect() |>
      mutate(age_sd = sqrt(age_m2 - age_mean^2))

    N <- format(round(summary_stats$people), big.mark = ",")

    age_row <- tibble(
      characteristic = c("Age, mean (SD)", "BMI, mean"),
      value = c(sprintf("%0.2f (%0.2f)", summary_stats$age_mean, summary_stats$age_sd),
                sprintf("%0.2f", summary_stats$bmi_mean)))

    pct_rows <- summary_stats |>
      select(starts_with("p_")) |>
      mutate(p_male = 1 - p_female) |>
      pivot_longer(everything(), names_to = "characteristic", values_to = "value") |>
      mutate(
        characteristic = str_to_title(str_remove(characteristic, "^p_")),
        value          = sprintf("%0.2f", value * 100))

    tbl <- bind_rows(age_row, pct_rows)
    colnames(tbl)[2] <- paste0("Global population, weighted (N=", N, "; sample n=",
                               format(summary_stats$sample_n, big.mark = ","), ")")

  fwrite(tbl, file = file.path(paths$out, paste0("summary_stats_y0fmt_", date, ".csv")))
}
