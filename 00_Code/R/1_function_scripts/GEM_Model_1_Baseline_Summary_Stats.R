##############
# Code to compute global summary statistics for year 0 (full population, no random sampling)
# Uses arrow to read over individual cohort files
##############

global_summary_stats <- function(config, date) {
  
  ds <- open_dataset(file.path(paths$working, paste0("baseline_yr0_segment_files_", date)), format = "parquet")
    # summarize without loading the full dataset into memory
  need_cols <- unique(c(
    "age", "female", "t2d", "bmi", "obesity",
    global_summary_conds
  )) 
   summary_stats <- ds |>
     select(all_of(need_cols)) |>
      summarise(
        people = n(),
        age_mean = mean(age, na.rm = T),
        age_sd   = sqrt(var(age, na.rm = T)),
        across(all_of(global_summary_conds),
               ~ mean(.x, na.rm = TRUE),
               .names = "p_{col}"),
        eligible_n = sum(
          if_else(
            (age >= 18 & t2d == 1 & bmi >= 27) | (age >= 12 & obesity == 1),
            1L, 0L
          ),
          na.rm = T
        ),
        p_treatment_eligible = eligible_n / people
      ) |>
      collect() 
    
    N <- format(summary_stats$people, big.mark = ",")
    
    age_row <- tibble(
      characteristic = "Age, mean (SD)",
      value = sprintf("%0.2f (%0.2f)",
                      summary_stats$age_mean,
                      summary_stats$age_sd))
    
    pct_rows <- summary_stats |> 
      select(starts_with("p_")) |> 
      mutate(p_male = 1 - p_female) |> 
      pivot_longer(
        everything(), names_to  = "characteristic", values_to = "value"
      ) |> 
      mutate(
        characteristic = str_to_title(str_remove(characteristic, "^p_")),
        value          = sprintf("%0.2f", value * 100)
        )
    
    tbl <- bind_rows(age_row, pct_rows)
    
    colnames(tbl)[2] <- paste0("Global Sample (N=", N, ")")

  # Save the summary table to CSV
  fwrite(tbl, file = file.path(paths$out, y0_summary_filename))
}
