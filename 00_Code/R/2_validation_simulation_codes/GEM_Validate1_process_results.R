# -------------------'
#' Generate validation tables and figures
#' Date: Aug. 2025
#' Author: Liz Mitchell
#' # -------------------'

clamp01  <- function(x, eps = 1e-12) pmin(pmax(x, eps), 1 - eps)

process_results <- function(vres_outfile, out_subdir ) {
  validation_output_comp_data <- read_csv(file.path(validation_out_dir, vres_outfile), show_col_types = F)
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  labels_map <- c(
    eskd_new      = "ESKD", ckd_new      = "CKD",
    cvd_new       = "CVD", t2d_new       = "T2D",
    ac_death  = "All-Cause Mortality", mace      = "MACE",  stroke = "Stroke")
   
  # Make regression table ---------
  
  ## label cleaning functions: 
  nice_outcome <- function(x) {
    lab <- labels_map[x]
    ifelse(!is.na(lab),
           lab,
           x %>%
             gsub("_", " ", ., fixed = TRUE) %>%
             stringr::str_to_title())
  }
  
  nice_measure <- function(x) case_when(
    x == "incidence_rate_1py"        ~ "1-year incidence",
    x == "incidence_rate_cumulative" ~ "Cumulative incidence",
    TRUE ~ x
  )
  
  ## prep output data
  reg_voutput_wide <- validation_output_comp_data %>%
    pivot_wider(names_from  = test_indicator,  values_from = estimate ) %>%
    mutate(
      outcome_label = nice_outcome(outcome),
      measure       = nice_measure(measure)
    ) %>%
    left_join(
      validation_data_input %>% select(study_id, cohort_id, futime),
      by = c("study_id","cohort_id")
    ) %>%
    mutate(
      T = if_else(str_detect(measure, regex("cumulative", ignore_case = TRUE)),
                  pmax(futime, 1), 1),
      ref_1y = if_else(
        str_detect(measure, regex("cumulative", ignore_case = TRUE)),
        -log1p(-clamp01(reference_study)) / T,  # cumulative risk -> avg annual rate
        reference_study                         # already 1-year rate
      ),
      sim_1y = if_else(
        str_detect(measure, regex("cumulative", ignore_case = TRUE)),
        -log1p(-clamp01(simulated)) / T,
        simulated
      )
    ) %>%
    arrange(outcome_label, measure)
  
  reg_voutput_long <- validation_output_comp_data |> 
    left_join(
      validation_data_input %>% select(study_id, cohort_id, futime),
      by = c("study_id","cohort_id")
    ) |> 
    mutate(
      outcome_label = nice_outcome(outcome),
      measure       = nice_measure(measure),
      T             = if_else(str_detect(measure, regex("cumulative", ignore_case = TRUE)),
                              pmax(futime, 1), 1),
      estimate_1yr   = if_else(
        str_detect(measure, regex("cumulative", ignore_case = TRUE)),
        -log1p(-clamp01(estimate)) / T,   # cumulative risk -> avg annual rate
        estimate                           # already 1-year rate
      ),
      test_indicator = factor(test_indicator, levels = c("reference_study", "simulated"))
    ) |> 
    arrange(outcome_label, measure, study_id, cohort_id)
  
  # outcome-level frequency weights
  w_map <- reg_voutput_wide |>
    count(outcome_label, name = "w_outcome")
  
  # attach to wide and long frames
  reg_voutput_wide <- reg_voutput_wide |>
    left_join(w_map, by = "outcome_label")
  
  reg_voutput_long <- reg_voutput_long |>
    left_join(w_map, by = "outcome_label")
  
  ## helper function for regression output 
  tidy_term <- function(fit, term_pattern, n_obs) {
    co <- coef(summary(fit))
    ci <- suppressMessages(as.data.frame(confint(fit)))
    rn <- rownames(co)
    idx <- grep(term_pattern, rn, ignore.case = TRUE, value = TRUE)
    if (length(idx) != 1) stop("Could not uniquely match term pattern: ", term_pattern,
                               " among: ", paste(rn, collapse = ", "))
    tibble(
      term = idx,
      est  = unname(co[idx, 1]),
      lcl  = ci[idx, 1],
      ucl  = ci[idx, 2],
      p    = unname(co[idx, 4]),
      n    = n_obs
    )
  }
  
  # --- (1) Pooled regression: estimate ~ test_indicator -------------------------
  # Intercept = Reference mean; Coef(test_indicator) = Simulated - Reference
  fit_pooled <- lm(estimate ~ test_indicator,
                    data = reg_voutput_long,
                    weights = w_outcome)
  
  pooled_table <- bind_rows(
    tidy_term(fit_pooled, "Intercept", nrow(reg_voutput_long)) %>%
      mutate(Label = "Reference mean"),
    tidy_term(fit_pooled, "test_indicatorsimulated", nrow(reg_voutput_long)) %>%
      mutate(Label = "Difference (Simulated − Reference)")
  ) %>%
    transmute(
      Term     = Label,
      Estimate = round(est, 3),
      `95% CI` = sprintf("[%.3f, %.3f]", lcl, ucl),
      `P value`= signif(p, 3),
      N        = n
    )
  # --- (2) Pooled calibration: Simulated ~ Reference ----------------------------
  
  fit_pooled_calib <- lm(simulated ~ reference_study,
                         data = reg_voutput_wide,
                         weights = w_outcome)
  s  <- summary(fit_pooled_calib)
  ci <- suppressMessages(confint(fit_pooled_calib))
  
  pooled_calibration <- tibble(
    Term     = c("Intercept", "Slope"),
    Estimate = round(coef(fit_pooled_calib), 3),
    `95% CI` = c(
      sprintf("[%.3f, %.3f]", ci["(Intercept)", 1], ci["(Intercept)", 2]),
      sprintf("[%.3f, %.3f]", ci["reference_study",   1], ci["reference_study", 2])
    ),
    `P value`= signif(coef(s)[, 4], 3),
    R2       = c(round(s$r.squared, 3), NA),
    N        = nrow(reg_voutput_wide)
  )
  
  # --- (3) By-outcome regression: estimate ~ test_indicator (stratified) --------
  by_outcome_table <- reg_voutput_long %>%
    group_by(Outcome = outcome_label) %>%
    group_modify(~{
      df <- .x
      if (length(unique(df$test_indicator)) < 2) {
        return(tibble(
          Term     = "Difference (Simulated − Reference)",
          Estimate = NA_real_, `95% CI` = NA_character_,
          `P value`= NA_real_, N = nrow(df)
        ))
      }
      fit <- lm(estimate ~ test_indicator, data = df)
      row <- tidy_term(fit, "test_indicatorsimulated", nrow(df))
      tibble(
        Term     = "Difference (Simulated − Reference)",
        Estimate = round(row$est, 3),
        `95% CI` = sprintf("[%.3f, %.3f]", row$lcl, row$ucl),
        `P value`= signif(row$p, 3),
        N        = row$n
      )
    }) %>%
    ungroup() %>%
    arrange(Outcome)
  
  # --- (4) Calibration (supplement): Simulated ~ Reference ----------------------
  
  calib_df <- reg_voutput_wide %>%
    group_by(outcome_label) %>%
    group_modify(~{
      df <- .x %>% filter(!is.na(reference_study), !is.na(simulated))
      if (nrow(df) < 3 || length(unique(df$reference_study)) < 2 ) {
        return(tibble(
          intercept = NA_real_, slope = NA_real_,
          intercept_l = NA_real_, intercept_u = NA_real_,
          slope_l = NA_real_, slope_u = NA_real_,
          r_squared = NA_real_, n = nrow(df)
        ))
      }
      fit <- lm(simulated ~ reference_study, data = df)
      s   <- summary(fit)
      ci  <- suppressMessages(as.data.frame(confint(fit)))
      tibble(
        intercept   = unname(coef(fit)[1]),
        slope       = unname(coef(fit)[2]),
        intercept_l = ci["(Intercept)", 1],
        intercept_u = ci["(Intercept)", 2],
        slope_l     = ci["reference", 1],
        slope_u     = ci["reference", 2],
        r_squared   = unname(s$r.squared),
        n           = nrow(df)
      )
    }) %>%
    ungroup()
  
  calibration_table <- calib_df %>%
    transmute(
      Outcome   = outcome_label,
      Intercept = round(intercept, 3),
      `Intercept 95% CI` = sprintf("[%.3f, %.3f]", intercept_l, intercept_u),
      Slope     = round(slope, 3),
      `Slope 95% CI` = sprintf("[%.3f, %.3f]", slope_l, slope_u),
      `R²`      = round(r_squared, 3),
      N         = n
    )
  
  
  # --- (5) Calibration excluding ESKD ----------------------
  reg_voutput_wide_noESKD <- reg_voutput_wide %>% filter(outcome_label != "ESKD")
  fit_pooled_calib_noESKD <- lm(simulated ~ reference_study,
                                data = reg_voutput_wide_noESKD,
                                weights = w_outcome)
  s_no  <- summary(fit_pooled_calib_noESKD)
  ci_no <- suppressMessages(confint(fit_pooled_calib_noESKD))
  
  pooled_calibration_noESKD <- tibble(
    Term     = c("Intercept", "Slope"),
    Estimate = round(coef(fit_pooled_calib_noESKD), 3),
    `95% CI` = c(
      sprintf("[%.3f, %.3f]", ci_no["(Intercept)", 1],    ci_no["(Intercept)", 2]),
      sprintf("[%.3f, %.3f]", ci_no["reference_study", 1], ci_no["reference_study", 2])
    ),
    `P value`= signif(coef(s_no)[, 4], 3),
    R2       = c(round(s_no$r.squared, 3), NA),
    N        = nrow(reg_voutput_wide_noESKD)
  )
  
  write_xlsx(
    list(
      "PooledCalibration"      = pooled_calibration,
      "PooledSimAccuracy"      = pooled_table,
      "Calibration" = calibration_table,
      "PooledCalibration-noESKD"      = pooled_calibration_noESKD,
      "ByOutcome"   = by_outcome_table
    ),
    path = file.path(out_subdir, "et5sup2_simulationaccuracy.xlsx")
  )
###########################'  
#Plot figure ---------
###########################'    
  # regression line:
  fit_pool <- lm(simulated ~ reference_study, data = reg_voutput_wide, weights =w_outcome )
    b0 <- unname(coef(fit_pool)[1])
    b1 <- unname(coef(fit_pool)[2])
    r2 <- summary(fit_pool)$r.squared
  fit_label <- sprintf("Fit line: y = %.3f + %.3f x (R² = %.3f)", b0, b1, r2)
  
  line_df <- data.frame(
    label     = c(fit_label, "45° line"),
    slope     = c(b1, 1),
    intercept = c(b0, 0)
  )
  
  p1 <- ggplot(reg_voutput_wide, aes(
    x = reference_study, y = simulated,
    color = outcome_label, shape = measure
  )) +
    geom_point(size = 3.2) +
    geom_abline(intercept = b0, slope = b1, 
                color = "darkgrey", linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, 
                linetype = "solid", linewidth = 0.4, color = "black") +
    geom_text_repel(
      aes(label = paste0("(", study_id, "-", cohort_id, ")")),
      size = 3, min.segment.length = 0, segment.color = NA,  max.overlaps = 20,
      box.padding = 0.25,
      point.padding = 0.15
    ) +
    labs(
      x = "Reference Estimate", y = "Simulated Estimate",
      title = "Simulated vs. Real-World Outcomes (by study-cohort)",
      color = "Outcome:", shape = "Measure:"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.title.position = "plot",
      axis.line = element_line(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.background = element_rect(colour = "grey55", fill = NA),
      legend.margin = margin(6, 6, 6, 6),
      plot.caption.position = "plot"
    ) +
    scale_color_tableau(palette = "Tableau 10") +
    coord_fixed(xlim = c(0, .5), ylim = c(0, .5), clip = "on")+
    annotate("segment", x=.40, xend=.47, y=b0+b1*.40, yend=b0+b1*.47,
             linetype="dashed", colour="grey35", linewidth=.35) +
    annotate("text", x=.47, y=b0+b1*.47, label=fit_label,
             hjust=1, vjust=.5, size=3.5, colour="grey25") +
    annotate("segment", x=.40, xend=.47, y=.40, yend=.47,
             colour="black", linewidth=.35) +
    annotate("text", x=.47, y=.47, label="45° line",
             hjust=1, vjust=.5, size=3.5)

  ggsave(
    filename = "validation_outcomes_comparison_plot.png",
    path     = out_subdir,
    width    = 12, height   = 8,  units    = "in", dpi      = 300
  )
  
  return(reg_voutput_wide)
}

###################
## death only ---------

process_results_death <- function(vres_outfile, out_subdir ) {
  validation_output_comp_data <- read_csv(file.path(validation_out_dir, vres_outfile), show_col_types = F)
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  labels_map <- c(
    eskd_new      = "ESKD", ckd_new      = "CKD",
    cvd_new       = "CVD", t2d_new       = "T2D",
    ac_death  = "All-Cause Mortality", mace      = "MACE",  stroke = "Stroke")
  
  # Make regression table ---------
  
  ## label cleaning functions: 
  nice_outcome <- function(x) {
    lab <- labels_map[x]
    ifelse(!is.na(lab),
           lab,
           x %>%
             gsub("_", " ", ., fixed = TRUE) %>%
             stringr::str_to_title())
  }
  
  nice_measure <- function(x) case_when(
    x == "incidence_rate_1py"        ~ "1-year incidence",
    x == "incidence_rate_cumulative" ~ "Cumulative incidence",
    TRUE ~ x
  )
  
  ## prep output data
  reg_voutput_wide <- validation_output_comp_data %>%
    filter(outcome == "ac_death") |> 
    pivot_wider(names_from  = test_indicator,  values_from = estimate ) %>%
    mutate(
      outcome_label = nice_outcome(outcome),
      measure       = nice_measure(measure)
    ) %>%
    left_join(
      validation_data_input %>% select(study_id, cohort_id, futime),
      by = c("study_id","cohort_id")
    ) %>%
    mutate(
      T = if_else(str_detect(measure, regex("cumulative", ignore_case = TRUE)),
                  pmax(futime, 1), 1),
      ref_1y = if_else(
        str_detect(measure, regex("cumulative", ignore_case = TRUE)),
        -log1p(-clamp01(reference_study)) / T,  # cumulative risk -> avg annual rate
        reference_study                         # already 1-year rate
      ),
      sim_1y = if_else(
        str_detect(measure, regex("cumulative", ignore_case = TRUE)),
        -log1p(-clamp01(simulated)) / T,
        simulated
      )
    ) %>%
    arrange(outcome_label, measure)
  
  reg_voutput_long <- validation_output_comp_data |> 
    filter(outcome == "ac_death") |> 
    left_join(
      validation_data_input %>% select(study_id, cohort_id, futime),
      by = c("study_id","cohort_id")
    ) |> 
    mutate(
      outcome_label = nice_outcome(outcome),
      measure       = nice_measure(measure),
      T             = if_else(str_detect(measure, regex("cumulative", ignore_case = TRUE)),
                              pmax(futime, 1), 1),
      estimate_1yr   = if_else(
        str_detect(measure, regex("cumulative", ignore_case = TRUE)),
        -log1p(-clamp01(estimate)) / T,   # cumulative risk -> avg annual rate
        estimate                           # already 1-year rate
      ),
      test_indicator = factor(test_indicator, levels = c("reference_study", "simulated"))
    ) |> 
    arrange(outcome_label, measure, study_id, cohort_id)
  
  # outcome-level frequency weights
  w_map <- reg_voutput_wide |>
    count(outcome_label, name = "w_outcome")
  
  # attach to wide and long frames
  reg_voutput_wide <- reg_voutput_wide |>
    left_join(w_map, by = "outcome_label")
  
  reg_voutput_long <- reg_voutput_long |>
    left_join(w_map, by = "outcome_label")
  
  ## helper function for regression output 
  tidy_term <- function(fit, term_pattern, n_obs) {
    co <- coef(summary(fit))
    ci <- suppressMessages(as.data.frame(confint(fit)))
    rn <- rownames(co)
    idx <- grep(term_pattern, rn, ignore.case = TRUE, value = TRUE)
    if (length(idx) != 1) stop("Could not uniquely match term pattern: ", term_pattern,
                               " among: ", paste(rn, collapse = ", "))
    tibble(
      term = idx,
      est  = unname(co[idx, 1]),
      lcl  = ci[idx, 1],
      ucl  = ci[idx, 2],
      p    = unname(co[idx, 4]),
      n    = n_obs
    )
  }
  
  # --- (1) Pooled regression: estimate ~ test_indicator -------------------------
  # Intercept = Reference mean; Coef(test_indicator) = Simulated - Reference
  fit_pooled <- lm(estimate ~ test_indicator,
                   data = reg_voutput_long,
                   weights = w_outcome)
  
  pooled_table <- bind_rows(
    tidy_term(fit_pooled, "Intercept", nrow(reg_voutput_long)) %>%
      mutate(Label = "Reference mean"),
    tidy_term(fit_pooled, "test_indicatorsimulated", nrow(reg_voutput_long)) %>%
      mutate(Label = "Difference (Simulated − Reference)")
  ) %>%
    transmute(
      Term     = Label,
      Estimate = round(est, 3),
      `95% CI` = sprintf("[%.3f, %.3f]", lcl, ucl),
      `P value`= signif(p, 3),
      N        = n
    )
  # --- (2) Pooled calibration: Simulated ~ Reference ----------------------------
  
  fit_pooled_calib <- lm(simulated ~ reference_study,
                         data = reg_voutput_wide,
                         weights = w_outcome)
  s  <- summary(fit_pooled_calib)
  ci <- suppressMessages(confint(fit_pooled_calib))
  
  pooled_calibration <- tibble(
    Term     = c("Intercept", "Slope"),
    Estimate = round(coef(fit_pooled_calib), 3),
    `95% CI` = c(
      sprintf("[%.3f, %.3f]", ci["(Intercept)", 1], ci["(Intercept)", 2]),
      sprintf("[%.3f, %.3f]", ci["reference_study",   1], ci["reference_study", 2])
    ),
    `P value`= signif(coef(s)[, 4], 3),
    R2       = c(round(s$r.squared, 3), NA),
    N        = nrow(reg_voutput_wide)
  )
  
  
  write_xlsx(
    list(
      "PooledCalibration"      = pooled_calibration,
      "PooledSimAccuracy"      = pooled_table
    ),
    path = file.path(out_subdir, "et5sup2_simulationaccuracy_deathonly.xlsx")
  )
  ###########################'  
  #Plot figure ---------
  ###########################'    
  # regression line:
  fit_pool <- lm(simulated ~ reference_study, data = reg_voutput_wide, weights =w_outcome )
  b0 <- unname(coef(fit_pool)[1])
  b1 <- unname(coef(fit_pool)[2])
  r2 <- summary(fit_pool)$r.squared
  fit_label <- sprintf("Fit line: y = %.3f + %.3f x (R² = %.3f)", b0, b1, r2)
  
  line_df <- data.frame(
    label     = c(fit_label, "45° line"),
    slope     = c(b1, 1),
    intercept = c(b0, 0)
  )
  
  p1 <- ggplot(reg_voutput_wide, aes(
    x = reference_study, y = simulated,
    color = outcome_label, shape = measure
  )) +
    geom_point(size = 3.2) +
    geom_abline(intercept = b0, slope = b1, 
                color = "darkgrey", linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, 
                linetype = "solid", linewidth = 0.4, color = "black") +
    geom_text_repel(
      aes(label = paste0("(", study_id, "-", cohort_id, ")")),
      size = 3, min.segment.length = 0, segment.color = NA,  max.overlaps = 20,
      box.padding = 0.25,
      point.padding = 0.15
    ) +
    labs(
      x = "Reference Estimate", y = "Simulated Estimate",
      title = "Simulated vs. Real-World Outcomes (by study-cohort)",
      color = "Outcome:", shape = "Measure:"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.title.position = "plot",
      axis.line = element_line(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.background = element_rect(colour = "grey55", fill = NA),
      legend.margin = margin(6, 6, 6, 6),
      plot.caption.position = "plot"
    ) +
    scale_color_tableau(palette = "Tableau 10") +
    coord_fixed(xlim = c(0, .5), ylim = c(0, .5), clip = "on")+
    annotate("segment", x=.40, xend=.47, y=b0+b1*.40, yend=b0+b1*.47,
             linetype="dashed", colour="grey35", linewidth=.35) +
    annotate("text", x=.47, y=b0+b1*.47, label=fit_label,
             hjust=1, vjust=.5, size=3.5, colour="grey25") +
    annotate("segment", x=.40, xend=.47, y=.40, yend=.47,
             colour="black", linewidth=.35) +
    annotate("text", x=.47, y=.47, label="45° line",
             hjust=1, vjust=.5, size=3.5)
  
  
  ggsave(
    filename = "validation_outcomes_comparison_plot.png",
    path     = out_subdir,
    width    = 12, height   = 8,  units    = "in", dpi      = 300
  )
  
  return(reg_voutput_wide)
}



## excluding MACE ---------

process_results_noMACE <- function(vres_outfile, out_subdir ) {
  validation_output_comp_data <- read_csv(file.path(validation_out_dir, vres_outfile), show_col_types = F)
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  labels_map <- c(
    eskd_new      = "ESKD", ckd_new      = "CKD",
    cvd_new       = "CVD", t2d_new       = "T2D",
    ac_death  = "All-Cause Mortality", mace      = "MACE",  stroke = "Stroke")
  
  # Make regression table ---------
  
  ## label cleaning functions: 
  nice_outcome <- function(x) {
    lab <- labels_map[x]
    ifelse(!is.na(lab),
           lab,
           x %>%
             gsub("_", " ", ., fixed = TRUE) %>%
             stringr::str_to_title())
  }
  
  nice_measure <- function(x) case_when(
    x == "incidence_rate_1py"        ~ "1-year incidence",
    x == "incidence_rate_cumulative" ~ "Cumulative incidence",
    TRUE ~ x
  )
  
  ## prep output data
  reg_voutput_wide <- validation_output_comp_data %>%
    filter(outcome != "mace") |> 
    pivot_wider(names_from  = test_indicator,  values_from = estimate ) %>%
    mutate(
      outcome_label = nice_outcome(outcome),
      measure       = nice_measure(measure)
    ) %>%
    left_join(
      validation_data_input %>% select(study_id, cohort_id, futime),
      by = c("study_id","cohort_id")
    ) %>%
    mutate(
      T = if_else(str_detect(measure, regex("cumulative", ignore_case = TRUE)),
                  pmax(futime, 1), 1),
      ref_1y = if_else(
        str_detect(measure, regex("cumulative", ignore_case = TRUE)),
        -log1p(-clamp01(reference_study)) / T,  # cumulative risk -> avg annual rate
        reference_study                         # already 1-year rate
      ),
      sim_1y = if_else(
        str_detect(measure, regex("cumulative", ignore_case = TRUE)),
        -log1p(-clamp01(simulated)) / T,
        simulated
      )
    ) %>%
    arrange(outcome_label, measure)
  
  reg_voutput_long <- validation_output_comp_data |> 
    filter(outcome != "mace") |> 
    left_join(
      validation_data_input %>% select(study_id, cohort_id, futime),
      by = c("study_id","cohort_id")
    ) |> 
    mutate(
      outcome_label = nice_outcome(outcome),
      measure       = nice_measure(measure),
      T             = if_else(str_detect(measure, regex("cumulative", ignore_case = TRUE)),
                              pmax(futime, 1), 1),
      estimate_1yr   = if_else(
        str_detect(measure, regex("cumulative", ignore_case = TRUE)),
        -log1p(-clamp01(estimate)) / T,   # cumulative risk -> avg annual rate
        estimate                           # already 1-year rate
      ),
      test_indicator = factor(test_indicator, levels = c("reference_study", "simulated"))
    ) |> 
    arrange(outcome_label, measure, study_id, cohort_id)
  
  # outcome-level frequency weights
  w_map <- reg_voutput_wide |>
    count(outcome_label, name = "w_outcome")
  
  # attach to wide and long frames
  reg_voutput_wide <- reg_voutput_wide |>
    left_join(w_map, by = "outcome_label")
  
  reg_voutput_long <- reg_voutput_long |>
    left_join(w_map, by = "outcome_label")
  
  ## helper function for regression output 
  tidy_term <- function(fit, term_pattern, n_obs) {
    co <- coef(summary(fit))
    ci <- suppressMessages(as.data.frame(confint(fit)))
    rn <- rownames(co)
    idx <- grep(term_pattern, rn, ignore.case = TRUE, value = TRUE)
    if (length(idx) != 1) stop("Could not uniquely match term pattern: ", term_pattern,
                               " among: ", paste(rn, collapse = ", "))
    tibble(
      term = idx,
      est  = unname(co[idx, 1]),
      lcl  = ci[idx, 1],
      ucl  = ci[idx, 2],
      p    = unname(co[idx, 4]),
      n    = n_obs
    )
  }
  
  # --- (1) Pooled regression: estimate ~ test_indicator -------------------------
  # Intercept = Reference mean; Coef(test_indicator) = Simulated - Reference
  fit_pooled <- lm(estimate ~ test_indicator,
                   data = reg_voutput_long,
                   weights = w_outcome)
  
  pooled_table <- bind_rows(
    tidy_term(fit_pooled, "Intercept", nrow(reg_voutput_long)) %>%
      mutate(Label = "Reference mean"),
    tidy_term(fit_pooled, "test_indicatorsimulated", nrow(reg_voutput_long)) %>%
      mutate(Label = "Difference (Simulated − Reference)")
  ) %>%
    transmute(
      Term     = Label,
      Estimate = round(est, 3),
      `95% CI` = sprintf("[%.3f, %.3f]", lcl, ucl),
      `P value`= signif(p, 3),
      N        = n
    )
  # --- (2) Pooled calibration: Simulated ~ Reference ----------------------------
  
  fit_pooled_calib <- lm(simulated ~ reference_study,
                         data = reg_voutput_wide,
                         weights = w_outcome)
  s  <- summary(fit_pooled_calib)
  ci <- suppressMessages(confint(fit_pooled_calib))
  
  pooled_calibration <- tibble(
    Term     = c("Intercept", "Slope"),
    Estimate = round(coef(fit_pooled_calib), 3),
    `95% CI` = c(
      sprintf("[%.3f, %.3f]", ci["(Intercept)", 1], ci["(Intercept)", 2]),
      sprintf("[%.3f, %.3f]", ci["reference_study",   1], ci["reference_study", 2])
    ),
    `P value`= signif(coef(s)[, 4], 3),
    R2       = c(round(s$r.squared, 3), NA),
    N        = nrow(reg_voutput_wide)
  )
  
  
  write_xlsx(
    list(
      "PooledCalibration"      = pooled_calibration,
      "PooledSimAccuracy"      = pooled_table
    ),
    path = file.path(out_subdir, "et5sup2_simulationaccuracy_noMACE.xlsx")
  )
  
  ###########################'  
  #Plot figure ---------
  ###########################'    
  # regression line:
  fit_pool <- lm(simulated ~ reference_study, data = reg_voutput_wide, weights =w_outcome )
  b0 <- unname(coef(fit_pool)[1])
  b1 <- unname(coef(fit_pool)[2])
  r2 <- summary(fit_pool)$r.squared
  fit_label <- sprintf("Fit line: y = %.3f + %.3f x (R² = %.3f)", b0, b1, r2)
  
  line_df <- data.frame(
    label     = c(fit_label, "45° line"),
    slope     = c(b1, 1),
    intercept = c(b0, 0)
  )
  
  p1 <- ggplot(reg_voutput_wide, aes(
    x = reference_study, y = simulated,
    color = outcome_label, shape = measure
  )) +
    geom_point(size = 3.2) +
    geom_abline(intercept = b0, slope = b1, 
                color = "darkgrey", linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, 
                linetype = "solid", linewidth = 0.4, color = "black") +
    geom_text_repel(
      aes(label = paste0("(", study_id, "-", cohort_id, ")")),
      size = 3, min.segment.length = 0, segment.color = NA,  max.overlaps = 20,
      box.padding = 0.25,
      point.padding = 0.15
    ) +
    labs(
      x = "Reference Estimate", y = "Simulated Estimate",
      title = "Simulated vs. Real-World Outcomes (by study-cohort)",
      color = "Outcome:", shape = "Measure:"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.title.position = "plot",
      axis.line = element_line(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.background = element_rect(colour = "grey55", fill = NA),
      legend.margin = margin(6, 6, 6, 6),
      plot.caption.position = "plot"
    ) +
    scale_color_tableau(palette = "Tableau 10") +
    coord_fixed(xlim = c(0, .5), ylim = c(0, .5), clip = "on")+
    annotate("segment", x=.40, xend=.47, y=b0+b1*.40, yend=b0+b1*.47,
             linetype="dashed", colour="grey35", linewidth=.35) +
    annotate("text", x=.47, y=b0+b1*.47, label=fit_label,
             hjust=1, vjust=.5, size=3.5, colour="grey25") +
    annotate("segment", x=.40, xend=.47, y=.40, yend=.47,
             colour="black", linewidth=.35) +
    annotate("text", x=.47, y=.47, label="45° line",
             hjust=1, vjust=.5, size=3.5)
  
  
  ggsave(
    filename = "validation_outcomes_comparison_plotnomace.png",
    path     = out_subdir,
    width    = 12, height   = 8,  units    = "in", dpi      = 300
  )
  
  return(reg_voutput_wide)
}