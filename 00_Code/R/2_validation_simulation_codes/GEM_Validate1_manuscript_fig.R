#######################
###' Function 2 to generate validation 2-panel figure for manuscript
###' Date: Sept. 2025
###' Author: Liz Mitchell
#######################

#' panel C: exclude = outcome labels left out of the pooled pre/post plot and fit lines
#'   (default MACE; c("MACE", "ESKD") for the version without ESKD); note_exclude = FALSE drops the
#'   "(without ...)" text and subtitle
#' panels B, C: fs = scale for text and point sizes; box = fit-line box position (x, y as fraction of axes, hjust, vjust);
#'   box_fs = extra scale for the fit-line box text; title = panel title (NULL = default); save = FALSE returns the plot only (used by make_calib_figure())
make_panel_figs <- function(vres_pre, vres_post=NULL, out_path_png, panel = NULL,
                           xlim = c(0, .5), ylim = c(0, .5), exclude = "MACE", note_exclude = TRUE,
                           fs = 1, box = c(0.3, 0.15, 0, 1), box_fs = 1, title = NULL, save = TRUE) {
  clamp01 <- function(x, eps = 1e-12) pmin(pmax(x, eps), 1 - eps)
  
  prep_wide <- function(df) {
    df |>
      pivot_wider(names_from = test_indicator, values_from = estimate) |>
      left_join(validation_data_input |> select(study_id, cohort_id, futime),
                by = c("study_id","cohort_id")) |>
      mutate(
        outcome_label = case_when(
          outcome == "eskd_new" ~ "ESKD",
          outcome == "ckd_new"  ~ "CKD",
          outcome == "cvd_new"  ~ "CVD",
          outcome == "t2d_new"  ~ "T2D",
          outcome == "ac_death" ~ "All-Cause Mortality",
          outcome == "mace"     ~ "MACE",
          outcome == "stroke"   ~ "Stroke",
          TRUE ~ outcome
        ),
        measure = case_when(
          measure == "incidence_rate_1py"        ~ "1-year incidence",
          measure == "incidence_rate_cumulative" ~ "Cumulative incidence",
          TRUE ~ measure
        )) |>
      arrange(outcome_label, measure)
  }
  
  
  # --------- Panel A: pooled pre-calibration results ----------
  if (is.null(vres_post)){
    dat_pre  <- read_csv(file.path(validation_out_dir, vres_pre),  show_col_types = FALSE)
    wide_pre  <- prep_wide(dat_pre)
    w_map <- wide_pre |>
      count(outcome_label, name = "w_outcome")
    wide_pre <-  wide_pre |>  left_join(w_map, by = "outcome_label")
    
    fit_pre <- lm(simulated ~ reference_study, data = wide_pre, weights = w_outcome )
    b0 <- unname(coef(fit_pre)[1]); b1 <- unname(coef(fit_pre)[2])
    r2 <- summary(fit_pre)$r.squared
    tbl_lab <- sprintf("Fit line:\nIntercept=%.2f, Slope=%.2f, R²=%.2f", b0, b1, r2)
    
    
    p <- ggplot(wide_pre, aes(x = reference_study, y = simulated,
                              color = outcome_label, shape = measure)) +
      geom_point(size = 3.5) +
      geom_abline(intercept = b0, slope = b1, color = "grey35",
                  linetype = "dashed", linewidth = 0.35) +
      geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.35) +
      geom_text_repel(
        aes(label = paste0("(", study_id, "-", cohort_id, ")")),
        size = 4.5, min.segment.length = 0, segment.color = NA,  max.overlaps = 20,
        box.padding = 0.25,
        point.padding = 0.15
      ) +
      annotate("label",
               x = xlim[2] - 0.5 * diff(xlim),  
               y = ylim[1] + 0.1 * diff(ylim), 
               hjust = 0, vjust = 1,
               label = tbl_lab,
               label.size = 0.3,size = 5.5,  fill = "white", color = "grey20"
      ) +
      labs(title = "A. Pooled Prediction Accuracy of All Five Outcomes\nBefore Calibration",
           x = "Reference Estimate", y = "Simulated Estimate",
           color = "Outcome:", shape = "Measure:") +
      theme_minimal(base_size = 16) +
        theme(       
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            legend.position = "bottom",
            legend.box = "vertical",
            legend.text = element_text(size = 15),  
            legend.title = element_text(size = 15),
            legend.box.background = element_rect(colour = "grey55", fill = NA),
            legend.spacing.y = unit(.3, "mm"),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 14),
            axis.line = element_line()) +
      scale_color_tableau(palette = "Tableau 10") +
      coord_fixed(xlim = xlim, ylim = ylim, clip = "on")
    
    ggsave(out_path_png, p, width = 9, height = 10, dpi = 300)
    return(invisible(p))
  } else if (panel == "B"){ 
  # --------- Panel B: pre- post- calibration results for mortality ----------
  dat_pre   <- read_csv(file.path(validation_out_dir, vres_pre),  show_col_types = FALSE)
  dat_post  <- read_csv(file.path(validation_out_dir, vres_post), show_col_types = FALSE)
  wide_pre  <- prep_wide(dat_pre)
  wide_post <- prep_wide(dat_post)
  
  death_pre  <- wide_pre  |> filter(outcome_label == "All-Cause Mortality") |>  
    select(study_id, cohort_id, measure, reference_study, simulated_pre = simulated)
  death_post <- wide_post |> filter(outcome_label == "All-Cause Mortality") |> 
    select(study_id, cohort_id, measure, simulated_post = simulated)
  
  death_combined <- death_pre |>
    left_join(death_post, by = c("study_id","cohort_id", "measure"))
  
  fit_pre  <- lm(simulated_pre   ~ reference_study, data = death_combined )
  fit_post <- lm(simulated_post  ~ reference_study, data = death_combined  )
  
  b0_pre  <- unname(coef(fit_pre)[1]);  b1_pre  <- unname(coef(fit_pre)[2]);  r2_pre  <- summary(fit_pre)$r.squared
  b0_post <- unname(coef(fit_post)[1]); b1_post <- unname(coef(fit_post)[2]); r2_post <- summary(fit_post)$r.squared
  
  tbl_lab <- sprintf(
    "Fit line:\nBefore: Intercept = %.2f, Slope = %.2f, R² = %.2f\nAfter:   Intercept = %.2f, Slope = %.2f, R² = %.2f",
    b0_pre, b1_pre, r2_pre, b0_post, b1_post, r2_post
  )
  
  death_long <- death_combined |>
    pivot_longer(c(simulated_pre, simulated_post),
                        names_to = "phase", values_to = "simulated") |>
   mutate( phase =recode(phase, simulated_pre = "Before", simulated_post = "After"))
  
  fit_df <- tibble(
    phase     = c("Before","After"),
    slope     = c(b1_pre, b1_post),
    intercept = c(b0_pre, b0_post)
  )
  
  p <- ggplot(death_long, aes(x = reference_study, y = simulated)) + 
    geom_point(aes(color = phase, shape = measure), size = 3.5 * fs) +
    # fit lines
    geom_abline(
      data = fit_df,
      aes(slope = slope, intercept = intercept, linetype = phase, color = phase),
      linewidth = 0.5,  show.legend = FALSE
    ) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "solid", linewidth = 0.4) +
    # legends
    scale_color_manual(
      values = c(Before = "#5A9AC5", After = "#174F78"),  
      name   = "Phase:"
    ) +
    scale_shape_manual(
      values = c("Cumulative incidence" = 17, "1-year incidence" = 16),
      name   = "Measure:" ) +
    scale_linetype_manual(
      values = c(Before = "longdash", After = "dotted"),
      name   = "Fit lines"
    ) +
    geom_text_repel(
      aes(label = paste0("(", study_id, "-", cohort_id, ")"), color = phase),
      size = 4.5 * fs, min.segment.length = 0, segment.color = NA,  max.overlaps = 20,
      box.padding = 0.25,
      point.padding = 0.15
    ) +
      # stats box (top-left inside panel)
      annotate("label",
               x = xlim[1] + box[1] * diff(xlim),
               y = ylim[1] + box[2] * diff(ylim),
               hjust = box[3], vjust = box[4],
               label = tbl_lab,
               label.size = 0.3,size = 5.5 * fs * box_fs,  fill = "white", color = "grey20") +
    labs(title = if (is.null(title)) "B. Prediction Accuracy of All-Cause Mortality\nBefore and After Calibration" else title,
         x = "Reference Estimate", y = "Simulated Estimate" )+
    theme_minimal(base_size = 16 * fs) +
      theme(       
          plot.title = element_text(hjust = 0.5, size = 16 * fs, face = "bold"),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.text = element_text(size = 15 * fs),  
          legend.title = element_text(size = 15 * fs),
          legend.box.background = element_rect(colour = "grey55", fill = NA),
          legend.spacing.y = unit(.3, "mm"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 14 * fs),
          axis.line = element_line()) +
    coord_fixed(xlim = xlim, ylim = ylim, clip = "on")
  
  if (save) ggsave(out_path_png, p, width = 9, height = 10, dpi = 300)
  return(invisible(p))
  } 
  # --------- Panel C: pre- post- pooled calibration results without excluded outcomes (default MACE) ----------
  dat_pre   <- read_csv(file.path(validation_out_dir, vres_pre),  show_col_types = FALSE)
  dat_post  <- read_csv(file.path(validation_out_dir, vres_post), show_col_types = FALSE)
  wide_pre  <- prep_wide(dat_pre)
  wide_post <- prep_wide(dat_post)
  
  ds_pre  <- wide_pre  |> filter(!outcome_label %in% exclude) |>  
    select(study_id, cohort_id,outcome_label, measure, reference_study, simulated_pre = simulated)
  ds_post <- wide_post |> filter(!outcome_label %in% exclude) |> 
    select(study_id, cohort_id, outcome_label, measure, simulated_post = simulated)
  
  ds_combined <- ds_pre |>
    left_join(ds_post, by = c("study_id","cohort_id", "measure", "outcome_label"))
  
  # outcome freq
  w_map <- ds_pre |> count(outcome_label, name = "w_outcome")
  ds_combined <-  ds_combined |>  left_join(w_map, by = "outcome_label")
  
  fit_pre  <- lm(simulated_pre   ~ reference_study, data = ds_combined, weights = w_outcome )
  fit_post <- lm(simulated_post  ~ reference_study, data = ds_combined, weights = w_outcome  )
  b0_pre  <- unname(coef(fit_pre)[1]);  b1_pre  <- unname(coef(fit_pre)[2]);  r2_pre  <- summary(fit_pre)$r.squared
  b0_post <- unname(coef(fit_post)[1]); b1_post <- unname(coef(fit_post)[2]); r2_post <- summary(fit_post)$r.squared
  
  tbl_lab <- sprintf(
    "Fit line%s:\nBefore: Intercept = %.2f, Slope = %.2f, R² = %.2f\nAfter:   Intercept = %.2f, Slope = %.2f, R² = %.2f",
    if (note_exclude) paste0(" (without ", paste(exclude, collapse = ", "), ")") else "", b0_pre, b1_pre, r2_pre, b0_post, b1_post, r2_post
  )
  
  ds_long <- ds_combined |>
    pivot_longer(c(simulated_pre, simulated_post),
                 names_to = "phase", values_to = "simulated") |>
    mutate( phase =recode(phase, simulated_pre = "Before", simulated_post = "After"))
  
  # colors match panel A: Tableau 10 over all outcomes (incl. excluded), in panel A order
  labs  <- levels(factor(wide_pre$outcome_label))
  colmap <- setNames(ggthemes::tableau_color_pal("Tableau 10")(length(labs)), labs)
  
  
  fit_df <- tibble(
    phase     = c("Before","After"),
    slope     = c(b1_pre, b1_post),
    intercept = c(b0_pre, b0_post)
  )
 
  p <- ggplot(ds_long, aes(x = reference_study, y = simulated)) + 
    geom_point(
      data = subset(ds_long, phase == "After"),
      aes(color = outcome_label, shape = measure),
      size = 3.5 * fs, stroke = 0.4
    ) +
    # fit lines
    geom_abline(
      data = fit_df,
      aes(slope = slope, intercept = intercept, linetype = phase),
      linewidth = 0.5,  show.legend = TRUE) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "solid", linewidth = 0.4) +
    # legends
    geom_text_repel(data = subset(ds_long, phase == "After"),
      aes(label = paste0("(", study_id, "-", cohort_id, ")"), color = outcome_label),
      size = 4.5 * fs, min.segment.length = 0, segment.color = NA,  max.overlaps = 20,
      box.padding = 0.25,
      point.padding = 0.15
    ) +
    scale_shape_manual(
      values = c("Cumulative incidence" = 17, "1-year incidence" = 16),
      name   = "Measure:" ) +
    scale_color_manual(values = colmap, name = "Outcome:") +
    scale_linetype_manual(values = c("Before" = "longdash", "After" = "dotted"),
                          name = "Phase:") +
    guides(
      shape    = guide_legend(order = 1, override.aes = list(linetype = 0, color = "grey30", fill = NA)),
      color    = guide_legend(order = 2, override.aes = list(linetype = 0)),
      linetype = guide_legend(order = 3, override.aes = list(color = "black"))
    ) +
    # stats box (top-left inside panel)
    annotate("label",
             x = xlim[1] + box[1] * diff(xlim),
             y = ylim[1] + box[2] * diff(ylim),
             hjust = box[3], vjust = box[4],
             label = tbl_lab,
             label.size = 0.3,size = 5.5 * fs * box_fs,  fill = "white", color = "grey20") +
    labs(title = if (is.null(title)) "C. Pooled Prediction Accuracy\nAfter Mortality Calibration" else title,
         subtitle = if (note_exclude) paste("Excluding", paste(exclude, collapse = " and ")) else NULL,
         x = "Reference Estimate", y = "Simulated Estimate" )+
    theme_minimal(base_size = 16 * fs) +
      theme(
          plot.title = element_text(hjust = 0.5, size = 16 * fs, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.text = element_text(size = 15 * fs),  
          legend.title = element_text(size = 15 * fs),
          legend.spacing.y = unit(.2, "mm"),
          legend.box.background = element_rect(colour = "grey55", fill = NA),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 14 * fs),
          axis.line = element_line()) +
    coord_fixed(xlim = xlim, ylim = ylim, clip = "on")
  
  if (save) ggsave(out_path_png, p, width = 9, height = 10.5, dpi = 300)
  invisible(p)
}

# Manuscript figure: 2 panels side by side, before vs after mortality calibration ----
#'  A: all-cause mortality (panel B layout); B: all outcomes except MACE (panel C layout; MACE not shown or named)
#'  sized for a double-column page (7.5 x 5.2 in, ~190 mm wide); text at about 7-8 pt; PNG 600 dpi + PDF
make_calib_figure <- function(vres_pre, vres_post, out_path, width = 7.5, height = 5.2,
                              xlim = c(0, .5), ylim = c(0, .5), fs = 0.5) {
  pa <- make_panel_figs(vres_pre, vres_post, NULL, panel = "B", xlim = xlim, ylim = ylim, fs = fs,
                        box = c(0.02, 0.98, 0, 1), box_fs = 0.85, save = FALSE,
                        title = "A. All-Cause Mortality\nBefore and After Calibration")
  pb <- make_panel_figs(vres_pre, vres_post, NULL, panel = "C", xlim = xlim, ylim = ylim, fs = fs,
                        exclude = "MACE", note_exclude = FALSE, box = c(0.02, 0.98, 0, 1), box_fs = 0.85, save = FALSE,
                        title = "B. All Outcomes\nBefore and After Mortality Calibration")
  p <- (pa | pb) & theme(plot.background = element_rect(fill = "white", colour = NA),
                         legend.box.spacing = unit(2, "mm"), legend.key.size = unit(3.5, "mm"),
                         legend.text = element_text(size = 13 * fs), legend.title = element_text(size = 13 * fs))
  out_path <- sub("\\.(png|pdf)$", "", out_path)
  ggsave(paste0(out_path, ".png"), p, width = width, height = height, dpi = 600, bg = "white")
  ggsave(paste0(out_path, ".pdf"), p, width = width, height = height,
         device = if (capabilities("cairo")) cairo_pdf else "pdf", bg = "white")
  invisible(p)
}

    
    
  
  ## to add hollow points for before:
# geom_point(
#   data = subset(ds_long, phase == "Before"),
#   aes(color = outcome_label, shape = measure),
#   size = 4.5 * fs, stroke = 0.9, fill = NA
# ) +
  ## use shapes 24 and 21 to use "fill" quality
  
  
  
  
  
  
  