# R/prep_main_inputs.R
# 1. reads UN, NCD-RisC, GBD data files in 01_Data folder -------------
#    writes cleaning log to paths$out/logs/1data_cleaning_log*.txt

prep_main_inputs <- function(paths, config) {
  # Set up logging------------
  dir.create(file.path(paths$out, "logs"), recursive = T, showWarnings=F)
  log_file <- file.path(paths$out, "logs", paste0("1data_cleaning_log",format(Sys.Date(), "%m.%d"), ".txt"))
  log_msg <- function(..., .sep = " ") {
    msg <- paste(..., sep = .sep)
    cat(msg, "\n", file = log_file, append = TRUE)
  }
  log_df <- function(x, header = NULL) {
    if (!is.null(header)) {log_msg(header)}
    capture.output(print(x), file = log_file, append = TRUE)
  }
  
  # Disability weights data=---------
  dw_file <- file.path(paths$data, paste0("IHME_GBD_disability_weights_", config$GBD_data_yr, ".xlsx"))
  if (!file.exists(dw_file)) { stop("disability weight file not found: ", dw_file)}
  log_msg("\nUsing disability weight file:", basename(dw_file))
  
  dw_data <- read_excel(dw_file) |>  
    select( analysis_indicator,
            healthstate_name,
            probability_weight,
            condition_logic,
            mean, lower, upper, `Application Notes`) |>  
    filter(analysis_indicator=="primary") |>  select(-analysis_indicator)
  
  write_csv(dw_data, file.path(paths$working, "condition_disability_wts.csv"))

  # Population data------------
  pop_file <- file.path(paths$data, paste0("un_pop_data_", config$UN_data_yr, ".csv"))
  if (!file.exists(pop_file)) { stop("UN population file not found: ", pop_file)}
  log_msg("\nUsing UN population file:", basename(pop_file))
  pop_data <- read_csv(pop_file, show_col_types = FALSE)
  
  pop_data <- pop_data |>
    rename_with(tolower) |> 
    rename(measure=variant) |> 
    pivot_wider(names_from = measure, values_from = value) |> 
    rename(ci_lb = `95% lower bound`, ci_ub = `95% upper bound`, median = Median,
           age_low = agestart, age_high = ageend) |> 
    mutate(
      female = if_else(sex == "Female", 1L, if_else(sex == "Male", 0L, NA_integer_)),
      pop_size = round(median),
      # Clean country names with accents/unicode issues:
      location_decoded = stri_unescape_unicode(location),         
      location_cleaned = stri_trans_general(location_decoded, "Latin-ASCII"), 
      location_cleaned = gsub("[[:punct:]]", "", location_cleaned), 
      location_cleaned = tolower(location_cleaned),
      location_cleaned = case_when(
        grepl("voire", location_cleaned, ignore.case = TRUE) ~ "cote divoire",
        grepl("cura.ao", location_cleaned, ignore.case = TRUE) ~ "curacao",  
        grepl("r.union", location_cleaned, ignore.case = TRUE) ~ "reunion", 
        grepl("t*kiye", location_cleaned, ignore.case = TRUE) ~ "turkiye", 
        grepl("barth*emy", location_cleaned, ignore.case = TRUE) ~ "saint barthelemy", 
        TRUE ~ location_cleaned
      )
    ) |> 
    filter(age_low < 100) |>
    mutate(age_high = as.double(age_high)) |> 
    select(-sex, -time, -median, -location_decoded, -location) |> 
    rename(location=location_cleaned) 
  
  # ISO mappings ----------------------------------------
  
  country_name_iso <- pop_data |> distinct(location, iso2, iso3) 
  country_groups = read.csv(file.path(paths$data, "who_world_bank_country_groups_iso.csv"))
  country_groups <- country_groups |> 
    mutate(
      # punctuation -> space ("Lower-middle-income" = "Lower middle income"; 6 rows in source file)
      world_bank_income_group = tolower(str_squish(gsub("[[:punct:]]", " ", world_bank_income_group))),
      who_region = tolower(who_region)
    ) |> 
    filter(!is.na(iso3), !grepl("^XX[1-9]$|^X1[0-3]$", iso3))
  wb_ok <- c("high income", "upper middle income", "lower middle income", "low income", "no wb income group")
  if (!all(country_groups$world_bank_income_group %in% wb_ok))
    stop("unknown World Bank income group: ", paste(setdiff(country_groups$world_bank_income_group, wb_ok), collapse = ", "))
  
  # hic/lmic grouping
  country_income_groups_2 <- country_groups |> 
    select(-who_region) |> 
    mutate(
      income_group = case_when(
        world_bank_income_group ==  "high income"  ~ "hic",
        world_bank_income_group == "no wb income group" ~ "not_assigned",
        .default = "lmic"
      )) |> 
    select(-world_bank_income_group)
  
  write_csv(country_name_iso, file.path(paths$working, "country_name_iso_mapping.csv"))
  write.csv(country_groups, file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv") )
  write.csv(country_income_groups_2, file.path(paths$working, "country_income_groups_2.csv") )
  
  ###############################'
  ###############################'
  
  #  Load data for multimorbidity prevalence, incidence, and death ---------------

  cc_fp <- file.path(paths$working, "pop_cc_cohort.csv")
  cc_ok <- file.exists(cc_fp) && {                     # rebuild if file predates BMI-category or adult-only columns
    h <- names(read.csv(cc_fp, nrows = 1))
    "prevalence_bmi_lt18_val" %in% h && !any(grepl("^incidence_(class[123]_)?obesity", h))
  }
  if (!cc_ok) { 
    # Helper function to clean GBD data
    gbd_location_clean <- function(data) {
      if ("location" %in% names(data)) {
        data <-  rename(data,  location_name =location)
      }
      remove_words_country_names <- c(
        "republic of nauru",
        "principality of monaco",
        "republic of niue",
        "republic of san marino",
        "republic of palau"
      )
      data |> 
        mutate(
          # clean location_name to match un population data source
          location_decoded = stri_unescape_unicode(location_name),         
          location_cleaned = stri_trans_general(location_decoded, "Latin-ASCII"), 
          location_cleaned = gsub("[[:punct:]]", "", location_cleaned), 
          location_cleaned = tolower(location_cleaned),
          location_cleaned = case_when(
            grepl("ivoire$", location_cleaned, ignore.case = TRUE) ~ "cote divoire",
            grepl("cura.ao", location_cleaned, ignore.case = TRUE) ~ "curacao", 
            grepl("r.union", location_cleaned, ignore.case = TRUE) ~ "reunion",  
            grepl("t.kiye", location_cleaned, ignore.case = TRUE) ~ "turkiye", 
            grepl("taiwan", location_cleaned) ~ "china taiwan province of china",
            grepl("lao", location_cleaned) ~ "lao peoples dem republic",
            grepl("democratic.*republic of korea", location_cleaned) ~ "dem peoples rep of korea",
            grepl("micronesia", location_cleaned) ~ "micronesia fed states of",
            grepl("viet nam", location_cleaned) ~ "viet nam",
            grepl("macedonia", location_cleaned) ~ "north macedonia",
            grepl("czech", location_cleaned) ~ "czechia",
            grepl("great britain and northern ireland", location_cleaned) ~ "united kingdom",
            grepl("bahamas", location_cleaned) ~ "bahamas",
            grepl("bolivia", location_cleaned) ~ "bolivia plurinational state of",
            grepl("iran", location_cleaned) ~ "iran islamic republic of",
            grepl("venezuela", location_cleaned) ~ "venezuela bolivarian republic of",
            grepl("palestine", location_cleaned) ~ "state of palestine",
            grepl("turkey", location_cleaned) ~ "turkiye",
            grepl("democratic republic of the congo", location_cleaned) ~ "dem rep of the congo",
            grepl("kingdom of eswatini", location_cleaned) ~ "eswatini",
            grepl("cabo verde", location_cleaned) ~ "cabo verde",
            grepl("republic of the gambia", location_cleaned) ~ "gambia",
            location_cleaned %in% remove_words_country_names ~ gsub(".*of\\s+", "", location_cleaned),
            TRUE ~ location_cleaned
          )
        ) |> 
        rename(location = location_cleaned) |> 
        select(-location_decoded, -location_name)
    }
    
    load_data <- function(data_source_identifier, data_year, measurement) {
      
        ### NCD-risc data load-------------------
      # clean prevalence rates and use to approximate incidence
      if (data_source_identifier == "obesity") {
  
        data_male   <- read_csv(file.path(paths$data, "NCD_RisC_Lancet_2024_BMI_male_age_specific_country.csv"), show_col_types = FALSE)
        data_female <- read_csv(file.path(paths$data, "NCD_RisC_Lancet_2024_BMI_female_age_specific_country.csv"), show_col_types = FALSE)
        # adults only: NCD-RisC 18-19 group represents the UN 15-19 band (only ages 18-19 are modeled)
        data <- bind_rows(data_male, data_female) |> 
          rename_with(tolower) |> 
          filter(year == data_year) |> 
          mutate(`age group` = if_else(`age group` == "18-19", "15-19", `age group`)) |> 
          select(year, sex, iso, `age group`, 
                 matches("\\(underweight\\)"),                 # full category set for continuous BMI distribution
                 matches("prevalence of bmi 18.5"),
                 matches("prevalence of bmi 20"),
                 matches("prevalence of bmi 25"),
                 matches("\\(obesity\\)"),
                 matches("prevalence of bmi 30"), 
                 matches("prevalence of bmi 35"),
                 matches("prevalence of bmi >=40")
          ) |> 
          mutate(`age group` = if_else(`age group` == "85plus", "85-89", `age group`)) |> 
          separate(`age group`, into = c("age_low", "age_high"), sep = "-", convert = TRUE) |>
          rename_with(~ case_when(
            str_detect(., "\\(underweight\\)") & str_detect(., "lower") ~ "prevalence_bmi_lt18_lower",
            str_detect(., "\\(underweight\\)") & str_detect(., "upper") ~ "prevalence_bmi_lt18_upper",
            str_detect(., "\\(underweight\\)") ~ "prevalence_bmi_lt18_val",
            str_detect(., "bmi 18.5") & str_detect(., "lower") ~ "prevalence_bmi_18_20_lower",
            str_detect(., "bmi 18.5") & str_detect(., "upper") ~ "prevalence_bmi_18_20_upper",
            str_detect(., "bmi 18.5") ~ "prevalence_bmi_18_20_val",
            str_detect(., "bmi 20 ") & str_detect(., "lower") ~ "prevalence_bmi_20_25_lower",
            str_detect(., "bmi 20 ") & str_detect(., "upper") ~ "prevalence_bmi_20_25_upper",
            str_detect(., "bmi 20 ") ~ "prevalence_bmi_20_25_val",
            str_detect(., "bmi 25 ") & str_detect(., "lower") ~ "prevalence_bmi_25_30_lower",
            str_detect(., "bmi 25 ") & str_detect(., "upper") ~ "prevalence_bmi_25_30_upper",
            str_detect(., "bmi 25 ") ~ "prevalence_bmi_25_30_val",
            (str_detect(., "bmi>=30")) & str_detect(., "lower") ~ "prevalence_obesity_lower",
            str_detect(., "bmi>=30") & str_detect(., "upper") ~ "prevalence_obesity_upper",
            str_detect(., "bmi>=30") ~ "prevalence_obesity_val",         # only applies if upper and lower not detected
            str_detect(., "bmi 30.*<35") & str_detect(., "lower") ~ "prevalence_class1_obesity_lower",
            str_detect(., "bmi 30.*<35") & str_detect(., "upper") ~ "prevalence_class1_obesity_upper",
            str_detect(., "bmi 30.*<35") ~ "prevalence_class1_obesity_val",
            str_detect(., "bmi 35.*<40") & str_detect(., "lower") ~ "prevalence_class2_obesity_lower",
            str_detect(., "bmi 35.*<40") & str_detect(., "upper") ~ "prevalence_class2_obesity_upper",
            str_detect(., "bmi 35.*<40") ~ "prevalence_class2_obesity_val",
            str_detect(., "bmi >=40") & str_detect(., "lower") ~ "prevalence_class3_obesity_lower",
            str_detect(., "bmi >=40") & str_detect(., "upper") ~ "prevalence_class3_obesity_upper",
            str_detect(., "bmi >=40") ~ "prevalence_class3_obesity_val",
            TRUE ~ .
          )) |> 
          rename(iso3 = iso) |> 
          mutate(female = if_else(sex == "Women", 1, 0)) |> 
          select(-sex, -year)
        
         ### GBD (prev/inc) data load-------------------  
      } else if (data_source_identifier == "gbd") {
        
        data_male   <- read_csv(file.path(paths$data, paste0("IHME-GBD_prev_inc_male_", data_year, ".csv")), show_col_types = FALSE)
        data_female <- read_csv(file.path(paths$data, paste0("IHME-GBD_prev_inc_female_",  data_year, ".csv")), show_col_types = FALSE)
        data <- bind_rows(data_male, data_female)
        
        conditions <- tibble(
          condition_code = c("t2d", "cvd", "ckd"),
          cause_name_pattern = c(
            "(?i)^Diabetes mellitus type 2$", 
            "(?i)^Cardiovascular diseases$", 
            "(?i)^Chronic kidney disease$"
          ))
        if (measurement == "incidence") { 
          outcomes <- tibble(
            outcome_code = c("stroke"), 
            outcome_cause_name_pattern = c(
              "Stroke"
            ))
        }
        data <- gbd_location_clean(data)
        data <- data |> 
          select(measure, location, sex, age, cause, metric, val, upper, lower) |> 
          filter(tolower(measure) == tolower(measurement) & metric == "Rate") |> 
          left_join(country_name_iso |> select(location, iso3), 
                    by = "location") |> 
          filter(!is.na(iso3)) |>        
          select(-location) |>
          mutate( condition_code = case_when(
            str_detect(cause, conditions$cause_name_pattern[1]) ~ conditions$condition_code[1],
            str_detect(cause, conditions$cause_name_pattern[2]) ~ conditions$condition_code[2],
            str_detect(cause, conditions$cause_name_pattern[3]) ~ conditions$condition_code[3],
            TRUE ~ NA_character_
          ))
        
        # Process additional incidence outcomes
        if (measurement == "incidence") {
          data <- data |>  
            mutate(
              condition_code = case_when(
                str_detect(cause, outcomes$outcome_cause_name_pattern[1]) ~ outcomes$outcome_code[1],
                TRUE ~ condition_code
              ))
        }
        
        data <- data |>
          filter(!is.na(condition_code)) |>
          select(-cause, -metric, -measure) |>
          mutate(
            across(c(val, lower, upper), ~ .x / 100000)    # Transform rate (given as per 100,000) to be per person
          ) |> 
          pivot_wider(
            names_from = condition_code,
            values_from = c(val, lower, upper),
            names_glue = paste0(measurement, "_{condition_code}_{.value}")
          ) |>
          mutate(
            female = if_else(sex == "Female", 1, 0),
            age = if_else(age == "95+ years", "95-99 years", age)          
          ) |>
          select(-sex) |>
          separate(age, into = c("age_low", "age_high"), sep = "-", convert = TRUE) |>
          mutate(
            age_high = as.numeric(str_remove_all(age_high, "\\s*years"))
          )
              
              ## GBD (death) data load-------------------  
      } else if (data_source_identifier == "gbd death") {
        data   <- read_csv(file.path(paths$data, paste0("IHME-GBD_mortality_by_cause_", data_year, ".csv")), show_col_types = FALSE)
           
        # Set var lists:
        conditions <- tibble(
          cause_code = c( "t2d", "cvd", "ckd_due_to_t2d", "ckd", "all_cause", "stroke"),
          cause_name_pattern = c(
            "(?i)^Diabetes mellitus type 2$", 
            "(?i)^Cardiovascular diseases$", 
            "Chronic kidney disease due to diabetes mellitus type 2$",
            "(?i)^Chronic kidney disease$",
            "All causes$",
            "Stroke$"
          ))
        data <- gbd_location_clean(data)
        data <- data |> 
          select(location, sex_name, age_name, cause_name, metric_name, val, upper, lower) |> 
          filter(metric_name == "Rate")  |>              #  (defined as deaths per 100,000 population by cause)
          left_join(country_name_iso |> select(location, iso3), by = "location") |> 
          filter(!is.na(iso3)) |>                     
          select(-location) |> 
          # Reshape estimates: 
          mutate(
            death_cause_code = case_when(
              str_detect(cause_name, conditions$cause_name_pattern[1]) ~ conditions$cause_code[1],
              str_detect(cause_name, conditions$cause_name_pattern[2]) ~ conditions$cause_code[2],
              str_detect(cause_name, conditions$cause_name_pattern[3]) ~ conditions$cause_code[3],
              str_detect(cause_name, conditions$cause_name_pattern[4]) ~ conditions$cause_code[4],
              str_detect(cause_name, conditions$cause_name_pattern[5]) ~ conditions$cause_code[5],
              str_detect(cause_name, conditions$cause_name_pattern[6]) ~ conditions$cause_code[6],
              TRUE ~ NA_character_
            )) |> 
          filter(!is.na(death_cause_code)) |>
          select(-cause_name, -metric_name) |>
          mutate(across(c(val, lower, upper), ~ .x / 100000))|>   # Transform rate (given as per 100,000) to be per person
          pivot_wider(
            names_from = death_cause_code,
            values_from = c(val, lower, upper),
            names_glue = paste0(measurement, "_{death_cause_code}_{.value}")
          )|>
          mutate(across(
            where(is.list),
            ~ map(.x, ~ {
              if (is.null(.x) || length(.x) == 0) {
                NA_real_
              } else if (length(.x) == 1) {
                .x
              } else {
                warning("Multiple values in list element; collapsing to the first value")
                .x[1] 
              }
            }) |> 
              unlist()
          )) |> 
          # Clean sex and age
          mutate(
            female   = if_else(sex_name == "Female", 1, 0),
            age_name = if_else(age_name == "95+ years", "95-99 years", age_name)    
          ) |>
          select(-sex_name) |>
          separate(age_name, into = c("age_low", "age_high"), sep = "-", convert = TRUE) |>
          mutate(age_high = as.numeric(str_remove_all(age_high, "\\s*years")))
          
        # ## GBD (life table) data load-------------------
        # not used, rates = probability of death in 5-yr age group not 1-yr rate
        # } else if (data_source_identifier == "gbd lt") {
        #   data   <- read_csv(file.path(paths$data, paste0("IHME-GBD_life_table_", data_year, ".csv")), show_col_types = FALSE)
        #   
        #   data <- gbd_location_clean(data)
        #   data <- data |>
        #     select(location, sex, age, val, upper, lower) |>
        #     left_join(country_name_iso |> select(location, iso3), by = "location") |> 
        #     filter(!is.na(iso3)) |>
        #     mutate(
        #       female = if_else(sex == "Female", 1, 0),
        #       age = if_else(age == "95+ years", "95-99 years", age)
        #     ) |>
        #     select( -location, -sex) |>
        #     separate(age, into = c("age_low", "age_high"), sep = "-", convert = TRUE) |>
        #     mutate(age_high = as.numeric(str_remove_all(age_high, "\\s*years")))|>
        #     rename(
        #       !!paste0(measurement, "_val")   := val,
        #       !!paste0(measurement, "_lower") := lower,
        #       !!paste0(measurement, "_upper") := upper
        #     )
          ## GBD (life expectancy) data load-------------------
        } else if (data_source_identifier == "gbd le") {
          data   <- read_csv(file.path(paths$data, paste0("IHME-GBD_life_expectancy_", data_year, ".csv")), show_col_types = FALSE)
          
          data <- gbd_location_clean(data)
          data <- data |>
            select(location, sex=sex_name, age=age_name, val, upper, lower) |> 
            left_join(country_name_iso |> select(location, iso3), by = "location") |> 
            filter(!is.na(iso3)) |>
            mutate(
              female = if_else(sex == "Female", 1, 0),
              age = if_else(age == "95+ years", "95-99 years", age)
            ) |>
            select( -location, -sex) |>
            separate(age, into = c("age_low", "age_high"), sep = "-", convert = TRUE) |>
            mutate(age_high = as.numeric(str_remove_all(age_high, "\\s*years")))|>
            rename(
              !!paste0(measurement, "_val")   := val,
              !!paste0(measurement, "_lower") := lower,
              !!paste0(measurement, "_upper") := upper
            )
        }
      return(data)
    }
    
      ## Incidence/prevalence working data construction---------------------------------------------
    
    obesity_data        <- load_data("obesity",    config$NCDRisC_data_yr, "prevalence")
    gbd_data_prevalence <- load_data("gbd",        config$GBD_data_yr, "prevalence")
    gbd_data_incidence  <- load_data("gbd",        config$GBD_data_yr, "incidence")
    gbd_data_mortality  <- load_data("gbd death",  config$GBD_data_yr, "death")
  
    # gbd_data_lifetable  <- load_data("gbd lt",     config$GBD_data_yr, "prob_death")  # caveat: = prob of dying in age interval not 1-yr rate 
    gbd_data_lifeexp  <- load_data("gbd le",     config$GBD_data_yr, "le_yrs")
    
    dfs <- list(obesity_data, gbd_data_prevalence, gbd_data_incidence, gbd_data_mortality, 
                # gbd_data_lifetable, 
                gbd_data_lifeexp)
    
    country_groups = read.csv(file.path(paths$working, "who_world_bank_country_groups_iso_cleaned.csv"))
    matching_vars = c("iso3", "female", "age_low", "age_high")
    
    # Bind data by cohort:
    rates_raw <- reduce(dfs, full_join, by = matching_vars)
   
    # Drop rate-row cohorts w/out population data
    unmatched_pop_in_rates <- anti_join(rates_raw, pop_data, by = matching_vars)
    min_age_unmatched <- if (nrow(unmatched_pop_in_rates) == 0) NA else
      suppressWarnings(min(unmatched_pop_in_rates$age_low, na.rm = TRUE))
    log_msg(sprintf(
      "Dropped %d/%d rate rows lacking UN pop data with min age: %s",
      nrow(unmatched_pop_in_rates), nrow(rates_raw),
      ifelse(is.na(min_age_unmatched), "NA", as.character(min_age_unmatched))
    ))
    rates2 <- semi_join(rates_raw, pop_data, by = matching_vars)  
    
    # Ensure all pop data cohorts appear in rate table
    cohort_template <- pop_data |>  select(all_of(matching_vars)) |>  distinct()
    rates3 <- cohort_template |> 
      left_join(rates2, by = matching_vars) |> 
      left_join(pop_data |> select(all_of(matching_vars), pop_size), by = matching_vars) |> 
      left_join(country_groups |>  select(iso3, who_region), by = "iso3") |> 
      mutate(age_high = if_else(is.na(age_high), age_low + 4, age_high)) |> 
      filter(age_high >= config$min_age_model) |>          # adult model: bands containing ages 18+
      arrange(iso3, female, age_low) 
    
    rate_vars <- names(rates3) |> grep(pattern =
                                         "^(death|incidence|prevalence|prob|le)_(.*)_(val|upper|lower)$", value = TRUE)
    
    rates3 <- rates3 |> 
    # 1. Within country-sex cohort-- nearest‐neighbor fill:
      group_by(iso3, female) |>
      mutate(across(all_of(rate_vars), ~{
        vals  <- .x; ages  <- age_low
        known <- which(!is.na(vals))
        if (length(known) == 0) return(vals)
        missing_idx <- which(is.na(vals))
        for (i in missing_idx) {
          # compute number of 5-year steps away
          d <- abs(ages[known] - ages[i]) / 5
          j <- which.min(d)
          # fill if within 3 bands (<=15 years)
          if (d[j] <= 3) {vals[i] <- vals[known[j]]}
        }
        vals
      }, .names = "{.col}")) |> 
      ungroup()
    
    # apply country inclusion setting:
    if (config$country_inclusion == "has_rate_data") {
      iso_all            <- sort(unique(pop_data$iso3))
      iso_with_rates     <- sort(unique(rates_raw$iso3))
      
      rates3 <- rates3 |> 
        filter(iso3 %in% iso_with_rates)
      }
    
    # 2. Impute other missing values by coalescing regional averages
    
    # Compute averages for imputation:
    ## global: 
    global_averages <- rates3 |> 
      group_by(age_low, age_high, female) |> 
      summarize(
        across(all_of(rate_vars),
               ~ mean(.x, na.rm =TRUE),
               .names = "{.col}_global_avg"), .groups = "drop")
    ## regional:  
    regional_averages <- rates3 |> 
      group_by(age_low, age_high, female, who_region) |> 
      summarize(
        across(all_of(rate_vars),
               ~ mean(.x, na.rm =TRUE),
               .names = "{.col}_region_avg"), .groups = "drop")
    
    ### Coalesce from different sources ----
    processed_rates_data <- rates3 |>   
      mutate(across(all_of(rate_vars), 
                    .names = "{.col}_raw", ~ .x)) |>
      left_join(global_averages, by = c("age_low", "age_high", "female")) |> 
      left_join(regional_averages, by = c("age_low", "age_high", "female", "who_region")) |>
      mutate(across(all_of(rate_vars),
                    ~ coalesce(.x, get(paste0(cur_column(), "_region_avg")), 
                               get(paste0(cur_column(), "_global_avg"))),
                    .names = "{.col}")) |> 
      # source used:
      mutate(across(all_of(rate_vars),
                    .fns = list(imputed_by = ~ case_when(
                      is.na(get(str_c(cur_column(), "_raw"))) & 
                        !is.na(get(str_c(cur_column(), "_region_avg"))) ~ "regional",
                      is.na(get(str_c(cur_column(), "_raw"))) & 
                        is.na(get(str_c(cur_column(), "_region_avg"))) &
                        !is.na(get(str_c(cur_column(), "_global_avg"))) ~ "global",
                      TRUE ~ NA_character_
                    )),
                    .names = "{.col}_{.fn}")) |>
      select(iso3, female, age_low, age_high, who_region, everything(), -ends_with("_region_avg"), -ends_with("_global_avg"), -ends_with("_raw"))
    
    # summary of imputation
    log_msg("For missing rate data cohorts: ")
    msg_lines <- processed_rates_data |>
      pivot_longer(ends_with("_imputed_by"),
                   names_to  = "flag_col",values_to = "imputed_by") |> 
      drop_na(imputed_by) |> 
      mutate(
        rate_var = str_remove(flag_col, "_imputed_by$"),
        category = case_when(
          str_starts(rate_var, "prevalence_") ~ "prevalence",
          str_starts(rate_var, "incidence_")  ~ "incidence",
          str_starts(rate_var, "death_")      ~ "death",
          TRUE                                ~ "other"
        ),
        imputed_by = factor(imputed_by, levels = c("regional", "global"))
      ) |>
      distinct(iso3, category, imputed_by) |>
      count(category, imputed_by, name = "country_count") |>
      complete(category, imputed_by, fill = list(country_count = 0L)) |>
      pivot_wider(
        names_from = imputed_by,
        values_from = country_count,
        values_fill = list(country_count = 0L)
      ) |>
      mutate(
        line = sprintf(
          "%s: %d country cohorts needed regional, %d needed global",
          str_to_title(category), regional, global
        )
      ) |>
      pull(line)
    
    log_msg(paste(msg_lines, collapse = "\n"))
    # unique countries w/ prev/inc imputation
    inc_prev_countries <- processed_rates_data  |> 
      pivot_longer(cols      = ends_with("_imputed_by"),
                   names_to  = "flag_col",
                   values_to = "imputed_by")  |> 
      drop_na(imputed_by) |>
      mutate(rate_var = str_remove(flag_col, "_imputed_by$")) |> 
      filter(str_starts(rate_var, "incidence_") |
               str_starts(rate_var, "prevalence_") |
               str_starts(rate_var, "death_"))|> 
      distinct(iso3) |> 
      pull(iso3)
    
    log_msg(length(inc_prev_countries), "Countries with any cohort level regional or global rate imputation:\n",
        paste(sort(inc_prev_countries), collapse = ", "), "\n")
    
    ### Imputation for unmatched UN pop data (iso3 & cohorts not in gbd/ncdr) ------------------
    
    # Identify rows of UN data that don't have matches in rates data 
    unmatched_popdata <- pop_data |> anti_join(rates_raw, by = matching_vars) |>  # cohort level unmatched age- sex- country- groups
      left_join(country_groups |> select(iso3, who_region), by = "iso3") 
    # Check for iso3 missing in rates and no who_region
    todrop_temp <- unmatched_popdata |> filter(is.na(who_region))
    if (nrow(todrop_temp) > 0) {
      dropped_iso3 <- unique(todrop_temp$iso3)
      log_msg( "Dropping", length(dropped_iso3),
           "countries not in rate data and no who_region:\n  ",
           paste(sort(dropped_iso3), collapse = ", "), "\n")
      unmatched_popdata <- unmatched_popdata |> filter(!is.na(who_region))
    }
    
    ### Combine and final clean merged data---------------------------------------------
    
    pop_cc_cohort <- processed_rates_data |> 
      rename_with(~ str_replace(., "^death_all_cause", "incidence_ac_death"), starts_with("death_all_cause")) |> 
      select(
        -any_of(c("ci_lb", "ci_ub", "iso2", "location", "who_region")),
        -ends_with("_imputed_by")
      )
    
    # debugging check:
    existing_rate_vars <- intersect(rate_vars, names(pop_cc_cohort))
    
    if (anyNA(pop_cc_cohort[, existing_rate_vars, drop = FALSE])) {
      log_msg("Warning: NA values found in rate_vars.\n")}
    
    # Save pop_cc_cohort:
    write.csv(pop_cc_cohort, file.path(paths$working,  "pop_cc_cohort.csv"), row.names = FALSE)
    log_msg("Main input file saved to working folder: pop_cc_cohort.csv")
    
    ## Summary missingness/imputation at country level -------
    
    iso_all            <- sort(unique(pop_data$iso3))
    iso_gbd            <- sort(unique(gbd_data_prevalence$iso3))
    iso_ncdr           <- sort(unique(obesity_data$iso3))
    missing_in_ncdr    <- setdiff(iso_gbd,  iso_ncdr)
    missing_in_gbd     <- setdiff(iso_ncdr, iso_gbd)
    iso_with_rates     <- sort(unique(rates_raw$iso3))
    iso_missing_rates  <- setdiff(iso_all, iso_with_rates)
    iso_missing_region <- setdiff(iso_all, country_groups$iso3)
    log_msg("### Country coverage summary ###")
    log_msg(sprintf("Total in UN pop data:          %3d countries", length(iso_all)))
    log_msg(sprintf("Total in GBD:                  %3d countries", length(iso_gbd)))
    log_msg(sprintf("Total in NCD-risc:             %3d countries", length(iso_ncdr)))
    log_msg("   In NCD-risc but not in GBD:")
    log_msg(paste(missing_in_gbd, collapse = ", "))
    log_msg("   In GBD but not NCD-risc:")
    log_msg(paste(missing_in_ncdr, collapse = ", "))
    log_msg(sprintf("Combined with any rates:       %3d countries", length(iso_with_rates)))
    log_msg(sprintf(
      "  Missing *all* rates:   %3d: %s",
      length(iso_missing_rates),
      paste(iso_missing_rates, collapse = ", ")
    ))
    log_msg(sprintf(
      "  Missing WHO region:    %3d: %s",
      length(iso_missing_region),
      paste(iso_missing_region, collapse = ", ")
    ))
    u1 <- n_distinct(pop_cc_cohort$iso3)
    u2 <- n_distinct(pop_data$iso3)
    log_msg("Countries in final / UN pop:", u1, "/", u2)
    
    # any duplicate (iso3,female,age_low)?
    dups <- pop_cc_cohort %>% 
      count(iso3, female, age_low) %>% 
      filter(n>1)
    
    if (nrow(dups)>0) {
      log_msg("Duplicates found in pop_cc_cohort:")
      log_df(dups)
    }
    return(invisible(pop_cc_cohort))
  } else { return(read.csv(file.path(paths$working,  "pop_cc_cohort.csv" )))}
}

# 2. generate weighted stratified baseline sample -------------
#    one file per country-sex-5yr age cohort; n = min(pop, max(n_min, ceiling(f x pop))); weight w = pop / n
cohort_sample_n <- function(pop, f, n_min) pmin(pop, pmax(n_min, ceiling(f * pop)))

initialize_baseline <- function(cohort_input_data, paths, config = config_inputs, output_dir = NULL){
  inc_prev_input_params <- read_excel(file.path(paths$data, "GEM_input_parameters.xlsx"), sheet = "baseline_dataset") |>
    select(parameter, value)
  country_income_groups_2 <- fread(file.path(paths$working, "country_income_groups_2.csv"), header = T)
   setkey(country_income_groups_2, iso3)
  param <- setNames(as.numeric(inc_prev_input_params$value), inc_prev_input_params$parameter)
  
  # Helper functions -------------
  # condition prevalence by t2d status
  assign_prev_by_t2d_status <- function(dt, cond, prevalence_field, target_t2d, cohort_data) {
    r <- cohort_data[["prevalence_t2d_val"]]
    p_total <- cohort_data[[prevalence_field]]
    # non-t2d probability so that overall prevalence remains p_total
    p_non <- if (r < 1) (p_total - r * target_t2d) / (1 - r) else p_total
    p_non <- pmin(1, pmax(0, p_non))
    dt[t2d == 0, (cond) := as.integer(runif(.N) <= p_non)]
    dt[t2d == 1, (cond) := as.integer(runif(.N) <= target_t2d)]
  }
  #---------------------------------------------'
  create_person_data <- function(start_pid, cohort_data) {
    age_low  <- as.numeric(cohort_data[["age_low"]][1])
    age_high <- as.numeric(cohort_data[["age_high"]][1])
    if (age_high < config$min_age_model) return(NULL)
    ages_ok  <- seq(max(age_low, config$min_age_model), age_high)          # integer ages in band >= min age
    pop_eff  <- cohort_data$pop_size * length(ages_ok) / length(age_low:age_high)   # uniform ages within band
    n        <- cohort_sample_n(round(pop_eff), config$sample_prop, config$min_cohort_n)
    if (n <= 0) return(NULL)
    
    # income group and probability for ESKD
    current_income_group <- country_income_groups_2[J(cohort_data$iso3[1]), income_group]
    if (is.na(current_income_group)) current_income_group <- "not_assigned"
    prob_eskd <- switch(current_income_group,
      lmic = param[["prob_eskd_lmic"]], hic = param[["prob_eskd_hic"]], param[["prob_eskd_na"]])
    #---------------------------'
    chunk_dt <- data.table(
      pid     = seq.int(from = start_pid, length.out = n),
      age     = ages_ok[sample.int(length(ages_ok), n, replace = TRUE)],
      iso3    = cohort_data$iso3,
      female  = cohort_data$female,
      w       = pop_eff / n,                       # survey weight
      stroke  = 0L, eskd = 0L, ac_death = 0L)
    
    ### t2d ----
    r_t2d <- cohort_data[["prevalence_t2d_val"]]
    chunk_dt[, t2d := as.integer(runif(.N) <= r_t2d)]
    
    ### BMI percentile (u) with obesity-t2d association ----
    #   P(obese | t2d) = proportion_t2d_obese; non-t2d set so overall obesity prevalence is preserved;
    #   u uniform within obese / non-obese ranges -> marginal BMI distribution matches NCD-RisC
    p_o <- min(max(1 - cohort_data$bmi_F4, 0), 1)                              # P(BMI >= 30)
    a   <- if (r_t2d > 0) min(param[["proportion_t2d_obese"]], p_o / r_t2d) else 0
    b   <- if (r_t2d < 1) min(max((p_o - r_t2d * a) / (1 - r_t2d), 0), 1) else 0
    chunk_dt[, ob_ := as.integer(runif(.N) < fifelse(t2d == 1, a, b))]
    chunk_dt[, bmi_u := fifelse(ob_ == 1, (1 - p_o) + runif(.N) * p_o, runif(.N) * (1 - p_o))][, ob_ := NULL]
    mp <- cohort_data[rep(1, n)]
    chunk_dt[, bmi := bmi_from_u(bmi_u, mp)]
    set_obesity_from_bmi(chunk_dt, seq_len(n))
    
    ### cvd and ckd --------
    assign_prev_by_t2d_status(chunk_dt, "cvd", "prevalence_cvd_val", param[["proportion_t2d_cvd"]], cohort_data)
    assign_prev_by_t2d_status(chunk_dt, "ckd", "prevalence_ckd_val", param[["proportion_t2d_ckd"]], cohort_data)
    # ESKD among CKD ------
    chunk_dt[ckd == 1, eskd := as.integer(runif(.N) < prob_eskd)]
    chunk_dt
  }
  #---------------------------------------------'
  # Set up ------------------------
  cohort_prev_temp <- as.data.table(cohort_input_data)[, !grep("^(incidence_|death_)", names(cohort_input_data)), with = FALSE
                                                       ][order(iso3, -pop_size)]
  add_bmi_dist_params(cohort_prev_temp)
  if (is.null(output_dir)) output_dir <- file.path(paths$local_temp, paste0("baseline_sample_f", config$sample_prop, "_nmin", config$min_cohort_n, subset_tag(config)))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  n_cohorts  <- nrow(cohort_prev_temp)
  existing   <- list.files(output_dir, pattern = "^cohort_[0-9]+\\.parquet$")
  idx <- setdiff(seq_len(n_cohorts), as.integer(sub("^cohort_([0-9]+)\\.parquet$", "\\1", existing)))
  if (length(idx) == 0L) {
    message("Already generated baseline sample files; not re-running")
  } else {
    message("Generating ", length(idx), " cohort sample files")
    n_i <- cohort_sample_n(round(cohort_prev_temp$pop_size), config$sample_prop, config$min_cohort_n)
    cohort_prev_temp[, start_pid_idx := shift(cumsum(as.numeric(n_i)), fill = 0) + 1]
    future_lapply(idx, function(i) {
      person_data <- create_person_data(cohort_prev_temp$start_pid_idx[i], cohort_prev_temp[i])
      if (!is.null(person_data)) write_parquet(person_data, file.path(output_dir, paste0("cohort_", i, ".parquet")))
      NULL
    }, future.seed = config$seed)
    gc(verbose = FALSE)
  }
  invisible(output_dir)
}
