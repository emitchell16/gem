# 00_Code/R/prep_main_inputs.R
#TODO add life expectancy table load and clean from GBD

#***********
#* 1. Reads UN, NCD-RisC, GBD data files in 01_Data folder -------------
#* writes log to paths$out/1data_cleaning_log*.txt
#***********

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
  if (!file.exists(dw_file)) { stop("disability weight file not found for year: ", config$GBD_data_yr)}
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
  if (!file.exists(pop_file)) { stop("UN population file not found for year: ", config$UN_data_yr)}
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
      world_bank_income_group = tolower(gsub("[[:punct:]]", "", world_bank_income_group)),
      who_region = tolower(who_region)
    ) |> 
    filter(!is.na(iso3), !grepl("^XX[1-9]$|^X1[0-3]$", iso3))
  
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

  if (!file.exists(file.path(paths$working,  "pop_cc_cohort.csv"))) { 
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
        data_adolescent <- read_csv(file.path(paths$data,"NCD_RisC_Lancet_2024_BMI_child_adolescent_country.csv"), show_col_types = FALSE) |> 
          rename(age = `Age group` ) |> 
          filter(age >= 10) |> 
          mutate(`Age group` =  case_when( 
            age >= 10 & age <= 14 ~ "10-14",
            age >= 15 & age <= 19 ~ "15-19",
            TRUE ~ NA_character_)) |>
          filter(!is.na(`Age group`)) |> 
          group_by(`ISO`, `Sex`, `Year`, `Age group`) |> 
          summarise(
            across(matches("prevalence", ignore.case = TRUE), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop"
          ) 
  
        data <- bind_rows(data_male, data_female, data_adolescent) |> 
          rename_with(tolower) |> 
          filter(year %in% c(data_year, data_year-1),`age group` != "18-19") |> 
          select(year, sex, iso, `age group`, 
                 matches("\\(obesity\\)"),
                 matches("prevalence of bmi 30"), 
                 matches("prevalence of bmi 35"),
                 matches("prevalence of bmi >=40")
          ) |> 
          mutate(`age group` = if_else(`age group` == "85plus", "85-89", `age group`)) |> 
          separate(`age group`, into = c("age_low", "age_high"), sep = "-", convert = TRUE) |>
          rename_with(~ case_when(
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
            # cols for adolescent obesity thresholds by SD (WHO standard)
            str_detect(., "> 2sd") & str_detect(., "lower") ~ "prevalence_obesity_adolsc_lower",
            str_detect(., "> 2sd") & str_detect(., "upper") ~ "prevalence_obesity_adolsc_upper",
            str_detect(., "> 2sd") ~ "prevalence_obesity_adolsc_val",
            TRUE ~ .
          )) |> 
          rename(iso3 = iso) |> 
          mutate(female = if_else(sex %in% c("Women", "Girls"), 1, 0)) |> 
          select(-sex) |>
          #  use prevalence rates from year n-1 -> n to approximate cohort-level incidence rates
          group_by(age_low, age_high, female, iso3) |> 
          arrange(year) |> 
          mutate(
            # Apply lagged subtraction for all relevant columns (_val, _lower, _upper)
            across(c(
                prevalence_obesity_val, prevalence_obesity_lower, prevalence_obesity_upper,
                prevalence_class1_obesity_val, prevalence_class1_obesity_lower, prevalence_class1_obesity_upper,
                prevalence_class2_obesity_val, prevalence_class2_obesity_lower, prevalence_class2_obesity_upper,
                prevalence_class3_obesity_val, prevalence_class3_obesity_lower, prevalence_class3_obesity_upper,
                prevalence_obesity_adolsc_val, prevalence_obesity_adolsc_lower, prevalence_obesity_adolsc_upper
              ),
              ~ . - lag(.),
              .names = "{gsub('prevalence_', 'incidence_', .col)}"
            )
          ) |> 
          ungroup() |> 
          filter( year == data_year ) |> 
          select(-year)
        
        # summary of proxy incidence rates for obesity
        summary_stats <- data |> 
          summarize(
            mean = mean(incidence_obesity_val, na.rm = TRUE),
            median = median(incidence_obesity_val, na.rm = TRUE),
            range = paste0(range(incidence_obesity_val, na.rm = TRUE), collapse = " to "),
            negative_count = sum(incidence_obesity_val < 0, na.rm = TRUE),
            total_count = sum(!is.na(incidence_obesity_val))
          )
        
        log_msg("Obesity proxy incidence summary:")
        log_df(summary_stats)
        
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
      mutate(age_high = if_else(is.na(age_high), age_low + 4, age_high), 
             # clean t2d 0s for age 10 which should be NA
             across(matches("t2d"), ~ ifelse(age_low == 10 & (.x == 0 | .x == 0.0), NA, .x)),
             prevalence_obesity_val = coalesce(
               prevalence_obesity_val,
               prevalence_obesity_adolsc_val
             ),
             prevalence_obesity_lower = coalesce(
               prevalence_obesity_lower,
               prevalence_obesity_adolsc_lower
             ),
             prevalence_obesity_upper = coalesce(
               prevalence_obesity_upper,
               prevalence_obesity_adolsc_upper
             )
             ) |> 
      select(-contains("adolsc_"))|>
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

#***********
#* 2. Generate person level dataset -------------
#***********
initialize_baseline <- function(cohort_input_data, paths, config = config_inputs){
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
    # Compute non-t2d probability so that overall prevalence remains p_total
    p_non <- if (r < 1) (p_total - r * target_t2d) / (1 - r) else p_total
    p_non <- pmin(1, pmax(0, p_non))
    
    # For non-t2d, assign based on p_non
    dt[t2d == 0, (cond) := as.integer(runif(.N) <= p_non)]
    # For t2d individuals, assign based on the target probability
    dt[t2d == 1, (cond) := as.integer(runif(.N) <= target_t2d)]
  }
  #---------------------------------------------'
  create_person_data <- function(start_pid, pop_size, cohort_data) {
    age_low  <- as.numeric(cohort_data[["age_low"]][1])
    age_high <- as.numeric(cohort_data[["age_high"]][1])
    
    # income group and probability for ESKD
    current_iso3 <- cohort_data$iso3[1]
    current_income_group <- country_income_groups_2[J(current_iso3), income_group]
    if (is.na(current_income_group)) current_income_group <- "not_assigned"
  
    prob_eskd_lmic <- param[["prob_eskd_lmic"]]  
    prob_eskd_hic  <- param[["prob_eskd_hic"]]
    prob_eskd_na   <- param[["prob_eskd_na"]]
    prob_eskd <- switch(current_income_group,
      lmic         = prob_eskd_lmic,
      hic          = prob_eskd_hic,
      not_assigned = prob_eskd_na, prob_eskd_na
    )
    #---------------------------'
    
    chunk_dt <- data.table(
      age     = round(runif(pop_size, age_low, age_high)),
      iso3    = rep(cohort_data$iso3, pop_size),
      female  = rep(cohort_data$female, pop_size),
      stroke  = 0L,
      eskd    = 0L,
      ac_death= 0L
    )
    
    # Filter out simulated age below min_age_model
    chunk_dt <- chunk_dt[age >= config$min_age_model] 
    pop_size <- nrow(chunk_dt) 
    chunk_dt[, pid := seq.int(from = start_pid, length.out = pop_size)]
    
    ## Add binary conditions -----
    
    ### t2d ----
    chunk_dt[, t2d := as.integer(runif(.N) <= cohort_data[["prevalence_t2d_val"]])]
    ###  obesity class assignment ----
    #' Computes mutually‑exclusive class‑1/2/3 probabilities conditional on t2d status and assigns 1 class pp
    prev_c1 <- cohort_data[["prevalence_class1_obesity_val"]]
    prev_c2 <- cohort_data[["prevalence_class2_obesity_val"]]
    prev_c3 <- cohort_data[["prevalence_class3_obesity_val"]]
    # overall obesity prevalence & class‐shares
    p_obs_total <- prev_c1 + prev_c2 + prev_c3
    w <- c(prev_c1, prev_c2, prev_c3) / p_obs_total
    r_t2d <- cohort_data[["prevalence_t2d_val"]]
    p_obs_t2d <- param[["proportion_t2d_obese"]]
    p_obs_non <- (p_obs_total - r_t2d * p_obs_t2d) / (1 - r_t2d)
    p_obs_non <- pmin(1, pmax(0, p_obs_non))
    # conditional class prevalence
    p_t2d <- p_obs_t2d * w
    p_non <- p_obs_non * w
    # cumulative thresholds
    thresh_t2d <- cumsum(p_t2d)
    thresh_non <- cumsum(p_non)
    # draw uniform random for each person
    chunk_dt[, u := runif(.N)]
    # assign mutually exclusive classes
    chunk_dt[t2d == 1, `:=`(
      class1_obesity = as.integer(u < thresh_t2d[1]),
      class2_obesity = as.integer(u >= thresh_t2d[1] & u < thresh_t2d[2]),
      class3_obesity = as.integer(u >= thresh_t2d[2] & u < thresh_t2d[3])
    )]
    chunk_dt[t2d == 0, `:=`(
      class1_obesity = as.integer(u < thresh_non[1]),
      class2_obesity = as.integer(u >= thresh_non[1] & u < thresh_non[2]),
      class3_obesity = as.integer(u >= thresh_non[2] & u < thresh_non[3])
    )]
    # binary obesity indicator
    chunk_dt[, obesity := as.integer(class1_obesity | class2_obesity | class3_obesity)]
    chunk_dt[, u := NULL]

    ### cvd and ckd --------
    assign_prev_by_t2d_status(chunk_dt, "cvd", "prevalence_cvd_val", param[["proportion_t2d_cvd"]], cohort_data)
    assign_prev_by_t2d_status(chunk_dt, "ckd", "prevalence_ckd_val", param[["proportion_t2d_ckd"]], cohort_data)
    
    ### BMI ------
    chunk_dt[, bmi := NA_real_] # bmi is background running variable for obese only
    # Assign a continuous BMI within class-specific range if obese
    chunk_dt[class1_obesity == 1, bmi := runif(.N, min = 30, max = 34.99)]
    chunk_dt[class2_obesity == 1, bmi := runif(.N, min = 35, max = 39.99)]
    chunk_dt[class3_obesity == 1, bmi := runif(.N, min = 40, max = 65)] # class 3 defined as >= 40. 65 is above 99th percentile bmi
    # overweight bmi for portion of non-obese with t2d
    non_obese_t2d_idx <- which(chunk_dt$obesity == 0 & chunk_dt$t2d == 1)
    if (length(non_obese_t2d_idx) > 0) {
      overweight_idx <- non_obese_t2d_idx[ 
        runif(length(non_obese_t2d_idx)) <= param[["proportion_t2d_overweight"]]
      ]
      if (length(overweight_idx) > 0) {
        # Assign a BMI uniformly drawn from [25, 30) for overweight individuals (based on WHO bmi category)
        chunk_dt[overweight_idx, bmi := runif(length(overweight_idx), min = 25, max = 29.99)]  # for glp1 treatment interested in 27+
      }
    }
    # ESKD among CKD ------
    chunk_dt[ckd == 1, eskd := as.integer(runif(.N) < prob_eskd)]
    return(chunk_dt)
  }
  #---------------------------------------------'
  # Set up ------------------------
  conditions <- c( "t2d", "cvd", "obesity", "class1_obesity", "class2_obesity", "class3_obesity", "ckd")
  cohort_prev_temp <- as.data.table(cohort_input_data)[, !grep("^(incidence_|death_)", names(cohort_input_data)), with = FALSE
                                                       ][order(iso3, -pop_size)]
  # file management
  output_dir <- file.path(paths$shared_working, directory_name)
  if (!dir.exists(output_dir)) { 
    dir.create(output_dir, recursive = TRUE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  
  scratch_base <- Sys.getenv("TMPDIR")
  if (scratch_base == "") scratch_base <- "/tmp"
  temp_dir <- file.path(scratch_base, paste0("temp_dir_", Sys.getpid()))
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)}
    
    # Process each cohort row ----
  
  n_cohorts <- nrow(cohort_prev_temp)
  n_existing <- length(list.files(output_dir, pattern = "^cohort_[0-9]+\\.parquet$", full.names = FALSE))
  if (n_existing >= n_cohorts) {
    message("Already generated baseline files; not re-running")
  } else {
    cohort_prev_temp[ , start_pid_idx := shift(cumsum(as.numeric(pop_size)), 
                                               fill = 0) + 1 ]
    existing_files <- list.files(output_dir, pattern="^cohort_[0-9]+\\.parquet$", full.names=FALSE)
    existing_dirs  <- list.files(output_dir, pattern="^cohort_[0-9]+$", full.names=FALSE)
    existing_dirs <- existing_dirs[file.info(file.path(output_dir, existing_dirs))$isdir %in% TRUE]
    file_ids <- as.integer(sub("^cohort_([0-9]+)\\.parquet$", "\\1", existing_files))
    dir_ids  <- as.integer(sub("^cohort_([0-9]+)$", "\\1", existing_dirs))
    existing_ids <- sort(unique(c(file_ids, dir_ids)))
    
    # run only missing cohorts
    idx <- setdiff(seq_len(n_cohorts), existing_ids)
    idx <- sort(idx)
    
    if (length(idx) == 0L) { message("No missing cohorts detected")
    } else {message("Generating ", length(idx), " missing cohorts. ")
      
      plan(multisession, workers = 5)
      future_lapply(idx, function(i){
        cohort_data <- cohort_prev_temp[i, ]
        pop_size    <- as.numeric(cohort_data[["pop_size"]][1])
        pid_start   <- cohort_data[["start_pid_idx"]]
        
        if (is.na(pid_start)) {
          stop("pid_start is NA on iteration ", i)
        }
        person_data <-  create_person_data(pid_start, pop_size, cohort_data)
        # chunk write if large cohort
        if (pop_size > 5e5) {
          cohort_temp_dir <- file.path(temp_dir, paste0("cohort_",i, "_chunks"))
          if (!dir.exists(cohort_temp_dir)) {
            dir.create(cohort_temp_dir, recursive = TRUE)
          }
          nrows <- nrow(person_data)
          chunk_size <- 1e5  
          for (row_start in seq(1, nrows, by = chunk_size)) {
            row_end <- min(nrows, row_start + chunk_size - 1)
            chunk_dt <- person_data[row_start:row_end, ]
            chunk_file <- file.path(cohort_temp_dir, paste0("part_", ceiling(row_start/chunk_size), ".parquet"))
            write_parquet(chunk_dt, chunk_file)
          }
          combined_data <- open_dataset(cohort_temp_dir, format = "parquet") |> collect()
          final_file <- file.path(output_dir, paste0("cohort_", i, ".parquet"))
          write_parquet(combined_data, final_file)
          unlink(cohort_temp_dir, recursive = TRUE)
        } else{ 
          final_file <- file.path(output_dir, paste0("cohort_", i, ".parquet"))
          write_parquet(person_data, final_file)
        }
        
        rm(person_data)
        #})
      }, future.seed = config$seed, future.scheduling = 1 )
    }
    
    # Clean up remaining temporary objects
    rm(cohort_prev_temp)
    unlink(temp_dir, recursive = TRUE)
    gc(verbose = FALSE)
    invisible(NULL)
  }
}

#***********
#* 3. Random sample of baseline dataset -------------
#***********
sample_baseline_data <- function() {
  short_date <- config$short_date
  files <- list.files(path = baseline_files_folder, pattern = "\\.parquet$", full.names = TRUE)
  dir.create(sample_baseline_files, recursive = TRUE, showWarnings = FALSE)
  dir.create(local_sample_files, recursive = TRUE, showWarnings = FALSE)
  existing_cohorts <- list.files(local_sample_files, pattern = 
                                   "^cohort_[0-9]+\\.parquet$", full.names = FALSE)
  dest_ids <- as.integer(sub("^cohort_([0-9]+)\\.parquet$", "\\1", existing_cohorts))
  remaining <- setdiff(seq_along(files), dest_ids)
  
  if (length(remaining) > nrow(cohort_input_data)) stop("alert, n files to run > n cohorts") 
  
  future_lapply(remaining, function(i) {
    out_file <- file.path(local_sample_files, paste0("cohort_", i, ".parquet"))
    if (file.exists(out_file)) return()
    df <- read_parquet(files[i], memory_map = F)
    current_sample_p <- if (nrow(df) < 100) sample_p * 2 else sample_p
    sampled_df <- df |> sample_frac(current_sample_p)
    
    #fix id
    if ("pid" %in% names(sampled_df)) sampled_df$pid <- NULL
    sampled_df$pid_local <- seq_len(nrow(sampled_df))
    sampled_df$pid <- paste0("c", i, "_", sampled_df$pid_local)
    sampled_df$pid_local <- NULL
    write_parquet(sampled_df, out_file)
  }, future.seed=config$seed )
  
  # copy over to onedrive:
  existing_cohorts <- list.files(sample_baseline_files, pattern = "^cohort_[0-9]+\\.parquet$", full.names = FALSE)
  dest_ids <- as.integer(sub("^cohort_([0-9]+)\\.parquet$", "\\1", existing_cohorts))
  src_all <- list.files(local_sample_files, full.names = TRUE, pattern = "\\.parquet$")
  
  src_ids <- as.integer(sub("^cohort_([0-9]+)\\.parquet$", "\\1", basename(src_all)))
  src <- src_all[!is.na(src_ids) & !(src_ids %in% dest_ids)]
  
  if (length(src) == 0) { message("No new local files to copy")
  } else {
    chunk_size <- 500
    idx <- split(seq_along(src), ceiling(seq_along(src) / chunk_size))
    ok_all <- rep(FALSE, length(src))
    for (g in idx) {
      ok <- file.copy(src[g], sample_baseline_files, overwrite = FALSE)
      if (length(ok) != length(g)) ok <- rep(FALSE, length(g))
      ok_all[g] <- ok
      Sys.sleep(0.001)
    }
    if (!all(ok_all)) {
      failed <- src[!ok_all]
      warning(sprintf("Copy failed for %d / %d files. Example: %s",
        length(failed), length(src), failed[1]))}
  }
    invisible(sample_baseline_files)
  }
  

