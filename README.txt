*This readme file was generated on 2025-11-04 and last updated on 2026-06-16 by Elizabeth Mitchell.*


GENERAL INFORMATION

Title of Package: Global Epidemiology Multimorbidity (GEM) Microsimulation Model

Author Information
Name: Elizabeth Mitchell
ORCID:0009-0005-7176-7465
Institution: Emory University
Address: Atlanta, GA
Email: estaton@emory.edu / ec.mitchell3@gmail.com

Co-Author Information
Name: Hui Shao
ORCID:0000-0002-4088-546   
Institution: Emory University
Address: Atlanta, GA
Email: hui.shao@emory.edu


- Date of data collection: January 2025 - May 2026
- Type of data: published estimates, cohort studies, epidemiological and population data
- Date of initial model development: 2025-2026
- Information about EM & HS funding sources that supported the project: T32HL130025 (NHLBI), P30DK111024 (NIDDK), R01DK133465 (NIDDK), and 1R01DK143456 (NIDDK)

*************************************************

SHARING/ACCESS INFORMATION

- Licenses/restrictions placed on the data: Input data from public sources
- Links to publications that cite or use the data: https://diabetesjournals.org/diabetes/article/74/Supplement_1/2043-LB/159453 
- Links to other publicly accessible locations of the data:
	KEY SOURCES:
		1. United Nations Population Division. Home Page | Data Portal [Internet]. Population Division Data Portal. Available from: 
		 	https://population.un.org/dataportal/home?df=f48f8e56-4fa5-4a22-bbd4-8b066cb8b1b9
		2. Global Burden of Disease (GBD) [Internet]. Available from: https://www.healthdata.org/research-analysis/gbd-data
			To download data from GBD: 
			Link: https://vizhub.healthdata.org/gbd-results/
		3. NCD-RisC. National Adult BMI > Data Download [Internet]. Available from: https://ncdrisc.org/data-downloads-adiposity.html
	
		* FOR OTHERS SEE FILE "GEM_documentation_workbook.xlsx"

*************************************************
		
DATA & FILE OVERVIEW

	 ------------------------
	|  	INPUT DATA 	 |
	 ------------------------

  - Location: [file path root]"\01_Data"
  - Necessary files:
	- who_world_bank_country_groups_iso.csv
			Country ISO3 code with WHO region, World Bank country income category 
	- un_pop_data_2024.csv
			Data download from UN Data Portal
	- NCD_RisC_Lancet_2024_BMI_male_age_specific_country.csv
			Data download from NCD-RisC by sex
	- NCD_RisC_Lancet_2024_BMI_female_age_specific_country.csv
			Data download from NCD-RisC by sex
 	- NCD_RisC_Lancet_2024_BMI_child_adolescent_country.csv
			Data download from NCD-RisC for children/adolescents
	- IHME-GBD_prev_inc_male_2023.csv
			Data download from GBD on disease prevalence and incidence rates by sex
	- IHME-GBD_prev_inc_female_2023.csv
			Data download from GBD on disease prevalence and incidence rates by sex
        - IHME-GBD_mortality_by_cause_2023.csv
			Data download from GBD on mortality by cause rates
	- IHME_GBD_disability_weights_2023.csv
			Data download from GBD on disability weights on scale from 0 (full health) to 1 (death)
	- IHME-GBD_life_expectancy_2023.csv
			Data download from GBD on life expectancy
	- IHME-GBD_life_table_2023.csv
			Data download from GBD on life table for years of life lost

	- GEM_input_parameters.xlsx
			Contains model input parameters used throughout the GEM microsimulation framework. Key worksheets include:
				baseline_dataset — baseline epidemiologic and demographic parameters used to initialize simulation cohorts
				semaglutide_treatment_effects — treatment efficacy parameters, including effects on BMI and disease outcomes, used in semaglutide access simulations
				*** RISK RATIO INPUTS ***
					obesity_rr_nejm — obesity-related risk estimates derived from published literature.
					obesity_rr_acdeath_byclass — obesity class-specific relative risks for all-cause mortality.
					obesity_rr_acdeath_bysex — sex-specific obesity relative risks for all-cause mortality.
					obesity_rr_acdeath_byage — age-specific obesity relative risks for all-cause mortality.
					obesity_rr_cvd — obesity relative risks for cardiovascular disease.
					obesity_rr_ckd — obesity relative risks for chronic kidney disease.
					t2d_rr — relative risks associated with type 2 diabetes.
					cvd_rr_stroke_bysex — sex-specific relative risks linking cardiovascular disease and stroke.
				*** COSTING INPUTS ***
					rwe_sema_costs — semaglutide pricing estimates derived from real-world evidence sources.
					us_health_costs — U.S. disease-specific healthcare cost estimates.
					alt_sema_costs — alternative semaglutide pricing scenarios used in sensitivity analyses.
					intl_health_costs — international disease-specific healthcare cost estimates.
					PPP conversion rate_LCUperID — purchasing power parity conversion factors used for international cost standardization.

	- "Country categorizations by iso3.xlsm"
			Contains country-level geographic and economic classifications used throughout the GEM model. 
			The primary worksheet is "manipulated", which includes standardized mappings for:

				iso3 — ISO 3166-1 alpha-3 country code
				who_region — World Health Organization region
				world_bank_income_group — World Bank income classification
				world_bank_income_group_lookup — World Bank income classification, different formatting from index(match())
	- OTHER COSTING DATA FILES:
		- ppp_LCU.xlsx
				Purchasing power parity (PPP) conversion factors used to convert local currency units into internationally comparable purchasing-power-adjusted values
		- WB_LCU_per_USD.xlsx
				World Bank exchange rate data providing local currency units (LCU) per U.S. dollar, used for currency conversions and cost standardization
		- WHO_GHED_CHE_percap_PPP_2024.xlsx
				World Health Organization (WHO) Global Health Expenditure Database estimates of per-capita current health expenditure (CHE) in PPP-adjusted international dollars for 2024, used to contextualize healthcare spending across countries

	- VALIDATION INPUT DATA:
		- GEM_data_validation_workbook.xlsx
				Contains reference study information and extracted data used in preliminary model validation

  - Note:  (1) Users can download and use different years of population and epidemiological data to model different years. 
	       It is recommended to use the same file naming conventions and specify  the "*_data_yr" variables in the config function.
	   
	 ------------------------
	|     MODEL CODE FILES	 |
	 ------------------------

  - Location: [file path root]"\00_Code" 
	Subfolders: 
		"\Python data vis"  ----- Python codes to generate figures including bar charts and circle plots of change over time by region or country
		"\R"  ------------------- R code files to set up environment, prepare input data, and run main model
		"\scripts"  ------------- R code that takes user input and runs functions defined in code in "\R" folder

*************************************************

METHODOLOGICAL INFORMATION 

Software Requirements:
	The primary analyses were conducted using:

		* R version: 4.5.2 (2025-10-31 ucrt)

	Key R package dependencies are documented in the analysis scripts.

	Data visualization and figure generation were conducted using:

		* Python version: 3.12.7 (64-bit, Windows)

	Key Python packages include:
		numpy
		pandas
		matplotlib
		Additional standard library modules used include: os, datetime, re

*************************************************