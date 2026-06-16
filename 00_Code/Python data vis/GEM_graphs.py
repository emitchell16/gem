# graphing code for GEM model
# date: 4/3/2026
# author: liz mitchell


# ********  SET UP ********
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime

# directories
data_folder     = r'C:/Users/Liz/OneDrive - Emory University/1_Research proj/GEM Model_shared/01_Data/'
working_folder  = r'C:/Users/Liz/OneDrive - Emory University/1_Research proj/GEM Model_shared/02_Working/'
py_codes_folder = r'C:/Users/Liz/OneDrive - Emory University/1_Research proj/GEM Model_shared/00_Code/Python data vis/'
out_folder      = r'C:/Users/Liz/OneDrive - Emory University/1_Research proj/Hui Shao Group/Global impact of GLP-1RAs/03_Out_326/'

# inputs
shortdate           = datetime.now().strftime("%m.%d")
who_region_with_pop = pd.read_excel(os.path.join(working_folder, "who_region_with_pop.xlsx"))
region_pop          = who_region_with_pop.set_index('who_region')['pop_size'].to_dict()

country_results = pd.read_excel(os.path.join(out_folder, "RES_access_scenarios", "100glp1_0.1sample_5yr_3.28",
                                              "sim_result_iso3_100glp1.xlsx"))
country_icer_results = pd.read_excel(os.path.join(out_folder, 
                                                  "country_level_icer_res_100_0.1.xlsx"))
country_icer_results["n"] *= 10 # rescale n from 10% sample
region_results = pd.read_excel(os.path.join(out_folder, 
                                            "region_level_uncertainty_int_res_100_0.1.xlsx"))

print(region_results.columns.tolist())

## DEBUGGING DALY AND COST * 100 in R: (remove after regenerating summary post-sim results)

DALY_COST_COLS_REGION_COUNTRY = [
    'daly_total_adj', 'yld_adj_total', 'yll_adj',
    'cost_total_intl_rwe', 'cost_total_id_sust_costx1', 'cost_total_id_sust_costx3',
    'cost_total_id_lowest_mp', 'cost_total_id_generic_entry', 'cost_conditions_total'
]
DALY_COST_COLS_ICER = [
    'incremental_cost_val', 'incremental_cost_lower', 'incremental_cost_upper',
    'cost_total_val', 'cost_total_lower', 'cost_total_upper',
    'yll_adj_val', 'yll_adj_lower', 'yll_adj_upper',
    'yld_adj_total_val', 'yld_adj_total_lower', 'yld_adj_total_upper',
    'daly_total_adj_val', 'daly_averted_val',
    'daly_total_adj_lower', 'daly_averted_lower',
    'daly_total_adj_upper', 'daly_averted_upper'
]
DALY_COST_ABS_DIFF_COLS = [
    'abs_diff_pct',
    'abs_diff_pct_lb',
    'abs_diff_pct_ub'
]

def divide_cols_by_100(df, cols, label=""):
    df = df.copy()
    found = [c for c in cols if c in df.columns]
    missing = [c for c in cols if c not in df.columns]
    for col in found:
        df[col] = pd.to_numeric(df[col], errors='coerce') / 100
    print(f"  [{label}] Divided by 100: {found}")
    if missing:
        print(f"  [{label}] WARNING - cols not found: {missing}")
    return df

# region: filter to target conditions then divide
target_conditions = set(DALY_COST_COLS_REGION_COUNTRY)
mask = region_results['condition'].isin(target_conditions)
region_results.loc[mask, [c for c in DALY_COST_COLS_REGION_COUNTRY if c in region_results.columns]] /= 100
print(f"  [region_results] Divided by 100 for conditions: {region_results.loc[mask, 'condition'].unique().tolist()}")

abs_cols_present = [c for c in DALY_COST_ABS_DIFF_COLS if c in region_results.columns]
region_results.loc[mask, abs_cols_present] /= 100

# country: same cols, all rows (wide format — cols are outcomes not conditions)
country_results      = divide_cols_by_100(country_results,      DALY_COST_COLS_REGION_COUNTRY, "country_results")
country_icer_results = divide_cols_by_100(country_icer_results, DALY_COST_COLS_ICER,           "country_icer_results")
## DEBUGGING end

iso_region_assignment_data = pd.read_csv(working_folder + 'who_world_bank_country_groups_iso_cleaned.csv')
country_results['who_region'] = country_results['iso3'].map(iso_region_assignment_data.set_index('iso3')['who_region'])
"""
# 1. bar chart (mortality, t2d, cvd, ckd, stroke) ----------------------------------------
"""
with open(os.path.join(py_codes_folder,  'GEM_bar_chart.py'), 'r') as file:
    exec(file.read())

create_grouped_barchart(
    region_results, 
    'ac_death', 't2d_new','cvd_new', 'ckd_new',
    max_rate = 11, 
    figsize=(9, 8),
    file_name = f"FIG_crm_outcome_groupbar_{shortdate}.png")
 
"""
create_grouped_barchart(
    region_results, 
    'daly_total_adj', #"yld_adj_total", "yll_adj",  these just show yll drives daly
    #max_rate = 11, 
    figsize=(9, 8),
    file_name = f"FIG_daly_outcome_groupbar_{shortdate}.pdf")
"""

# 2. circle plot (obesity) ----------------------------------------
with open(os.path.join(py_codes_folder,  'GEM_circle_plot.py'), 'r') as file:
    exec(file.read())

create_circleplot_1fig(
    data= country_results,
    comparison_data = region_results, # this has _lb and _ub
    outcomes=['obesity'],
    y_min=0,
    y_max=0.69, # prevalence rate max in plot
    figsize=(6,8),
    file_name =  f"FIG_cirlceplot_obesity_{shortdate}.png" # can do pdf or png by changing extension and format = in file
)

"""  
# 3. icer / cost-daly plots ----------------------------------------

with open(os.path.join(py_codes_folder,  'GEM_cost_effect_plot.py'), 'r') as file:
    exec(file.read())

#plot_ce_plane(["intl_rwe", "id_sust_costx1"], "ce_plane_rwe_costbased.png")
# plot_ce_plane(["intl_rwe", "id_sust_costx1"], "ce_plane_rwe_costbased_v2.png", share_y_ax=False)
plot_ce_plane(["id_sust_costx1", "id_sust_costx3", "intl_rwe" ],
               "ce_plane_rwe_costbased_3panelV2.png", 
              share_y_ax=False)
              
              ## RERUN LAST: 5.7.26
""" 
plot_country_icer(
    scenarios  = ["intl_rwe", "id_sust_costx1"],
    out_file   = "FIG_iso3_icer2.png",
)

plot_incomegroup_icer(
    scenarios = ["intl_rwe", "id_sust_costx1", "id_sust_costx3", "tc_zero"],
    out_file  = "FIG_incomelvl_icer2.png",
    y_cap     = 60,
)
"""

# code graveyard:
## tested versions: 
# plot_ce_plane(["id_sust_costx1", "id_sust_costx3"], "ce_plane_costbased_2.png")  # note: 3x cost based prices not *that* different 
# plot_ce_plane(["tc_zero"],                    "ce_plane_zerocost.png")
