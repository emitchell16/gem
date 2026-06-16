# circle plot function

import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
import re

def darken(color, factor=0.6):
    """Darken a matplotlib color by multiplying RGB by factor < 1."""
    rgb = mcolors.to_rgb(color)
    return tuple(factor * c for c in rgb)

def fmt_val(x, sign=False):
    """Format to 1 decimal place, or 1 sig fig if |x| < 0.1."""
    if abs(x) >= 0.1 or x == 0:
        result = f"{x:.1f}"
    else:
        result = f"{x:.1g}"
    if sign and x > 0:
        result = '+' + result
    return result

def create_circleplot_1fig(
    data: pd.DataFrame, 
    comparison_data: pd.DataFrame,
    outcomes: list[str],                   
    y_min: float = 0, 
    y_max: float = 0.6, 
    base_color: str = 'darkslateblue', 
    intervention_color: str = 'skyblue',   
    figsize: tuple = (12, 6),
    file_name: str = 'circle_plot.png'
):
    """
    Plot base and intervention simulation results for each iso3 by region in circle plot panels by outcome
    """
    # Create color map for each region
    unique_regions = list(data.sort_values(by=['who_region'])['who_region'].unique())
    cmap_region = plt.colormaps['tab10']
    region_colors = {region: cmap_region(i/len(unique_regions)) for i, region in enumerate(unique_regions)}

    # Create subplots 1 outcome version
    fig = plt.figure(figsize=figsize)
    gs  = fig.add_gridspec(nrows=3, ncols=1,  height_ratios=[2, 9, 2])  # top plot, small footer
    header_ax = fig.add_subplot(gs[0]); header_ax.axis('off')
    ax        = fig.add_subplot(gs[1], projection='polar')  
    footer_ax = fig.add_subplot(gs[2]); footer_ax.axis('off')
    gs.update(hspace=0.3) # adjust spacing between fig, legend 

    # Consistent font sizes
    body_fsize = 9

    for outcome in  outcomes:
        # Order iso3 by outcome and region
        rate_varname = "rate__"+ outcome 
        
        global_order = list(
            data
            .sort_values(by=['who_region', rate_varname], ascending=[True, False])
            ['iso3'].unique())
        
        subset = (
            data
            .loc[data['sim'].isin(['5yr_basecase','5yr_100glp1'])]
            .copy()
            )
        subset['iso3'] = pd.Categorical(
            subset['iso3'],
            categories=global_order,
            ordered=True
        )
        subset.sort_values(['iso3'], inplace=True)
        # extract final order and plot
        iso3_order = subset['iso3'].cat.categories
              
        wide = (
            subset
            .pivot(index='iso3', columns='sim', values=rate_varname)
            .reindex(iso3_order))
        
        # build a map ISO3 → WHO region
        iso3_to_region = (
            subset[['iso3', 'who_region']]
            .drop_duplicates()
            .set_index('iso3')['who_region']
            .to_dict())

        # add it to wide
        wide['region'] = [iso3_to_region[i] for i in wide.index]
        iso3_order = list(wide.index)
        n = len(iso3_order)
        angles = np.linspace(0, 2*np.pi, n, endpoint=False)
        width  = 2*np.pi/n * 0.9
       
        # Region formatting
        for region, group_df in wide.groupby('region'):
            # Get the positions (indices in sorted order) for this region
            pos = [wide.index.get_loc(idx_val) for idx_val in group_df.index]
            # Compute the angular span
            region_start = angles[min(pos)]
            region_end   = angles[max(pos)] + width  # right edge of last bar
            
            # Create a set of angles for the polygon
            # If the region does not wrap around 2π:
            if region_end <= 2*np.pi:
                theta = np.linspace(region_start, region_end, 100)
                # Outer edge at y_max and inner edge at y_min.
                theta_full = np.concatenate([theta, theta[::-1]])
                r_full = np.concatenate([np.full_like(theta, y_max), np.full_like(theta, y_min)])
                ax.fill(theta_full, r_full, color=region_colors[region], alpha=0.15, zorder=0)
                # For label placement, use midpoint:
                region_center = (region_start + region_end) / 2.
            else:
                # If wrapping around, split into two segments.
                theta1 = np.linspace(region_start, 2*np.pi, 50)
                theta2 = np.linspace(0, region_end - 2*np.pi, 50)
                theta_full = np.concatenate([theta1, theta2[::-1]])
                r_full = np.concatenate([np.full_like(theta1, y_max), np.full_like(theta2, y_min)])
                ax.fill(theta_full, r_full, color=region_colors[region], alpha=0.15, zorder=0)
                # For label, choose a midpoint near the split.
                region_center = (region_start + (region_end - 2*np.pi)) / 2.
                if region_center < 0:
                    region_center += 2*np.pi
             
            # Region labels from comparison df: AC, RC, Prevalences
            region_row = (
                comparison_data[
                    (comparison_data['who_region']
                     .str.lower()
                     == region) &
                    (comparison_data['condition'] 
                     .str.lower()
                     .str.replace(r'[\s_]+', '', regex=True)
                     == re.sub(r'[\s_]+', '', outcome.lower()))
                ]
                .reset_index(drop=True)
            )
            
            y5_rate_sq = region_row.loc[0, 'outcome_value_end_base']
            y5_rate_int = region_row.loc[0, 'outcome_value_end_int']
            rel_dp = region_row.loc[0, 'rel_diff_pct']
            abs_dp = region_row.loc[0, 'abs_diff_pct']
            region_rate_label = f"{y5_rate_sq:.1f}% status quo vs. {y5_rate_int:.1f}% intervention;"
            rc_label = fmt_val(rel_dp, sign=True) + '%'
            ac_label = fmt_val(abs_dp, sign=True) + 'pp'
                    
            # uncertainty intervals
            rel_low = region_row.loc[0, 'rel_diff_pct_lb']
            rel_high = region_row.loc[0, 'rel_diff_pct_ub']
            abs_low = region_row.loc[0, 'abs_diff_pct_lb']
            abs_high = region_row.loc[0, 'abs_diff_pct_ub']
            rc_si_label = fmt_val(rel_low, sign=True) + '% to ' + fmt_val(rel_high, sign=True) + '%'
            ac_si_label = fmt_val(abs_low, sign=True) + ' to ' + fmt_val(abs_high, sign=True)

            impact_line = f"{region_rate_label}\nRD: {rc_label} ({rc_si_label});\nAD: {ac_label} ({ac_si_label})"
        
            # orient labels
            span = (region_end - region_start) % (2*np.pi)
            region_center = (region_start + span/2.0) % (2*np.pi)
            medium_pad = {"africa"}
            big_pad = {"eastern mediterranean", "western pacific", "europe"}
            r_pad_base = 0.59 if region.lower() in big_pad else 0.4 if region.lower() in medium_pad else 0.3
            r_pad_overrides = {
                "western pacific":  0.6,
                "europe":           0.45,
                "americas":         0.3,
            }
            r_pad    = r_pad_overrides.get(region.lower(), r_pad_base)
            r_anchor = y_max + r_pad
            offsets = {
                "europe":   -0.25,
                "americas": 0.1,
                "western pacific":  0.2
            }
            region_center += offsets.get(region.lower(), 0)
            header = r"$\mathbf{" + region.title().replace('-', r'\text{-}').replace(' ', r'\ ') + r"}$"
            block  = header + "\n" + impact_line

            ax.text(
                region_center, r_anchor, block,
                ha='center', va='center',
                fontsize=body_fsize, linespacing=1.05,
                color=darken(region_colors[region], 0.8),
                clip_on=False,
                # bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="lightgray", lw=1, alpha=0.8)
            )

        # Plot formatting    
        ax.bar(angles, wide['5yr_basecase'], width=width, color=base_color, zorder=2)
        ax.bar(angles, wide['5yr_100glp1'],  width=width, color=intervention_color, zorder=2)
        
        # Label each angle with its corresponding iso3 code (color-coded by region)
        label_colors = [region_colors[r] for r in wide['region']] 
        
        ax.set_xticks([])

        for iso, theta, region in zip(
            iso3_order,
            angles + width/2,           # center of each bar
            wide['region']):
            r_label = y_max + 0.02      # push just outside the circle; tweak as needed
            rot = np.degrees(theta)

            # For labels on the left/bottom side, flip so text isn't upside down
            if 90 < rot < 270:
                rot = rot + 180
                ha = 'right'
            else:
                ha = 'left'

            ax.text(
                theta,
                r_label,
                iso,
                rotation=rot,               
                rotation_mode='anchor',
                ha=ha,
                va='center',
                fontsize=3.8,
                color=region_colors[region],
                clip_on=False,
            )
        
        # Set radial limits and title
        ax.set_ylim(y_min, y_max)
        ax.yaxis.set_major_locator(mticker.MultipleLocator(0.1))
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0, decimals=0))
        ax.yaxis.grid(True, which='major', color='whitesmoke', linestyle='-', linewidth=0.45, zorder= 1) 
        ax.xaxis.grid(True, which='major', color='whitesmoke', linestyle='-', linewidth=0.45, zorder= 1) 
  
        ax.set_rlabel_position(45)  # degrees: 0=right, 90=top, 180=left, 270=bottom
        ax.set_axisbelow(True)
        for label in ax.get_yticklabels():
            label.set_fontsize(7)
            label.set_color('#333333')
            label.set_bbox(dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7))
       
    # Legend:
    base_patch = mpatches.Patch(facecolor=base_color, label='Prevalence Under Status Quo')
    int_patch  = mpatches.Patch(facecolor=intervention_color, label='Prevalence Under Universal Access')
    footer_ax.legend(handles = [base_patch, int_patch],  
            title = "Scenario:",
            loc='center', 
            bbox_to_anchor=(0.5, 0.03), 
            ncol=2,
            fontsize=body_fsize,
            title_fontsize=body_fsize, 
            frameon=True,
            prop={'weight': 'normal'})
    leg = footer_ax.get_legend()
    if leg and leg.get_title():
        leg.get_title().set_fontweight('bold')

    # Global change labels 
    global_row = (
        comparison_data[
            (comparison_data['who_region'] == 'global') &
            (comparison_data['condition'] 
                .str.lower()
                .str.replace(r'[\s_]+', '', regex=True)
                == re.sub(r'[\s_]+', '', outcome.lower()))
        ]
        .reset_index(drop=True)
    )
    y5_rate_sq_g  = global_row.loc[0, 'outcome_value_end_base'] #TODO fix to actual var names
    y5_rate_int_g = global_row.loc[0, 'outcome_value_end_int']
    rel_g = global_row.loc[0, 'rel_diff_pct']
    abs_g = global_row.loc[0, 'abs_diff_pct']
    global_rate_label = f"{y5_rate_sq_g:.1f}% status quo vs. {y5_rate_int_g:.1f}% intervention;"
 
    rel_low  = global_row.loc[0, 'rel_diff_pct_lb']
    rel_high = global_row.loc[0, 'rel_diff_pct_ub']
    abs_low  = global_row.loc[0, 'abs_diff_pct_lb']
    abs_high = global_row.loc[0, 'abs_diff_pct_ub']
    # format
    rc_si_label = fmt_val(rel_low, sign=True) + '% to ' + fmt_val(rel_high, sign=True) + '%'
    ac_si_label = fmt_val(abs_low, sign=True) + ' to ' + fmt_val(abs_high, sign=True)

    global_body = f"{global_rate_label}\nRD: {fmt_val(rel_g, sign=True)}% ({rc_si_label});\nAD: {fmt_val(abs_g, sign=True)}pp ({ac_si_label})"

    # placement of global change label  
    ta = TextArea(
        r"$\mathbf{Global}$" + "\n" + global_body,
        textprops=dict(ha="center", fontsize=body_fsize)
    )
    ab = AnnotationBbox(
        ta, (0.5, 0.98),                         
        xycoords=header_ax.transAxes,
        box_alignment=(0.5, 1.0), 
        frameon=True,
        bboxprops=dict(boxstyle="square,pad=0.25", fc="white", ec="lightgray", lw=1),
        zorder=5)
    
    header_ax.add_artist(ab)

    fig.suptitle(
    "Figure 1. Comparison of Simulated 5-Year Obesity Prevalence "
    "With and Without Universal Semaglutide Access",
    fontsize=11,          
    fontweight="bold",
    #fontfamily="serif",
    y=0.98                )
    
    #plt.tight_layout()
    #plt.savefig(f"{out_folder}/{file_name}", dpi=320, bbox_inches='tight', format = 'pdf')
    plt.savefig(f"{out_folder}/{file_name}", dpi=320, bbox_inches='tight', format = 'png')
    plt.show()