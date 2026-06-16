# bar chart function

import math
import matplotlib.colors as mcolors
from matplotlib.ticker import PercentFormatter
from matplotlib.patches import Rectangle, Patch
import matplotlib.ticker as mticker

def darken(color, factor=0.6):
    rgb = mcolors.to_rgb(color)
    return tuple(factor * c for c in rgb)

def _fmt_diff(val, sign=False):
    """1 decimal place; 2 decimal for small values that would round to 0.0; <0.01 for near-zero."""
    fmt_sign = '+' if sign else ''
    abs_val = abs(val)
    if abs_val < 0.01:
        return ("<0.01" if val >= 0 else ">-0.01")
    elif abs_val < 0.1:
        return f"{val:{fmt_sign}.2f}" 
    else: 
        return f"{val:{fmt_sign}.1f}"

def create_grouped_barchart(
    comparison_data: pd.DataFrame,  
    *outcomes,
    max_rate: float = 20, 
    basecase_color: str = 'darkslateblue',
    intervention_color: str = 'skyblue', 
    figsize: tuple = (6, 6),
    file_name: str = 'groupbar.png'
):

    panel_outcomes = [o for o in outcomes if o is not None]

    label_map = {
        "ac_death":        "All-Cause Mortality",
        "obesity":         "Obesity",
        "class1_obesity":  "Obesity class 1",
        "class2_obesity":  "Obesity class 2",
        "class3_obesity":  "Obesity class 3",
        "eskd_new":        "ESKD",
        "cvd_new":         "CVD",
        "ckd_new":         "CKD",
        "stroke_new":      "Stroke",
        "t2d_new":         "T2D",
        "daly":            "DALYs (millions)",
        "yld":             "YLDs (millions)",
        "yll":             "YLLs (millions)",
    }
    # Outcomes whose values should be divided by 1,000,000 for display
    scale_1m = {"daly", "yld", "yll"}
  
    body_fsize = 8
    heading_fsize = 9

    # Sort and color regions
    comparison_data['who_region'] = (
        comparison_data['who_region'].astype(str).str.lower()
    )

    # build the region order using lowercased names
    all_regions = sorted(comparison_data['who_region'].unique().tolist())
    non_global  = [r for r in all_regions if r != 'global']
    regions     = ['global'] + non_global


    # Create a color map for each region & global category
    cmap = plt.get_cmap('tab10')
    palette = cmap(np.linspace(0, 1, len(regions)))

    region_colors = {r: palette[i] for i, r in enumerate(non_global)}
    region_colors['global'] = palette[len(non_global)]
    region_colors = {k: darken(v, 0.8) for k, v in region_colors.items()}

    # layout
    n_panels = len(panel_outcomes)
    ncols = max(1, int(math.ceil(math.sqrt(n_panels))))
    nrows = int(math.ceil(n_panels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes_flat = list(axes.flat)
    

    # plot each panel
    for i, (ax, outcome) in enumerate(zip(axes_flat, panel_outcomes)):
        col_idx = i % ncols
        # match outcome input list to condition values in data
        if outcome == "ac_death": 
            cond_name = "all-cause mortality"
        elif outcome.endswith("_new"):
            cond_name = outcome[:-4]
        else:  
            cond_name = outcome
        
        # Subset to this outcome 
        mask = comparison_data['condition'].str.lower() == outcome.lower()
        df_out = (
            comparison_data
            .loc[mask]
            .set_index('who_region')
            .reindex(regions)
        )

        base_vals = df_out['outcome_value_end_base'].to_numpy()
        int_vals  = df_out['outcome_value_end_int'].to_numpy()
        diff_pct_rel   = df_out['rel_diff_pct'].to_numpy()
        diff_pct_abs  = df_out['abs_diff_pct'].to_numpy()

        outcome_lc = outcome.lower()
        is_count = any(k in outcome_lc for k in ("daly", "yld", "yll"))
        if any(k in outcome_lc for k in scale_1m):
            base_vals = base_vals / 1_000_000
            int_vals  = int_vals / 1_000_000

        # bars
        bar_width = 0.36
        bar_gap = 0.04
        bar_offset = (bar_width + bar_gap)/2
        y = np.arange(len(regions))
        ax.barh(y - bar_offset, base_vals,   bar_width, color=basecase_color)
        ax.barh(y + bar_offset, int_vals,    bar_width, color=intervention_color)
        ax.invert_yaxis() # reorder regions

        # annotate change
        for yi in range(len(regions)):
            region = regions[yi]  
            abs_dp = diff_pct_abs[yi]
            rel_dp = diff_pct_rel[yi]
            rel_str = _fmt_diff(rel_dp) + "%"
            rel_low,  rel_high  = df_out.loc[region, ['rel_diff_pct_lb','rel_diff_pct_ub']]
            abs_low,  abs_high  = df_out.loc[region, ['abs_diff_pct_lb','abs_diff_pct_ub']]
            rc_si_label = f"{_fmt_diff(rel_low, sign=True)}% to {_fmt_diff(rel_high, sign=True)}%"
            if is_count:
                abs_dp_s  = abs_dp / 1_000_000
                abs_low_s = abs_low / 1_000_000
                abs_high_s = abs_high / 1_000_000
                abs_str = _fmt_diff(abs_dp_s) + "M"
                ac_si_label = f"{_fmt_diff(abs_low_s, sign=True)}M to {_fmt_diff(abs_high_s, sign=True)}M"
            else:
                abs_str = _fmt_diff(abs_dp) + " pp"
                ac_si_label = f"{_fmt_diff(abs_low, sign=True)} to {_fmt_diff(abs_high, sign=True)}"

            label_text =(
                f" RD: {rel_str} ({rc_si_label});\n "
                f"AD: {abs_str} ({ac_si_label})"
            )

            bar_max = max(base_vals[yi], int_vals[yi])
            x_pos = bar_max + (bar_max * 0.02 if is_count else 0.02)
            y_pos = y[yi]
            ax.text(x_pos, y_pos, label_text, ha='left', va='center', fontsize=body_fsize-1)
                
        # styling
        panel_title = next((v for k, v in label_map.items() if k in outcome_lc),
                           outcome.replace("_"," ").title())
        ax.set_title(panel_title, fontsize=heading_fsize, fontweight='bold')

        if is_count:
            ax.set_xlabel(panel_title, fontsize=body_fsize)
            x_max = max(np.nanmax(base_vals), np.nanmax(int_vals))
            ax.set_xlim(0, x_max * 1.6)
            ax.xaxis.set_major_formatter(mticker.FuncFormatter(
                lambda x, pos: f"{x:,.0f}"))
        elif outcome == 'eskd_new':
            ax.set_xlabel("5-Year Incidence", fontsize=body_fsize)
            ax.set_xlim(0, 0.25)
            ax.xaxis.set_major_formatter(mticker.FuncFormatter(
                lambda x, pos: f"{x:.2f}%"))
        else:
            ax.set_xlabel("5-Year Incidence", fontsize=body_fsize)
            ax.set_xlim(0, max_rate)
            ax.xaxis.set_major_formatter(mticker.FuncFormatter(
                lambda x, pos: f"{x:.1f}%"))
        ax.tick_params(axis='x', labelsize=body_fsize )
        # set and color region labels
        wrapped_labels = []
        for r in regions:
            title = r.title()
            if " " in title:
                first, rest = title.split(" ", 1)
                wrapped_labels.append(f"{first}\n{rest}")
            else:
                wrapped_labels.append(title)
        ax.set_yticks(y)
        ax.set_yticklabels(wrapped_labels, fontsize=body_fsize, fontweight='bold')

        # color each label to match its region
        for yi, lbl in enumerate(ax.get_yticklabels()):
            lbl.set_color(region_colors[regions[yi]])

        # Remove y-axis labels for non-leftmost columns
        if col_idx != 0:
            ax.set_yticklabels([])
            ax.set_ylabel(None)
            ax.tick_params(axis="y", length=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(True)
        ax.spines["bottom"].set_visible(True)

    # Legend:
    base_patch = Patch(color=basecase_color, label='Cumulative Incidence Under Status Quo')
    int_patch  = Patch(color=intervention_color, label='Cumulative Incidence Under Universal Access')
    leg = fig.legend(handles = [base_patch, int_patch],  
                     title = "Scenario:", loc="lower center",
                     bbox_to_anchor=(0.5, 0.04), ncol=2,
                     fontsize=body_fsize, frameon=True, 
                     prop={'size': body_fsize, 'weight': 'normal'})
    if leg.get_title():
        leg.get_title().set_fontweight('bold')
                    
    axes_flat = axes.flat  
    for ax in list(axes_flat)[n_panels:]:
        ax.axis('off')
        
    fig.subplots_adjust(
        left=0.03, right=0.98,
        top=0.88, bottom=0.2,
        hspace=0.35, wspace=0.2
        )
    any_burden = any(any(k in o.lower() for k in ("daly","yld","yll")) for o in panel_outcomes)
    suptitle_text = (
        "Figure 2. Comparison of Simulated 5-Year Disease Burden\n"
        "With and Without Universal Semaglutide Access by Region"
        if any_burden else
        "Figure 2. Comparison of Simulated 5-Year Mortality and Cardio-Renal-Metabolic\n"
        "Condition Incidence With and Without Universal Semaglutide Access by Region"
    )
    fig.suptitle(
        suptitle_text,
        fontsize=heading_fsize+.5,          
        fontweight="bold",
        #fontfamily="serif",
        y=0.98 , x=0.4    )
    st = fig._suptitle 

    output_file = os.path.join(out_folder, file_name)
    plt.savefig(output_file, dpi=320, format = 'png', bbox_inches="tight",
    pad_inches=0.1, bbox_extra_artists=[leg, st] if st is not None else [leg])
    plt.show()
