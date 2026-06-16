# Cost effectiveness plots
# Translated from R — harmonized with Wong colorblind palette used in GEM_graphs.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# 3.1. CE plane plot:

all_scenario_labels = {
    "id_sust_costx1": "1. Low Cost-Based Prices",
    "id_sust_costx3": "2. Reduced Market Prices",
    "intl_rwe":       "3. Current Market Prices",
    "tc_zero":        "Zero Treatment Cost",
}
# Income group custom order and display labels (low → high)
_ig_order  = ["low income", "lower middle income", "upper middle income", "high income"]
ig_colors  = dict(zip(_ig_order, ["#009E73", "#AA490C", "#E69D00CF", "#0072B2"]))
ig_labels  = {
    "low income":          "Low Income",
    "lower middle income": "Lower Middle Income",
    "upper middle income": "Upper Middle Income",
    "high income":         "High Income",
}

def _add_labels_no_overlap(ax, x_arr, y_arr, sz_arr, labels, fontsize=8, offset=8):
    """Place iso3 labels offset from dots with consistent connector lines."""
    fig = ax.get_figure()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_bbox = ax.get_window_extent(renderer=renderer)
    placed_bboxes = []
    arrow_threshold = offset -5

    candidates = [
        ('center', 'bottom', ( 0,       offset)),
        ('center', 'top',    ( 0,      -offset)),
        ('left',   'center', ( offset,  0)),
        ('right',  'center', (-offset,  0)),
        ('left',   'bottom', ( offset,  offset)),
        ('right',  'bottom', (-offset,  offset)),
        ('left',   'top',    ( offset, -offset)),
        ('right',  'top',    (-offset, -offset)),
    ]

    for i in sorted(range(len(labels)), key=lambda i: -sz_arr[i]):
        for ha, va, (xt, yt) in candidates:
            dist = (xt**2 + yt**2) ** 0.5
            use_arrow = dist > arrow_threshold

            t = ax.annotate(
                labels[i],
                xy=(x_arr[i], y_arr[i]),
                xytext=(xt, yt),
                textcoords='offset points',
                fontsize=fontsize, ha=ha, va=va,
                annotation_clip=True,
                arrowprops=dict(
                    arrowstyle='-',
                    color='grey',
                    lw=0.5,
                    shrinkA=0,
                    shrinkB=3,
                ) if use_arrow else None
            )
            bb = t.get_tightbbox(renderer)
            if bb is None:
                t.remove()
                continue
            outside = (not ax_bbox.contains(bb.x0, bb.y0) or
                       not ax_bbox.contains(bb.x1, bb.y1))
            overlaps = any(bb.overlaps(pb) for pb in placed_bboxes)
            if outside or overlaps:
                t.remove()
            else:
                placed_bboxes.append(bb)
                break

def plot_ce_plane(scenarios, out_file, share_y_ax=True):
    df_ce = country_icer_results[country_icer_results["cost_scenario"].isin(scenarios)].copy()

    # drop countries with no WB income group; report exclusions
    no_ig = df_ce.loc[df_ce["income_group"] == "no wb income group", "iso3"].unique()
    if len(no_ig):
        print(f"[{out_file}] Dropping {len(no_ig)} iso3 with no WB income group: {sorted(no_ig)}")
    df_ce = df_ce[df_ce["income_group"] != "no wb income group"].reset_index(drop=True)

    n_vals = df_ce["n"].values
    df_ce["_dot_size"] = 20 + 180 * (n_vals - n_vals.min()) / (n_vals.max() - n_vals.min() + 1e-9)

    # preserve custom low→high order, only include groups present in data
    income_groups = [ig for ig in _ig_order if ig in df_ce["income_group"].unique()]
    n_panels = len(scenarios)
    fig, axes = plt.subplots(1, n_panels, figsize=(2.5 * n_panels + 1, 5.0), squeeze=False)
    axes = axes[0]

    # pass 1: scatter and collect label data
    ax_label_data = {}
    for ax, scenario in zip(axes, scenarios):
        sub = df_ce[df_ce["cost_scenario"] == scenario]
        all_x, all_y, all_sz, all_iso = [], [], [], []

        for ig, grp in sub.groupby("income_group"):
            x   = (grp["daly_averted_val"] / 1e6).values
            y   = (grp["incremental_cost_val"] / 1e9).values
            sz  = grp["_dot_size"].values
            iso = grp["iso3"].values

            ax.scatter(x, y, s=sz, color=ig_colors.get(ig, "grey"), alpha=0.85, zorder=3)
            all_x.extend(x); all_y.extend(y)
            all_sz.extend(sz); all_iso.extend(iso)

        ax.axhline(0, color='k', lw=0.6)
        ax.axvline(0, color='k', lw=0.6)
        ax.set_title(all_scenario_labels[scenario], fontsize=9)
        ax.set_xlabel("DALYs Averted (Million)", fontsize=9)
        ax.set_ylabel("Incremental Cost (Billion Int$)", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.yaxis.set_major_formatter(plt.matplotlib.ticker.FuncFormatter(lambda x, _: f'{x:,.0f}'))
        ax_label_data[scenario] = (ax, all_x, all_y, all_sz, all_iso)

    # pass 2: y-axis limits — lock autoscale immediately so WTP lines can't expand range
    if share_y_ax and n_panels > 1:
        y_min = min(a.get_ylim()[0] for a in axes)
        y_max = max(a.get_ylim()[1] for a in axes)
        for ax in axes:
            ax.set_ylim(y_min, y_max)
    elif not share_y_ax:
        # cost-based panels: fixed range [-20, 20] B; report clipped countries
        cost_based = {"id_sust_costx1", "id_sust_costx3", "tc_zero"}
        cb_ymin, cb_ymax = -1, 1
        for ax, sc in zip(axes, scenarios):
            if sc in cost_based:
                clipped = df_ce[df_ce["cost_scenario"] == sc].copy()
                clipped = clipped[
                    (clipped["incremental_cost_val"] / 1e9 < cb_ymin) |
                    (clipped["incremental_cost_val"] / 1e9 > cb_ymax)
                ]
                if not clipped.empty:
                    print(f"[{out_file} | {sc}] {len(clipped)} countries clipped by y=[{cb_ymin},{cb_ymax}]B: "
                          f"{sorted(clipped['iso3'].tolist())}")
                ax.set_ylim(cb_ymin, cb_ymax)
                ax.set_yticks([cb_ymin, -0.5, 0, 0.5, cb_ymax])
                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:g}'))
                ax.set_autoscale_on(False) #TODO update here when change min/max
    for ax in axes:
        ax.set_autoscale_on(False)

    # WTP threshold lines — (value $/DALY, linestyle, color)
    # slope units: WTP / 1000 converts ($/DALY) × (M DALYs) → B USD
    cost_based_set = {"id_sust_costx1", "id_sust_costx3", "tc_zero"}
    wtp_rw         = [(50_000, '--', "#7E7EFC"), (5_000, ':', '#999999')]
    wtp_costbased  = [(50_000, '--', "#7E7EFC"), (5_000,  ':',  '#999999'), (1_000,  '--', '#17becf')]

    all_wtp_lines = []  # union for legend
    for ax, sc in zip(axes, scenarios):
        lines = wtp_costbased if sc in cost_based_set else wtp_rw
        x0, x1 = ax.get_xlim()
        for wtp, ls, col in lines:
            slope = wtp / 1_000
            ax.plot([x0, x1], [slope * x0, slope * x1],
                    color=col, lw=1, linestyle=ls, zorder=1)
        for entry in lines:
            if entry not in all_wtp_lines:
                all_wtp_lines.append(entry)

    # pass 3: place iso3 labels
    for scenario, (ax, all_x, all_y, all_sz, all_iso) in ax_label_data.items():
        _add_labels_no_overlap(ax, all_x, all_y, all_sz, all_iso)

    # two-row legend with subheadings; no frame
    ig_handles = [plt.Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=ig_colors[ig], markersize=7, label=ig_labels[ig])
                  for ig in income_groups]
    wtp_handles = [plt.Line2D([0], [0], color=col, lw=1.2, linestyle=ls,
                               label=f'${wtp:,}/DALY')
                   for wtp, ls, col in sorted(all_wtp_lines, key=lambda x: -x[0])]
    leg1 = fig.legend(handles=ig_handles, loc="lower center",
                      bbox_to_anchor=(0.5, 0.08), ncol=len(income_groups),
                      fontsize=9, frameon=False,
                      title="Country/Territory Income Group", title_fontsize=9)
    leg2 = fig.legend(handles=wtp_handles, loc="lower center",
                      bbox_to_anchor=(0.5, 0.01), ncol=len(wtp_handles),
                      fontsize=9, frameon=False,
                      title="Willingness to Pay Threshold", title_fontsize=9)
    
    for leg in (leg1, leg2):
        leg.get_title().set_fontweight('bold')
    fig.add_artist(leg1)

    plt.suptitle("Cost-Effectiveness Plane by Cost Scenario*", fontsize=9)
    plt.tight_layout(rect=[0, 0.18, 1, 0.97])
    plt.savefig(os.path.join(out_folder, out_file), dpi=300, bbox_inches='tight')
    plt.show()

# tab10 palette for scenarios (distinct from Wong income-group colors)
scenario_colors = {
    "intl_rwe":       "#ff7f0e", 
    "id_sust_costx1": "#1f77b4",  
    "id_sust_costx3": "#2ca02c",  # tab:green
    "tc_zero":        "#9467bd",  # tab:purple
}

_ig_order = ["low income", "lower middle income", "upper middle income", "high income"]
ig_display = {
    "low income":          "Low Income",
    "lower middle income": "Lower Middle Income",
    "upper middle income": "Upper Middle Income",
    "high income":         "High Income",
}

def _prep_icer_df(scenarios):
    """Filter country_icer_results, drop no-WB-group, add display columns."""
    df = country_icer_results[country_icer_results["cost_scenario"].isin(scenarios)].copy()
    no_ig = df.loc[df["income_group"] == "no wb income group", "iso3"].unique()
    if len(no_ig):
        print(f"ICER plot: dropping {len(no_ig)} iso3 with no WB income group: {sorted(no_ig)}")
    df = df[df["income_group"] != "no wb income group"].copy()
    df["icer_k"] = df["icer_val"] / 1e3
    df["ig_disp"] = pd.Categorical(
        df["income_group"].map(ig_display),
        categories=[ig_display[ig] for ig in _ig_order if ig in df["income_group"].unique()],
        ordered=True,
    )
    df["sc_disp"] = df["cost_scenario"].map(all_scenario_labels)
    return df


# ── plot 1: country-level ICER, faceted by cost scenario ───────────────────
def plot_country_icer(scenarios, out_file, icer_k_floor=-20, outlier_thresh=200,
                      label_icer_k=25, label_n=5e8):
    df = _prep_icer_df(scenarios)
    df = df[df["icer_k"] > icer_k_floor].copy()

    # report outliers
    outliers = df[df["icer_k"] > outlier_thresh]
    if not outliers.empty:
        print(f"\n[{out_file}] Outliers (ICER > {outlier_thresh}k):")
        print(outliers.groupby(["sc_disp", "ig_disp"])
                      .agg(n_outliers=("icer_k", "count"), max_icer=("icer_k", "max"))
                      .to_string())
    df = df[df["icer_k"] <= outlier_thresh].copy()

    present_scenarios = [s for s in scenarios if s in df["cost_scenario"].unique()]
    n_panels = len(present_scenarios)
    fig, axes = plt.subplots(1, n_panels, figsize=(4.5 * n_panels, 5), squeeze=False)
    axes = axes[0]

    ig_cats = [ig_display[ig] for ig in _ig_order if ig in df["income_group"].unique()]
    ig_pos  = {ig: i for i, ig in enumerate(ig_cats)}
    np.random.seed(42)

    for ax, sc in zip(axes, present_scenarios):
        sub = df[df["cost_scenario"] == sc].copy()
        sub = sub.dropna(subset=["ig_disp"])
        color = scenario_colors.get(sc, "#888888")

        # dot size scaled by n
        n_vals = sub["n"].values
        sz = 20 + 150 * (n_vals - n_vals.min()) / (n_vals.max() - n_vals.min() + 1e-9)

        x_pos = sub["ig_disp"].map(ig_pos).values.astype(float)
        x_pos += np.random.uniform(-0.15, 0.15, len(x_pos))

        ax.scatter(x_pos, sub["icer_k"].values, s=sz,
                   color=color, alpha=0.7, zorder=3)

        # label notable countries
        mask = (sub["icer_k"] > label_icer_k) | (sub["n"] > label_n)
        for _, row in sub[mask].iterrows():
            xi = ig_pos[row["ig_disp"]] + np.random.uniform(-0.1, 0.1)
            ax.annotate(row["iso3"], (xi, row["icer_k"]),
                        fontsize=6, ha="center", va="bottom",
                        xytext=(0, 3), textcoords="offset points")

        # reference lines
        for yref in [1, 50, 100]:
            ax.axhline(yref, linestyle="--", color="grey", alpha=0.5, lw=0.8)
        ax.axhline(0, color="k", lw=0.8, alpha=0.8)

        ax.set_xticks(range(len(ig_cats)))
        ax.set_xticklabels(ig_cats, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("ICER\n(1,000 Int$ per DALY averted)", fontsize=8)
        ax.set_xlabel(None)
        ax.set_title(all_scenario_labels[sc], fontsize=9)
        ax.tick_params(labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.suptitle("Country-Level ICER by Income Group", fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, out_file), dpi=300, bbox_inches="tight")
    plt.show()


# ── plot 2: income-group summary ICER, dodged by scenario ──────────────────
def plot_incomegroup_icer(scenarios, out_file, y_cap=60):
    df = _prep_icer_df(scenarios)

    # aggregate to income-group level
    agg = (df.groupby(["income_group", "ig_disp", "cost_scenario", "sc_disp"], observed=True)
             .agg(
                 incremental_cost_val   =("incremental_cost_val",   "sum"),
                 incremental_cost_lower =("incremental_cost_lower",  "sum"),
                 incremental_cost_upper =("incremental_cost_upper",  "sum"),
                 daly_averted_val       =("daly_averted_val",        "sum"),
                 daly_averted_lower     =("daly_averted_lower",      "sum"),
                 daly_averted_upper     =("daly_averted_upper",      "sum"),
             )
             .reset_index())
    agg["icer_k"]       = agg["incremental_cost_val"]   / agg["daly_averted_val"]   / 1e3
    agg["icer_k_lower"] = agg["incremental_cost_lower"] / agg["daly_averted_lower"] / 1e3
    agg["icer_k_upper"] = agg["incremental_cost_upper"] / agg["daly_averted_upper"] / 1e3

    ig_cats = [ig_display[ig] for ig in _ig_order if ig in agg["income_group"].unique()]
    ig_pos  = {ig: i for i, ig in enumerate(ig_cats)}
    present_scenarios = [s for s in scenarios if s in agg["cost_scenario"].unique()]
    n_sc = len(present_scenarios)
    dodge_offsets = np.linspace(-0.2, 0.2, n_sc)

    fig, ax = plt.subplots(figsize=(7, 5))

    for sc, dx in zip(present_scenarios, dodge_offsets):
        sub = agg[agg["cost_scenario"] == sc].copy()
        sub = sub.dropna(subset=["ig_disp"])
        color = scenario_colors.get(sc, "#888888")
        x = sub["ig_disp"].map(ig_pos).values.astype(float) + dx

        ax.scatter(x, sub["icer_k"].values, color=color, s=60, alpha=0.85,
                   zorder=4, label=all_scenario_labels[sc])

    # reference lines and annotations
    for yref, lbl in [(1, "1k Int$/DALY"), (5, "5k Int$/DALY"), (50, "50k Int$/DALY")]:
        ax.axhline(yref, linestyle="--", color="grey", alpha=0.5, lw=0.8)
        ax.annotate(lbl, xy=(len(ig_cats) - 0.5, yref),
                    fontsize=7, color="grey", va="bottom", ha="right")
    ax.axhline(0, color="k", lw=0.8, alpha=0.8)

    # vertical grid lines between groups
    for xi in np.arange(0.5, len(ig_cats) - 0.5, 1):
        ax.axvline(xi, color="grey", lw=0.4, alpha=0.4, zorder=1)

    ax.set_xticks(range(len(ig_cats)))
    ax.set_xticklabels(ig_cats, rotation=30, ha="right", fontsize=8)
    ax.set_ylim(bottom=agg["icer_k"].min() - 2, top=y_cap)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:g}'))
    ax.set_ylabel("ICER\n(1,000 Int$ per DALY averted)", fontsize=9)
    ax.set_xlabel(None)
    ax.set_title("Income Group Summary ICER by Cost Scenario", fontsize=10)
    ax.tick_params(labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(title="Semaglutide Cost Scenario", fontsize=7, title_fontsize=8,
              frameon=False, loc="upper left")

    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, out_file), dpi=300, bbox_inches="tight")
    plt.show()