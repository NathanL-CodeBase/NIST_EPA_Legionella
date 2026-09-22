#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle Emission Rate Variant Figures
========================================

Reads particle_analysis_summary.xlsx (all_results sheet) produced by
scripts/particle_decay_analysis.py and builds two families of figures. It
does not recompute p, beta, or E. Run particle_decay_analysis.py first.

1. Per-bin variant-comparison figures (Bokeh, interactive): for each of the
   12 particle-size bins and both beta_other/E_mean, plots one point per
   event per variant over time, with a whisker showing the metric's std, for
   three comparison pairs:
       entry_outside      E(t)2 vs. E(t)1 -- the dual air-change-rate bracket
       outside_adjusted   E(t)3 vs. E(t)1 -- ratio-corrected vs. raw/fleet
                           concentration series (pre-cutover events only,
                           since E(t)3 is NaN after the room-correction cutover)
       outside_bathroom   E(t)4 vs. E(t)1 -- bedroom+bathroom vs. bedroom-only
                           control volume (beta_other is identical for this
                           pair since it has no volume term -- see
                           src/particle_calculations.calculate_other_process_rate)
   A vertical line at 2026-06-03 marks the room-concentration correction
   cutover (src/particle_room_correction.py): before that date the primary
   concentration series is raw C_bed1(t); from that date it is the
   position-weighted MODULAIR-PM fleet average C_room(t). Informational only
   -- it does not split the plotted points into groups.

2. Summary/boxplot/comparison figures (matplotlib, static; moved here from
   scripts/particle_decay_analysis.py): bar-chart summaries of penetration,
   deposition, and emission by bin size; fixed-water-temp-axis and
   continuous-metric-axis boxplots of E_total, beta_other, E_mean, and p; and
   the water-temp/spray-pattern/shower-head/mannequin/door/fan comparison
   boxplot families. All outside/entry (E1/E2) only, matching today's
   reported values; FLOW_RATE_MIN/FLOW_RATE_MAX restrict these figures (not
   the Excel workbook or the per-event pm_decay figures) to events with a
   measured flow rate of 4.1-5.6 LPM.

Style matches the shared MODULAIR Bokeh figures (1600x800, 12pt, no title,
click-to-hide legend, hover enabled) via src.plot_style.style_moduair_figure.

Output Files:
    output/plots/particle/beta_other_entry_outside_bin{N}.html     (12 figures)
    output/plots/particle/emission_entry_outside_bin{N}.html       (12 figures)
    output/plots/particle/beta_other_outside_adjusted_bin{N}.html  (12 figures)
    output/plots/particle/emission_outside_adjusted_bin{N}.html    (12 figures)
    output/plots/particle/beta_other_outside_bathroom_bin{N}.html  (12 figures)
    output/plots/particle/emission_outside_bathroom_bin{N}.html    (12 figures)
    output/plots/penetration_summary.png, deposition_summary_{outside,entry}.png,
        emission_summary_{outside,entry}.png
    output/plots/emission_etotal_boxplot_{bin0-2,bin3-6,bin7-11}_{outside,entry}.png,
        other_process_rate_boxplot_..., emission_rate_boxplot_...,
        penetration_factor_boxplot_{bin0-2,bin3-6,bin7-11}.png
    output/plots/emission_etotal_by_{bedroom_rh,bedroom_temp,acr,beta,p,showerhead}_
        boxplot_{bin0-2,bin3-6,bin7-11}_{outside,entry}.png
    output/plots/{spray_pattern,head_type,mannequin,door_position,fan_status}_
        {emission_etotal,other_process_rate}_boxplot_{bin0-2,bin3-6,bin7-11}_{outside,entry}.png,
        ..._penetration_factor_boxplot_{bin0-2,bin3-6,bin7-11}.png

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-18
Update log:
    2026-09-18  Initial version (entry vs. outside comparison figures only).
    2026-09-22  Renamed from particle_beta_emission_source_figures.py. Added
        the outside_adjusted (E1 vs. E3) and outside_bathroom (E1 vs. E4)
        comparison figure pairs. Folded in particle_decay_analysis.py's
        _generate_summary_plots (summary bar charts, fixed/metric-axis
        boxplots, and the 5 condition-comparison families), which now reads
        all_results from this script's already-loaded DataFrame instead of
        an in-memory results_df from the pipeline script.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from bokeh.layouts import column as bokeh_column
from bokeh.models import ColumnDataSource, DatetimeTickFormatter, Div, HoverTool, Label, Span, Whisker
from bokeh.plotting import figure, output_file, save
from scipy.stats import ttest_rel

# Add project root to path for src/ imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import src.sig_figs as sf  # noqa: E402
from src.data_paths import get_data_root  # noqa: E402
from src.particle_calculations import PARTICLE_BINS  # noqa: E402
from src.particle_room_correction import ROOM_CUTOVER  # noqa: E402
from src.plot_style import MODUAIR_FIGURE_WIDTH, MODUAIR_TEXT_PT, SENSOR_COLORS, style_moduair_figure  # noqa: E402

# =============================================================================
# Configuration: Per-Bin Variant Comparison Figures
# =============================================================================

SUMMARY_XLSX_NAME = "particle_analysis_summary.xlsx"
SUMMARY_SHEET = "all_results"

FIGURE_SUBDIR = ("plots", "particle")

SHOWER_ON_COL = "shower_on"

# Lambda-bracketed outside/entry sheets in particle_analysis_summary.xlsx stay
# on this tuple (matches scripts/particle_decay_analysis.py's _LAMBDA_SOURCES).
_LAMBDA_SOURCES = ("entry", "outside")

# (variant_a, variant_b, filename pair stem). Each pair reuses the same
# make_figure/compute_variant_change/make_summary_div machinery -- only the
# two column names being compared change.
COMPARISON_PAIRS = [
    ("entry", "outside", "entry_outside"),
    ("outside", "adjusted", "outside_adjusted"),
    ("outside", "bathroom", "outside_bathroom"),
]

VARIANT_LEGEND = {
    "outside": "Outside (E1)",
    "entry": "Entry (E2)",
    "adjusted": "Adjusted (E3)",
    "bathroom": "Bathroom (E4)",
}
VARIANT_COLORS = {
    "entry": SENSOR_COLORS[0],
    "outside": SENSOR_COLORS[1],
    "adjusted": SENSOR_COLORS[2],
    "bathroom": SENSOR_COLORS[3],
}

# One metric spec per figure kind: value/std/r2 column templates ({bin} and
# {variant} filled in per figure/variant), axis label, hover value format,
# unit string for the summary Div, and output filename template ({pair} and
# {bin} filled in per comparison pair).
METRICS = {
    "beta": {
        "value_col": "bin{bin}_{variant}_beta_other (h-1)",
        "std_col": "bin{bin}_{variant}_beta_other_std (h-1)",
        "r2_col": "bin{bin}_{variant}_beta_other_r_squared",
        "y_label": "Other process rate β",
        "y_units": "h⁻¹",
        "value_fmt": "0.000",
        "unit": " h⁻¹",
        "filename": "beta_other_{pair}_bin{bin}.html",
    },
    "E": {
        "value_col": "bin{bin}_{variant}_E_mean (#/min)",
        "std_col": "bin{bin}_{variant}_E_std (#/min)",
        "r2_col": "bin{bin}_{variant}_E_r_squared",
        "y_label": "Emission rate E",
        "y_units": "#/min",
        "value_fmt": "0.00e+00",
        "unit": " #/min",
        "filename": "emission_{pair}_bin{bin}.html",
    },
}

LEGEND_LOCATION = "top_right"

# Cutover-region labels either side of the vertical Span (left, right).
CUTOVER_LABELS = ("Raw C_bed1 / C_adjusted room", "C_room")
CUTOVER_LABEL_COLOR = "gray"

# =============================================================================
# Configuration: Summary/Boxplot/Comparison Figures
# =============================================================================

# Flow-rate filter for summary plots.  Events whose measured shower head flow
# rate falls outside this range are excluded from all boxplot and comparison
# figures.  The Excel results file (and the per-event pm_decay figures)
# retain all events regardless of flow rate.
FLOW_RATE_MIN: float = 4.1  # LPM (inclusive)
FLOW_RATE_MAX: float = 5.6  # LPM (inclusive)


# =============================================================================
# Data Loading
# =============================================================================


def load_summary(summary_path: Path) -> pd.DataFrame:
    """
    Load the all_results sheet and parse the shower_on column.

    Parameters
    ----------
    summary_path : Path
        Path to particle_analysis_summary.xlsx from particle_decay_analysis.py.

    Returns
    -------
    pd.DataFrame
        all_results rows with shower_on coerced to datetime. Rows without a
        valid shower_on time are dropped, since shower_on is the figure x
        position.
    """
    if not summary_path.exists():
        raise FileNotFoundError(
            f"Summary workbook not found: {summary_path}\n"
            "Run scripts/particle_decay_analysis.py first to generate it."
        )

    df = pd.read_excel(summary_path, sheet_name=SUMMARY_SHEET, engine="openpyxl")

    if SHOWER_ON_COL not in df.columns:
        raise KeyError(
            f"'{SHOWER_ON_COL}' column missing from the {SUMMARY_SHEET} sheet. "
            "Regenerate it with the updated particle_decay_analysis.py."
        )

    df[SHOWER_ON_COL] = pd.to_datetime(df[SHOWER_ON_COL], errors="coerce")
    n_before = len(df)
    df = df.dropna(subset=[SHOWER_ON_COL]).reset_index(drop=True)
    n_dropped = n_before - len(df)
    if n_dropped:
        print(f"  Dropped {n_dropped} row(s) with no shower_on time.")
    return df


_UNIT_SUFFIX_RE = re.compile(r" \([^)]*\)$")


def _strip_unit_suffixes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Undo scripts.particle_decay_analysis._save_results's column_rename (e.g.
    "bin3_outside_beta_other (h-1)" -> "bin3_outside_beta_other"), so
    Excel-loaded columns match the plain names src.plot_particle's summary and
    boxplot functions expect. The per-bin variant comparison figures
    (make_figure) read the suffixed Excel names directly and don't call this.
    """
    return df.rename(columns={c: _UNIT_SUFFIX_RE.sub("", c) for c in df.columns})


# =============================================================================
# Paired Change
# =============================================================================


def compute_variant_change(df: pd.DataFrame, col_a: str, col_b: str) -> dict:
    """
    Compute paired (b - a) change across events for one metric/bin/variant pair.

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary rows.
    col_a, col_b : str
        Column names holding the two variants' values.

    Returns
    -------
    dict
        n (paired event count), mean_diff, std_diff, t_stat, p_value. t_stat
        and p_value are NaN when n < 2 (ttest_rel requires at least 2 pairs).
        NaN rows in either column (e.g. E3 on post-cutover events) are
        excluded before pairing.
    """
    col_a_vals = df[col_a]
    col_b_vals = df[col_b]
    valid = df[col_a_vals.notna() & col_b_vals.notna()]

    a = valid[col_a].to_numpy(dtype=float)
    b = valid[col_b].to_numpy(dtype=float)
    diff = b - a
    n = len(valid)

    result = {
        "n": n,
        "mean_diff": diff.mean() if n else np.nan,
        "std_diff": diff.std() if n else np.nan,
        "t_stat": np.nan,
        "p_value": np.nan,
    }
    if n >= 2:
        result["t_stat"], result["p_value"] = ttest_rel(b, a)
    return result


def _format_p_value(p_value: float) -> str:
    """Format a p-value for display, using '< 0.001' below that threshold."""
    if not np.isfinite(p_value):
        return "n/a"
    if p_value < 0.001:
        return "p < 0.001"
    return f"p = {sf.fmt_fig(p_value)}"


def make_summary_div(
    df: pd.DataFrame, col_a: str, label_a: str, col_b: str, label_b: str, unit: str
) -> Div:
    """
    Build a Div summarizing the paired (b - a) change for one metric/bin/pair.

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary rows.
    col_a, label_a, col_b, label_b : str
        Column names and display labels for the two variants being compared.
    unit : str
        Unit string appended to the reported values.

    Returns
    -------
    Div
        Bokeh Div with the paired mean change, std, n, and t-test result.
    """
    change = compute_variant_change(df, col_a, col_b)

    if change["n"] == 0:
        text = f"Paired {label_b} minus {label_a} change: no events valid in both variants."
    else:
        mean_str = sf.fmt_fig(change["mean_diff"])
        std_str = sf.fmt_fig(change["std_diff"]) if change["n"] > 1 else "n/a"
        p_str = _format_p_value(change["p_value"])
        t_str = sf.fmt_fig(change["t_stat"]) if np.isfinite(change["t_stat"]) else "n/a"
        text = (
            f"Paired {label_a}→{label_b} change (n = {change['n']} events): "
            f"mean Δ ({label_b} − {label_a}) = {mean_str} ± {std_str}{unit} "
            f"(paired t-test: t = {t_str}, {p_str})"
        )

    return Div(
        text=text,
        width=MODUAIR_FIGURE_WIDTH,
        styles={"font-size": MODUAIR_TEXT_PT},
    )


# =============================================================================
# Plotting: Per-Bin Variant Comparison
# =============================================================================


def make_figure(
    df: pd.DataFrame,
    bin_index: int,
    metric_key: str,
    variant_a: str,
    variant_b: str,
    pair_stem: str,
    output_path: Path,
) -> None:
    """
    Build and save one bin's variant-comparison figure for one metric (beta or E).

    Each event contributes one point per variant at the event shower_on time,
    with a whisker showing the metric's std. A vertical line at ROOM_CUTOVER
    (2026-06-03) marks where the primary concentration series switches from
    raw C_bed1(t) to the fleet-averaged C_room(t); it does not split the
    plotted points.

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary rows.
    bin_index : int
        Particle-size bin index (0-11).
    metric_key : str
        Key into METRICS ("beta" or "E").
    variant_a, variant_b : str
        Keys into VARIANT_LEGEND/VARIANT_COLORS for the two variants compared.
    pair_stem : str
        Filename stem for this comparison pair (e.g. "entry_outside").
    output_path : Path
        Destination HTML path.
    """
    spec = METRICS[metric_key]
    bin_name = PARTICLE_BINS[bin_index]["name"]

    output_file(
        str(output_path), title=f"{metric_key} {VARIANT_LEGEND[variant_a]} vs {VARIANT_LEGEND[variant_b]}, bin {bin_index}"
    )

    fig = figure(
        x_axis_type="datetime",
        x_axis_label="Event shower start (date and time)",
        y_axis_label=f"{spec['y_label']} ({spec['y_units']}), bin {bin_index} ({bin_name} µm)",
        tools="pan,box_zoom,wheel_zoom,reset,save",
    )

    hover = HoverTool(
        tooltips=[
            ("Variant", "@variant"),
            ("Shower on", "@time_str"),
            ("Value", f"@value{{{spec['value_fmt']}}}"),
            ("Std", f"@std{{{spec['value_fmt']}}}"),
            ("R²", "@r2{0.000}"),
        ]
    )
    fig.add_tools(hover)

    cols_by_variant = {}
    for variant in (variant_a, variant_b):
        value_col = spec["value_col"].format(bin=bin_index, variant=variant)
        std_col = spec["std_col"].format(bin=bin_index, variant=variant)
        r2_col = spec["r2_col"].format(bin=bin_index, variant=variant)
        cols_by_variant[variant] = value_col

        if value_col not in df.columns:
            print(f"    [WARN] Bin {bin_index} ({metric_key}, {variant}): column not in results.")
            continue

        sub = df[df[value_col].notna()].copy()
        if sub.empty:
            print(f"    [WARN] Bin {bin_index} ({metric_key}, {variant}): no valid values.")
            continue

        mean_vals = sub[value_col].astype(float)
        std_vals = sub[std_col].astype(float).fillna(0.0) if std_col in sub.columns else pd.Series(0.0, index=sub.index)
        r2_vals = sub[r2_col].astype(float) if r2_col in sub.columns else pd.Series(np.nan, index=sub.index)
        color = VARIANT_COLORS[variant]

        source_data = ColumnDataSource(
            data={
                "time": sub[SHOWER_ON_COL],
                "value": mean_vals,
                "std": std_vals,
                "upper": mean_vals + std_vals,
                "lower": mean_vals - std_vals,
                "r2": r2_vals,
                "variant": [VARIANT_LEGEND[variant]] * len(sub),
                "time_str": sub[SHOWER_ON_COL].dt.strftime("%Y-%m-%d %H:%M:%S"),
            }
        )
        fig.scatter(
            "time",
            "value",
            source=source_data,
            marker="circle",
            size=7,
            color=color,
            alpha=0.8,
            legend_label=VARIANT_LEGEND[variant],
        )
        whisker = Whisker(
            base="time",
            upper="upper",
            lower="lower",
            source=source_data,
            line_color=color,
            line_alpha=0.6,
        )
        whisker.upper_head.line_color = color
        whisker.lower_head.line_color = color
        fig.add_layout(whisker)

    # Shared MODULAIR style: 1600x800, 12pt, no title, click-to-hide legend.
    style_moduair_figure(fig, legend_title="Variant", legend_location=LEGEND_LOCATION)

    fig.xaxis.formatter = DatetimeTickFormatter(
        days="%Y-%m-%d", hours="%m-%d %H:%M", minutes="%H:%M"
    )

    # Vertical line at the room-correction cutover, with region labels either side.
    fig.add_layout(
        Span(
            location=ROOM_CUTOVER,
            dimension="height",
            line_color=CUTOVER_LABEL_COLOR,
            line_dash="dashed",
            line_width=1.5,
        )
    )
    label_kwargs = dict(
        y=fig.height - 40,
        y_units="screen",
        text_font_size=MODUAIR_TEXT_PT,
        text_color=CUTOVER_LABEL_COLOR,
    )
    fig.add_layout(
        Label(x=ROOM_CUTOVER, x_offset=-10, text=CUTOVER_LABELS[0], text_align="right", **label_kwargs)
    )
    fig.add_layout(
        Label(x=ROOM_CUTOVER, x_offset=10, text=CUTOVER_LABELS[1], text_align="left", **label_kwargs)
    )

    col_a = cols_by_variant.get(variant_a)
    col_b = cols_by_variant.get(variant_b)
    if col_a and col_b and col_a in df.columns and col_b in df.columns:
        summary_div = make_summary_div(
            df, col_a, VARIANT_LEGEND[variant_a], col_b, VARIANT_LEGEND[variant_b], spec["unit"]
        )
        save(bokeh_column(fig, summary_div))
    else:
        save(fig)


# =============================================================================
# Plotting: Summary/Boxplot/Comparison Figures
# =============================================================================


def generate_summary_plots(results_df: pd.DataFrame, output_dir: Path) -> None:
    """
    Generate the summary bar charts and boxplot families (outside/entry, E1/E2
    only). Ported from scripts/particle_decay_analysis.py's
    _generate_summary_plots; reads all_results already loaded by main() rather
    than an in-memory results_df from the pipeline script.

    Applies the FLOW_RATE_MIN/FLOW_RATE_MAX filter before any figure is drawn
    (Excel and the per-event pm_decay figures are unaffected).
    """
    # Undo the Excel unit-suffix rename (e.g. "bin3_outside_beta_other (h-1)"
    # -> "bin3_outside_beta_other") -- src.plot_particle expects plain names.
    results_df = _strip_unit_suffixes(results_df)

    if "flow_rate" in results_df.columns and results_df["flow_rate"].notna().any():
        plot_df = results_df[
            results_df["flow_rate"].between(FLOW_RATE_MIN, FLOW_RATE_MAX)
        ].copy()
    else:
        plot_df = results_df.copy()

    if plot_df.empty:
        print("\nSkipping summary plot generation - no results within flow rate range.")
        return

    print("\nGenerating summary/boxplot/comparison figures...")
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(exist_ok=True)

    try:
        from src.plot_particle import (
            plot_deposition_rate_boxplot,
            plot_deposition_summary,
            plot_door_comparison_boxplots,
            plot_emission_boxplot,
            plot_emission_etotal_by_metric_boxplot,
            plot_emission_etotal_by_showerhead_boxplot,
            plot_emission_rate_boxplot,
            plot_emission_summary,
            plot_fan_comparison_boxplots,
            plot_mannequin_comparison_boxplots,
            plot_penetration_factor_boxplot,
            plot_penetration_summary,
            plot_shower_head_comparison_boxplots,
            plot_spray_pattern_comparison_boxplots,
        )
    except ImportError:
        print("  Warning: plot_particle module not found. Skipping plots.")
        return

    results_df = plot_df

    # Load Bedroom_Conditions RH data once for n= / RH= boxplot annotations.
    # This is best-effort: if the summary file is unavailable the annotation
    # will show n= only (no RH line).  The Bedroom_Conditions sheet has columns
    # "shower_on" (datetime) and "rh_mean (%)" which are renamed to the
    # "datetime" / "RH_bedroom" convention expected by the boxplot helpers.
    rh_data = None
    try:
        import pandas as _pd

        from src.data_paths import get_common_file

        _summary_path = get_common_file("rh_temp_wind_summary")
        _bc = _pd.read_excel(_summary_path, sheet_name="Bedroom_Conditions")
        _bc = _bc.rename(columns={"shower_on": "datetime", "rh_mean (%)": "RH_bedroom"})
        rh_data = _bc[["datetime", "RH_bedroom"]].copy()
        rh_data["datetime"] = _pd.to_datetime(rh_data["datetime"])
        print("  Loaded Bedroom_Conditions RH for boxplot annotations.")
    except Exception as _rh_err:
        print(f"  Note: Could not load Bedroom_Conditions RH data (n= only annotations): {_rh_err}")

    # Bar-chart summary plots (no RH annotation needed). Penetration is
    # lambda-independent (one file); deposition and emission are
    # lambda-dependent and each produce an _outside/_entry pair.
    for plot_func, filename, has_source in [
        (plot_penetration_summary, "penetration_summary.png", False),
        (plot_deposition_summary, "deposition_summary.png", True),
        (plot_emission_summary, "emission_summary.png", True),
    ]:
        try:
            plot_func(results_df, PARTICLE_BINS, plot_dir / filename)
            if has_source:
                stem, suffix = filename.rsplit(".", 1)
                print(f"  Generated: {stem}_{{outside,entry}}.{suffix}")
            else:
                print(f"  Generated: {filename}")
        except Exception as e:
            print(f"  Error generating {filename}: {e}")

    # ── emission E_total vs. continuous metric ──────────────────────────────
    # Add per-event computed columns to a local copy so results_df is not mutated.
    _df7 = results_df.copy()
    _bin_nums = list(PARTICLE_BINS.keys())

    # Average beta across all bins per event, per lambda source (lambda-dependent).
    # Average p across all bins per event, shared (lambda-independent).
    for _source in _LAMBDA_SOURCES:
        _df7[f"avg_beta_{_source}"] = _df7[
            [
                f"bin{b}_{_source}_beta_other"
                for b in _bin_nums
                if f"bin{b}_{_source}_beta_other" in _df7.columns
            ]
        ].mean(axis=1)
    _df7["avg_p"] = _df7[
        [f"bin{b}_p_mean" for b in _bin_nums if f"bin{b}_p_mean" in _df7.columns]
    ].mean(axis=1)

    # Merge bedroom RH and temperature from Bedroom_Conditions sheet (best-effort)
    try:
        import pandas as _pd7

        from src.data_paths import get_common_file as _gcf7

        _bc7 = _pd7.read_excel(_gcf7("rh_temp_wind_summary"), sheet_name="Bedroom_Conditions")
        _bc7 = _bc7[["event_number", "rh_mean (%)", "temp_mean (degC)"]].rename(
            columns={"rh_mean (%)": "bedroom_rh", "temp_mean (degC)": "bedroom_temp"}
        )
        _df7 = _df7.merge(_bc7, on="event_number", how="left")
        print("  Merged Bedroom_Conditions for metric-axis figures.")
    except Exception as _e7:
        _df7["bedroom_rh"] = np.nan
        _df7["bedroom_temp"] = np.nan
        print(f"  Note: Bedroom_Conditions not merged for metric-axis figures: {_e7}")

    # ── Boxplot x-range configuration ────────────────────────────────────────
    # Fixed water-temperature axis: x_range=(xmin, xmax, xtick_step) in °C
    _fixed_axis_boxplots = [
        (plot_emission_boxplot, "emission_etotal_boxplot.png", (5, 55, 5)),
        (plot_deposition_rate_boxplot, "other_process_rate_boxplot.png", (5, 55, 5)),
        (plot_emission_rate_boxplot, "emission_rate_boxplot.png", (5, 55, 5)),
        (plot_penetration_factor_boxplot, "penetration_factor_boxplot.png", (5, 55, 5)),
    ]
    for plot_func, filename, x_range in _fixed_axis_boxplots:
        try:
            plot_func(
                results_df,
                PARTICLE_BINS,
                plot_dir / filename,
                rh_data=rh_data,
                x_range=x_range,
            )
            print(f"  Generated: {filename}")
        except Exception as e:
            print(f"  Error generating {filename}: {e}")

    # Continuous metric axis: x_range=(xmin, xmax, step) in metric units.
    # E_total (the y-axis) is lambda-dependent, so every entry is drawn once
    # per source. Source-independent x-axis metrics (bedroom_rh, bedroom_temp,
    # avg_p) reuse the same metric_col for both source calls; lambda-dependent
    # x-axis metrics (air change rate, avg beta) pair each source's own column
    # with that source's E_total rather than crossing sources.
    _metric_axes_shared = [
        # (metric_col, metric_label, filename, x_range=(xmin, xmax, step))
        (
            "bedroom_rh",
            "Bedroom RH (%)",
            "emission_etotal_by_bedroom_rh_boxplot.png",
            (23, 43, 1),
        ),
        (
            "bedroom_temp",
            "Bedroom Temperature (°C)",
            "emission_etotal_by_bedroom_temp_boxplot.png",
            (14.9, 18.2, 0.1),
        ),
        (
            "avg_p",
            "Avg. Penetration Factor p",
            "emission_etotal_by_p_boxplot.png",
            (0.48, 0.8, 0.02),
        ),
    ]
    for metric_col, metric_label, filename, x_range in _metric_axes_shared:
        for _source in _LAMBDA_SOURCES:
            try:
                plot_emission_etotal_by_metric_boxplot(
                    _df7,
                    PARTICLE_BINS,
                    plot_dir / filename,
                    metric_col=metric_col,
                    metric_label=metric_label,
                    source=_source,
                    rh_data=rh_data,
                    x_range=x_range,
                )
                print(f"  Generated: {filename} ({_source} source, bin0-2 and bin3-6)")
            except Exception as e:
                print(f"  Error generating {filename} ({_source} source): {e}")

    # Lambda-dependent x-axis metrics: each source's own column paired with
    # that source's E_total (e.g. lambda_outside x-axis with outside E_total).
    _metric_axes_per_source = [
        # (metric_col_template, metric_label, filename, x_range=(xmin, xmax, step))
        (
            "lambda_{source}",
            "Air Change Rate λ (h⁻¹)",
            "emission_etotal_by_acr_boxplot.png",
            (0.75, 1.65, 0.05),
        ),
        (
            "avg_beta_{source}",
            "Avg. Other Process Rate β (h⁻¹)",
            "emission_etotal_by_beta_boxplot.png",
            (-0.35, 0.35, 0.05),
        ),
    ]
    for metric_col_template, metric_label, filename, x_range in _metric_axes_per_source:
        for _source in _LAMBDA_SOURCES:
            try:
                plot_emission_etotal_by_metric_boxplot(
                    _df7,
                    PARTICLE_BINS,
                    plot_dir / filename,
                    metric_col=metric_col_template.format(source=_source),
                    metric_label=metric_label,
                    source=_source,
                    rh_data=rh_data,
                    x_range=x_range,
                )
                print(f"  Generated: {filename} ({_source} source, bin0-2 and bin3-6)")
            except Exception as e:
                print(f"  Error generating {filename} ({_source} source): {e}")

    # Shower head type comparison: W53 (base) vs. W52pw (Pepco)
    for _source in _LAMBDA_SOURCES:
        try:
            _sh_filename = "emission_etotal_by_showerhead_boxplot.png"
            plot_emission_etotal_by_showerhead_boxplot(
                results_df,
                PARTICLE_BINS,
                plot_dir / _sh_filename,
                source=_source,
                rh_data=rh_data,
            )
            print(f"  Generated: {_sh_filename} ({_source} source, bin0-2 and bin3-6)")
        except Exception as e:
            print(f"  Error generating emission_etotal_by_showerhead_boxplot ({_source}): {e}")

    # ── Condition-comparison boxplots (5 families x 3 metrics x 3 bin groups,
    #    x2 sources for the two lambda-dependent metrics = 15 files/family) ──
    _comparison_families = [
        (
            plot_spray_pattern_comparison_boxplots,
            "spray_pattern",
            "spray pattern",
        ),
        (
            plot_shower_head_comparison_boxplots,
            "head_type",
            "shower head type",
        ),
        (
            plot_mannequin_comparison_boxplots,
            "mannequin",
            "mannequin presence",
        ),
        (
            plot_door_comparison_boxplots,
            "door_position",
            "door position",
        ),
        (
            plot_fan_comparison_boxplots,
            "fan_status",
            "fan status",
        ),
    ]
    for _cmp_func, _cmp_stem, _cmp_label in _comparison_families:
        try:
            _cmp_func(results_df, PARTICLE_BINS, plot_dir, rh_data=rh_data)
            print(
                f"  Generated: {_cmp_stem}_{{emission_etotal,other_process_rate}}"
                f"_boxplot_{{bin0-2,bin3-6,bin7-11}}_{{outside,entry}}.png, "
                f"{_cmp_stem}_penetration_factor_boxplot_{{bin0-2,bin3-6,bin7-11}}.png"
                f"  ({_cmp_label} comparison)"
            )
        except Exception as e:
            print(f"  Error generating {_cmp_stem} comparison figures: {e}")

    print(f"  Plots saved to: {plot_dir}")


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Particle emission-rate variant figures: per-bin entry/outside/"
        "adjusted/bathroom comparisons plus the summary and boxplot figure suite."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help=(
            "Directory holding particle_analysis_summary.xlsx and where figures are "
            "written under plots/ and plots/particle/ (default: data_root/output, "
            "matching particle_decay_analysis.py)"
        ),
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else get_data_root() / "output"

    print("\n" + "=" * 70)
    print("Particle Emission Rate Variant Figures")
    print("=" * 70)

    summary_path = output_dir / SUMMARY_XLSX_NAME
    print(f"Reading: {summary_path}")
    df = load_summary(summary_path)
    print(f"  {len(df)} event row(s) with a shower_on time")

    figure_dir = output_dir.joinpath(*FIGURE_SUBDIR)
    figure_dir.mkdir(parents=True, exist_ok=True)

    print("\nBuilding per-bin variant comparison figures...")
    for bin_index in PARTICLE_BINS:
        for metric_key in METRICS:
            for variant_a, variant_b, pair_stem in COMPARISON_PAIRS:
                spec = METRICS[metric_key]
                output_path = figure_dir / spec["filename"].format(pair=pair_stem, bin=bin_index)
                make_figure(df, bin_index, metric_key, variant_a, variant_b, pair_stem, output_path)
                print(f"  Saved {output_path.name}")

    generate_summary_plots(df, output_dir)

    print("\n" + "=" * 70)
    print("Done")
    print("=" * 70)


if __name__ == "__main__":
    main()
