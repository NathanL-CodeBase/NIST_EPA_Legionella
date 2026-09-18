#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle Beta/Emission by Source Concentration Figures
========================================================

Builds two interactive Bokeh figures per particle-size bin (24 total): the
other-process rate (beta_other) and the mean emission rate (E_mean), each
plotted per event over time with an entry-source point and an outside-source
point (the dual air-change-rate bounding from
scripts/particle_decay_analysis.py), matching the layout of
scripts/co2_lambda_source_figure.py. A vertical line at 2026-06-03 marks the
room-concentration correction cutover (src/particle_room_correction.py):
before that date the indoor concentration feeding beta/E is a per-event
ratio-corrected MOD-PM-00195 reading ("C_adjusted room"); from that date it is
the position-weighted MODULAIR-PM fleet average ("C_room"). The line is
informational only -- it does not split the plotted points into groups the
way the entry/outside color does.

This script reads particle_analysis_summary.xlsx (all_results sheet) produced
by scripts/particle_decay_analysis.py. It does not recompute beta or E. Run
particle_decay_analysis.py first.

Style matches the shared MODULAIR Bokeh figures (1600x800, 12pt, no title,
click-to-hide legend, hover enabled) via src.plot_style.style_moduair_figure.

Output Files:
    output/plots/particle/beta_other_entry_outside_bin{N}.html   (12 figures)
    output/plots/particle/emission_entry_outside_bin{N}.html     (12 figures)

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-18
Update log:
    2026-09-18  Initial version.
"""

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
# Configuration
# =============================================================================

SUMMARY_XLSX_NAME = "particle_analysis_summary.xlsx"
SUMMARY_SHEET = "all_results"

FIGURE_SUBDIR = ("plots", "particle")

SHOWER_ON_COL = "shower_on"

_LAMBDA_SOURCES = ("entry", "outside")

# One metric spec per figure kind: value/std/r2 column templates ({bin} and
# {source} filled in per figure/source), axis label, hover value format, unit
# string for the summary Div, and output filename template.
METRICS = {
    "beta": {
        "value_col": "bin{bin}_{source}_beta_other (h-1)",
        "std_col": "bin{bin}_{source}_beta_other_std (h-1)",
        "r2_col": "bin{bin}_{source}_beta_other_r_squared",
        "y_label": "Other process rate β",
        "y_units": "h⁻¹",
        "value_fmt": "0.000",
        "unit": " h⁻¹",
        "filename": "beta_other_entry_outside_bin{bin}.html",
    },
    "E": {
        "value_col": "bin{bin}_{source}_E_mean (#/min)",
        "std_col": "bin{bin}_{source}_E_std (#/min)",
        "r2_col": "bin{bin}_{source}_E_r_squared",
        "y_label": "Emission rate E",
        "y_units": "#/min",
        "value_fmt": "0.00e+00",
        "unit": " #/min",
        "filename": "emission_entry_outside_bin{bin}.html",
    },
}

SOURCE_LEGEND = {"entry": "Entry", "outside": "Outside"}
SOURCE_COLORS = {"entry": SENSOR_COLORS[0], "outside": SENSOR_COLORS[1]}

LEGEND_LOCATION = "top_right"

# Cutover-region labels either side of the vertical Span (left, right).
CUTOVER_LABELS = ("C_adjusted room", "C_room")
CUTOVER_LABEL_COLOR = "gray"


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


# =============================================================================
# Paired Change
# =============================================================================


def compute_source_change(df: pd.DataFrame, entry_col: str, outside_col: str) -> dict:
    """
    Compute paired outside minus entry change across events for one metric/bin.

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary rows.
    entry_col, outside_col : str
        Column names holding the entry and outside values.

    Returns
    -------
    dict
        n (paired event count), mean_diff, std_diff, t_stat, p_value. t_stat
        and p_value are NaN when n < 2 (ttest_rel requires at least 2 pairs).
    """
    entry = df[entry_col]
    outside = df[outside_col]
    valid = df[entry.notna() & outside.notna()]

    e = valid[entry_col].to_numpy(dtype=float)
    o = valid[outside_col].to_numpy(dtype=float)
    diff = o - e
    n = len(valid)

    result = {
        "n": n,
        "mean_diff": diff.mean() if n else np.nan,
        "std_diff": diff.std() if n else np.nan,
        "t_stat": np.nan,
        "p_value": np.nan,
    }
    if n >= 2:
        result["t_stat"], result["p_value"] = ttest_rel(o, e)
    return result


def _format_p_value(p_value: float) -> str:
    """Format a p-value for display, using '< 0.001' below that threshold."""
    if not np.isfinite(p_value):
        return "n/a"
    if p_value < 0.001:
        return "p < 0.001"
    return f"p = {sf.fmt_fig(p_value)}"


def make_summary_div(df: pd.DataFrame, entry_col: str, outside_col: str, unit: str) -> Div:
    """
    Build a Div summarizing the paired outside minus entry change for one metric/bin.

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary rows.
    entry_col, outside_col : str
        Column names holding the entry and outside values.
    unit : str
        Unit string appended to the reported values.

    Returns
    -------
    Div
        Bokeh Div with the paired mean change, std, n, and t-test result.
    """
    change = compute_source_change(df, entry_col, outside_col)

    if change["n"] == 0:
        text = "Paired outside minus entry change: no events valid in both sources."
    else:
        mean_str = sf.fmt_fig(change["mean_diff"])
        std_str = sf.fmt_fig(change["std_diff"]) if change["n"] > 1 else "n/a"
        p_str = _format_p_value(change["p_value"])
        t_str = sf.fmt_fig(change["t_stat"]) if np.isfinite(change["t_stat"]) else "n/a"
        text = (
            f"Paired outside→entry change (n = {change['n']} events): "
            f"mean Δ (outside − entry) = {mean_str} ± {std_str}{unit} "
            f"(paired t-test: t = {t_str}, {p_str})"
        )

    return Div(
        text=text,
        width=MODUAIR_FIGURE_WIDTH,
        styles={"font-size": MODUAIR_TEXT_PT},
    )


# =============================================================================
# Plotting
# =============================================================================


def make_figure(df: pd.DataFrame, bin_index: int, metric_key: str, output_path: Path) -> None:
    """
    Build and save one bin's entry vs outside figure for one metric (beta or E).

    Each event contributes one entry point and one outside point at the event
    shower_on time, with a whisker showing the metric's std. A vertical line at
    ROOM_CUTOVER (2026-06-03) marks where the indoor concentration feeding
    beta/E switches from the per-event ratio-corrected MOD-PM-00195 reading to
    the fleet-averaged C_room; it does not split the plotted points.

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary rows.
    bin_index : int
        Particle-size bin index (0-11).
    metric_key : str
        Key into METRICS ("beta" or "E").
    output_path : Path
        Destination HTML path.
    """
    spec = METRICS[metric_key]
    bin_name = PARTICLE_BINS[bin_index]["name"]

    output_file(str(output_path), title=f"{metric_key} by source, bin {bin_index}")

    fig = figure(
        x_axis_type="datetime",
        x_axis_label="Event shower start (date and time)",
        y_axis_label=f"{spec['y_label']} ({spec['y_units']}), bin {bin_index} ({bin_name} µm)",
        tools="pan,box_zoom,wheel_zoom,reset,save",
    )

    hover = HoverTool(
        tooltips=[
            ("Source", "@source"),
            ("Shower on", "@time_str"),
            ("Value", f"@value{{{spec['value_fmt']}}}"),
            ("Std", f"@std{{{spec['value_fmt']}}}"),
            ("R²", "@r2{0.000}"),
        ]
    )
    fig.add_tools(hover)

    cols_by_source = {}
    for source in _LAMBDA_SOURCES:
        value_col = spec["value_col"].format(bin=bin_index, source=source)
        std_col = spec["std_col"].format(bin=bin_index, source=source)
        r2_col = spec["r2_col"].format(bin=bin_index, source=source)
        cols_by_source[source] = (value_col, std_col)

        sub = df[df[value_col].notna()].copy()
        if sub.empty:
            print(f"    [WARN] Bin {bin_index} ({metric_key}, {source}): no valid values.")
            continue

        mean_vals = sub[value_col].astype(float)
        std_vals = sub[std_col].astype(float).fillna(0.0)
        r2_vals = sub[r2_col].astype(float) if r2_col in sub.columns else pd.Series(np.nan, index=sub.index)
        color = SOURCE_COLORS[source]

        source_data = ColumnDataSource(
            data={
                "time": sub[SHOWER_ON_COL],
                "value": mean_vals,
                "std": std_vals,
                "upper": mean_vals + std_vals,
                "lower": mean_vals - std_vals,
                "r2": r2_vals,
                "source": [SOURCE_LEGEND[source]] * len(sub),
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
            legend_label=SOURCE_LEGEND[source],
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
    style_moduair_figure(fig, legend_title="Source", legend_location=LEGEND_LOCATION)

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

    entry_col, _ = cols_by_source["entry"]
    outside_col, _ = cols_by_source["outside"]
    summary_div = make_summary_div(df, entry_col, outside_col, spec["unit"])
    save(bokeh_column(fig, summary_div))


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Particle beta/emission by source concentration figures (entry vs outside), "
        "one pair of figures per particle-size bin, with the room-correction cutover marked."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help=(
            "Directory holding particle_analysis_summary.xlsx and where figures are "
            "written under plots/particle/ (default: data_root/output, matching "
            "particle_decay_analysis.py)"
        ),
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else get_data_root() / "output"

    print("\n" + "=" * 70)
    print("Particle Beta/Emission by Source Concentration Figures")
    print("=" * 70)

    summary_path = output_dir / SUMMARY_XLSX_NAME
    print(f"Reading: {summary_path}")
    df = load_summary(summary_path)
    print(f"  {len(df)} event row(s) with a shower_on time")

    figure_dir = output_dir.joinpath(*FIGURE_SUBDIR)
    figure_dir.mkdir(parents=True, exist_ok=True)

    print("\nBuilding figures...")
    for bin_index in PARTICLE_BINS:
        for metric_key, spec in METRICS.items():
            output_path = figure_dir / spec["filename"].format(bin=bin_index)
            make_figure(df, bin_index, metric_key, output_path)
            print(f"  Saved {output_path.name}")

    print("\n" + "=" * 70)
    print("Done")
    print("=" * 70)


if __name__ == "__main__":
    main()
