#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CO2 Air-Change Rate by Source Concentration Figure
===================================================

Builds one interactive Bokeh figure of the per-event air-change rate (lambda)
computed from two source concentrations: the entry zone and the outdoor
reference. Each event contributes two points at its shower_on time, one lambda
from the entry source and one from the outside source, matching the pre/post
per-event layout in scripts/hobo_onset_decay_figures.py. Color encodes source.
Whiskers show the regression standard error (lambda_*_std). A summary Div below
the figure reports the paired outside minus entry change across events with a
valid fit in both sources, with a paired t-test (scipy.stats.ttest_rel).

This script reads co2_lambda_summary.csv produced by
scripts/co2_decay_analysis.py. It does not recompute lambda. Lambda values that
failed the R2 gate in that script are stored as NaN and are skipped here, so a
plotted point means its fit passed the R2 threshold. Run co2_decay_analysis.py
first, then this script.

Style matches the shared MODULAIR Bokeh figures (1600x800, 12pt, no title,
click-to-hide legend, hover enabled) via src.plot_style.style_moduair_figure.

Output Files:
    output/plots/co2/co2_lambda_entry_outside.html

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-15
Update log:
    2026-09-15  Initial version.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from bokeh.layouts import column as bokeh_column
from bokeh.models import ColumnDataSource, DatetimeTickFormatter, Div, HoverTool, Whisker
from bokeh.plotting import figure, output_file, save
from scipy.stats import ttest_rel

# Add project root to path for src/ imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import src.sig_figs as sf  # noqa: E402
from src.data_paths import get_data_root  # noqa: E402
from src.plot_style import SENSOR_COLORS, style_moduair_figure  # noqa: E402
from src.plot_style import MODUAIR_FIGURE_WIDTH, MODUAIR_TEXT_PT  # noqa: E402

# =============================================================================
# Configuration
# =============================================================================

# Input summary produced by scripts/co2_decay_analysis.py.
SUMMARY_CSV_NAME = "co2_lambda_summary.csv"

# Output figure, placed alongside the other CO2 plots under output/plots/.
FIGURE_SUBDIR = ("plots", "co2")
FIGURE_NAME = "co2_lambda_entry_outside.html"

# Value column formatting and axis label for lambda (h^-1).
VALUE_FMT = "0.000"
Y_AXIS_LABEL = "Air-change rate \u03bb (h\u207b\u00b9)"
LEGEND_LOCATION = "top_right"

# Source concentration methods to plot. Column names match the renamed columns
# written by co2_decay_analysis.py (units appended on lambda mean/std; the
# r_squared columns are written without a unit suffix). Color encodes source.
SOURCES = {
    "entry": {
        "mean_col": "lambda_entry_mean (h-1)",
        "std_col": "lambda_entry_std (h-1)",
        "r2_col": "lambda_entry_r_squared",
        "legend": "Entry",
        "color": SENSOR_COLORS[0],
    },
    "outside": {
        "mean_col": "lambda_outside_mean (h-1)",
        "std_col": "lambda_outside_std (h-1)",
        "r2_col": "lambda_outside_r_squared",
        "legend": "Outside",
        "color": SENSOR_COLORS[1],
    },
}

# Column holding the event shower_on time, used as the x position for both points.
SHOWER_ON_COL = "shower_on"


# =============================================================================
# Data Loading
# =============================================================================


def load_summary(summary_path: Path) -> pd.DataFrame:
    """
    Load the CO2 lambda summary CSV and parse the shower_on column.

    Parameters
    ----------
    summary_path : Path
        Path to co2_lambda_summary.csv from co2_decay_analysis.py.

    Returns
    -------
    pd.DataFrame
        Summary rows with shower_on coerced to datetime. Rows without a valid
        shower_on time are dropped, since shower_on is the figure x position.
    """
    if not summary_path.exists():
        raise FileNotFoundError(
            f"Summary CSV not found: {summary_path}\n"
            "Run scripts/co2_decay_analysis.py first to generate it."
        )

    # co2_decay_analysis.py writes the per-event summary with a UTF-8 BOM.
    df = pd.read_csv(summary_path, encoding="utf-8-sig")

    required = {SHOWER_ON_COL}
    for spec in SOURCES.values():
        required.update({spec["mean_col"], spec["std_col"]})
    missing = required - set(df.columns)
    if missing:
        raise KeyError(
            f"Summary CSV is missing expected columns: {sorted(missing)}. "
            "Regenerate it with the updated co2_decay_analysis.py."
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


def compute_source_change(df: pd.DataFrame) -> dict:
    """
    Compute paired outside minus entry change in lambda across events.

    Events are paired within a row: each event carries both an entry lambda and
    an outside lambda. Only events with a valid (non-NaN) lambda for both sources
    are included, which means both fits passed the R2 gate upstream. A paired
    t-test checks whether the mean outside minus entry difference differs from
    zero.

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary rows.

    Returns
    -------
    dict
        n (paired event count), mean_diff, std_diff, t_stat, p_value. t_stat and
        p_value are NaN when n < 2 (ttest_rel requires at least 2 pairs).
    """
    entry = df[SOURCES["entry"]["mean_col"]]
    outside = df[SOURCES["outside"]["mean_col"]]
    valid = df[entry.notna() & outside.notna()]

    e = valid[SOURCES["entry"]["mean_col"]].to_numpy(dtype=float)
    o = valid[SOURCES["outside"]["mean_col"]].to_numpy(dtype=float)
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


def make_summary_div(df: pd.DataFrame) -> Div:
    """
    Build a Div summarizing the paired outside minus entry lambda change.

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary rows.

    Returns
    -------
    Div
        Bokeh Div with the paired mean change, std, n, and t-test result.
    """
    change = compute_source_change(df)
    unit = " h\u207b\u00b9"

    if change["n"] == 0:
        text = "Paired outside minus entry change: no events valid in both sources."
    else:
        mean_str = sf.fmt_fig(change["mean_diff"])
        std_str = sf.fmt_fig(change["std_diff"]) if change["n"] > 1 else "n/a"
        p_str = _format_p_value(change["p_value"])
        t_str = sf.fmt_fig(change["t_stat"]) if np.isfinite(change["t_stat"]) else "n/a"
        text = (
            f"Paired outside\u2192entry change (n = {change['n']} events): "
            f"mean \u0394 (outside \u2212 entry) = {mean_str} \u00b1 {std_str}{unit} "
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


def make_figure(df: pd.DataFrame, output_path: Path) -> None:
    """
    Build and save the entry vs outside lambda figure.

    Each event contributes one entry point and one outside point at the event
    shower_on time, with a whisker showing lambda plus or minus its regression
    standard error. R2 and the fit standard error are in the hover. Color encodes
    source with a two-entry legend (Entry, Outside).

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary rows.
    output_path : Path
        Destination HTML path.
    """
    output_file(str(output_path), title="CO2 air-change rate by source")

    fig = figure(
        x_axis_type="datetime",
        x_axis_label="Event shower start (date and time)",
        y_axis_label=Y_AXIS_LABEL,
        tools="pan,box_zoom,wheel_zoom,reset,save",
    )

    hover = HoverTool(
        tooltips=[
            ("Source", "@source"),
            ("Shower on", "@time_str"),
            ("\u03bb", f"@value{{{VALUE_FMT}}}"),
            ("Std err", f"@std{{{VALUE_FMT}}}"),
            ("R\u00b2", "@r2{0.000}"),
        ]
    )
    fig.add_tools(hover)

    for spec in SOURCES.values():
        mean_col = spec["mean_col"]
        std_col = spec["std_col"]
        r2_col = spec["r2_col"]

        sub = df[df[mean_col].notna()].copy()
        if sub.empty:
            print(f"  [WARN] No valid {spec['legend']} lambda values to plot.")
            continue

        mean_vals = sub[mean_col].astype(float)
        std_vals = sub[std_col].astype(float).fillna(0.0)
        r2_vals = (
            sub[r2_col].astype(float)
            if r2_col in sub.columns
            else pd.Series(np.nan, index=sub.index)
        )
        color = spec["color"]

        source = ColumnDataSource(
            data={
                "time": sub[SHOWER_ON_COL],
                "value": mean_vals,
                "std": std_vals,
                "upper": mean_vals + std_vals,
                "lower": mean_vals - std_vals,
                "r2": r2_vals,
                "source": [spec["legend"]] * len(sub),
                "time_str": sub[SHOWER_ON_COL].dt.strftime("%Y-%m-%d %H:%M:%S"),
            }
        )
        fig.scatter(
            "time",
            "value",
            source=source,
            marker="circle",
            size=7,
            color=color,
            alpha=0.8,
            legend_label=spec["legend"],
        )
        whisker = Whisker(
            base="time",
            upper="upper",
            lower="lower",
            source=source,
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

    summary_div = make_summary_div(df)
    save(bokeh_column(fig, summary_div))
    print(f"  Saved {output_path}")


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CO2 air-change rate by source concentration (entry vs outside)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help=(
            "Directory holding co2_lambda_summary.csv and where the figure is "
            "written under plots/co2/ (default: data_root/output, matching "
            "co2_decay_analysis.py)"
        ),
    )
    args = parser.parse_args()

    output_dir = (
        Path(args.output_dir) if args.output_dir else get_data_root() / "output"
    )

    print("\n" + "=" * 70)
    print("CO2 Air-Change Rate by Source Concentration Figure")
    print("=" * 70)

    summary_path = output_dir / SUMMARY_CSV_NAME
    print(f"Reading: {summary_path}")
    df = load_summary(summary_path)
    print(f"  {len(df)} event row(s) with a shower_on time")

    n_entry = int(df[SOURCES["entry"]["mean_col"]].notna().sum())
    n_outside = int(df[SOURCES["outside"]["mean_col"]].notna().sum())
    print(f"  Valid Entry lambda: {n_entry}   Valid Outside lambda: {n_outside}")
    if n_entry == 0 and n_outside == 0:
        print("ERROR: no valid entry or outside lambda values to plot.")
        sys.exit(1)

    figure_dir = output_dir.joinpath(*FIGURE_SUBDIR)
    figure_dir.mkdir(parents=True, exist_ok=True)
    output_path = figure_dir / FIGURE_NAME

    print("\nBuilding figure...")
    make_figure(df, output_path)

    print("\n" + "=" * 70)
    print("Done")
    print("=" * 70)


if __name__ == "__main__":
    main()
