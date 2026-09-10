#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HOBO UX100 Onset/Decay Time Series Figures
==========================================

Builds two interactive Bokeh time series figures for HOBO UX100 temperature
and relative humidity in the EPA Legionella project. Data are restricted to
onset and decay windows around each shower event, defined exactly as in
scripts/moduair_cave_ratio.py:

    onset: shower_on           to shower_on  + 60 min   (inclusive)
    decay: shower_off + 60 min to shower_off + 120 min  (inclusive)

Only HOBO measurements inside an onset or decay window of any event are kept.
Onset points are drawn as circles, decay points as squares, one color per
sensor. Figures match the shared MODULAIR Bokeh style (1600x800, 12pt, no
title, click-to-hide legend, hover enabled) via src.plot_style.style_moduair_figure.

Experimental period: 2026-01-15 00:00:00 through 2026-07-16 23:59:59 (events
ended 2026-07-16).

Sensor display labels follow the report section, which intentionally differs
from the loader room config for MB_E:

    MB_C    -> Bedroom East
    MB_Bed  -> Bedroom Central
    MB_E    -> Bedroom West
    MB_F    -> Bedroom North
    MB_Bath -> Bathroom

Output Files:
    output/plots/hobo/hobo_temp_onset_decay.html
    output/plots/hobo/hobo_rh_onset_decay.html

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-09
Update log:
    2026-09-09  Initial version.
    2026-09-10  Add pre/post per-event figures. Each event shows one pre point
                and one post point: the pooled mean across all five sensors' raw
                points in the 30-min pre and 2-hr post windows (matching
                rh_temp_other_analysis.py), with pooled std whiskers. Color
                encodes window (pre vs post); markers at window midpoints.
    2026-09-10  Add a summary Div below each pre/post figure with the paired
                pre->post change: mean difference +/- std across events with
                data in both windows, and a paired t-test (scipy.stats.ttest_rel).
"""

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
from src.data_paths import get_common_file  # noqa: E402
from src.env_data_loader import (  # noqa: E402
    identify_shower_events,
    load_hobo_data,
    load_shower_log,
)
from src.plot_style import COLORS, SENSOR_COLORS, style_moduair_figure  # noqa: E402
from src.plot_style import MODUAIR_FIGURE_WIDTH, MODUAIR_TEXT_PT  # noqa: E402

# =============================================================================
# Configuration
# =============================================================================

# Experimental period. Events ended 2026-07-16; end_date bounds the last decay
# window (shower_off + 120 min) well within end of day.
PERIOD_START = pd.Timestamp("2026-01-15 00:00:00")
PERIOD_END = pd.Timestamp("2026-07-16 23:59:59")

# Sensor location keys (env_data_loader) mapped to report display labels.
# NOTE: MB_E is labeled "Bedroom West" per the report section, which differs
# from the loader room config (doorway); this is intentional for that figure.
SENSOR_LABELS = {
    "MB_C": "Bedroom East",
    "MB_Bed": "Bedroom Central",
    "MB_E": "Bedroom West",
    "MB_F": "Bedroom North",
    "MB_Bath": "Bathroom",
}

# One color per sensor, drawn from the shared project palette in draw order.
SENSOR_ORDER = list(SENSOR_LABELS.keys())
SENSOR_COLOR = {sn: SENSOR_COLORS[i] for i, sn in enumerate(SENSOR_ORDER)}

# Onset/decay windows relative to shower_on/shower_off, inclusive on both ends
# (matches SCATTER_WINDOWS in scripts/moduair_cave_ratio.py).
WINDOWS = {
    "onset": {
        "anchor": "shower_on",
        "offset_start": pd.Timedelta(minutes=0),
        "offset_end": pd.Timedelta(minutes=60),
        "marker": "circle",
    },
    "decay": {
        "anchor": "shower_off",
        "offset_start": pd.Timedelta(minutes=60),
        "offset_end": pd.Timedelta(minutes=120),
        "marker": "square",
    },
}

# Pre/post windows for the per-event pooled-mean figures. These match the windows
# used in scripts/rh_temp_other_analysis.py (PRE_SHOWER_MINUTES = 30 before
# shower_on, POST_SHOWER_HOURS = 2 after shower_off). For each event, the raw
# points from all five sensors in a window are pooled; the marker is the mean of
# that pool and the whisker is its std. One pre point and one post point per
# event, placed at the window midpoint. Color encodes window (pre vs post).
PRE_SHOWER_MINUTES = 30
POST_SHOWER_HOURS = 2
PRE_POST_WINDOWS = {
    "pre": {
        "start_anchor": "shower_on",
        "start_offset": pd.Timedelta(minutes=-PRE_SHOWER_MINUTES),
        "end_anchor": "shower_on",
        "end_offset": pd.Timedelta(0),
        "color": COLORS["pre_shower"],
        "legend": "Pre-shower",
    },
    "post": {
        "start_anchor": "shower_off",
        "start_offset": pd.Timedelta(0),
        "end_anchor": "shower_off",
        "end_offset": pd.Timedelta(hours=POST_SHOWER_HOURS),
        "color": COLORS["post_shower"],
        "legend": "Post-shower",
    },
}

# Output directory: data_root/output/plots/hobo/ from data_config.json
FIGURE_DIR = get_common_file("output_folder") / "plots" / "hobo"

# Per-figure y-axis label, value column in the loader output, and legend corner.
FIGURES = {
    "temp": {
        "column": "temp_c",
        "y_axis_label": "Temperature (\u00b0C)",
        "legend_location": "top_right",
        "value_fmt": "0.00",
        "filename": "hobo_temp_onset_decay.html",
        "title_tag": "HOBO temperature onset/decay",
    },
    "rh": {
        "column": "rh_pct",
        "y_axis_label": "Relative humidity (%)",
        "legend_location": "top_left",
        "value_fmt": "0.0",
        "filename": "hobo_rh_onset_decay.html",
        "title_tag": "HOBO relative humidity onset/decay",
    },
}

# Per-figure config for the pre/post per-event mean figures.
PRE_POST_FIGURES = {
    "temp": {
        "column": "temp_c",
        "y_axis_label": "Temperature (\u00b0C)",
        "legend_location": "top_right",
        "value_fmt": "0.00",
        "filename": "hobo_temp_pre_post.html",
        "title_tag": "HOBO temperature pre/post event means",
    },
    "rh": {
        "column": "rh_pct",
        "y_axis_label": "Relative humidity (%)",
        "legend_location": "top_left",
        "value_fmt": "0.0",
        "filename": "hobo_rh_pre_post.html",
        "title_tag": "HOBO relative humidity pre/post event means",
    },
}


# =============================================================================
# Event Windows
# =============================================================================


def load_events_in_period(start: pd.Timestamp, end: pd.Timestamp) -> list:
    """
    Load shower events whose shower_on falls within the experimental period.

    Parameters
    ----------
    start, end : pd.Timestamp
        Inclusive period bounds applied to each event's shower_on time.

    Returns
    -------
    list
        Event dicts (from identify_shower_events) with shower_on/shower_off
        inside the period.
    """
    shower_log = load_shower_log()
    events = identify_shower_events(shower_log)
    return [
        ev
        for ev in events
        if pd.notna(ev["shower_on"])
        and start <= ev["shower_on"] <= end
        and pd.notna(ev["shower_off"])
    ]


def window_mask(index: pd.DatetimeIndex, events: list, window: str) -> np.ndarray:
    """
    Boolean mask selecting timestamps inside any event's onset or decay window.

    Windows are inclusive on both ends and unioned across events, matching
    build_scatter_window_mask in scripts/moduair_cave_ratio.py.

    Parameters
    ----------
    index : pd.DatetimeIndex
        Timestamps to test.
    events : list
        Event dicts with shower_on and shower_off keys.
    window : str
        Key into WINDOWS ("onset" or "decay").

    Returns
    -------
    np.ndarray
        Boolean mask aligned to index.
    """
    spec = WINDOWS[window]
    mask = np.zeros(len(index), dtype=bool)
    for ev in events:
        anchor = ev[spec["anchor"]]
        start = anchor + spec["offset_start"]
        end = anchor + spec["offset_end"]
        mask |= (index >= start) & (index <= end)
    return mask


# =============================================================================
# Data Assembly
# =============================================================================


def build_windowed_frames(events: list) -> dict:
    """
    Load each sensor and split its readings into onset and decay windowed frames.

    Parameters
    ----------
    events : list
        Shower event dicts within the experimental period.

    Returns
    -------
    dict
        Nested mapping {sensor_key: {window: DataFrame}} where each DataFrame
        has datetime, temp_c, and rh_pct columns for points inside that window.
        Sensors or windows with no surviving points map to an empty DataFrame.
    """
    frames = {}
    for sn in SENSOR_ORDER:
        df = load_hobo_data(sn, PERIOD_START, PERIOD_END)
        frames[sn] = {}
        if df.empty:
            print(f"  [WARN] No HOBO data for {sn} in period; skipping sensor.")
            for window in WINDOWS:
                frames[sn][window] = pd.DataFrame()
            continue

        idx = pd.DatetimeIndex(df["datetime"])
        for window in WINDOWS:
            mask = window_mask(idx, events, window)
            frames[sn][window] = df.loc[mask].reset_index(drop=True)
            n = len(frames[sn][window])
            print(f"  {sn} {window}: {n} points")
    return frames


def build_pre_post_stats(events: list) -> dict:
    """
    Compute per-event pooled pre and post window stats across all five sensors.

    For each event and window, the raw points from all five sensors that fall in
    the window are pooled into a single set; the marker value is the mean of that
    pool and the whisker is its std. The pre window is the 30 minutes before
    shower_on and the post window is the 2 hours after shower_off, matching
    PRE_SHOWER_MINUTES and POST_SHOWER_HOURS in scripts/rh_temp_other_analysis.py.
    Windows are inclusive on both ends. The marker is placed at the window
    midpoint. Events with no points from any sensor in a window are counted and
    skipped for that window.

    Parameters
    ----------
    events : list
        Shower event dicts with shower_on and shower_off.

    Returns
    -------
    dict
        Mapping {window: DataFrame} where each DataFrame has midpoint, shower_on,
        temp_c_mean, temp_c_std, rh_pct_mean, rh_pct_std, n_points (pooled point
        count), and n_sensors (sensors contributing) columns, one row per event
        with pooled data in that window.
    """
    # Load each sensor once; reuse across all events and both windows.
    sensor_frames = {}
    for sn in SENSOR_ORDER:
        df = load_hobo_data(sn, PERIOD_START, PERIOD_END)
        if df.empty:
            print(f"  [WARN] No HOBO data for {sn} in period; excluded from pool.")
            continue
        sensor_frames[sn] = df.assign(_times=pd.DatetimeIndex(df["datetime"]))

    stats = {}
    for window, spec in PRE_POST_WINDOWS.items():
        rows = []
        n_empty = 0
        for ev in events:
            start = ev[spec["start_anchor"]] + spec["start_offset"]
            end = ev[spec["end_anchor"]] + spec["end_offset"]

            pooled = []
            n_sensors = 0
            for df in sensor_frames.values():
                mask = (df["_times"] >= start) & (df["_times"] <= end)
                sub = df.loc[mask, ["temp_c", "rh_pct"]]
                if sub.empty:
                    continue
                pooled.append(sub)
                n_sensors += 1

            if not pooled:
                n_empty += 1
                continue

            pool = pd.concat(pooled, ignore_index=True)
            midpoint = start + (end - start) / 2
            rows.append(
                {
                    "midpoint": midpoint,
                    "shower_on": ev["shower_on"],
                    "temp_c_mean": pool["temp_c"].mean(),
                    "temp_c_std": pool["temp_c"].std(),
                    "rh_pct_mean": pool["rh_pct"].mean(),
                    "rh_pct_std": pool["rh_pct"].std(),
                    "n_points": len(pool),
                    "n_sensors": n_sensors,
                }
            )
        stats[window] = pd.DataFrame(rows)
        print(
            f"  {window}: {len(rows)} events with pooled data"
            + (f", {n_empty} events with no points (skipped)" if n_empty else "")
        )
    return stats


def compute_pre_post_change(stats: dict, column: str) -> dict:
    """
    Compute paired pre->post change statistics for one variable.

    Events are paired on shower_on, which is present in both the pre and post
    window stats. For each paired event, the change is the post pooled mean
    minus the pre pooled mean. A paired t-test checks whether the mean change
    differs from zero.

    Parameters
    ----------
    stats : dict
        Output of build_pre_post_stats.
    column : str
        Value column to summarize ("temp_c" or "rh_pct").

    Returns
    -------
    dict
        n (paired event count), mean_diff, std_diff, t_stat, p_value.
        t_stat and p_value are NaN when n < 2 (ttest_rel requires at least 2
        paired observations).
    """
    mean_col = f"{column}_mean"
    pre = stats["pre"][["shower_on", mean_col]].rename(columns={mean_col: "pre_mean"})
    post = stats["post"][["shower_on", mean_col]].rename(columns={mean_col: "post_mean"})
    merged = pre.merge(post, on="shower_on", how="inner")

    diff = merged["post_mean"] - merged["pre_mean"]
    n = len(merged)
    result = {
        "n": n,
        "mean_diff": diff.mean() if n else np.nan,
        "std_diff": diff.std() if n else np.nan,
        "t_stat": np.nan,
        "p_value": np.nan,
    }
    if n >= 2:
        result["t_stat"], result["p_value"] = ttest_rel(merged["post_mean"], merged["pre_mean"])
    return result


def _format_p_value(p_value: float) -> str:
    """Format a p-value for display, using '< 0.001' below that threshold."""
    if not np.isfinite(p_value):
        return "n/a"
    if p_value < 0.001:
        return "p < 0.001"
    return f"p = {sf.fmt_fig(p_value)}"


def make_pre_post_summary_div(stats: dict, kind: str) -> Div:
    """
    Build a Div summarizing paired pre->post change for one variable.

    Parameters
    ----------
    stats : dict
        Output of build_pre_post_stats.
    kind : str
        Key into PRE_POST_FIGURES ("temp" or "rh").

    Returns
    -------
    Div
        Bokeh Div with the paired mean change, std, n, and t-test result.
    """
    cfg = PRE_POST_FIGURES[kind]
    change = compute_pre_post_change(stats, cfg["column"])
    unit = "°C" if kind == "temp" else " pct RH"

    if change["n"] == 0:
        text = "Paired pre/post change: no events with data in both windows."
    else:
        mean_str = sf.fmt_fig(change["mean_diff"])
        std_str = sf.fmt_fig(change["std_diff"]) if change["n"] > 1 else "n/a"
        p_str = _format_p_value(change["p_value"])
        t_str = sf.fmt_fig(change["t_stat"]) if np.isfinite(change["t_stat"]) else "n/a"
        text = (
            f"Paired pre→post change (n = {change['n']} events): "
            f"mean Δ = {mean_str} ± {std_str}{unit} "
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


def make_figure(frames: dict, kind: str) -> None:
    """
    Build and save one onset/decay Bokeh figure (temperature or RH).

    Parameters
    ----------
    frames : dict
        Output of build_windowed_frames.
    kind : str
        Key into FIGURES ("temp" or "rh").
    """
    cfg = FIGURES[kind]
    column = cfg["column"]

    output_path = FIGURE_DIR / cfg["filename"]
    output_file(str(output_path), title=cfg["title_tag"])

    fig = figure(
        x_axis_type="datetime",
        x_axis_label="Date and time",
        y_axis_label=cfg["y_axis_label"],
        tools="pan,box_zoom,wheel_zoom,reset,save",
    )

    hover = HoverTool(
        tooltips=[
            ("Sensor", "@sensor"),
            ("Time", "@time_str"),
            ("Value", f"@value{{{cfg['value_fmt']}}}"),
        ]
    )
    fig.add_tools(hover)

    for sn in SENSOR_ORDER:
        color = SENSOR_COLOR[sn]
        label = SENSOR_LABELS[sn]
        for window, spec in WINDOWS.items():
            wdf = frames[sn][window]
            if wdf.empty:
                continue
            source = ColumnDataSource(
                data={
                    "time": wdf["datetime"],
                    "value": wdf[column],
                    "sensor": [label] * len(wdf),
                    "window": [window] * len(wdf),
                    "time_str": wdf["datetime"].dt.strftime("%Y-%m-%d %H:%M:%S"),
                }
            )
            # Both windows share one legend entry per sensor so click-to-hide
            # toggles onset and decay together. Marker encodes window type:
            # circle = onset, square = decay.
            fig.scatter(
                "time",
                "value",
                source=source,
                marker=spec["marker"],
                size=6,
                color=color,
                alpha=0.7,
                legend_label=label,
            )

    # Shared MODULAIR style: 1600x800, 12pt, no title, click-to-hide legend.
    style_moduair_figure(fig, legend_title="Sensor", legend_location=cfg["legend_location"])

    fig.xaxis.formatter = DatetimeTickFormatter(
        days="%Y-%m-%d", hours="%m-%d %H:%M", minutes="%H:%M"
    )

    save(fig)
    print(f"  Saved {output_path}")


def make_pre_post_figure(stats: dict, kind: str) -> None:
    """
    Build and save one pre/post per-event pooled-mean figure (temperature or RH).

    Each event contributes one pre point and one post point, each the pooled mean
    across all five sensors' raw points in that window, at the window midpoint,
    with a whisker showing pooled mean +/- pooled std. Std, pooled point count,
    and contributing-sensor count are in the hover. Color encodes window: pre and
    post are two colors, with a two-entry legend (Pre-shower, Post-shower).

    Parameters
    ----------
    stats : dict
        Output of build_pre_post_stats.
    kind : str
        Key into PRE_POST_FIGURES ("temp" or "rh").
    """
    cfg = PRE_POST_FIGURES[kind]
    column = cfg["column"]
    mean_col = f"{column}_mean"
    std_col = f"{column}_std"

    output_path = FIGURE_DIR / cfg["filename"]
    output_file(str(output_path), title=cfg["title_tag"])

    fig = figure(
        x_axis_type="datetime",
        x_axis_label="Date and time",
        y_axis_label=cfg["y_axis_label"],
        tools="pan,box_zoom,wheel_zoom,reset,save",
    )

    hover = HoverTool(
        tooltips=[
            ("Window", "@window"),
            ("Time", "@time_str"),
            ("Mean", f"@value{{{cfg['value_fmt']}}}"),
            ("Std", f"@std{{{cfg['value_fmt']}}}"),
            ("Pooled points", "@n_points"),
            ("Sensors", "@n_sensors"),
        ]
    )
    fig.add_tools(hover)

    for window, spec in PRE_POST_WINDOWS.items():
        wdf = stats[window]
        if wdf.empty:
            continue
        color = spec["color"]
        mean_vals = wdf[mean_col]
        std_vals = wdf[std_col].fillna(0.0)
        source = ColumnDataSource(
            data={
                "time": wdf["midpoint"],
                "value": mean_vals,
                "std": std_vals,
                "upper": mean_vals + std_vals,
                "lower": mean_vals - std_vals,
                "window": [spec["legend"]] * len(wdf),
                "n_points": wdf["n_points"],
                "n_sensors": wdf["n_sensors"],
                "time_str": wdf["midpoint"].dt.strftime("%Y-%m-%d %H:%M:%S"),
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
    style_moduair_figure(fig, legend_title="Window", legend_location=cfg["legend_location"])

    fig.xaxis.formatter = DatetimeTickFormatter(
        days="%Y-%m-%d", hours="%m-%d %H:%M", minutes="%H:%M"
    )

    summary_div = make_pre_post_summary_div(stats, kind)
    save(bokeh_column(fig, summary_div))
    print(f"  Saved {output_path}")


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    print("\n" + "=" * 70)
    print("HOBO UX100 Onset/Decay Time Series Figures")
    print("=" * 70)
    print(f"Period: {PERIOD_START} to {PERIOD_END}")
    print(f"Sensors: {', '.join(f'{k} ({v})' for k, v in SENSOR_LABELS.items())}")

    print("\nLoading shower events...")
    events = load_events_in_period(PERIOD_START, PERIOD_END)
    print(f"  {len(events)} shower events in period")
    if not events:
        print("ERROR: no shower events in the period; nothing to plot.")
        sys.exit(1)

    print("\nBuilding windowed sensor frames...")
    frames = build_windowed_frames(events)

    total_points = sum(len(frames[sn][w]) for sn in SENSOR_ORDER for w in WINDOWS)
    if total_points == 0:
        print("ERROR: no HOBO points fell inside any onset/decay window.")
        sys.exit(1)

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    print("\nBuilding onset/decay figures...")
    make_figure(frames, "temp")
    make_figure(frames, "rh")

    print("\nBuilding pre/post per-event mean frames...")
    stats = build_pre_post_stats(events)

    pre_post_points = sum(len(stats[w]) for w in PRE_POST_WINDOWS)
    if pre_post_points == 0:
        print("  [WARN] No events with data in any pre/post window; skipping pre/post figures.")
    else:
        print("\nBuilding pre/post figures...")
        make_pre_post_figure(stats, "temp")
        make_pre_post_figure(stats, "rh")

    print("\n" + "=" * 70)
    print("Done")
    print("=" * 70)


if __name__ == "__main__":
    main()
