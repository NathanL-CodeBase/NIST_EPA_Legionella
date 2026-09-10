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
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource, DatetimeTickFormatter, HoverTool
from bokeh.models.annotations import Whisker
from bokeh.plotting import figure, output_file, save

# Add project root to path for src/ imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.env_data_loader import (  # noqa: E402
    identify_shower_events,
    load_hobo_data,
    load_shower_log,
)
from src.data_paths import get_data_root  # noqa: E402
from src.plot_style import SENSOR_COLORS, style_moduair_figure  # noqa: E402

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

# Pre/post windows for the per-event mean figures. These match the windows used
# in scripts/rh_temp_other_analysis.py (PRE_SHOWER_MINUTES = 30 before shower_on,
# POST_SHOWER_HOURS = 2 after shower_off). Each event contributes one mean with
# std per window, placed at the window midpoint. Marker convention parallels the
# onset/decay figures: circle = pre-shower, square = post-shower.
PRE_SHOWER_MINUTES = 30
POST_SHOWER_HOURS = 2
PRE_POST_WINDOWS = {
    "pre": {
        "start_anchor": "shower_on",
        "start_offset": pd.Timedelta(minutes=-PRE_SHOWER_MINUTES),
        "end_anchor": "shower_on",
        "end_offset": pd.Timedelta(0),
        "marker": "circle",
    },
    "post": {
        "start_anchor": "shower_off",
        "start_offset": pd.Timedelta(0),
        "end_anchor": "shower_off",
        "end_offset": pd.Timedelta(hours=POST_SHOWER_HOURS),
        "marker": "square",
    },
}

# Output directory: output/plots/hobo/ under data_root from data_config.json.
FIGURE_DIR = Path(get_data_root()) / "output" / "plots" / "hobo"

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
    Compute per-event pre and post window mean, std, and midpoint per sensor.

    For each sensor and each event, the pre window is the 30 minutes before
    shower_on and the post window is the 2 hours after shower_off, matching
    PRE_SHOWER_MINUTES and POST_SHOWER_HOURS in scripts/rh_temp_other_analysis.py.
    Windows are inclusive on both ends. The marker for an event is placed at the
    window midpoint. Events with no sensor points in a window are counted and
    skipped for that sensor and window.

    Parameters
    ----------
    events : list
        Shower event dicts with shower_on and shower_off.

    Returns
    -------
    dict
        Nested mapping {sensor_key: {window: DataFrame}} where each DataFrame has
        midpoint, temp_c_mean, temp_c_std, rh_pct_mean, rh_pct_std, shower_on,
        and n_points columns, one row per event with data in that window.
    """
    stats = {}
    for sn in SENSOR_ORDER:
        df = load_hobo_data(sn, PERIOD_START, PERIOD_END)
        stats[sn] = {}
        if df.empty:
            print(f"  [WARN] No HOBO data for {sn} in period; skipping sensor.")
            for window in PRE_POST_WINDOWS:
                stats[sn][window] = pd.DataFrame()
            continue

        times = pd.DatetimeIndex(df["datetime"])
        for window, spec in PRE_POST_WINDOWS.items():
            rows = []
            n_empty = 0
            for ev in events:
                start = ev[spec["start_anchor"]] + spec["start_offset"]
                end = ev[spec["end_anchor"]] + spec["end_offset"]
                mask = (times >= start) & (times <= end)
                sub = df.loc[mask]
                if sub.empty:
                    n_empty += 1
                    continue
                midpoint = start + (end - start) / 2
                rows.append(
                    {
                        "midpoint": midpoint,
                        "shower_on": ev["shower_on"],
                        "temp_c_mean": sub["temp_c"].mean(),
                        "temp_c_std": sub["temp_c"].std(),
                        "rh_pct_mean": sub["rh_pct"].mean(),
                        "rh_pct_std": sub["rh_pct"].std(),
                        "n_points": len(sub),
                    }
                )
            stats[sn][window] = pd.DataFrame(rows)
            print(
                f"  {sn} {window}: {len(rows)} events with data"
                + (f", {n_empty} events with no points (skipped)" if n_empty else "")
            )
    return stats


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
    Build and save one pre/post per-event mean figure (temperature or RH).

    Each event contributes one circle (pre-shower mean) and one square
    (post-shower mean) at the respective window midpoint, with a whisker showing
    mean +/- std. Std is also in the hover. Five sensors are distinguished by
    color; circle = pre and square = post parallel the onset/decay figures.

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
            ("Sensor", "@sensor"),
            ("Window", "@window"),
            ("Time", "@time_str"),
            ("Mean", f"@value{{{cfg['value_fmt']}}}"),
            ("Std", f"@std{{{cfg['value_fmt']}}}"),
        ]
    )
    fig.add_tools(hover)

    for sn in SENSOR_ORDER:
        color = SENSOR_COLOR[sn]
        label = SENSOR_LABELS[sn]
        for window, spec in PRE_POST_WINDOWS.items():
            wdf = stats[sn][window]
            if wdf.empty:
                continue
            mean_vals = wdf[mean_col]
            std_vals = wdf[std_col].fillna(0.0)
            source = ColumnDataSource(
                data={
                    "time": wdf["midpoint"],
                    "value": mean_vals,
                    "std": std_vals,
                    "upper": mean_vals + std_vals,
                    "lower": mean_vals - std_vals,
                    "sensor": [label] * len(wdf),
                    "window": [window] * len(wdf),
                    "time_str": wdf["midpoint"].dt.strftime("%Y-%m-%d %H:%M:%S"),
                }
            )
            # One legend entry per sensor; circle = pre, square = post share it
            # so click-to-hide toggles both windows together.
            fig.scatter(
                "time",
                "value",
                source=source,
                marker=spec["marker"],
                size=7,
                color=color,
                alpha=0.8,
                legend_label=label,
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
    style_moduair_figure(fig, legend_title="Sensor", legend_location=cfg["legend_location"])

    fig.xaxis.formatter = DatetimeTickFormatter(
        days="%Y-%m-%d", hours="%m-%d %H:%M", minutes="%H:%M"
    )

    save(fig)
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

    pre_post_points = sum(len(stats[sn][w]) for sn in SENSOR_ORDER for w in PRE_POST_WINDOWS)
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
