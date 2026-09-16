#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle Analysis — Event Plots and Summary Bar Charts
=======================================================

Individual event concentration/emission figures and cross-event bar charts for
the EPA Legionella particle analysis.  Boxplot and comparison functions live in
companion modules that are re-exported here for backward-compatible imports:

  - plot_particle_boxplots.py: temperature-axis boxplots and shower-head boxplot
  - plot_comparison.py: categorical comparison boxplots (spray pattern, head type,
    mannequin, door position, fan status)

Functions defined here:
    - plot_particle_decay_event: Four-panel individual event figure
      (concentration top panel + three emission panels by bin group)
    - plot_penetration_summary: Bar chart of penetration factors by size bin
    - plot_deposition_summary: Bar chart of other process rates by size bin
    - plot_emission_summary: Bar chart of emission rates by size bin
    - plot_size_distribution_summary: Multi-panel summary of all three metrics

Re-exported from plot_particle_boxplots:
    - plot_emission_boxplot
    - plot_deposition_rate_boxplot
    - plot_emission_rate_boxplot
    - plot_penetration_factor_boxplot
    - plot_emission_etotal_by_metric_boxplot
    - plot_emission_etotal_by_showerhead_boxplot

Re-exported from plot_comparison:
    - plot_spray_pattern_comparison_boxplots
    - plot_shower_head_comparison_boxplots
    - plot_mannequin_comparison_boxplots
    - plot_door_comparison_boxplots
    - plot_fan_comparison_boxplots

Plot Features:
    - Four-panel event plots: concentration time series (top) + three per-step
      emission panels (bins 0–2 small, 3–6 medium, 7–11 large)
    - Color-coded particle size bins; all bins plotted as solid lines regardless
      of beta validity
    - Optional outdoor PM overlay: dotted lines at 55% alpha for events in
      OUTDOOR_PM_EVENTS (e.g., event 77)
    - Shaded deposition analysis window (2 hr post-shower)
    - Shower ON/OFF dotted markers distinct from fit/predicted dashed lines
    - Decay R² values listed in top-panel text box alongside lambda and valid-bin count
    - Log-scale concentration axis for wide dynamic range
    - Emission panel x-axis matches concentration panel (shared time axis)
    - Emission panel y-axis clipped to 2nd–98th percentile of E_per_step data

Output Files (from this module):
    - {test_name}_particle_decay.png: Individual event concentration + emission figure
    - penetration_summary.png: Bar chart of penetration factors across all events
    - deposition_summary.png: Bar chart of other process (beta) rates across all events
    - emission_summary.png: Bar chart of emission rates across all events
    - size_distribution_summary.png: Multi-panel summary of all metrics

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Date: 2026
"""

from datetime import timedelta
from pathlib import Path
from typing import Dict, Optional, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import src.sig_figs as sf
from src.event_manager import sort_config_keys_by_water_temp
from src.plot_style import (
    COLORS,
    FONT_SIZE_LABEL,
    FONT_SIZE_LEGEND,
    FONT_SIZE_TICK,
    FONT_SIZE_TITLE,
    LINE_WIDTH_DATA,
    LINE_WIDTH_FIT,
    SENSOR_COLORS,
    TITLE_FONTWEIGHT,
    WINDOW_ALPHA,
    add_shaded_window,
    add_shower_off_marker,
    add_shower_on_marker,
    apply_style,
    create_figure,
    format_datetime_axis,
    format_test_name_for_title,
    format_title,
    get_config_color,
    save_figure,
)


def plot_particle_decay_event(
    particle_data: pd.DataFrame,
    event: Dict,
    particle_bins: Dict,
    result: Dict,
    output_path: Path,
    event_number: int,
    test_name: Optional[str] = None,
    show_outdoor: bool = False,
) -> None:
    """
    Plot particle concentration decay for a single event showing all bins.

    Creates a two-panel figure matching the CO2 analytical plot style:

    Top panel (ax1) — concentration time series:
      - Measured indoor particle concentrations for all bins
      - Two continuous predicted Ct curves per valid bin (emission + decay
        phases), one per air-change-rate source: outside (dashed) and entry
        (dotted), bracketing the measured curve instead of a single blend
      - Shower ON/OFF markers and shaded deposition window
      - Optional outdoor PM overlay (when show_outdoor=True): plots
        opc_binN_outside columns using the same per-bin colours with
        reduced alpha and dotted lines

    Bottom panels (ax2–ax4) — per-step emission rates by bin group:
      - Per-step E_t values as faint lines for each valid bin, both sources
      - E_mean as a horizontal line spanning shower_on to peak_time, dashed
        for the outside source and dotted for the entry source
      - Emission R² annotation in each panel, per source

    Parameters:
        particle_data: DataFrame with particle concentrations
        event: Event timing dictionary
        particle_bins: Dictionary of all particle bin information
        result: Analysis results for this event (all bins)
        output_path: Path to save the figure
        event_number: Event number for title
        test_name: Test name for title (e.g., "0114_HW_Morning_R01")
        show_outdoor: If True, overlay outdoor PM concentration on the top
            panel using opc_binN_outside columns (same colours, dotted lines).
    """
    apply_style()

    # Extract data for plotting window (1 hour before shower to 1 hour after deposition end)
    plot_start = event["shower_on"] - timedelta(hours=1)
    plot_end = event["deposition_end"] + timedelta(hours=1)

    mask = (particle_data["datetime"] >= plot_start) & (particle_data["datetime"] <= plot_end)
    plot_data = particle_data[mask].copy()

    if plot_data.empty:
        print(f"    Warning: No data for event {event_number}")
        return

    # Four-panel figure: concentration (top) + three emission panels by bin group
    fig, axes_array = create_figure(nrows=4, ncols=1, figsize=(12, 14), height_ratios=[2, 1, 1, 1])
    ax1, ax2, ax3, ax4 = cast(np.ndarray, axes_array)

    lambda_outside = result.get("lambda_outside", np.nan)
    lambda_entry = result.get("lambda_entry", np.nan)
    # Outside/entry line styles used throughout this figure: dashed brackets
    # the measured curve from above, dotted from below (or vice versa),
    # same bin color for both.
    _SOURCE_LINESTYLE = {"outside": "--", "entry": ":"}

    # =========================================================================
    # Top panel: Particle concentrations with decay predictions
    # =========================================================================
    for bin_num, bin_info in particle_bins.items():
        col_inside = f"{bin_info['column']}_inside"
        color = SENSOR_COLORS[bin_num % len(SENSOR_COLORS)]

        if col_inside in plot_data.columns:
            # Plot raw instrument data — always solid regardless of beta validity
            ax1.plot(
                plot_data["datetime"],
                plot_data[col_inside],
                label=f"Bin {bin_num} ({bin_info['name']} µm)",
                color=color,
                linewidth=LINE_WIDTH_DATA,
                linestyle="-",
                alpha=0.9,
            )

            # Plot continuous predicted Ct: emission phase then decay phase
            # as two segments of the same simulation (they connect at peak_time).
            # Each bin carries two predictions, one per lambda source (outside,
            # entry), bracketing the measured curve rather than a single blend.
            for source, linestyle in _SOURCE_LINESTYLE.items():
                prefix = f"bin{bin_num}_{source}"
                emission_dts = result.get(f"{prefix}_emission_datetimes", [])
                emission_pred = result.get(f"{prefix}_emission_predicted", [])
                decay_dts = result.get(f"{prefix}_decay_datetimes", [])
                decay_pred = result.get(f"{prefix}_decay_predicted", [])

                if len(emission_dts) > 0 and len(emission_pred) > 0:
                    ep_arr = np.array(emission_pred, dtype=float)
                    ep_arr = np.where(ep_arr > 0, ep_arr, np.nan)
                    if np.any(np.isfinite(ep_arr)):
                        ax1.plot(
                            pd.to_datetime(emission_dts),
                            ep_arr,
                            color=color,
                            linewidth=LINE_WIDTH_FIT,
                            linestyle=linestyle,
                            alpha=0.8,
                        )
                if len(decay_dts) > 0 and len(decay_pred) > 0:
                    dp_arr = np.array(decay_pred, dtype=float)
                    dp_arr = np.where(dp_arr > 0, dp_arr, np.nan)
                    if np.any(np.isfinite(dp_arr)):
                        ax1.plot(
                            pd.to_datetime(decay_dts),
                            dp_arr,
                            color=color,
                            linewidth=LINE_WIDTH_FIT,
                            linestyle=linestyle,
                            alpha=0.8,
                        )

    # Add legend entries for the two predicted Ct lines (one per source)
    has_predictions = any(
        len(result.get(f"bin{bn}_{source}_emission_predicted", [])) > 0
        or len(result.get(f"bin{bn}_{source}_decay_predicted", [])) > 0
        for bn in particle_bins.keys()
        for source in _SOURCE_LINESTYLE
    )
    if has_predictions:
        for source, linestyle in _SOURCE_LINESTYLE.items():
            ax1.plot(
                [],
                [],
                color="gray",
                linestyle=linestyle,
                linewidth=LINE_WIDTH_FIT,
                label=f"Predicted Ct ({source})",
            )

    # Optional outdoor PM concentration overlay
    if show_outdoor:
        outdoor_added = False
        for bin_num, bin_info in particle_bins.items():
            col_outside = f"{bin_info['column']}_outside"
            color = SENSOR_COLORS[bin_num % len(SENSOR_COLORS)]
            if col_outside in plot_data.columns:
                ax1.plot(
                    plot_data["datetime"],
                    plot_data[col_outside],
                    color=color,
                    linewidth=LINE_WIDTH_DATA * 0.8,
                    linestyle=":",
                    alpha=0.55,
                )
                outdoor_added = True
        if outdoor_added:
            ax1.plot(
                [],
                [],
                color="gray",
                linestyle=":",
                linewidth=LINE_WIDTH_DATA * 0.8,
                alpha=0.55,
                label="Outdoor Concentration",
            )

    # Add shaded window for deposition analysis period
    add_shaded_window(
        ax1,
        event["shower_off"],
        event["deposition_end"],
        color=COLORS["post_shower"],
        label="Deposition window (2 hr)",
        alpha=WINDOW_ALPHA,
    )

    # Add shower ON/OFF markers
    add_shower_on_marker(ax1, event["shower_on"], label="Shower ON")
    add_shower_off_marker(ax1, event["shower_off"], label="Shower OFF")

    # Axis formatting
    ax1.set_ylabel("Particle Concentration (#/cm³)", fontsize=FONT_SIZE_LABEL)

    # Use consistent title formatting
    if test_name:
        formatted_name = format_test_name_for_title(test_name)
        title = f"Event {event_number:02d} - {formatted_name}: PM Decay"
    else:
        title = format_title(
            "Particle Decay - All Size Bins",
            event_number=event_number,
            event_datetime=event["shower_on"],
        )
    ax1.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight=TITLE_FONTWEIGHT)

    # Count valid bins (valid in at least one source) and build a per-bin,
    # per-source decay R² summary for the text box
    valid_bins = 0
    decay_r2_lines = []
    for bin_num in particle_bins.keys():
        source_r2 = {}
        bin_valid = False
        for source in _SOURCE_LINESTYLE:
            beta_val = result.get(f"bin{bin_num}_{source}_beta_other", np.nan)
            if not np.isnan(beta_val):
                bin_valid = True
                r2_val = result.get(f"bin{bin_num}_{source}_beta_other_r_squared", np.nan)
                source_r2[source] = (
                    sf.fmt_fig(r2_val, fallback=".3f") if not np.isnan(r2_val) else "N/A"
                )
        if bin_valid:
            valid_bins += 1
            r2_text = ", ".join(f"{s}={v}" for s, v in source_r2.items())
            decay_r2_lines.append(f" B{bin_num}: R²({r2_text})")

    # Build text box content: lambda (both sources), valid-bin count, and
    # per-bin decay R²
    textstr = (
        f"λ_outside = {sf.fmt_fig(lambda_outside, fallback='.4f')} h⁻¹\n"
        f"λ_entry = {sf.fmt_fig(lambda_entry, fallback='.4f')} h⁻¹\n"
    )
    if decay_r2_lines:
        textstr += "\nDecay R²:\n" + "\n".join(decay_r2_lines)

    props = dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="gray")
    ax1.text(
        0.02,
        0.98,
        textstr,
        transform=ax1.transAxes,
        fontsize=FONT_SIZE_LEGEND,
        verticalalignment="top",
        bbox=props,
    )

    ax1.legend(
        loc="upper right",
        fontsize=FONT_SIZE_LEGEND - 1,
        framealpha=0.9,
        ncol=1,
    )
    # Compute y limits from all positive plotted values (measured + predicted)
    _all_pos_y = []
    for _bn, _bi in particle_bins.items():
        _col = f"{_bi['column']}_inside"
        if _col in plot_data.columns:
            _vals = plot_data[_col].dropna().values
            _all_pos_y.extend(_vals[_vals > 0].tolist())
        for _source in _SOURCE_LINESTYLE:
            for _key in (
                f"bin{_bn}_{_source}_emission_predicted",
                f"bin{_bn}_{_source}_decay_predicted",
            ):
                _all_pos_y.extend(
                    [
                        v
                        for v in result.get(_key, [])
                        if isinstance(v, (int, float)) and v > 0 and not np.isnan(v)
                    ]
                )
    ax1.set_yscale("log")
    if _all_pos_y:
        _ylo = min(_all_pos_y) * 0.3
        _yhi = max(_all_pos_y) * 3.0
        ax1.set_ylim(bottom=max(_ylo, 1e-6), top=_yhi)
    else:
        ax1.set_ylim(bottom=1e-3)
    ax1.grid(True, alpha=0.3, which="both")
    ax1.tick_params(labelsize=FONT_SIZE_TICK)
    format_datetime_axis(ax1)
    ax1.set_xlim(plot_start, plot_end)

    # =========================================================================
    # Emission panels: per-step emission rates split by bin group
    # ax2 = Bins 0–2  (small particles)
    # ax3 = Bins 3–6  (medium particles)
    # ax4 = Bins 7–11 (large particles)
    # Each panel has its own y-axis scale so all groups are readable.
    # =========================================================================
    bins_small = [bn for bn in particle_bins.keys() if bn <= 2]
    bins_medium = [bn for bn in particle_bins.keys() if 3 <= bn <= 6]
    bins_large = [bn for bn in particle_bins.keys() if bn >= 7]

    # Shared x-axis limits: 10 min before shower_on to 10 min after latest peak_time
    peak_times = []
    for bn in particle_bins.keys():
        pt = result.get(f"bin{bn}_peak_time", None)
        if pt is not None:
            peak_times.append(pd.Timestamp(pt))
    ax_e_left = event["shower_on"] - timedelta(minutes=10)
    ax_e_right = (
        max(peak_times) + timedelta(minutes=10)
        if peak_times
        else event["shower_off"] + timedelta(minutes=10)
    )

    def _populate_emission_panel(ax, bin_group, panel_label):
        """Plot E_per_step scatter and E_mean line (per source) for a subset of bins."""
        has_data = False
        r2_lines = []
        e_steps_panel = []
        e_means_panel = []

        for bin_num in bin_group:
            bin_info = particle_bins[bin_num]
            color = SENSOR_COLORS[bin_num % len(SENSOR_COLORS)]
            peak_time = result.get(f"bin{bin_num}_peak_time", None)
            source_r2 = {}

            for source, linestyle in _SOURCE_LINESTYLE.items():
                prefix = f"bin{bin_num}_{source}"
                E_times = result.get(f"{prefix}_E_times", [])
                E_per_step = result.get(f"{prefix}_E_per_step", [])
                E_mean_val = result.get(f"{prefix}_E_mean", np.nan)
                E_r2_val = result.get(f"{prefix}_E_r_squared", np.nan)

                if len(E_times) > 0 and len(E_per_step) > 0:
                    has_data = True
                    e_steps_panel.extend(E_per_step)
                    # Per-step E_t as a faint line (all values, including negative)
                    ax.plot(
                        pd.to_datetime(E_times),
                        np.array(E_per_step),
                        color=color,
                        linewidth=0.8,
                        alpha=0.35,
                    )

                # E_mean as a horizontal line spanning shower_on -> peak_time,
                # dashed for outside / dotted for entry (matches top panel)
                if not np.isnan(E_mean_val) and peak_time is not None:
                    has_data = True
                    e_means_panel.append(E_mean_val)
                    ax.hlines(
                        E_mean_val,
                        event["shower_on"],
                        pd.Timestamp(peak_time),
                        color=color,
                        linewidth=LINE_WIDTH_FIT,
                        linestyle=linestyle,
                        alpha=0.9,
                    )
                    source_r2[source] = (
                        sf.fmt_fig(E_r2_val, fallback=".3f") if not np.isnan(E_r2_val) else "N/A"
                    )

            if source_r2:
                r2_text = ", ".join(f"{s}={v}" for s, v in source_r2.items())
                r2_lines.append(f"B{bin_num}: R²({r2_text})")

        add_shower_on_marker(ax, event["shower_on"])
        add_shower_off_marker(ax, event["shower_off"])

        ax.axhline(0, color="gray", linewidth=0.8, linestyle=":", alpha=0.6)
        ax.set_ylabel(f"E (#/cm³·min)\n{panel_label}", fontsize=FONT_SIZE_LABEL - 1)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=False)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=FONT_SIZE_TICK)
        ax.set_xlim(ax_e_left, ax_e_right)
        format_datetime_axis(ax, interval_minutes=5)

        # Per-panel percentile y-limits so E_mean dashed lines are never clipped
        if e_steps_panel:
            e_arr = np.array([v for v in e_steps_panel if not np.isnan(v)], dtype=float)
            if len(e_arr) >= 4:
                p2 = float(np.percentile(e_arr, 2))
                p98 = float(np.percentile(e_arr, 98))
                if e_means_panel:
                    p2 = min(p2, min(e_means_panel))
                    p98 = max(p98, max(e_means_panel))
                margin = max((p98 - p2) * 0.15, abs(p98) * 0.05, 1e-6)
                ax.set_ylim(p2 - margin, p98 + margin)

        if r2_lines:
            r2_text = "Emission R²:\n" + "\n".join(r2_lines)
            props_r2 = dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="gray")
            ax.text(
                0.02,
                0.98,
                r2_text,
                transform=ax.transAxes,
                fontsize=FONT_SIZE_LEGEND - 1,
                verticalalignment="top",
                bbox=props_r2,
            )

        if not has_data:
            ax.text(
                0.5,
                0.5,
                "No emission data available",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=FONT_SIZE_LABEL,
                color="gray",
            )

    _populate_emission_panel(ax2, bins_small, "Bins 0–2")
    _populate_emission_panel(ax3, bins_medium, "Bins 3–6")
    _populate_emission_panel(ax4, bins_large, "Bins 7–11")

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close(fig)


# Configuration for the three per-bin summary bar charts.
# Keys: col_template, ylabel, title, label_fmt, label_offset, ylim, log_scale,
#       has_source (bool — True for metrics computed once per lambda source;
#       col_template includes a {source} placeholder for those entries).
_BAR_CHART_CONFIG = {
    "penetration": dict(
        col_template="bin{n}_p_mean",
        ylabel="Penetration Factor (p)",
        title="Penetration Factor by Particle Size\n(Mean ± Std Dev)",
        label_fmt=".3f",
        label_offset=0.02,
        ylim=(0, 1.1),
        log_scale=False,
        has_source=False,
    ),
    "deposition": dict(
        col_template="bin{n}_{source}_beta_other",
        ylabel="Other Process Rate β (h⁻¹)",
        title="Other Process Rate by Particle Size\n(Mean ± Std Dev)",
        label_fmt=".2f",
        label_offset=0.1,
        ylim=None,
        log_scale=False,
        has_source=True,
    ),
    "emission": dict(
        col_template="bin{n}_{source}_E_mean",
        ylabel="Emission Rate E (#/min)",
        title="Shower Emission Rate by Particle Size\n(Mean ± Std Dev)",
        label_fmt=".1e",
        label_offset=0,
        ylim=None,
        log_scale="auto",
        has_source=True,
    ),
}


def _plot_summary_bar_chart(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    output_path: Path,
    cfg: dict,
) -> None:
    """
    Shared implementation for per-bin summary bar charts.

    Handles all three metric types (penetration, other process, emission) via a
    configuration dict.  If config_key is present, creates one subplot per
    configuration; otherwise draws a single panel. Lambda-dependent metrics
    (``cfg["has_source"]``) are drawn twice, once per air-change-rate source,
    with an ``_outside``/``_entry`` suffix appended to the output filename;
    source-independent metrics (penetration factor) are drawn once.

    Parameters:
        results_df: DataFrame with analysis results.
        particle_bins: Dictionary of particle bin information.
        output_path: Path to save the figure.
        cfg: Configuration dict from _BAR_CHART_CONFIG.
    """
    apply_style()

    bin_nums = list(particle_bins.keys())
    bin_labels = [particle_bins[i]["name"] for i in bin_nums]

    has_config = "config_key" in results_df.columns
    if has_config:
        config_keys = sort_config_keys_by_water_temp(
            list(results_df["config_key"].dropna().unique())
        )
        n_configs = len(config_keys)
    else:
        config_keys = ["All"]
        n_configs = 1

    sources = ("outside", "entry") if cfg.get("has_source", True) else (None,)
    for source in sources:
        _plot_summary_bar_chart_one_source(
            results_df, particle_bins, output_path, cfg, config_keys, n_configs, source
        )


def _plot_summary_bar_chart_one_source(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    output_path: Path,
    cfg: dict,
    config_keys: list,
    n_configs: int,
    source: "Optional[str]",
) -> None:
    """Draw one summary bar chart (all configs/subplots) for a single lambda source."""
    bin_nums = list(particle_bins.keys())
    bin_labels = [particle_bins[i]["name"] for i in bin_nums]
    has_config = "config_key" in results_df.columns
    source_note = f" ({source} source)" if source else ""

    # Scale figure width with bin count (at least 1.4 in per bin, min 14 in wide)
    fig_width = max(14, len(bin_nums) * 1.4)
    if n_configs > 1:
        fig, axes = plt.subplots(n_configs, 1, figsize=(fig_width, 6 * n_configs), squeeze=False)
        axes = axes.flatten()
    else:
        fig, _ax = create_figure(figsize=(fig_width, 6))
        if isinstance(_ax, list):
            _ax = _ax[0]
        axes = [_ax]

    for ax_idx, config_key in enumerate(config_keys):
        ax = axes[ax_idx]

        if has_config and config_key != "All":
            config_df = results_df[results_df["config_key"] == config_key]
            color = get_config_color(config_key)
        else:
            config_df = results_df
            color = SENSOR_COLORS[0]

        x = np.arange(len(bin_nums))
        means = []
        stds = []
        for bin_num in bin_nums:
            col = cfg["col_template"].format(n=bin_num, source=source)
            valid_values = config_df[col].dropna() if col in config_df.columns else pd.Series([], dtype=float)
            means.append(valid_values.mean() if len(valid_values) > 0 else np.nan)
            stds.append(valid_values.std() if len(valid_values) > 0 else np.nan)

        bars = ax.bar(
            x,
            means,
            yerr=stds,
            capsize=5,
            color=color,
            alpha=0.7,
            edgecolor="black",
            linewidth=0.5,
        )

        # Value labels above bars
        for bar, mean, std in zip(bars, means, stds):
            if mean > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + std + cfg["label_offset"],
                    sf.fmt_fig(mean, fallback=cfg["label_fmt"]),
                    ha="center",
                    va="bottom",
                    fontsize=FONT_SIZE_TICK - 2,
                )

        ax.set_xlabel("Particle Size Bin (µm)", fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel(cfg["ylabel"], fontsize=FONT_SIZE_LABEL)
        if n_configs > 1:
            ax.set_title(
                f"Configuration: {config_key} (n={len(config_df)}){source_note}",
                fontsize=FONT_SIZE_TITLE,
                fontweight=TITLE_FONTWEIGHT,
            )
        else:
            ax.set_title(
                cfg["title"] + source_note, fontsize=FONT_SIZE_TITLE, fontweight=TITLE_FONTWEIGHT
            )

        ax.set_xticks(x)
        ax.set_xticklabels(bin_labels, rotation=45, ha="right", fontsize=FONT_SIZE_TICK - 1)

        if cfg["ylim"] is not None:
            ax.set_ylim(*cfg["ylim"])
        elif cfg["log_scale"] == "auto":
            valid_means = [m for m in means if m > 0]
            if valid_means and max(valid_means) / min(valid_means) > 100:
                ax.set_yscale("log")

        ax.grid(True, alpha=0.3, axis="y")
        ax.tick_params(labelsize=FONT_SIZE_TICK)

    if n_configs > 1:
        fig.suptitle(
            cfg["title"] + source_note,
            fontsize=FONT_SIZE_TITLE + 2,
            fontweight=TITLE_FONTWEIGHT,
        )

    source_output = (
        output_path.parent / f"{output_path.stem}_{source}{output_path.suffix}"
        if source
        else output_path
    )
    plt.tight_layout(pad=2.0)
    save_figure(fig, source_output)
    plt.close(fig)


def plot_penetration_summary(
    results_df: pd.DataFrame, particle_bins: Dict, output_path: Path
) -> None:
    """Bar chart of penetration factors across all bins (mean ± std per bin).

    Parameters:
        results_df: DataFrame with analysis results.
        particle_bins: Dictionary of particle bin information.
        output_path: Path to save the figure.
    """
    _plot_summary_bar_chart(
        results_df, particle_bins, output_path, _BAR_CHART_CONFIG["penetration"]
    )


def plot_deposition_summary(
    results_df: pd.DataFrame, particle_bins: Dict, output_path: Path
) -> None:
    """Bar chart of other process rates across all bins (mean ± std per bin).

    Parameters:
        results_df: DataFrame with analysis results.
        particle_bins: Dictionary of particle bin information.
        output_path: Path to save the figure.
    """
    _plot_summary_bar_chart(results_df, particle_bins, output_path, _BAR_CHART_CONFIG["deposition"])


def plot_emission_summary(results_df: pd.DataFrame, particle_bins: Dict, output_path: Path) -> None:
    """Bar chart of emission rates across all bins (mean ± std per bin).

    Parameters:
        results_df: DataFrame with analysis results.
        particle_bins: Dictionary of particle bin information.
        output_path: Path to save the figure.
    """
    _plot_summary_bar_chart(results_df, particle_bins, output_path, _BAR_CHART_CONFIG["emission"])


def plot_size_distribution_summary(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    output_path: Path,
    source: str,
) -> None:
    """
    Create multi-panel figure showing all three metrics vs particle size.

    Panel (a) penetration factor is lambda-independent and identical across
    calls. Panels (b) other process rate and (c) emission rate are
    lambda-dependent, so the caller invokes this once per air-change-rate
    source; the output filename gets an ``_outside``/``_entry`` suffix.

    Parameters:
        results_df: DataFrame with analysis results
        particle_bins: Dictionary of particle bin information
        output_path: Path to save the figure (suffixed with ``_{source}``)
        source: Air-change-rate source for panels (b)/(c): "outside" or "entry"
    """
    apply_style()

    bin_nums = list(particle_bins.keys())
    bin_centers = [(particle_bins[i]["min"] + particle_bins[i]["max"]) / 2 for i in bin_nums]

    fig, axes_temp = create_figure(nrows=1, ncols=3, figsize=(15, 5))
    axes = cast(np.ndarray, axes_temp)

    # Panel 1: Penetration factor (lambda-independent, shared across sources)
    p_means = []
    p_stds = []
    for bin_num in bin_nums:
        col = f"bin{bin_num}_p_mean"
        valid_values = results_df[col].dropna()
        p_means.append(valid_values.mean() if len(valid_values) > 0 else np.nan)
        p_stds.append(valid_values.std() if len(valid_values) > 0 else np.nan)

    axes[0].errorbar(
        bin_centers,
        p_means,
        yerr=p_stds,
        marker="o",
        markersize=8,
        capsize=5,
        color=SENSOR_COLORS[0],
        linewidth=LINE_WIDTH_DATA,
    )
    axes[0].set_xlabel("Particle Size (µm)", fontsize=FONT_SIZE_LABEL)
    axes[0].set_ylabel("Penetration Factor (p)", fontsize=FONT_SIZE_LABEL)
    axes[0].set_title(
        "(a) Penetration Factor", fontsize=FONT_SIZE_TITLE, fontweight=TITLE_FONTWEIGHT
    )
    axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim(0, 1.1)

    # Panel 2: Other Process Rate (this source)
    beta_means = []
    beta_stds = []
    for bin_num in bin_nums:
        col = f"bin{bin_num}_{source}_beta_other"
        valid_values = results_df[col].dropna() if col in results_df.columns else pd.Series([], dtype=float)
        beta_means.append(valid_values.mean() if len(valid_values) > 0 else np.nan)
        beta_stds.append(valid_values.std() if len(valid_values) > 0 else np.nan)

    axes[1].errorbar(
        bin_centers,
        beta_means,
        yerr=beta_stds,
        marker="s",
        markersize=8,
        capsize=5,
        color=SENSOR_COLORS[1],
        linewidth=LINE_WIDTH_DATA,
    )
    axes[1].set_xlabel("Particle Size (µm)", fontsize=FONT_SIZE_LABEL)
    axes[1].set_ylabel("Other Process Rate β (h⁻¹)", fontsize=FONT_SIZE_LABEL)
    axes[1].set_title(
        f"(b) Other Process Rate ({source})", fontsize=FONT_SIZE_TITLE, fontweight=TITLE_FONTWEIGHT
    )
    axes[1].grid(True, alpha=0.3)

    # Panel 3: Emission rate (this source)
    E_means = []
    E_stds = []
    for bin_num in bin_nums:
        col = f"bin{bin_num}_{source}_E_mean"
        valid_values = results_df[col].dropna() if col in results_df.columns else pd.Series([], dtype=float)
        E_means.append(valid_values.mean() if len(valid_values) > 0 else np.nan)
        E_stds.append(valid_values.std() if len(valid_values) > 0 else np.nan)

    axes[2].errorbar(
        bin_centers,
        E_means,
        yerr=E_stds,
        marker="^",
        markersize=8,
        capsize=5,
        color=SENSOR_COLORS[2],
        linewidth=LINE_WIDTH_DATA,
    )
    axes[2].set_xlabel("Particle Size (µm)", fontsize=FONT_SIZE_LABEL)
    axes[2].set_ylabel("Emission Rate E (#/min)", fontsize=FONT_SIZE_LABEL)
    axes[2].set_title(
        f"(c) Emission Rate ({source})", fontsize=FONT_SIZE_TITLE, fontweight=TITLE_FONTWEIGHT
    )
    axes[2].grid(True, alpha=0.3)

    # Apply log scale if needed
    if max([e for e in E_means if not np.isnan(e)] + [0]) > 0:
        valid_E = [e for e in E_means if not np.isnan(e) and e > 0]
        if valid_E and max(valid_E) / min(valid_E) > 100:
            axes[2].set_yscale("log")

    for ax in axes:
        ax.tick_params(labelsize=FONT_SIZE_TICK)

    source_output = output_path.parent / f"{output_path.stem}_{source}{output_path.suffix}"
    plt.tight_layout()
    save_figure(fig, source_output)
    plt.close(fig)


# =============================================================================
# Re-exports from companion modules for backward-compatible imports
# =============================================================================
# All imports from src.plot_particle continue to work after the module split.

from src.plot_particle_boxplots import (  # noqa: E402
    plot_deposition_rate_boxplot,
    plot_emission_boxplot,
    plot_emission_etotal_by_metric_boxplot,
    plot_emission_etotal_by_showerhead_boxplot,
    plot_emission_rate_boxplot,
    plot_penetration_factor_boxplot,
)
from src.plot_comparison import (  # noqa: E402
    plot_door_comparison_boxplots,
    plot_fan_comparison_boxplots,
    plot_mannequin_comparison_boxplots,
    plot_shower_head_comparison_boxplots,
    plot_spray_pattern_comparison_boxplots,
)

__all__ = [
    # Event plots and bar charts (defined here)
    "plot_particle_decay_event",
    "plot_penetration_summary",
    "plot_deposition_summary",
    "plot_emission_summary",
    "plot_size_distribution_summary",
    # Temperature-axis boxplots (re-exported from plot_particle_boxplots)
    "plot_emission_boxplot",
    "plot_deposition_rate_boxplot",
    "plot_emission_rate_boxplot",
    "plot_penetration_factor_boxplot",
    "plot_emission_etotal_by_metric_boxplot",
    "plot_emission_etotal_by_showerhead_boxplot",
    # Condition-comparison boxplots (re-exported from plot_comparison)
    "plot_spray_pattern_comparison_boxplots",
    "plot_shower_head_comparison_boxplots",
    "plot_mannequin_comparison_boxplots",
    "plot_door_comparison_boxplots",
    "plot_fan_comparison_boxplots",
]
