#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle Analysis — Categorical Condition-Comparison Boxplot Functions
=======================================================================

Five families of boxplots comparing primary experimental conditions (spray pattern,
shower head type, mannequin presence, door position, bath fan status). Each family
produces fifteen figures: three metrics (E_total, penetration factor, other process
rate) × three particle size bin groups (Bin 0–2, 3–6, 7–11), with E_total and other
process rate additionally split by air-change-rate source (outside, entry) since
both are lambda-dependent; penetration factor is not and stays a single file per
bin group.

All comparison figures filter to a ±2 °C temperature window (default 36–42 °C)
so that temperature is not a confounding variable.

Functions:
    - plot_spray_pattern_comparison_boxplots: Wide / Narrow / Mid Pepco spray patterns
    - plot_shower_head_comparison_boxplots: Standard, Pepco variants, FilterWand, Used
    - plot_mannequin_comparison_boxplots: with vs. without mannequin
    - plot_door_comparison_boxplots: door open vs. door closed (fan on, Pepco Wide)
    - plot_fan_comparison_boxplots: fan off vs. fan on (door open, Pepco Wide)
    - plot_condition_dose_comparison_boxplot: generic inhaled-dose-only panel for
      any of the above families, reusing their group_defs/title/x_label/temp_filter
      constants (SPRAY_PATTERN_*, HEAD_TYPE_*, MANNEQUIN_*, BATH_DOOR_*,
      BEDROOM_DOOR_*, FAN_*)

Output naming:
    - {stem}_{metric}_boxplot_{bin0-2,bin3-6,bin7-11}.png (penetration_factor)
    - {stem}_{metric}_boxplot_{bin0-2,bin3-6,bin7-11}_{outside,entry}.png
      (emission_etotal, other_process_rate)
    - {stem}_inhaled_dose_boxplot_{bin0-2,bin3-6,bin7-11}_{quantaq,croom}.png
      (plot_condition_dose_comparison_boxplot)

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Date: 2026
"""

from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

import src.sig_figs as sf
from src.plot_particle_boxplots import (
    _TEMP_BOXPLOT_CONFIG,
    _extract_config_temp,
    _get_rh_at_shower_on,
    _write_boxplot_companion_md,
)
from src.plot_style import (
    BOXPLOT_CONFIG,
    FONT_SIZE_ANNOTATION,
    FONT_SIZE_LABEL,
    FONT_SIZE_LEGEND,
    FONT_SIZE_TICK,
    FONT_SIZE_TITLE,
    SENSOR_COLORS,
    TITLE_FONTWEIGHT,
    apply_style,
    create_figure,
    save_figure,
)

# GroupDef: (group_key, tick_label, filter_fn)
#   group_key  — unique identifier used as dict key
#   tick_label — text shown on x-axis (may contain \n for wrapping)
#   filter_fn  — callable(config_key: str) -> bool; True means the row belongs here
_GroupDef = Tuple[str, str, Callable[[str], bool]]


def _draw_categorical_comparison_boxplot(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    output_path: Path,
    cfg: dict,
    group_defs: "List[_GroupDef]",
    title_base: str,
    x_label: str,
    rh_data: "Optional[pd.DataFrame]" = None,
    temp_filter: "Optional[Tuple[float, float]]" = None,
) -> None:
    """
    Shared implementation for all five condition-comparison boxplot families.

    Draws three output files (bin0-2, bin3-6, bin7-11) with a categorical
    x-axis whose groups are defined by *group_defs*.  Each group may pool
    multiple config_key values selected by an arbitrary filter function.

    Parameters:
        results_df: Full particle analysis results DataFrame.
        particle_bins: Dict of bin metadata.
        output_path: Base path; ``_bin0-2`` / ``_bin3-6`` / ``_bin7-11``
            are appended to the stem.
        cfg: One of the ``_TEMP_BOXPLOT_CONFIG`` entries (supplies col_template,
            ylabel, title_note, hline).
        group_defs: Ordered list of (group_key, tick_label, filter_fn) tuples.
        title_base: Figure title prefix (e.g. "Spray Pattern Effect on Emission").
        x_label: x-axis label string.
        rh_data: Optional Bedroom_Conditions RH DataFrame for n=/RH= annotations.
        temp_filter: Optional (min_temp, max_temp) in °C; rows whose config_key
            encodes a temperature outside this range are excluded before grouping.
    """
    apply_style()

    if results_df.empty or "config_key" not in results_df.columns:
        return

    # Optional temperature pre-filter (±2 °C rule)
    if temp_filter is not None:
        t_lo, t_hi = temp_filter

        def _in_range(ck: object) -> bool:
            if not isinstance(ck, str) or not ck:
                return False
            t = _extract_config_temp(ck)
            return t is not None and t_lo <= t <= t_hi

        work_df = results_df[results_df["config_key"].apply(_in_range)].copy()
    else:
        work_df = results_df.copy()

    if work_df.empty:
        return

    # Assign each row to the first matching group; NaN/non-string keys are skipped.
    def _assign(ck: object) -> "Optional[str]":
        if not isinstance(ck, str) or not ck:
            return None
        for gk, _, fn in group_defs:
            if fn(ck):
                return gk
        return None

    work_df["_cmp_group"] = work_df["config_key"].apply(_assign)
    work_df = work_df[work_df["_cmp_group"].notna()].copy()
    if work_df.empty:
        return

    x_pos = {gk: i for i, (gk, _, _) in enumerate(group_defs)}
    tick_label_map = {gk: lbl for gk, lbl, _ in group_defs}
    n_pos = len(group_defs)

    all_bin_nums = list(particle_bins.keys())
    # Fixed-width boxes scaled to unit categorical spacing
    all_bin_widths = np.linspace(0.10, 0.38, len(all_bin_nums))

    bin_groups = [
        ([b for b in all_bin_nums if b <= 2], "bin0-2"),
        ([b for b in all_bin_nums if 3 <= b <= 6], "bin3-6"),
        ([b for b in all_bin_nums if b >= 7], "bin7-11"),
    ]

    # Build companion .md once (event membership is identical across bin groups)
    md_groups = []
    for gk, tick_lbl, _ in group_defs:
        gdf = work_df[work_df["_cmp_group"] == gk]
        header = f"{tick_lbl.replace(chr(10), ' ')} | n={len(gdf)}"
        events_list = []
        for _, row in gdf.sort_values("event_number").iterrows():
            events_list.append(
                (
                    row.get("event_number"),
                    row.get("test_name", ""),
                    row.get("config_key", gk),
                )
            )
        md_groups.append({"header": header, "events": events_list})
    _write_boxplot_companion_md(output_path, title_base, md_groups)

    # Source-dependent metrics (has_source=True) are drawn once per entry in
    # cfg["sources"] (default ("outside", "entry"), the air-change-rate pair;
    # e.g. inhaled_dose uses ("quantaq", "croom") instead), bounding the
    # result; source-independent metrics (e.g. penetration factor) are drawn
    # once with no source suffix.
    if cfg.get("has_source", True):
        sources = cfg.get("sources", ("outside", "entry"))
    else:
        sources = (None,)

    for source in sources:
        source_suffix = f"_{source}" if source else ""
        source_note = f" ({source} source)" if source else ""

        for group_bins, group_label in bin_groups:
            if not group_bins:
                continue

            value_cols = [cfg["col_template"].format(n=b, source=source) for b in group_bins]

            # Per-group annotation stats (keyed by categorical x position)
            annot_stats: dict = {}
            for gk, tick_lbl, _ in group_defs:
                gdf = work_df[work_df["_cmp_group"] == gk]
                if gdf.empty:
                    continue
                all_vals: list = []
                for col in value_cols:
                    if col in gdf.columns:
                        all_vals.extend(gdf[col].dropna().tolist())
                annot_stats[x_pos[gk]] = {
                    "n": len(gdf),
                    "max_val": float(np.max(all_vals)) if all_vals else np.nan,
                    "shower_on": gdf["shower_on"].copy()
                    if "shower_on" in gdf.columns
                    else pd.Series([], dtype="object"),
                    "tick_label": tick_lbl,
                }

            fig, ax = create_figure(figsize=BOXPLOT_CONFIG["figsize"])
            if isinstance(ax, list):
                ax = ax[0]

            for bin_num in group_bins:
                global_idx = all_bin_nums.index(bin_num)
                color = SENSOR_COLORS[global_idx % len(SENSOR_COLORS)]
                col = cfg["col_template"].format(n=bin_num, source=source)
                if col not in work_df.columns:
                    continue

                box_width = float(all_bin_widths[global_idx])
                positions: list = []
                data: list = []
                for gk, _, _ in group_defs:
                    gdf = work_df[work_df["_cmp_group"] == gk]
                    vals = gdf[col].dropna().values if not gdf.empty else np.array([])
                    if len(vals) > 0:
                        positions.append(x_pos[gk])
                        data.append(vals)

                if not data:
                    continue

                bp = ax.boxplot(
                    data,
                    positions=positions,
                    widths=box_width,
                    patch_artist=True,
                    showfliers=True,
                    flierprops=dict(
                        marker=BOXPLOT_CONFIG["flier_marker"],
                        markersize=BOXPLOT_CONFIG["flier_markersize"],
                        alpha=BOXPLOT_CONFIG["flier_alpha"],
                        color=color,
                    ),
                )
                for patch in bp["boxes"]:
                    patch.set_facecolor(color)
                    patch.set_alpha(BOXPLOT_CONFIG["box_alpha"])
                for element in ("whiskers", "caps"):
                    for line in bp[element]:
                        line.set_color(color)
                        line.set_alpha(BOXPLOT_CONFIG["box_alpha"])
                for med in bp["medians"]:
                    med.set_color(BOXPLOT_CONFIG["median_color"])
                    med.set_linewidth(BOXPLOT_CONFIG["median_linewidth"])

            # Annotations above each group
            if annot_stats:
                y_min, y_max = ax.get_ylim()
                y_rng = y_max - y_min
                ax.set_ylim(y_min, y_max + 0.25 * y_rng)
                y_min, y_max = ax.get_ylim()
                offset = 0.02 * (y_max - y_min)
                for xp, stats in sorted(annot_stats.items()):
                    max_val = stats.get("max_val", np.nan)
                    if np.isnan(max_val):
                        continue
                    lines = [stats["tick_label"], f"n={stats['n']}"]
                    if rh_data is not None:
                        avg_rh = _get_rh_at_shower_on(stats["shower_on"], rh_data)
                        if not np.isnan(avg_rh):
                            lines.append(f"RH={sf.fmt_fig(avg_rh, fallback='.0f')}%")
                    ax.text(
                        xp,
                        max_val + offset,
                        "\n".join(lines),
                        ha="center",
                        va="bottom",
                        fontsize=FONT_SIZE_ANNOTATION,
                        color="black",
                    )

            all_xp = [x_pos[gk] for gk, _, _ in group_defs]
            all_tl = [tick_label_map[gk] for gk, _, _ in group_defs]
            ax.set_xticks(all_xp)
            ax.set_xticklabels(all_tl, fontsize=FONT_SIZE_TICK - 1, ha="center")
            ax.set_xlim(-0.6, n_pos - 0.4)

            ax.set_xlabel(x_label, fontsize=FONT_SIZE_LABEL)
            ax.set_ylabel(cfg["ylabel"], fontsize=FONT_SIZE_LABEL)
            bin_lbl = group_label.replace("-", "–").replace("bin", "Bin ")
            ax.set_title(
                f"{title_base} — {bin_lbl}{source_note}\n{cfg['title_note']}",
                fontsize=FONT_SIZE_TITLE,
                fontweight=TITLE_FONTWEIGHT,
            )
            if cfg.get("hline") is not None:
                ax.axhline(cfg["hline"], color="gray", linewidth=0.8, linestyle=":", alpha=0.6)
            ax.grid(True, alpha=0.3, axis="y")
            ax.tick_params(labelsize=FONT_SIZE_TICK)

            legend_elements = [
                Patch(
                    facecolor=SENSOR_COLORS[all_bin_nums.index(b) % len(SENSOR_COLORS)],
                    alpha=0.7,
                    label=f"Bin {b} ({particle_bins[b]['name']} µm)",
                )
                for b in group_bins
            ]
            ax.legend(
                handles=legend_elements,
                loc="upper right",
                fontsize=FONT_SIZE_LEGEND - 1,
                ncol=1,
            )

            out = (
                output_path.parent
                / f"{output_path.stem}_{group_label}{source_suffix}{output_path.suffix}"
            )
            plt.tight_layout()
            save_figure(fig, out)
            plt.close(fig)


# ---------------------------------------------------------------------------
# Metric keys used by all five comparison families:
#   emission_etotal    → bin{n}_E_total
#   penetration_factor → bin{n}_p_mean
#   deposition_rate    → bin{n}_beta_other_raw_mean
# ---------------------------------------------------------------------------
_COMPARISON_METRICS: "List[Tuple[str, str]]" = [
    ("emission_etotal", "emission_etotal"),
    ("penetration_factor", "penetration_factor"),
    ("deposition_rate", "other_process_rate"),
]


def _run_comparison_family(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    plot_dir: Path,
    stem_prefix: str,
    group_defs: "List[_GroupDef]",
    title_base: str,
    x_label: str,
    rh_data: "Optional[pd.DataFrame]" = None,
    temp_filter: "Optional[Tuple[float, float]]" = None,
) -> None:
    """Loop over the three comparison metrics and call the shared drawing function."""
    for cfg_key, file_metric in _COMPARISON_METRICS:
        cfg = _TEMP_BOXPLOT_CONFIG[cfg_key]
        base_path = plot_dir / f"{stem_prefix}_{file_metric}_boxplot.png"
        _draw_categorical_comparison_boxplot(
            results_df,
            particle_bins,
            base_path,
            cfg,
            group_defs,
            title_base=f"{title_base} — {cfg['title_metric']}",
            x_label=x_label,
            rh_data=rh_data,
            temp_filter=temp_filter,
        )


# ---------------------------------------------------------------------------
# 1. Spray Pattern comparison  (W36–42 °C, Pepco head, no mannequin, bath door open, fan off)
# ---------------------------------------------------------------------------

# Group definitions, title, x-label, and temperature filter are module-level
# constants (rather than local to plot_spray_pattern_comparison_boxplots) so
# plot_condition_dose_comparison_boxplot can reuse the same condition filters
# for the inhaled-dose panel without duplicating them.
SPRAY_PATTERN_GROUP_DEFS: "List[_GroupDef]" = [
    (
        "Wide",
        "Pepco\nWide",
        lambda ck: "_Pepco_Wide_BathDoorOpen_" in ck and "_FanOff" in ck and "_Mannequin" not in ck,
    ),
    (
        "Narrow",
        "Pepco\nNarrow",
        lambda ck: (
            "_Pepco_Narrow_BathDoorOpen_" in ck and "_FanOff" in ck and "_Mannequin" not in ck
        ),
    ),
]
SPRAY_PATTERN_TITLE = "Spray Pattern Effect"
SPRAY_PATTERN_XLABEL = "Spray Pattern"
SPRAY_PATTERN_TEMP_FILTER = (36.0, 42.0)


def plot_spray_pattern_comparison_boxplots(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    plot_dir: Path,
    rh_data: "Optional[pd.DataFrame]" = None,
) -> None:
    """Boxplots comparing Wide / Narrow Pepco spray patterns.

    Fixed conditions: Pepco head, W36–42 °C, no mannequin, bath door open, fan off.
    Produces fifteen figures: 3 metrics × 3 bin groups, with the two
    lambda-dependent metrics (E_total, other process rate) split outside/entry.

    Parameters:
        results_df: Full particle analysis results DataFrame.
        particle_bins: Dict of bin metadata.
        plot_dir: Directory for output files.
        rh_data: Optional RH DataFrame for annotations.
    """
    _run_comparison_family(
        results_df,
        particle_bins,
        plot_dir,
        stem_prefix="spray_pattern",
        group_defs=SPRAY_PATTERN_GROUP_DEFS,
        title_base=SPRAY_PATTERN_TITLE,
        x_label=SPRAY_PATTERN_XLABEL,
        rh_data=rh_data,
        temp_filter=SPRAY_PATTERN_TEMP_FILTER,
    )


# ---------------------------------------------------------------------------
# 2. Shower Head comparison  (W36–42 °C, bath door open, fan off)
# ---------------------------------------------------------------------------

HEAD_TYPE_GROUP_DEFS: "List[_GroupDef]" = [
    (
        "Standard",
        "Standard",
        lambda ck: (
            "_Pepco" not in ck
            and "_FilterWand" not in ck
            and "_Used" not in ck
            and "_Mannequin" not in ck
            and "_BathDoorOpen_" in ck
            and "_FanOff" in ck
        ),
    ),
    (
        "Pepco_Wide",
        "Pepco\nWide",
        lambda ck: "_Pepco_Wide_BathDoorOpen_" in ck and "_FanOff" in ck and "_Mannequin" not in ck,
    ),
    (
        "Pepco_Narrow",
        "Pepco\nNarrow",
        lambda ck: (
            "_Pepco_Narrow_BathDoorOpen_" in ck and "_FanOff" in ck and "_Mannequin" not in ck
        ),
    ),
    (
        "FilterWand",
        "Filter\nWand",
        lambda ck: "_FilterWand_" in ck and "_BathDoorOpen_" in ck and "_FanOff" in ck,
    ),
    (
        "Used_rainfall",
        "Used\nRainfall",
        lambda ck: "_Used_rainfall_BathDoorOpen_" in ck and "_FanOff" in ck,
    ),
    (
        "Used_12Nozzle",
        "Used\n12Nozzle",
        lambda ck: "_Used_12Nozzle_BathDoorOpen_" in ck and "_FanOff" in ck,
    ),
    (
        "Used_SingleWide",
        "Used\nSingleWide",
        lambda ck: (
            (
                "_Used_SingleWide_BathDoorOpen_" in ck
                or "_Used_SingleWide1_BathDoorOpen_" in ck
                or "_Used_SingleWide2_BathDoorOpen_" in ck
            )
            and "_FanOff" in ck
        ),
    ),
    (
        "Used_SingleNarrow",
        "Used\nSingleNarrow",
        lambda ck: "_Used_SingleNarrow_BathDoorOpen_" in ck and "_FanOff" in ck,
    ),
]
HEAD_TYPE_TITLE = "Shower Head Type Effect"
HEAD_TYPE_XLABEL = "Shower Head Configuration"
HEAD_TYPE_TEMP_FILTER = (36.0, 42.0)


def plot_shower_head_comparison_boxplots(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    plot_dir: Path,
    rh_data: "Optional[pd.DataFrame]" = None,
) -> None:
    """Boxplots comparing Standard, Pepco (Wide/Narrow), FilterWand, and Used heads.

    Fixed conditions: W36–42 °C, bath door open, fan off.  Standard-head events
    at W37 satisfy the ±2 °C rule relative to W38/W40.
    Produces fifteen figures: 3 metrics × 3 bin groups, with the two
    lambda-dependent metrics (E_total, other process rate) split outside/entry.

    Parameters:
        results_df: Full particle analysis results DataFrame.
        particle_bins: Dict of bin metadata.
        plot_dir: Directory for output files.
        rh_data: Optional RH DataFrame for annotations.
    """
    _run_comparison_family(
        results_df,
        particle_bins,
        plot_dir,
        stem_prefix="head_type",
        group_defs=HEAD_TYPE_GROUP_DEFS,
        title_base=HEAD_TYPE_TITLE,
        x_label=HEAD_TYPE_XLABEL,
        rh_data=rh_data,
        temp_filter=HEAD_TYPE_TEMP_FILTER,
    )


# ---------------------------------------------------------------------------
# 3. Mannequin comparison  (W36–42 °C, Pepco head, bath door open, fan off)
# ---------------------------------------------------------------------------

MANNEQUIN_GROUP_DEFS: "List[_GroupDef]" = [
    (
        "No_Mannequin_Narrow",
        "No Mannequin\nNarrow",
        lambda ck: (
            ("_Pepco_Narrow_" in ck or "_Used_SingleNarrow_" in ck)
            and "_Mannequin" not in ck
            and "_BathDoorOpen_" in ck
            and "_FanOff" in ck
        ),
    ),
    (
        "No_Mannequin_Wide",
        "No Mannequin\nWide",
        lambda ck: (
            "_Pepco_Wide_" in ck
            and "_Mannequin" not in ck
            and "_BathDoorOpen_" in ck
            and "_FanOff" in ck
        ),
    ),
    (
        "Mannequin_Narrow",
        "With Mannequin\nNarrow",
        lambda ck: (
            ("_Pepco_Narrow_" in ck or "_Used_SingleNarrow_" in ck)
            and "_Mannequin_BathDoorOpen_" in ck
            and "_FanOff" in ck
        ),
    ),
    (
        "Mannequin_Wide",
        "With Mannequin\nWide",
        lambda ck: "_Pepco_Wide_" in ck and "_Mannequin_BathDoorOpen_" in ck and "_FanOff" in ck,
    ),
]
MANNEQUIN_TITLE = "Mannequin Presence Effect by Head Geometry"
MANNEQUIN_XLABEL = "Mannequin / Head Geometry"
MANNEQUIN_TEMP_FILTER = (36.0, 42.0)


def plot_mannequin_comparison_boxplots(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    plot_dir: Path,
    rh_data: "Optional[pd.DataFrame]" = None,
) -> None:
    """Boxplots comparing with-mannequin vs. without-mannequin conditions.

    Fixed conditions: Pepco head, W36–42 °C, bath door open, fan off.
    Produces fifteen figures: 3 metrics × 3 bin groups, with the two
    lambda-dependent metrics (E_total, other process rate) split outside/entry.

    Parameters:
        results_df: Full particle analysis results DataFrame.
        particle_bins: Dict of bin metadata.
        plot_dir: Directory for output files.
        rh_data: Optional RH DataFrame for annotations.
    """
    _run_comparison_family(
        results_df,
        particle_bins,
        plot_dir,
        stem_prefix="mannequin",
        group_defs=MANNEQUIN_GROUP_DEFS,
        title_base=MANNEQUIN_TITLE,
        x_label=MANNEQUIN_XLABEL,
        rh_data=rh_data,
        temp_filter=MANNEQUIN_TEMP_FILTER,
    )


# ---------------------------------------------------------------------------
# 4. Door Position comparisons
#    4a. Bath door: Open vs. Closed  (W36–42 °C, Pepco Wide, fan on)
#    4b. Bedroom door: Closed vs. Open vs. Ajar  (W36–42 °C, Used SingleWide,
#        bath door open, fan off)
# ---------------------------------------------------------------------------

BATH_DOOR_GROUP_DEFS: "List[_GroupDef]" = [
    (
        "BathDoorOpen_FanOn",
        "Bath Door Open\n(Fan On)",
        lambda ck: "_Pepco_Wide_BathDoorOpen_" in ck and "_FanOn" in ck,
    ),
    (
        "BathDoorClosed_FanOn",
        "Bath Door Closed\n(Fan On)",
        lambda ck: "_Pepco_Wide_BathDoorClosed_" in ck and "_FanOn" in ck,
    ),
]
BATH_DOOR_TITLE = "Bath Door Position Effect"
BATH_DOOR_XLABEL = "Bath Door Position"
BATH_DOOR_TEMP_FILTER = (36.0, 42.0)

# Bedroom door — Pepco Narrow head, W36–42 °C, fan off, no mannequin.
# Note: the Ajar group (May 4–7) coincided with the bath door also being Ajar,
# so bath door position is not perfectly constant across the three groups.
BEDROOM_DOOR_GROUP_DEFS: "List[_GroupDef]" = [
    (
        "BdrmDoorClosed",
        "Bedroom Door\nClosed",
        lambda ck: (
            "_Pepco_Narrow_BathDoorOpen_BdrmDoorClosed_FanOff" in ck and "_Mannequin" not in ck
        ),
    ),
    (
        "BdrmDoorOpen",
        "Bedroom Door\nOpen",
        lambda ck: "_Pepco_Narrow_BathDoorOpen_BdrmDoorOpen_FanOff" in ck,
    ),
    (
        "BdrmDoorAjar",
        "Bedroom Door\nAjar",
        lambda ck: "_Pepco_Narrow_BathDoorAjar_BdrmDoorAjar_FanOff" in ck,
    ),
]
BEDROOM_DOOR_TITLE = "Bedroom Door Position Effect"
BEDROOM_DOOR_XLABEL = "Bedroom Door Position"
BEDROOM_DOOR_TEMP_FILTER = (36.0, 42.0)


def plot_door_comparison_boxplots(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    plot_dir: Path,
    rh_data: "Optional[pd.DataFrame]" = None,
) -> None:
    """Boxplots comparing bath door and bedroom door position conditions.

    4a — Bath door (bath_door_position): Open vs. Closed.
         Fixed conditions: Pepco Wide, W36–42 °C, fan on (door-closed period
         coincided with fan-on operation, so fan status is held constant).
         Output stem prefix: ``bath_door_position``.

    4b — Bedroom door (bedroom_door_position): Closed vs. Open vs. Ajar.
         Fixed conditions: Pepco Narrow head, W36–42 °C, fan off, no mannequin.
         Note: the Ajar group (May 4–7) had bath door Ajar; Closed and Open groups
         had bath door Open — bath door position is not held constant across groups.
         Output stem prefix: ``bedroom_door_position``.

    Each sub-comparison produces fifteen figures: 3 metrics × 3 bin groups,
    with the two lambda-dependent metrics (E_total, other process rate) split
    outside/entry.

    Parameters:
        results_df: Full particle analysis results DataFrame.
        particle_bins: Dict of bin metadata.
        plot_dir: Directory for output files.
        rh_data: Optional RH DataFrame for annotations.
    """
    _run_comparison_family(
        results_df,
        particle_bins,
        plot_dir,
        stem_prefix="bath_door_position",
        group_defs=BATH_DOOR_GROUP_DEFS,
        title_base=BATH_DOOR_TITLE,
        x_label=BATH_DOOR_XLABEL,
        rh_data=rh_data,
        temp_filter=BATH_DOOR_TEMP_FILTER,
    )
    _run_comparison_family(
        results_df,
        particle_bins,
        plot_dir,
        stem_prefix="bedroom_door_position",
        group_defs=BEDROOM_DOOR_GROUP_DEFS,
        title_base=BEDROOM_DOOR_TITLE,
        x_label=BEDROOM_DOOR_XLABEL,
        rh_data=rh_data,
        temp_filter=BEDROOM_DOOR_TEMP_FILTER,
    )


# ---------------------------------------------------------------------------
# 5. Fan Status comparison  (W36–42 °C, Pepco Wide, bath door open, no mannequin)
# ---------------------------------------------------------------------------

FAN_GROUP_DEFS: "List[_GroupDef]" = [
    (
        "FanOff",
        "Fan Off",
        lambda ck: "_Pepco_Wide_BathDoorOpen_" in ck and "_FanOff" in ck and "_Mannequin" not in ck,
    ),
    (
        "FanOn",
        "Fan On",
        lambda ck: "_Pepco_Wide_BathDoorOpen_" in ck and "_FanOn" in ck,
    ),
]
FAN_TITLE = "Bath Fan Effect"
FAN_XLABEL = "Bath Fan Status"
FAN_TEMP_FILTER = (36.0, 42.0)


def plot_fan_comparison_boxplots(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    plot_dir: Path,
    rh_data: "Optional[pd.DataFrame]" = None,
) -> None:
    """Boxplots comparing fan-off vs. fan-on conditions.

    Fixed conditions: Pepco Wide, W36–42 °C, bath door open, no mannequin.
    Produces fifteen figures: 3 metrics × 3 bin groups, with the two
    lambda-dependent metrics (E_total, other process rate) split outside/entry.

    Parameters:
        results_df: Full particle analysis results DataFrame.
        particle_bins: Dict of bin metadata.
        plot_dir: Directory for output files.
        rh_data: Optional RH DataFrame for annotations.
    """
    _run_comparison_family(
        results_df,
        particle_bins,
        plot_dir,
        stem_prefix="fan_status",
        group_defs=FAN_GROUP_DEFS,
        title_base=FAN_TITLE,
        x_label=FAN_XLABEL,
        rh_data=rh_data,
        temp_filter=FAN_TEMP_FILTER,
    )


# ---------------------------------------------------------------------------
# Inhaled-dose comparison panel (generic) — reuses the same group_defs/title/
# x_label/temp_filter constants as the corresponding emission_etotal family
# above, but draws only the inhaled_dose metric (quantaq/croom sources; see
# src/inhalation_dose.py) instead of looping over all three
# _COMPARISON_METRICS. Callers pass one family's constants per call.
# ---------------------------------------------------------------------------


def plot_condition_dose_comparison_boxplot(
    results_df: pd.DataFrame,
    particle_bins: Dict,
    plot_dir: Path,
    stem_prefix: str,
    group_defs: "List[_GroupDef]",
    title_base: str,
    x_label: str,
    rh_data: "Optional[pd.DataFrame]" = None,
    temp_filter: "Optional[Tuple[float, float]]" = None,
) -> None:
    """Inhaled-dose version of one condition-comparison family.

    Produces six figures: 3 bin groups × 2 dose sources (quantaq, croom).
    Unlike the five named ``plot_*_comparison_boxplots`` functions, this
    draws only the inhaled_dose metric, not the full emission_etotal /
    other_process_rate / penetration_factor trio.

    Parameters:
        results_df: Particle analysis results DataFrame merged with
            bin{n}_{quantaq,croom}_inhaled_dose columns (see
            scripts/particle_inhalation_dose_boxplots.py).
        particle_bins: Dict of bin metadata.
        plot_dir: Directory for output files.
        stem_prefix: Output filename stem (e.g. "spray_pattern").
        group_defs: Same (group_key, tick_label, filter_fn) list used by the
            corresponding emission_etotal family function.
        title_base: Figure title prefix (e.g. "Spray Pattern Effect").
        x_label: x-axis label string.
        rh_data: Optional Bedroom_Conditions RH DataFrame for annotations.
        temp_filter: Optional (min_temp, max_temp) in °C.
    """
    cfg = _TEMP_BOXPLOT_CONFIG["inhaled_dose"]
    base_path = plot_dir / f"{stem_prefix}_inhaled_dose_boxplot.png"
    _draw_categorical_comparison_boxplot(
        results_df,
        particle_bins,
        base_path,
        cfg,
        group_defs,
        title_base=f"{title_base} — {cfg['title_metric']}",
        x_label=x_label,
        rh_data=rh_data,
        temp_filter=temp_filter,
    )
