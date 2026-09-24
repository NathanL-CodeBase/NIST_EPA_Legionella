#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle Inhalation Dose Boxplots
====================================

Reproduces the emission_etotal boxplot suite from
scripts/particle_emission_variant_figures.py with cumulative inhaled dose
(src/inhalation_dose.py) as the plotted metric in place of E_total. It does
not recompute the dose; run scripts/particle_inhalation_dose_analysis.py
first to generate output/inhaled_dose_summary.xlsx.

Wherever the original figures split by air-change-rate source (outside,
entry), these figures split by the two independent inhaled-dose
concentration sources instead: quantaq (raw MOD-PM-00195, every event) and
croom (fleet-averaged room concentration, valid only 2026-06-03 through
2026-07-16 -- boxes for water-temperature groups outside that window are
empty or sparse for the croom series). The two dose sources are never
blended.

The ACR and avg-beta metric-axis variants have no natural dose-source
counterpart on their x-axis (ACR/beta come from the outside/entry
air-change-rate box model, not the dose calculation), so both axes are
crossed: each x-axis air-change-rate source (outside, entry) is paired with
each y-axis dose source (quantaq, croom) independently, giving four files
per bin group instead of two.

Data sources (joined on event_number):
    - output/inhaled_dose_summary.xlsx (inhaled_dose_by_bin sheet) --
      bin{n}_{quantaq,croom}_inhaled_dose, config_key, water_temp, shower_on.
    - output/particle_analysis_summary.xlsx (all_results sheet) -- flow_rate
      (for the same FLOW_RATE_MIN/MAX filter as
      particle_emission_variant_figures.py), lambda_outside, lambda_entry,
      and per-bin beta_other/p_mean (averaged here into avg_beta_outside,
      avg_beta_entry, avg_p) for the metric-axis figures.
    - Bedroom_Conditions sheet of rh_temp_wind_summary.xlsx -- bedroom_rh,
      bedroom_temp for the metric-axis figures, and RH for the n=/RH=
      boxplot annotations (best-effort; see
      particle_emission_variant_figures.py).

Output Files:
    output/plots/inhaled_dose_boxplot_{bin0-2,bin3-6,bin7-11}_{quantaq,croom}.png
    output/plots/inhaled_dose_by_{bedroom_rh,bedroom_temp}_boxplot_
        {bin0-2,bin3-6,bin7-11}_{quantaq,croom}.png
    output/plots/inhaled_dose_by_acr_{outside,entry}_boxplot_
        {bin0-2,bin3-6,bin7-11}_{quantaq,croom}.png
    output/plots/inhaled_dose_by_beta_{outside,entry}_boxplot_
        {bin0-2,bin3-6,bin7-11}_{quantaq,croom}.png
    output/plots/inhaled_dose_by_showerhead_boxplot_{bin0-2,bin3-6,bin7-11}_
        {quantaq,croom}.png
    output/plots/{spray_pattern,head_type,mannequin,bath_door_position,
        bedroom_door_position,fan_status}_inhaled_dose_boxplot_
        {bin0-2,bin3-6,bin7-11}_{quantaq,croom}.png

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-24
Update log:
    2026-09-24  Initial version.
"""

import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Add project root to path for src/ imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_paths import get_common_file, get_data_root  # noqa: E402
from src.particle_calculations import PARTICLE_BINS  # noqa: E402
from src.plot_comparison import (  # noqa: E402
    BATH_DOOR_GROUP_DEFS,
    BATH_DOOR_TEMP_FILTER,
    BATH_DOOR_TITLE,
    BATH_DOOR_XLABEL,
    BEDROOM_DOOR_GROUP_DEFS,
    BEDROOM_DOOR_TEMP_FILTER,
    BEDROOM_DOOR_TITLE,
    BEDROOM_DOOR_XLABEL,
    FAN_GROUP_DEFS,
    FAN_TEMP_FILTER,
    FAN_TITLE,
    FAN_XLABEL,
    HEAD_TYPE_GROUP_DEFS,
    HEAD_TYPE_TEMP_FILTER,
    HEAD_TYPE_TITLE,
    HEAD_TYPE_XLABEL,
    MANNEQUIN_GROUP_DEFS,
    MANNEQUIN_TEMP_FILTER,
    MANNEQUIN_TITLE,
    MANNEQUIN_XLABEL,
    SPRAY_PATTERN_GROUP_DEFS,
    SPRAY_PATTERN_TEMP_FILTER,
    SPRAY_PATTERN_TITLE,
    SPRAY_PATTERN_XLABEL,
    plot_condition_dose_comparison_boxplot,
)
from src.plot_particle_boxplots import (  # noqa: E402
    plot_emission_etotal_by_metric_boxplot,
    plot_emission_etotal_by_showerhead_boxplot,
    plot_inhaled_dose_boxplot,
)

# =============================================================================
# Configuration
# =============================================================================

DOSE_XLSX_NAME = "inhaled_dose_summary.xlsx"
DOSE_SHEET = "inhaled_dose_by_bin"

DECAY_XLSX_NAME = "particle_analysis_summary.xlsx"
DECAY_SHEET = "all_results"

# Same water-temperature-sweep flow-rate filter as
# particle_emission_variant_figures.py, applied to these boxplot/comparison
# figures only (not to inhaled_dose_summary.xlsx itself).
FLOW_RATE_MIN: float = 4.1  # LPM (inclusive)
FLOW_RATE_MAX: float = 5.6  # LPM (inclusive)

DOSE_VALUE_LABEL = "Cumulative Inhaled Dose (#)"
DOSE_METRIC_TITLE = "Cumulative Inhaled Dose"
DOSE_COL_TEMPLATE = "bin{n}_{source}_inhaled_dose"
DOSE_SOURCES = ("quantaq", "croom")
ACH_SOURCES = ("outside", "entry")

_UNIT_SUFFIX_RE = re.compile(r" \([^)]*\)$")


def _strip_unit_suffixes(df: pd.DataFrame) -> pd.DataFrame:
    """Undo the Excel unit-suffix rename (e.g. "bin3_quantaq_inhaled_dose (#)"
    -> "bin3_quantaq_inhaled_dose"), matching
    particle_emission_variant_figures.py's helper of the same purpose.
    """
    return df.rename(columns={c: _UNIT_SUFFIX_RE.sub("", c) for c in df.columns})


# =============================================================================
# Data Loading
# =============================================================================


def _load_dose_summary(output_dir: Path) -> pd.DataFrame:
    """Load inhaled_dose_summary.xlsx and strip Excel unit suffixes."""
    path = output_dir / DOSE_XLSX_NAME
    if not path.exists():
        raise FileNotFoundError(
            f"Inhaled dose summary not found: {path}\n"
            "Run scripts/particle_inhalation_dose_analysis.py first to generate it."
        )
    df = pd.read_excel(path, sheet_name=DOSE_SHEET, engine="openpyxl")
    df["shower_on"] = pd.to_datetime(df["shower_on"], errors="coerce")
    return _strip_unit_suffixes(df)


def _load_decay_metrics(output_dir: Path) -> pd.DataFrame:
    """
    Load flow_rate and the ACR/beta/p metrics needed for the metric-axis
    figures from particle_analysis_summary.xlsx (scripts/particle_decay_analysis.py).

    Returns
    -------
    pd.DataFrame
        One row per event: event_number, flow_rate, lambda_outside,
        lambda_entry, avg_beta_outside, avg_beta_entry, avg_p. avg_beta_* and
        avg_p are the mean of the per-bin beta_other / p_mean columns
        (matching particle_emission_variant_figures.py's own computation).
    """
    path = output_dir / DECAY_XLSX_NAME
    if not path.exists():
        raise FileNotFoundError(
            f"Particle analysis summary not found: {path}\n"
            "Run scripts/particle_decay_analysis.py first to generate it."
        )
    df = pd.read_excel(path, sheet_name=DECAY_SHEET, engine="openpyxl")
    df = _strip_unit_suffixes(df)

    bin_nums = list(PARTICLE_BINS.keys())
    for source in ACH_SOURCES:
        cols = [f"bin{b}_{source}_beta_other" for b in bin_nums if f"bin{b}_{source}_beta_other" in df.columns]
        df[f"avg_beta_{source}"] = df[cols].mean(axis=1) if cols else np.nan
    p_cols = [f"bin{b}_p_mean" for b in bin_nums if f"bin{b}_p_mean" in df.columns]
    df["avg_p"] = df[p_cols].mean(axis=1) if p_cols else np.nan

    keep_cols = [
        "event_number",
        "flow_rate",
        "lambda_outside",
        "lambda_entry",
        "avg_beta_outside",
        "avg_beta_entry",
        "avg_p",
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    return df[keep_cols].copy()


def _load_bedroom_conditions() -> "tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]":
    """
    Load bedroom RH/temperature for the metric-axis x-axis and the RH data
    used for n=/RH= boxplot annotations (same source, two shapes).

    Returns
    -------
    tuple of (metrics_df, rh_data)
        metrics_df: event_number, bedroom_rh, bedroom_temp (or None on failure).
        rh_data: 'datetime', 'RH_bedroom' for _get_rh_at_shower_on (or None).
    """
    try:
        path = get_common_file("rh_temp_wind_summary")
        bc = pd.read_excel(path, sheet_name="Bedroom_Conditions", engine="openpyxl")
        metrics_df = bc[["event_number", "rh_mean (%)", "temp_mean (degC)"]].rename(
            columns={"rh_mean (%)": "bedroom_rh", "temp_mean (degC)": "bedroom_temp"}
        )
        rh_data = bc[["shower_on", "rh_mean (%)"]].rename(
            columns={"shower_on": "datetime", "rh_mean (%)": "RH_bedroom"}
        )
        rh_data["datetime"] = pd.to_datetime(rh_data["datetime"])
        print("  Loaded Bedroom_Conditions RH/temperature.")
        return metrics_df, rh_data
    except Exception as e:
        print(f"  Note: Could not load Bedroom_Conditions (metric-axis x-axis and RH annotations skipped): {e}")
        return None, None


def build_dose_plot_df(output_dir: Path) -> "tuple[pd.DataFrame, Optional[pd.DataFrame]]":
    """
    Assemble the merged DataFrame used by every figure in this script.

    Parameters
    ----------
    output_dir : Path
        Directory holding inhaled_dose_summary.xlsx and particle_analysis_summary.xlsx.

    Returns
    -------
    tuple of (plot_df, rh_data)
        plot_df: inhaled dose columns plus flow_rate, lambda_outside/entry,
            avg_beta_outside/entry, avg_p, bedroom_rh, bedroom_temp -- filtered
            to FLOW_RATE_MIN/MAX (rows with no flow_rate reading are kept).
        rh_data: Bedroom_Conditions RH DataFrame for boxplot annotations, or
            None if unavailable.
    """
    print("Loading inhaled dose summary...")
    dose_df = _load_dose_summary(output_dir)
    print(f"  {len(dose_df)} event row(s)")

    print("Loading flow rate and ACR/beta/p metrics...")
    decay_df = _load_decay_metrics(output_dir)
    dose_df = dose_df.merge(decay_df, on="event_number", how="left")

    print("Loading bedroom RH/temperature...")
    bedroom_df, rh_data = _load_bedroom_conditions()
    if bedroom_df is not None:
        dose_df = dose_df.merge(bedroom_df, on="event_number", how="left")

    if "flow_rate" in dose_df.columns and dose_df["flow_rate"].notna().any():
        n_before = len(dose_df)
        dose_df = dose_df[
            dose_df["flow_rate"].between(FLOW_RATE_MIN, FLOW_RATE_MAX) | dose_df["flow_rate"].isna()
        ].copy()
        print(f"  Flow rate filter ({FLOW_RATE_MIN}-{FLOW_RATE_MAX} LPM): {n_before} -> {len(dose_df)} rows")

    return dose_df, rh_data


# =============================================================================
# Figure Generation
# =============================================================================


def generate_dose_boxplots(plot_df: pd.DataFrame, rh_data: "Optional[pd.DataFrame]", plot_dir: Path) -> None:
    """Generate the full inhaled-dose boxplot suite into *plot_dir*."""
    plot_dir.mkdir(parents=True, exist_ok=True)

    # ── Fixed water-temperature axis (3 bin groups x 2 dose sources) ────────
    try:
        plot_inhaled_dose_boxplot(
            plot_df,
            PARTICLE_BINS,
            plot_dir / "inhaled_dose_boxplot.png",
            rh_data=rh_data,
            x_range=(5, 55, 5),
        )
        print("  Generated: inhaled_dose_boxplot_{bin0-2,bin3-6,bin7-11}_{quantaq,croom}.png")
    except Exception as e:
        print(f"  Error generating inhaled_dose_boxplot: {e}")

    # ── Metric axis: bedroom RH, bedroom temperature (dose-source-only) ─────
    _metric_axes_shared = [
        ("bedroom_rh", "Bedroom RH (%)", "inhaled_dose_by_bedroom_rh_boxplot.png", (23, 43, 1)),
        ("bedroom_temp", "Bedroom Temperature (°C)", "inhaled_dose_by_bedroom_temp_boxplot.png", (14.9, 18.2, 0.1)),
    ]
    for metric_col, metric_label, filename, x_range in _metric_axes_shared:
        for dose_source in DOSE_SOURCES:
            try:
                plot_emission_etotal_by_metric_boxplot(
                    plot_df,
                    PARTICLE_BINS,
                    plot_dir / filename,
                    metric_col=metric_col,
                    metric_label=metric_label,
                    source=dose_source,
                    rh_data=rh_data,
                    x_range=x_range,
                    value_col_template=DOSE_COL_TEMPLATE,
                    value_label=DOSE_VALUE_LABEL,
                    metric_title=DOSE_METRIC_TITLE,
                )
                print(f"  Generated: {filename} ({dose_source} dose source)")
            except Exception as e:
                print(f"  Error generating {filename} ({dose_source} dose source): {e}")

    # ── Metric axis: ACR, avg beta -- crossed with both dose sources, since ─
    # dose has no air-change-rate source of its own (see module docstring).
    _metric_axes_per_ach_source = [
        ("lambda_{source}", "Air Change Rate λ (h⁻¹)", "inhaled_dose_by_acr_{source}_boxplot.png", (0.75, 1.65, 0.05)),
        ("avg_beta_{source}", "Avg. Other Process Rate β (h⁻¹)", "inhaled_dose_by_beta_{source}_boxplot.png", (-0.35, 0.35, 0.05)),
    ]
    for metric_col_template, metric_label, filename_template, x_range in _metric_axes_per_ach_source:
        for ach_source in ACH_SOURCES:
            metric_col = metric_col_template.format(source=ach_source)
            filename = filename_template.format(source=ach_source)
            for dose_source in DOSE_SOURCES:
                try:
                    plot_emission_etotal_by_metric_boxplot(
                        plot_df,
                        PARTICLE_BINS,
                        plot_dir / filename,
                        metric_col=metric_col,
                        metric_label=f"{metric_label} ({ach_source} source)",
                        source=dose_source,
                        rh_data=rh_data,
                        x_range=x_range,
                        value_col_template=DOSE_COL_TEMPLATE,
                        value_label=DOSE_VALUE_LABEL,
                        metric_title=DOSE_METRIC_TITLE,
                    )
                    print(f"  Generated: {filename} ({ach_source} ACH source, {dose_source} dose source)")
                except Exception as e:
                    print(f"  Error generating {filename} ({ach_source}/{dose_source}): {e}")

    # ── Shower head type comparison ──────────────────────────────────────────
    for dose_source in DOSE_SOURCES:
        try:
            plot_emission_etotal_by_showerhead_boxplot(
                plot_df,
                PARTICLE_BINS,
                plot_dir / "inhaled_dose_by_showerhead_boxplot.png",
                source=dose_source,
                rh_data=rh_data,
                value_col_template=DOSE_COL_TEMPLATE,
                value_label=DOSE_VALUE_LABEL,
                metric_title=DOSE_METRIC_TITLE,
            )
            print(f"  Generated: inhaled_dose_by_showerhead_boxplot.png ({dose_source} dose source)")
        except Exception as e:
            print(f"  Error generating inhaled_dose_by_showerhead_boxplot ({dose_source}): {e}")

    # ── Condition-comparison families (inhaled-dose panel only) ─────────────
    _dose_comparison_families = [
        ("spray_pattern", SPRAY_PATTERN_GROUP_DEFS, SPRAY_PATTERN_TITLE, SPRAY_PATTERN_XLABEL, SPRAY_PATTERN_TEMP_FILTER),
        ("head_type", HEAD_TYPE_GROUP_DEFS, HEAD_TYPE_TITLE, HEAD_TYPE_XLABEL, HEAD_TYPE_TEMP_FILTER),
        ("mannequin", MANNEQUIN_GROUP_DEFS, MANNEQUIN_TITLE, MANNEQUIN_XLABEL, MANNEQUIN_TEMP_FILTER),
        ("bath_door_position", BATH_DOOR_GROUP_DEFS, BATH_DOOR_TITLE, BATH_DOOR_XLABEL, BATH_DOOR_TEMP_FILTER),
        ("bedroom_door_position", BEDROOM_DOOR_GROUP_DEFS, BEDROOM_DOOR_TITLE, BEDROOM_DOOR_XLABEL, BEDROOM_DOOR_TEMP_FILTER),
        ("fan_status", FAN_GROUP_DEFS, FAN_TITLE, FAN_XLABEL, FAN_TEMP_FILTER),
    ]
    for stem_prefix, group_defs, title_base, x_label, temp_filter in _dose_comparison_families:
        try:
            plot_condition_dose_comparison_boxplot(
                plot_df,
                PARTICLE_BINS,
                plot_dir,
                stem_prefix,
                group_defs,
                title_base,
                x_label,
                rh_data=rh_data,
                temp_filter=temp_filter,
            )
            print(f"  Generated: {stem_prefix}_inhaled_dose_boxplot_{{bin0-2,bin3-6,bin7-11}}_{{quantaq,croom}}.png")
        except Exception as e:
            print(f"  Error generating {stem_prefix} inhaled-dose comparison figures: {e}")

    print(f"  Plots saved to: {plot_dir}")


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Inhaled-dose boxplots: reproduces the emission_etotal boxplot suite "
        "from particle_emission_variant_figures.py with cumulative inhaled dose "
        "(quantaq/croom sources) in place of E_total."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Directory holding inhaled_dose_summary.xlsx and particle_analysis_summary.xlsx, "
        "and where figures are written under plots/ (default: data_root/output).",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else get_data_root() / "output"

    print("\n" + "=" * 70)
    print("Particle Inhalation Dose Boxplots")
    print("=" * 70)

    plot_df, rh_data = build_dose_plot_df(output_dir)

    plot_dir = output_dir / "plots"
    print("\nGenerating inhaled-dose boxplot suite...")
    generate_dose_boxplots(plot_df, rh_data, plot_dir)

    print("\n" + "=" * 70)
    print("Done")
    print("=" * 70)


if __name__ == "__main__":
    main()
