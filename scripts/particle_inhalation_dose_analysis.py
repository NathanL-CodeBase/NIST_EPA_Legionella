#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle Inhalation Dose Analysis
====================================

Computes and plots the cumulative excess particle count inhaled per shower
event, per particle-size bin, over a 2 h 10 min window starting at shower-on
(see src/inhalation_dose.py for the calculation). This is a standalone
exposure metric, independent of the penetration/deposition/emission-rate
pipeline in scripts/particle_decay_analysis.py -- it needs only the indoor
particle concentration series and the event registry.

Two concentration sources are computed and reported side by side (not
blended into one series):
    - QuantAQ-inside: raw MOD-PM-00195 reading, for every event.
    - C_room: the position-weighted MODULAIR-PM fleet average (see
      src.particle_room_correction.build_croom_data), which only exists
      2026-06-03 through 2026-07-16 -- events outside that window get NaN.
Requires the unified event registry (scripts/event_registry.py);
duration-excluded and otherwise-excluded events are skipped, matching
scripts/particle_decay_analysis.py's convention.

Output Files:
    output/inhaled_dose_summary.xlsx (sheet: inhaled_dose_by_bin) -- per-event
        cumulative inhaled dose (#), per bin, for both sources
        (bin{n}_quantaq_inhaled_dose, bin{n}_croom_inhaled_dose) plus
        n_minutes_covered_quantaq / n_minutes_covered_croom.
    output/plots/particle/inhaled_dose_bin{N}.html (12 figures) -- Bokeh
        time-series: x = shower_on, y = cumulative inhaled dose, one point per
        event per source (QuantAQ-inside vs. C_room). Style matches
        scripts/particle_emission_variant_figures.py
        (src.plot_style.style_moduair_figure: 1600x800, 12pt, no title,
        click-to-hide legend, hover enabled), with a vertical dashed line at
        ROOM_CUTOVER marking the start of C_room availability.

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-23
Update log:
    2026-09-23  Initial version (single blended concentration series, plus
        water-temperature boxplots).
    2026-09-23  Replaced the blended series with two independent sources
        (QuantAQ-inside for all events, C_room for the fleet window only)
        plotted side by side on each bin figure; removed the boxplots.
    2026-09-24  Changed C_room availability indicator from shaded span with label
        to a vertical dashed line at ROOM_CUTOVER, matching
        particle_emission_variant_figures.py style.
"""

import sys
import warnings
from pathlib import Path
from typing import Optional

import pandas as pd
from bokeh.models import ColumnDataSource, DatetimeTickFormatter, HoverTool, Span
from bokeh.plotting import figure, output_file, save

warnings.filterwarnings("ignore")

# Add project root to path for src/ imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import src.sig_figs as sf  # noqa: E402
from src.data_paths import get_data_root  # noqa: E402
from src.event_manager import is_event_excluded  # noqa: E402
from src.inhalation_dose import compute_cumulative_dose, resample_to_1min  # noqa: E402
from src.particle_calculations import PARTICLE_BINS  # noqa: E402
from src.particle_data_loader import get_events_from_registry, load_quantaq_data  # noqa: E402
from src.particle_room_correction import ROOM_CUTOVER, build_croom_data  # noqa: E402
from src.plot_style import SENSOR_COLORS, style_moduair_figure  # noqa: E402

# =============================================================================
# Configuration
# =============================================================================

SUMMARY_XLSX_NAME = "inhaled_dose_summary.xlsx"
SUMMARY_SHEET = "inhaled_dose_by_bin"

FIGURE_SUBDIR = ("plots", "particle")

# (label, legend text, color) for the two concentration sources plotted on
# every bin figure.
SOURCES = [
    ("quantaq", "QuantAQ-inside (MOD-PM-00195)", SENSOR_COLORS[0]),
    ("croom", "C_room (fleet average)", SENSOR_COLORS[1]),
]


# =============================================================================
# Event Loading
# =============================================================================


def _load_valid_events(output_dir: Path) -> list:
    """
    Load events from the unified event registry, skipping excluded events.

    Parameters
    ----------
    output_dir : Path
        Output directory containing the event registry file.

    Returns
    -------
    list of dict
        Event dicts with duration-excluded (no event_number) and otherwise
        excluded (is_event_excluded or registry is_excluded) events removed.
    """
    events, _co2_results, used_registry = get_events_from_registry(output_dir)
    if not used_registry:
        raise FileNotFoundError(
            "Event registry not found in "
            f"{output_dir}. Run scripts/event_registry.py first -- "
            "this analysis requires the unified registry for event numbering."
        )

    valid_events = []
    for event in events:
        is_excluded_flag, _reason = is_event_excluded(event["shower_on"])
        if not is_excluded_flag:
            is_excluded_flag = event.get("is_excluded", False)
        if is_excluded_flag:
            continue
        valid_events.append(event)

    print(f"  Loaded {len(events)} registry events; {len(valid_events)} valid (non-excluded)")
    return valid_events


# =============================================================================
# Dose Calculation
# =============================================================================


def _compute_dose(events: list) -> pd.DataFrame:
    """
    Compute cumulative inhaled dose for both concentration sources and merge
    them into one per-event table.

    Parameters
    ----------
    events : list of dict
        Valid (non-excluded) event dicts.

    Returns
    -------
    pd.DataFrame
        One row per event: event_number, test_name, config_key, water_temp,
        shower_on, bin{n}_quantaq_inhaled_dose, bin{n}_croom_inhaled_dose,
        n_minutes_covered_quantaq, n_minutes_covered_croom.
    """
    print("\nLoading QuantAQ-inside (raw) concentration...")
    quantaq_conc = resample_to_1min(load_quantaq_data("inside"))
    quantaq_dose = compute_cumulative_dose(quantaq_conc, events, label="quantaq")

    print("\nLoading C_room (fleet average) concentration...")
    croom_conc = resample_to_1min(build_croom_data())
    croom_dose = compute_cumulative_dose(croom_conc, events, label="croom")

    shared_id_cols = ["test_name", "config_key", "water_temp", "shower_on"]
    dose_df = quantaq_dose.merge(
        croom_dose.drop(columns=shared_id_cols),
        on="event_number",
        how="left",
    )
    return dose_df


# =============================================================================
# Excel Output
# =============================================================================


def _save_results(dose_df: pd.DataFrame, output_dir: Path) -> None:
    """
    Write the per-event, per-bin cumulative inhaled dose to Excel.

    Parameters
    ----------
    dose_df : pd.DataFrame
        Output of _compute_dose.
    output_dir : Path
        Directory to write inhaled_dose_summary.xlsx into.
    """
    output_file_path = output_dir / SUMMARY_XLSX_NAME

    if dose_df.empty:
        print(f"\nNo results to save - skipping {output_file_path}")
        return

    column_rename = {}
    for n in PARTICLE_BINS:
        for label, _legend, _color in SOURCES:
            column_rename[f"bin{n}_{label}_inhaled_dose"] = f"bin{n}_{label}_inhaled_dose (#)"

    export_df = dose_df.rename(columns=column_rename)
    export_df = sf.apply_sig_figs_to_df(export_df)

    with pd.ExcelWriter(output_file_path, engine="openpyxl") as writer:
        export_df.to_excel(writer, sheet_name=SUMMARY_SHEET, index=False)

    print(f"\nResults saved to: {output_file_path}")


# =============================================================================
# Bokeh Per-Bin Time Series
# =============================================================================


def _make_bin_figure(dose_df: pd.DataFrame, bin_index: int, output_path: Path) -> None:
    """
    Build and save one bin's cumulative-inhaled-dose time-series figure.

    Plots both concentration sources (QuantAQ-inside, C_room) as separate
    point series, one point per event per source at its shower_on time.
    Style matches scripts/particle_emission_variant_figures.py's per-bin
    figures.

    Parameters
    ----------
    dose_df : pd.DataFrame
        Output of _compute_dose.
    bin_index : int
        Particle-size bin index (0-11).
    output_path : Path
        Destination HTML path.
    """
    bin_name = PARTICLE_BINS[bin_index]["name"]

    output_file(str(output_path), title=f"Cumulative inhaled dose, bin {bin_index}")

    fig = figure(
        x_axis_type="datetime",
        x_axis_label="Event shower start (date and time)",
        y_axis_label=f"Cumulative inhaled dose (#), bin {bin_index} ({bin_name} µm)",
        tools="pan,box_zoom,wheel_zoom,reset,save",
    )

    hover = HoverTool(
        tooltips=[
            ("Source", "@source"),
            ("Event", "@event_number"),
            ("Test", "@test_name"),
            ("Config", "@config_key"),
            ("Shower on", "@time_str"),
            ("Dose (#)", "@value{0.00e+00}"),
            ("Minutes covered", "@n_minutes_covered"),
        ]
    )
    fig.add_tools(hover)

    fig.add_layout(
        Span(
            location=ROOM_CUTOVER,
            dimension="height",
            line_color="gray",
            line_dash="dashed",
            line_width=1.5,
        )
    )

    for label, legend_text, color in SOURCES:
        value_col = f"bin{bin_index}_{label}_inhaled_dose"
        coverage_col = f"n_minutes_covered_{label}"

        sub = dose_df[dose_df[value_col].notna()].copy()
        if sub.empty:
            print(f"    [WARN] Bin {bin_index} ({label}): no valid values.")
            continue

        source_data = ColumnDataSource(
            data={
                "time": sub["shower_on"],
                "value": sub[value_col],
                "event_number": sub["event_number"],
                "test_name": sub["test_name"],
                "config_key": sub["config_key"],
                "n_minutes_covered": sub[coverage_col],
                "source": [legend_text] * len(sub),
                "time_str": sub["shower_on"].dt.strftime("%Y-%m-%d %H:%M:%S"),
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
            legend_label=legend_text,
        )

    style_moduair_figure(fig, legend_title="Source", legend_location="top_right")

    fig.xaxis.formatter = DatetimeTickFormatter(days="%Y-%m-%d", hours="%m-%d %H:%M", minutes="%H:%M")

    save(fig)


# =============================================================================
# Main Pipeline
# =============================================================================


def run_inhalation_dose_analysis(
    output_dir: Optional[Path] = None,
    apply_sig_figs: bool = True,
) -> pd.DataFrame:
    """
    Run the complete cumulative inhaled-dose analysis.

    Parameters
    ----------
    output_dir : Path, optional
        Output directory (defaults to data_root/output).
    apply_sig_figs : bool
        If True (default), round output data to SIG_FIGS_DATA significant
        figures. Pass False (via --no-sig-figs) to preserve full precision.

    Returns
    -------
    pd.DataFrame
        Per-event, per-bin cumulative inhaled dose results for both sources.
    """
    sf.set_enabled(apply_sig_figs)
    print("=" * 80)
    print("Particle Inhalation Dose Analysis")
    print("=" * 80)

    if output_dir is None:
        output_dir = get_data_root() / "output"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\nLoading events from registry...")
    events = _load_valid_events(output_dir)
    if not events:
        print("No valid events found -- nothing to do.")
        return pd.DataFrame()

    print("\nComputing cumulative inhaled dose...")
    dose_df = _compute_dose(events)
    print(f"  Computed dose for {len(dose_df)} event(s)")

    _save_results(dose_df, output_dir)

    figure_dir = output_dir.joinpath(*FIGURE_SUBDIR)
    figure_dir.mkdir(parents=True, exist_ok=True)

    print("\nBuilding per-bin time-series figures...")
    for bin_index in PARTICLE_BINS:
        output_path = figure_dir / f"inhaled_dose_bin{bin_index}.html"
        _make_bin_figure(dose_df, bin_index, output_path)
        print(f"  Saved {output_path.name}")

    print("\n" + "=" * 80)
    print("Done")
    print("=" * 80)

    return dose_df


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Cumulative inhaled particle dose analysis: per-event, per-bin "
        "excess particle count inhaled over 2 h 10 min from shower-on."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for results (default: data_root/output)",
    )
    parser.add_argument(
        "--no-sig-figs",
        action="store_true",
        help="Disable significant figure rounding on output data "
        "(default: 3 sig figs)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else None

    run_inhalation_dose_analysis(
        output_dir=output_dir,
        apply_sig_figs=not args.no_sig_figs,
    )


if __name__ == "__main__":
    main()
