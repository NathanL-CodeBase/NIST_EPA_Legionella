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

Concentration source: the primary "inside" series from
src.particle_data_loader.load_and_merge_quantaq_data (default
inside_builder=build_raw_inside_data) -- raw QuantAQ-inside (MOD-PM-00195)
for all events, except 2026-06-03 through 2026-07-16 where it is the
position-weighted MODULAIR-PM fleet average C_room(t) (see
src/particle_room_correction.py). Requires the unified event registry
(scripts/event_registry.py); duration-excluded and otherwise-excluded events
are skipped, matching scripts/particle_decay_analysis.py's convention.

Output Files:
    output/inhaled_dose_summary.xlsx (sheet: inhaled_dose_by_bin) -- per-event,
        per-bin cumulative inhaled dose (#) and n_minutes_covered.
    output/plots/particle/inhaled_dose_bin{N}.html (12 figures) -- Bokeh
        time-series: x = shower_on, y = cumulative inhaled dose, one point per
        event. Style matches scripts/particle_emission_variant_figures.py
        (src.plot_style.style_moduair_figure: 1600x800, 12pt, no title,
        click-to-hide legend, hover enabled), including the ROOM_CUTOVER
        reference line.
    output/plots/inhaled_dose_boxplot_{bin0-2,bin3-6,bin7-11}.png (3 figures)
        -- matplotlib box-and-whisker by water temperature (base W## events
        only), via src.plot_particle_boxplots.plot_inhaled_dose_boxplot.

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-23
Update log:
    2026-09-23  Initial version.
"""

import sys
import warnings
from pathlib import Path
from typing import Optional

import pandas as pd
from bokeh.models import ColumnDataSource, DatetimeTickFormatter, HoverTool, Label, Span
from bokeh.plotting import figure, output_file, save

warnings.filterwarnings("ignore")

# Add project root to path for src/ imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import src.sig_figs as sf  # noqa: E402
from src.data_paths import get_data_root  # noqa: E402
from src.event_manager import is_event_excluded  # noqa: E402
from src.inhalation_dose import compute_cumulative_dose  # noqa: E402
from src.particle_calculations import PARTICLE_BINS  # noqa: E402
from src.particle_data_loader import get_events_from_registry, load_and_merge_quantaq_data  # noqa: E402
from src.particle_room_correction import ROOM_CUTOVER  # noqa: E402
from src.plot_particle_boxplots import plot_inhaled_dose_boxplot  # noqa: E402
from src.plot_style import MODUAIR_TEXT_PT, SENSOR_COLORS, style_moduair_figure  # noqa: E402

# =============================================================================
# Configuration
# =============================================================================

SUMMARY_XLSX_NAME = "inhaled_dose_summary.xlsx"
SUMMARY_SHEET = "inhaled_dose_by_bin"

FIGURE_SUBDIR = ("plots", "particle")
BOXPLOT_XRANGE = (5, 55, 5)  # (xmin, xmax, xtick_step) in degrees C

CUTOVER_LABELS = ("Raw C_bed1 / C_adjusted room", "C_room")
CUTOVER_LABEL_COLOR = "gray"


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
# Excel Output
# =============================================================================


def _save_results(dose_df: pd.DataFrame, output_dir: Path) -> None:
    """
    Write the per-event, per-bin cumulative inhaled dose to Excel.

    Parameters
    ----------
    dose_df : pd.DataFrame
        Output of src.inhalation_dose.compute_cumulative_dose.
    output_dir : Path
        Directory to write inhaled_dose_summary.xlsx into.
    """
    output_file_path = output_dir / SUMMARY_XLSX_NAME

    if dose_df.empty:
        print(f"\nNo results to save - skipping {output_file_path}")
        return

    column_rename = {f"bin{n}_inhaled_dose": f"bin{n}_inhaled_dose (#)" for n in PARTICLE_BINS}
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

    One point per event at its shower_on time. Style matches
    scripts/particle_emission_variant_figures.py's per-bin figures.

    Parameters
    ----------
    dose_df : pd.DataFrame
        Output of src.inhalation_dose.compute_cumulative_dose.
    bin_index : int
        Particle-size bin index (0-11).
    output_path : Path
        Destination HTML path.
    """
    bin_name = PARTICLE_BINS[bin_index]["name"]
    value_col = f"bin{bin_index}_inhaled_dose"

    output_file(str(output_path), title=f"Cumulative inhaled dose, bin {bin_index}")

    fig = figure(
        x_axis_type="datetime",
        x_axis_label="Event shower start (date and time)",
        y_axis_label=f"Cumulative inhaled dose (#), bin {bin_index} ({bin_name} µm)",
        tools="pan,box_zoom,wheel_zoom,reset,save",
    )

    hover = HoverTool(
        tooltips=[
            ("Event", "@event_number"),
            ("Test", "@test_name"),
            ("Config", "@config_key"),
            ("Shower on", "@time_str"),
            ("Dose (#)", "@value{0.00e+00}"),
            ("Minutes covered", "@n_minutes_covered"),
        ]
    )
    fig.add_tools(hover)

    sub = dose_df[dose_df[value_col].notna()].copy()
    if not sub.empty:
        source_data = ColumnDataSource(
            data={
                "time": sub["shower_on"],
                "value": sub[value_col],
                "event_number": sub["event_number"],
                "test_name": sub["test_name"],
                "config_key": sub["config_key"],
                "n_minutes_covered": sub["n_minutes_covered"],
                "time_str": sub["shower_on"].dt.strftime("%Y-%m-%d %H:%M:%S"),
            }
        )
        fig.scatter(
            "time",
            "value",
            source=source_data,
            marker="circle",
            size=7,
            color=SENSOR_COLORS[0],
            alpha=0.8,
            legend_label="Cumulative inhaled dose",
        )

    style_moduair_figure(fig, legend_title="Series", legend_location="top_right")

    fig.xaxis.formatter = DatetimeTickFormatter(days="%Y-%m-%d", hours="%m-%d %H:%M", minutes="%H:%M")

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
    fig.add_layout(Label(x=ROOM_CUTOVER, x_offset=-10, text=CUTOVER_LABELS[0], text_align="right", **label_kwargs))
    fig.add_layout(Label(x=ROOM_CUTOVER, x_offset=10, text=CUTOVER_LABELS[1], text_align="left", **label_kwargs))

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
        Per-event, per-bin cumulative inhaled dose results.
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

    particle_data = load_and_merge_quantaq_data(events)

    print("\nComputing cumulative inhaled dose...")
    dose_df = compute_cumulative_dose(particle_data, events)
    print(f"  Computed dose for {len(dose_df)} event(s)")

    _save_results(dose_df, output_dir)

    figure_dir = output_dir.joinpath(*FIGURE_SUBDIR)
    figure_dir.mkdir(parents=True, exist_ok=True)

    print("\nBuilding per-bin time-series figures...")
    for bin_index in PARTICLE_BINS:
        output_path = figure_dir / f"inhaled_dose_bin{bin_index}.html"
        _make_bin_figure(dose_df, bin_index, output_path)
        print(f"  Saved {output_path.name}")

    print("\nBuilding water-temperature boxplots...")
    rh_data = None
    try:
        from src.data_paths import get_common_file

        bc = pd.read_excel(get_common_file("rh_temp_wind_summary"), sheet_name="Bedroom_Conditions")
        bc = bc.rename(columns={"shower_on": "datetime", "rh_mean (%)": "RH_bedroom"})
        rh_data = bc[["datetime", "RH_bedroom"]].copy()
        rh_data["datetime"] = pd.to_datetime(rh_data["datetime"])
        print("  Loaded Bedroom_Conditions RH for boxplot annotations.")
    except Exception as rh_err:
        print(f"  Note: Could not load Bedroom_Conditions RH data (n= only annotations): {rh_err}")

    try:
        plot_inhaled_dose_boxplot(
            dose_df,
            PARTICLE_BINS,
            output_dir / "plots" / "inhaled_dose_boxplot.png",
            rh_data=rh_data,
            x_range=BOXPLOT_XRANGE,
        )
        print("  Generated: inhaled_dose_boxplot_{bin0-2,bin3-6,bin7-11}.png")
    except Exception as e:
        print(f"  Error generating inhaled_dose_boxplot: {e}")

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
