#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cumulative Inhaled Particle Dose
==================================

Computes, per shower event and per particle-size bin, the cumulative excess
particle count inhaled over a 2 h 10 min window starting at shower-on. This is
a standalone exposure metric: it does not depend on the penetration/
deposition/emission-rate calculations in src/particle_calculations.py, only on
an indoor particle concentration series.

Methodology
-----------
For each event and each of the 12 OPC bins:
    1. baseline = mean indoor concentration over the 15 minutes immediately
       before shower-on (BASELINE_WINDOW_MIN), i.e. [shower_on - 15 min,
       shower_on).
    2. For each minute t in [shower_on, shower_on + 130 min] (DOSE_WINDOW_MIN
       = 2 h 10 min):
           excess(t) = max(C(t) - baseline, 0)   [#/cm3]
           dose(t)   = excess(t) * BREATH_VOLUME_M3 * CM3_PER_M3   [#]
       Excess is clipped to zero so only above-baseline concentration
       contributes; minutes below baseline do not offset the running total.
    3. cumulative_dose = sum of dose(t) over the window (NaN minutes skipped).

Two concentration series, computed as two independent, non-blended dose
values (see scripts/particle_inhalation_dose_analysis.py):
    - QuantAQ-inside (MOD-PM-00195), raw and unmodified, for every event.
    - C_room(t), the position-weighted MODULAIR-PM fleet average (see
      src.particle_room_correction.build_croom_data), which only exists
      2026-06-03 through 2026-07-16; events outside that window get NaN.
These are NOT merged into one series -- both are reported side by side so
the raw single-sensor reading and the fleet average can be compared directly
during the co-location period.

Key Functions:
    - resample_to_1min: Regularize a particle-bin DataFrame to a 1-minute
      grid, interpolating gaps up to 5 minutes.
    - compute_cumulative_dose: Per-event, per-bin cumulative inhaled dose for
      one concentration series.

Input Files:
    - None directly; consumes DataFrames from src.particle_data_loader and
      src.particle_room_correction, and the events list from
      src.particle_data_loader.get_events_from_registry.

Output Files:
    - None; returns DataFrames consumed by
      scripts/particle_inhalation_dose_analysis.py.

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-23
Update log:
    2026-09-23  Initial version: single blended (raw pre-cutover, C_room
        fleet-period) concentration series per event.
    2026-09-23  Replaced the blended series with two independent, labeled
        series (QuantAQ-inside for all events, C_room for the fleet window
        only) reported side by side rather than substituted into one column;
        added resample_to_1min so the caller can prepare either source the
        same way.
"""

from typing import Dict, List

import numpy as np
import pandas as pd

from src.particle_calculations import CM3_PER_M3, PARTICLE_BINS

# =============================================================================
# Configuration
# =============================================================================

# Average tidal/minute breath volume used for the dose calculation (m^3).
BREATH_VOLUME_M3 = 0.006

# Minutes before shower-on averaged to obtain the pre-shower baseline
# concentration for each bin.
BASELINE_WINDOW_MIN = 15

# Minutes after shower-on over which the cumulative dose is summed
# (2 h 10 min total).
DOSE_WINDOW_MIN = 130


# =============================================================================
# Concentration Preparation
# =============================================================================


def resample_to_1min(df: pd.DataFrame) -> pd.DataFrame:
    """
    Regularize a particle-bin DataFrame to a 1-minute grid.

    Matches the resampling step in
    src.particle_data_loader.load_and_merge_quantaq_data (mean-resample to
    1 minute, then linearly interpolate gaps up to 5 minutes), so both the
    QuantAQ-inside and C_room series are prepared identically before the dose
    calculation.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with a "datetime" column and particle-bin columns.

    Returns
    -------
    pd.DataFrame
        1-minute-resolution DataFrame with the same columns.
    """
    resampled = df.set_index("datetime").sort_index().resample("1min").mean()
    resampled = resampled.interpolate(method="linear", limit=5)
    return resampled.reset_index()


# =============================================================================
# Cumulative Dose Calculation
# =============================================================================


def compute_cumulative_dose(
    concentration_df: pd.DataFrame,
    events: List[Dict],
    label: str,
    bin_col_template: str = "opc_bin{n}",
) -> pd.DataFrame:
    """
    Compute cumulative excess inhaled particle dose per event, per bin, for
    one concentration series.

    Parameters
    ----------
    concentration_df : pd.DataFrame
        1-minute-resolution particle data (see resample_to_1min) with a
        "datetime" column and bin_col_template-formatted columns.
    events : list of dict
        Event dicts (from src.particle_data_loader.get_events_from_registry),
        each with at least event_number, test_name, config_key, water_temp,
        and shower_on.
    label : str
        Short name for this concentration source (e.g. "quantaq", "croom"),
        used to name the output columns.
    bin_col_template : str
        Format template for the input bin columns, filled with n=0..11.

    Returns
    -------
    pd.DataFrame
        One row per event with columns: event_number, test_name, config_key,
        water_temp, shower_on, bin0_{label}_inhaled_dose ..
        bin11_{label}_inhaled_dose (particle count, #), and
        n_minutes_covered_{label} (number of non-NaN minutes found in the
        dose window, out of DOSE_WINDOW_MIN + 1 expected). A bin is NaN for
        an event when either the pre-shower baseline or the entire dose
        window has no data (e.g. every non-fleet-period event when label is
        "croom").
    """
    concentration = concentration_df.set_index("datetime").sort_index()
    coverage_col = f"n_minutes_covered_{label}"

    rows = []
    for event in events:
        shower_on = event["shower_on"]
        baseline_start = shower_on - pd.Timedelta(minutes=BASELINE_WINDOW_MIN)
        window_end = shower_on + pd.Timedelta(minutes=DOSE_WINDOW_MIN)

        row = {
            "event_number": event.get("event_number"),
            "test_name": event.get("test_name"),
            "config_key": event.get("config_key", ""),
            "water_temp": event.get("water_temp", ""),
            "shower_on": shower_on,
        }

        n_minutes_covered = 0
        for bin_num in PARTICLE_BINS:
            col = bin_col_template.format(n=bin_num)
            dose_col = f"bin{bin_num}_{label}_inhaled_dose"

            if col not in concentration.columns:
                row[dose_col] = np.nan
                continue

            baseline_series = concentration.loc[baseline_start:shower_on, col]
            baseline_series = baseline_series[baseline_series.index < shower_on]
            baseline = baseline_series.mean()

            window_series = concentration.loc[shower_on:window_end, col]

            if pd.isna(baseline) or window_series.notna().sum() == 0:
                row[dose_col] = np.nan
                continue

            excess = (window_series - baseline).clip(lower=0.0)
            dose_per_min = excess * BREATH_VOLUME_M3 * CM3_PER_M3
            row[dose_col] = dose_per_min.sum(skipna=True)
            n_minutes_covered = max(n_minutes_covered, int(window_series.notna().sum()))

        row[coverage_col] = n_minutes_covered
        rows.append(row)

    return pd.DataFrame(rows)
