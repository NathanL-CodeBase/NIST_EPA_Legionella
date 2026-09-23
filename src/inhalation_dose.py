#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cumulative Inhaled Particle Dose
==================================

Computes, per shower event and per particle-size bin, the cumulative excess
particle count inhaled over a 2 h 10 min window starting at shower-on. This is
a standalone exposure metric: it does not depend on the penetration/
deposition/emission-rate calculations in src/particle_calculations.py, only on
the indoor particle concentration series itself.

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

Concentration source: the caller supplies the 1-minute merged "inside" series
from src.particle_data_loader.load_and_merge_quantaq_data (default
inside_builder=build_raw_inside_data), i.e. raw QuantAQ-inside (MOD-PM-00195)
except 2026-06-03 through 2026-07-16, where it is the position-weighted
MODULAIR-PM fleet average C_room(t) (see src/particle_room_correction.py).
After 2026-07-16 the series reverts to raw QuantAQ-inside.

Key Functions:
    - compute_cumulative_dose: Per-event, per-bin cumulative inhaled dose.

Input Files:
    - None directly; consumes the merged DataFrame from
      src.particle_data_loader.load_and_merge_quantaq_data and the events list
      from src.particle_data_loader.get_events_from_registry.

Output Files:
    - None; returns a DataFrame consumed by
      scripts/particle_inhalation_dose_analysis.py.

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-23
Update log:
    2026-09-23  Initial version.
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
# Cumulative Dose Calculation
# =============================================================================


def compute_cumulative_dose(merged_df: pd.DataFrame, events: List[Dict]) -> pd.DataFrame:
    """
    Compute cumulative excess inhaled particle dose per event, per bin.

    Parameters
    ----------
    merged_df : pd.DataFrame
        1-minute merged particle data with a "datetime" column and
        "opc_bin{n}_inside" columns (from
        src.particle_data_loader.load_and_merge_quantaq_data).
    events : list of dict
        Event dicts (from src.particle_data_loader.get_events_from_registry),
        each with at least event_number, test_name, config_key, water_temp,
        and shower_on.

    Returns
    -------
    pd.DataFrame
        One row per event with columns: event_number, test_name, config_key,
        water_temp, shower_on, bin0_inhaled_dose .. bin11_inhaled_dose
        (particle count, #), and n_minutes_covered (number of non-NaN minutes
        found in the dose window, out of DOSE_WINDOW_MIN + 1 expected). A bin
        is NaN for an event when either the pre-shower baseline or the entire
        dose window has no data.
    """
    merged = merged_df.set_index("datetime").sort_index()

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
            col = f"opc_bin{bin_num}_inside"
            dose_col = f"bin{bin_num}_inhaled_dose"

            if col not in merged.columns:
                row[dose_col] = np.nan
                continue

            baseline_series = merged.loc[baseline_start:shower_on, col]
            baseline_series = baseline_series[baseline_series.index < shower_on]
            baseline = baseline_series.mean()

            window_series = merged.loc[shower_on:window_end, col]

            if pd.isna(baseline) or window_series.notna().sum() == 0:
                row[dose_col] = np.nan
                continue

            excess = (window_series - baseline).clip(lower=0.0)
            dose_per_min = excess * BREATH_VOLUME_M3 * CM3_PER_M3
            row[dose_col] = dose_per_min.sum(skipna=True)
            n_minutes_covered = max(n_minutes_covered, int(window_series.notna().sum()))

        row["n_minutes_covered"] = n_minutes_covered
        rows.append(row)

    return pd.DataFrame(rows)
