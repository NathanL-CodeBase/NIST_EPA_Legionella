#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle Room-Concentration Correction
=======================================

Builds the corrected indoor ("inside") particle concentration series consumed
by scripts/particle_decay_analysis.py, replacing the single-sensor MOD-PM-00195
reading with a position-weighted MODULAIR-PM fleet average where the fleet is
available, and with a fleet-derived correction factor where it is not.

Background
----------
From 2026-06-03 through 2026-07-16, more than eight additional MODULAIR-PM
sensors were co-located with MOD-PM-00195 (see scripts/moduair_cave_ratio.py).
Before that period, MOD-PM-00195 was the only indoor sensor, so its reading
carries a position-specific bias relative to well-mixed room air that the
fleet comparison can correct for.

Two regimes:
    - 2026-06-03 to 2026-07-16 (fleet period): the "inside" concentration is
      replaced with C_room(t), the position-weighted fleet average defined in
      moduair_cave_ratio.py, at every timestamp where more than 8 fleet
      sensors report (moduair_cave_ratio.MIN_QUANTS). Timestamps with 8 or
      fewer reporting sensors are left as NaN; there is no fallback to the raw
      sensor for those gaps.
    - Before 2026-06-03: MOD-PM-00195 is corrected only inside each shower
      event's analysis window (shower_on - 15 min through deposition_end, the
      same window moduair_cave_ratio.py aligns events on -- the window that
      covers the penetration-factor calculation's beta_other, E, and
      peak-time windows), using:
          C_adjusted_room(t) = C_bed1(t) * C_room,average(minute) / C_bed1,average(minute)
      where the averaged curves come from moduair_cave_ratio.py's
      group_average_curve, computed once over the fleet period for whichever
      water-temperature bucket the event belongs to:
          water_temp < 38 C        -> "W24" bucket (n=2 events)
          38 <= water_temp <= 41 C -> "W38-W41" bucket
          water_temp > 41 C        -> "W49" bucket (n=2 events)
      These are the only three buckets moduair_cave_ratio.py computes a ratio
      curve for. Outside any event's window (e.g. the wider before/after
      windows the penetration-factor calculation uses), MOD-PM-00195 is left
      uncorrected.

Key Functions:
    - get_water_temp_bucket: Map a water-temp code to one of the three ratio
      buckets.
    - build_corrected_inside_data: Build the full corrected "inside" DataFrame
      consumed by src.particle_data_loader.load_and_merge_quantaq_data.

Input Files:
    - QuantAQ MOD-PM-00195 processed CSVs (via
      src.particle_data_loader.load_quantaq_data)
    - MODULAIR-PM fleet chunk files (via src.moduair_loader.load_fleet_bins)
    - event_log.csv (via scripts.moduair_cave_ratio.load_group_events)

Output Files:
    - None; returns a DataFrame consumed by src.particle_data_loader.

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-18
Update log:
    2026-09-18  Initial version.
"""

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from src.event_manager import get_water_temp_sort_key
from src.particle_calculations import DEPOSITION_WINDOW_HOURS, PARTICLE_BINS

# =============================================================================
# Configuration
# =============================================================================

# Fleet co-location window, matching scripts.moduair_cave_ratio DEFAULT_START/END.
ROOM_CUTOVER = pd.Timestamp("2026-06-03 00:00:00")
FLEET_END = pd.Timestamp("2026-07-16 23:59:59")


def get_water_temp_bucket(water_temp_code: str) -> Optional[str]:
    """
    Map a water-temperature code to the fleet ratio-curve bucket it should use.

    Only three buckets have a computed ratio curve (the water-temp groups in
    scripts.moduair_cave_ratio.WATER_TEMP_GROUPS): "W38-W41", "W24", "W49".
    Every pre-2026-06-03 event is assigned to whichever bucket is closest in
    water temperature.

    Parameters
    ----------
    water_temp_code : str
        Registry water_temp code, e.g. "W48", "W38", "W52pw".

    Returns
    -------
    str or None
        "W24" (< 38 C), "W38-W41" (38-41 C inclusive), "W49" (> 41 C), or None
        if no numeric temperature could be extracted from the code.
    """
    temp = get_water_temp_sort_key(str(water_temp_code) if water_temp_code else "")
    if not np.isfinite(temp):
        return None
    if temp < 38:
        return "W24"
    if temp <= 41:
        return "W38-W41"
    return "W49"


# =============================================================================
# Ratio-Curve Application
# =============================================================================


def _ratio_at_minutes(curve: pd.Series, minutes: np.ndarray) -> np.ndarray:
    """
    Look up a per-minute ratio curve at arbitrary integer minute offsets.

    moduair_cave_ratio.py's group_average_curve only has entries at minutes
    where at least one event contributed a valid (>8-sensor) reading, so
    interior gaps are interpolated. Minutes requested outside the curve's own
    range (a small tail past a slightly longer deposition window) reuse the
    nearest available ratio rather than propagating NaN into the corrected
    concentration.

    Parameters
    ----------
    curve : pd.Series
        Ratio indexed by integer minute offset from shower_on.
    minutes : np.ndarray
        Integer minute offsets to evaluate the curve at.

    Returns
    -------
    np.ndarray
        Ratio value at each requested minute.
    """
    full_range = range(int(curve.index.min()), int(curve.index.max()) + 1)
    dense = curve.reindex(full_range).interpolate(limit_direction="both")
    return dense.reindex(minutes).ffill().bfill().to_numpy()


# =============================================================================
# Main Correction
# =============================================================================


def build_corrected_inside_data(events: List[Dict]) -> pd.DataFrame:
    """
    Build the corrected "inside" particle concentration series.

    Parameters
    ----------
    events : list of dict
        Event dicts as produced by
        particle_data_loader.get_events_from_registry (or the
        process_events_with_management fallback), each with at least
        shower_on, shower_off, deposition_end, and water_temp.

    Returns
    -------
    pd.DataFrame
        Columns datetime, opc_bin0..opc_bin11, same shape as
        particle_data_loader.load_quantaq_data("inside"), with the
        fleet-period and pre-cutover corrections applied.
    """
    # Local imports: particle_data_loader imports this module, so importing it
    # back at module scope would create a cycle.
    from scripts.moduair_cave_ratio import (
        FLEET_SNS,
        MIN_QUANTS,
        PRE_SHOWER_LEAD,
        TARGET_SN,
        build_position_totals,
        compute_room_frame,
        group_average_curve,
        load_group_events,
    )
    from src.moduair_loader import list_available_sensors, load_fleet_bins
    from src.particle_data_loader import load_quantaq_data

    bin_cols = [f"opc_bin{n}" for n in PARTICLE_BINS]

    raw_inside = load_quantaq_data("inside").set_index("datetime").sort_index()
    corrected = raw_inside.copy()

    print("\nBuilding room-corrected inside concentration...")
    print("  Loading MODULAIR-PM fleet data for room-concentration correction...")
    available = set(list_available_sensors("raw"))
    wanted = [sn for sn in (FLEET_SNS + [TARGET_SN]) if sn in available]
    fleet = load_fleet_bins(wanted, start=ROOM_CUTOVER, end=FLEET_END)
    if TARGET_SN not in fleet:
        print(
            "  [WARN] MODULAIR-PM fleet data unavailable; inside concentration "
            "left as raw MOD-PM-00195 for the full record."
        )
        return corrected.reset_index()

    totals = build_position_totals(fleet)
    groups = load_group_events(ROOM_CUTOVER, FLEET_END)

    fleet_index = pd.date_range(ROOM_CUTOVER, FLEET_END, freq="1min", name="datetime")
    fleet_frame = pd.DataFrame(index=fleet_index, columns=bin_cols, dtype=float)

    # Pre-cutover events: resolve each event's window and water-temp bucket once,
    # shared across all 12 bins below.
    pre_event_info = []
    for event in events:
        shower_on = event.get("shower_on")
        if shower_on is None or pd.isna(shower_on) or shower_on >= ROOM_CUTOVER:
            continue
        bucket = get_water_temp_bucket(event.get("water_temp", ""))
        if bucket is None:
            continue
        deposition_end = event.get("deposition_end")
        if pd.isna(deposition_end):
            deposition_end = event.get("shower_off", shower_on) + pd.Timedelta(
                hours=DEPOSITION_WINDOW_HOURS
            )
        window_start = shower_on - PRE_SHOWER_LEAD
        win_mask = (corrected.index >= window_start) & (corrected.index <= deposition_end)
        idx = corrected.index[win_mask]
        if len(idx) == 0:
            continue
        minute_offsets = ((idx - shower_on).total_seconds() / 60).round().astype(int)
        pre_event_info.append((idx, minute_offsets, bucket))

    print(f"  Fleet-period replacement window: {len(fleet_index)} minutes")
    print(f"  Pre-cutover events to correct: {len(pre_event_info)}")

    for bin_num in PARTICLE_BINS:
        bin_col = f"opc_bin{bin_num}"
        if bin_col not in corrected.columns:
            continue

        cave = compute_room_frame(totals, bin_col)
        if cave.empty:
            print(f"  [WARN] Bin {bin_num}: no fleet data; leaving raw MOD-PM-00195.")
            continue

        # Fleet period: full continuous replacement, masked to the >8-sensor filter.
        c_room = cave["C_room"].where(cave["n_sensors"] > MIN_QUANTS)
        fleet_frame[bin_col] = c_room.reindex(fleet_index)

        # Pre-cutover: one ratio curve per bucket, built from this bin's fleet data.
        bucket_curves: Dict[str, pd.Series] = {}
        for bucket, bucket_events in groups.items():
            c_bed1_curve = group_average_curve(cave, bucket_events, "C_bed1")
            c_room_curve = group_average_curve(cave, bucket_events, "C_room")
            if c_bed1_curve.empty or c_room_curve.empty:
                continue
            common = c_bed1_curve.index.intersection(c_room_curve.index)
            safe_bed1 = c_bed1_curve.loc[common].where(c_bed1_curve.loc[common] != 0)
            bucket_curves[bucket] = (c_room_curve.loc[common] / safe_bed1).sort_index()

        for idx, minute_offsets, bucket in pre_event_info:
            curve = bucket_curves.get(bucket)
            if curve is None or curve.empty:
                continue
            ratio_vals = _ratio_at_minutes(curve, minute_offsets)
            raw_vals = corrected.loc[idx, bin_col].to_numpy(dtype=float)
            corrected.loc[idx, bin_col] = raw_vals * ratio_vals

    pre_part = corrected[corrected.index < ROOM_CUTOVER]
    post_fleet_part = corrected[corrected.index > FLEET_END]
    result = pd.concat([pre_part, fleet_frame, post_fleet_part]).sort_index()
    return result.reset_index().rename(columns={"index": "datetime"})
