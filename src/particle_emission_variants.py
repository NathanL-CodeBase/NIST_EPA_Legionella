#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle Emission Rate Variants (E(t)1-E(t)4)
================================================

Per-event, per-bin analysis producing four emission-rate variants from the
combinations below (the report's Table 7 / Equation 10), rather than
committing to one set of inputs. Used by scripts/particle_decay_analysis.py;
not run directly.

    Variant   Lambda / beta source   Concentration series         Volume
    -------   ---------------------  ----------------------------  -----------
    E(t)1     lambda_outside          primary (raw C_bed1(t) /      36.1086 m3
              (bin{n}_outside_*)      C_room(t))                    (bedroom)
    E(t)2     lambda_entry            primary (same series as E1)   36.1086 m3
              (bin{n}_entry_*)
    E(t)3     lambda_outside          adjusted (ratio-corrected      36.1086 m3
              (bin{n}_adjusted_*)     C_adjusted_room(t))
    E(t)4     lambda_outside          primary (same series as E1)   54.5868 m3
              (bin{n}_bathroom_*)                                   (+bathroom)

The "primary" and "adjusted" concentration series are built by
src.particle_room_correction.build_raw_inside_data and
build_corrected_inside_data respectively, and merged into two parallel
particle_data DataFrames by src.particle_data_loader.load_and_merge_quantaq_data
(inside_builder parameter) before this module is called.

E(t)3 is computed only for pre-cutover (single-monitor) events -- before
src.particle_room_correction.ROOM_CUTOVER (2026-06-03) -- since
C_adjusted_room(t) is a correction for single-monitor position bias that
doesn't apply once the fleet is in place. Post-cutover bins are NaN with
skip_reason noting the restriction.

E(t)4 does not refit beta_other: calculate_other_process_rate has no volume
term (confirmed in src/particle_calculations.py), so E(t)4 reuses E(t)1's
already-fitted p, beta_other, and peak_time, and only reruns
calculate_emission_rate/calculate_ct_prediction with the larger volume.

Key Functions:
    - analyze_event_all_bins: Compute all four variants, all twelve bins, for
      one shower event.

Output Files:
    None -- returns a dict consumed by scripts/particle_decay_analysis.py.

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-22
Update log:
    2026-09-22  Split out of scripts/particle_decay_analysis.py's
        analyze_event_all_bins (outside/entry only) and extended with the
        adjusted (E3) and bathroom (E4) variants.
"""

from datetime import timedelta
from typing import Dict

import numpy as np
import pandas as pd

from src.particle_calculations import (
    BEDROOM_BATHROOM_VOLUME_M3,
    BEDROOM_VOLUME_M3,
    PARTICLE_BINS,
    calculate_ct_prediction,
    calculate_emission_rate,
    calculate_other_process_rate,
    calculate_penetration_factor,
    find_peak_time,
)
from src.particle_room_correction import ROOM_CUTOVER

# =============================================================================
# Variant Bookkeeping
# =============================================================================

# The two CO2-bounded, primary-series variants (E1, E2). "adjusted" (E3) and
# "bathroom" (E4) are handled separately below -- they aren't independent
# lambda choices in the same sense (E4 reuses E1's lambda/beta/p; E3 uses a
# different concentration series, not a different lambda).
_LAMBDA_SOURCES = ("outside", "entry")

# Per-source/variant result keys reset to NaN/empty placeholders whenever a
# bin can't be analyzed for that variant (missing lambda, invalid beta, etc.).
# Listed once here so every early-exit branch below fills the same set.
_SOURCE_SCALAR_NAN_KEYS = (
    "beta_other",
    "beta_other_raw_mean",
    "beta_other_std",
    "beta_other_r_squared",
    "beta_step",
    "c_steady_state",
    "E_mean",
    "E_std",
    "E_total",
    "E_r_squared",
    "peak_measured",
    "peak_predicted",
    "deposition_end_measured",
    "deposition_end_predicted",
)
_SOURCE_LIST_KEYS = (
    "ct_datetimes",
    "ct_predicted",
    "emission_datetimes",
    "emission_predicted",
    "decay_datetimes",
    "decay_predicted",
    "E_times",
    "E_per_step",
)


def _fill_source_nan(results: Dict, bin_num: int, source: str, skip_reason: str) -> None:
    """Fill NaN/empty placeholders for one bin/variant pair on early exit."""
    prefix = f"bin{bin_num}_{source}"
    for key in _SOURCE_SCALAR_NAN_KEYS:
        results[f"{prefix}_{key}"] = np.nan
    for key in _SOURCE_LIST_KEYS:
        results[f"{prefix}_{key}"] = []
    results[f"{prefix}_skip_reason"] = skip_reason


def _fill_adjusted_nan(results: Dict, bin_num: int, reason: str) -> None:
    """Fill NaN placeholders for the adjusted (E3) variant, including its own
    independent p/peak_time fields (unlike bathroom, adjusted uses a
    different concentration series so it can't share the primary p/peak_time)."""
    results[f"bin{bin_num}_adjusted_p_mean"] = np.nan
    results[f"bin{bin_num}_adjusted_p_std"] = np.nan
    results[f"bin{bin_num}_adjusted_peak_time"] = None
    _fill_source_nan(results, bin_num, "adjusted", reason)


# =============================================================================
# Per-Variant Computation
# =============================================================================


def _store_peak_deposition_comparison(
    results: Dict,
    particle_data: pd.DataFrame,
    event: Dict,
    bin_num: int,
    prefix: str,
    peak_time,
    ct_result: Dict,
) -> None:
    """Store measured-vs-predicted concentration at peak_time and
    deposition_end for one bin/variant, given its Ct prediction result."""
    col_inside = f"{PARTICLE_BINS[bin_num]['column']}_inside"
    emission_pred = ct_result.get("emission_predicted", [])
    decay_pred = ct_result.get("decay_predicted", [])

    # Use argmin() + to_numpy() to avoid pandas Scalar typing ambiguity.
    col_values = particle_data[col_inside].to_numpy(dtype=float, na_value=np.nan)

    if peak_time is not None and col_inside in particle_data.columns:
        time_diffs = (particle_data["datetime"] - pd.Timestamp(peak_time)).abs()
        peak_row = int(time_diffs.argmin())
        results[f"{prefix}_peak_measured"] = float(col_values[peak_row])
    else:
        results[f"{prefix}_peak_measured"] = np.nan
    results[f"{prefix}_peak_predicted"] = (
        float(emission_pred[-1]) if len(emission_pred) > 0 else np.nan
    )

    depo_end = event["deposition_end"]
    if col_inside in particle_data.columns:
        time_diffs_end = (particle_data["datetime"] - pd.Timestamp(depo_end)).abs()
        end_row = int(time_diffs_end.argmin())
        results[f"{prefix}_deposition_end_measured"] = float(col_values[end_row])
    else:
        results[f"{prefix}_deposition_end_measured"] = np.nan
    results[f"{prefix}_deposition_end_predicted"] = (
        float(decay_pred[-1]) if len(decay_pred) > 0 else np.nan
    )


def _compute_variant(
    results: Dict,
    particle_data: pd.DataFrame,
    event: Dict,
    bin_num: int,
    prefix: str,
    p_mean: float,
    peak_time,
    lambda_val: float,
    beta_window_start,
    beta_window_end,
    volume_m3: float,
) -> float:
    """
    Fit beta, emission rate, and Ct prediction for one independent bin/variant
    combination (its own concentration series, lambda, p, and peak_time), and
    store every result key under `prefix`. Shared by the outside, entry, and
    adjusted variants. The bathroom variant does not call this -- it reuses
    the outside variant's beta_val rather than refitting (see
    _compute_bathroom_variant).

    Returns
    -------
    float
        The fitted beta (np.nan if the bin was skipped for this variant).
    """
    beta_result = calculate_other_process_rate(
        particle_data, beta_window_start, beta_window_end, bin_num, p_mean, lambda_val
    )

    results[f"{prefix}_beta_other"] = beta_result.get("beta", np.nan)
    results[f"{prefix}_beta_other_raw_mean"] = beta_result.get("beta_raw_mean", np.nan)
    results[f"{prefix}_beta_other_std"] = beta_result.get("beta_std", np.nan)
    results[f"{prefix}_beta_other_r_squared"] = beta_result.get("beta_r_squared", np.nan)
    results[f"{prefix}_beta_step"] = beta_result.get("beta_step", np.nan)
    results[f"{prefix}_c_steady_state"] = beta_result.get("c_steady_state", np.nan)

    if np.isnan(beta_result.get("beta", np.nan)):
        results[f"{prefix}_E_mean"] = np.nan
        results[f"{prefix}_E_std"] = np.nan
        results[f"{prefix}_E_total"] = np.nan
        results[f"{prefix}_E_r_squared"] = np.nan
        results[f"{prefix}_skip_reason"] = beta_result.get("skip_reason", "Unknown")
        results[f"{prefix}_peak_measured"] = np.nan
        results[f"{prefix}_peak_predicted"] = np.nan
        results[f"{prefix}_deposition_end_measured"] = np.nan
        results[f"{prefix}_deposition_end_predicted"] = np.nan
        for key in _SOURCE_LIST_KEYS:
            results[f"{prefix}_{key}"] = []
        return np.nan

    beta_val = beta_result["beta"]

    E_result = calculate_emission_rate(
        particle_data,
        event["shower_on"],
        peak_time,
        bin_num,
        p_mean,
        lambda_val,
        beta_val,
        volume_m3=volume_m3,
    )
    results[f"{prefix}_E_mean"] = E_result.get("E_mean", np.nan)
    results[f"{prefix}_E_std"] = E_result.get("E_std", np.nan)
    results[f"{prefix}_E_total"] = E_result.get("E_total", np.nan)
    results[f"{prefix}_skip_reason"] = E_result.get("skip_reason", None)
    results[f"{prefix}_E_times"] = E_result.get("E_times", [])
    results[f"{prefix}_E_per_step"] = E_result.get("E_per_step", [])

    E_mean_val = E_result.get("E_mean", np.nan)
    effective_E = E_mean_val if not np.isnan(E_mean_val) else 0.0
    ct_result = calculate_ct_prediction(
        particle_data,
        event["shower_on"],
        event["shower_off"],
        event["deposition_end"],
        bin_num,
        p_mean,
        lambda_val,
        beta_val,
        effective_E,
        peak_time,
        volume_m3=volume_m3,
    )
    results[f"{prefix}_ct_datetimes"] = ct_result.get("datetimes", [])
    results[f"{prefix}_ct_predicted"] = ct_result.get("predicted_ct", [])
    results[f"{prefix}_emission_datetimes"] = ct_result.get("emission_datetimes", [])
    results[f"{prefix}_emission_predicted"] = ct_result.get("emission_predicted", [])
    results[f"{prefix}_decay_datetimes"] = ct_result.get("decay_datetimes", [])
    results[f"{prefix}_decay_predicted"] = ct_result.get("decay_predicted", [])
    results[f"{prefix}_E_r_squared"] = ct_result.get("E_r_squared", np.nan)

    _store_peak_deposition_comparison(
        results, particle_data, event, bin_num, prefix, peak_time, ct_result
    )

    return beta_val


def _compute_bathroom_variant(
    results: Dict,
    particle_data: pd.DataFrame,
    event: Dict,
    bin_num: int,
    p_mean: float,
    peak_time,
    lambda_outside: float,
    beta_outside_val: float,
) -> None:
    """
    Compute E(t)4 (bathroom-included volume) for one bin, reusing the outside
    variant's already-fitted p, beta_other, and peak_time -- beta_other has no
    volume term (src/particle_calculations.calculate_other_process_rate), so
    refitting it would just reproduce the outside value.
    """
    prefix = f"bin{bin_num}_bathroom"
    outside_prefix = f"bin{bin_num}_outside"

    if np.isnan(beta_outside_val):
        _fill_source_nan(
            results,
            bin_num,
            "bathroom",
            results.get(f"{outside_prefix}_skip_reason", "Outside variant unavailable"),
        )
        return

    results[f"{prefix}_beta_other"] = results[f"{outside_prefix}_beta_other"]
    results[f"{prefix}_beta_other_raw_mean"] = results[f"{outside_prefix}_beta_other_raw_mean"]
    results[f"{prefix}_beta_other_std"] = results[f"{outside_prefix}_beta_other_std"]
    results[f"{prefix}_beta_other_r_squared"] = results[f"{outside_prefix}_beta_other_r_squared"]
    results[f"{prefix}_beta_step"] = results[f"{outside_prefix}_beta_step"]
    results[f"{prefix}_c_steady_state"] = results[f"{outside_prefix}_c_steady_state"]

    E_result = calculate_emission_rate(
        particle_data,
        event["shower_on"],
        peak_time,
        bin_num,
        p_mean,
        lambda_outside,
        beta_outside_val,
        volume_m3=BEDROOM_BATHROOM_VOLUME_M3,
    )
    results[f"{prefix}_E_mean"] = E_result.get("E_mean", np.nan)
    results[f"{prefix}_E_std"] = E_result.get("E_std", np.nan)
    results[f"{prefix}_E_total"] = E_result.get("E_total", np.nan)
    results[f"{prefix}_skip_reason"] = E_result.get("skip_reason", None)
    results[f"{prefix}_E_times"] = E_result.get("E_times", [])
    results[f"{prefix}_E_per_step"] = E_result.get("E_per_step", [])

    E_mean_val = E_result.get("E_mean", np.nan)
    effective_E = E_mean_val if not np.isnan(E_mean_val) else 0.0
    ct_result = calculate_ct_prediction(
        particle_data,
        event["shower_on"],
        event["shower_off"],
        event["deposition_end"],
        bin_num,
        p_mean,
        lambda_outside,
        beta_outside_val,
        effective_E,
        peak_time,
        volume_m3=BEDROOM_BATHROOM_VOLUME_M3,
    )
    results[f"{prefix}_ct_datetimes"] = ct_result.get("datetimes", [])
    results[f"{prefix}_ct_predicted"] = ct_result.get("predicted_ct", [])
    results[f"{prefix}_emission_datetimes"] = ct_result.get("emission_datetimes", [])
    results[f"{prefix}_emission_predicted"] = ct_result.get("emission_predicted", [])
    results[f"{prefix}_decay_datetimes"] = ct_result.get("decay_datetimes", [])
    results[f"{prefix}_decay_predicted"] = ct_result.get("decay_predicted", [])
    results[f"{prefix}_E_r_squared"] = ct_result.get("E_r_squared", np.nan)

    _store_peak_deposition_comparison(
        results, particle_data, event, bin_num, prefix, peak_time, ct_result
    )


# =============================================================================
# Per-Event Orchestration
# =============================================================================


def analyze_event_all_bins(
    particle_data: pd.DataFrame,
    particle_data_adjusted: pd.DataFrame,
    event: Dict,
    lambda_outside: float,
    lambda_entry: float,
) -> Dict:
    """
    Analyze all particle bins for a single shower event, producing all four
    emission-rate variants (E1-E4; see module docstring for the parameter
    table).

    For each bin:
        1. Primary series: penetration factor (p) and peak time, shared by
           the outside, entry, and bathroom variants.
        2. Outside and entry variants (E1, E2): independent beta/E/Ct fits on
           the primary series, bounded by lambda_outside and lambda_entry.
        3. Bathroom variant (E4): reuses the outside variant's p, beta, and
           peak_time; refits only E and Ct with the larger volume.
        4. Adjusted variant (E3): independent beta/E/Ct fit on the adjusted
           (ratio-corrected) concentration series, using lambda_outside; only
           for events before ROOM_CUTOVER (single-monitor events).

    Per-bin result keys:
        bin{n}_p_mean, bin{n}_p_std, bin{n}_peak_time  (primary series,
            shared by outside/entry/bathroom)
        bin{n}_adjusted_p_mean, ..._p_std, ..._peak_time  (adjusted series,
            E3 only)
        bin{n}_{variant}_beta_other, ..._beta_other_raw_mean, ..._beta_other_std,
            ..._beta_other_r_squared, ..._beta_step, ..._c_steady_state
        bin{n}_{variant}_E_mean, ..._E_std, ..._E_total, ..._E_r_squared
        bin{n}_{variant}_emission_datetimes, ..._emission_predicted  (shower->peak)
        bin{n}_{variant}_decay_datetimes, ..._decay_predicted  (peak->deposition_end)
        bin{n}_{variant}_ct_datetimes, ..._ct_predicted  (full window)
        bin{n}_{variant}_skip_reason
        where {variant} is "outside", "entry", "adjusted", or "bathroom"

    Parameters:
        particle_data (pd.DataFrame): Primary concentration series (raw
            C_bed1(t) pre-cutover, fleet C_room(t)); feeds outside/entry/bathroom.
        particle_data_adjusted (pd.DataFrame): Ratio-corrected concentration
            series (C_adjusted_room(t) pre-cutover, fleet C_room(t)); feeds
            the adjusted (E3) variant only.
        event (Dict): Event timing information
        lambda_outside (float): Air change rate from the outside CO2 source (h-1)
        lambda_entry (float): Air change rate from the entry-zone CO2 source (h-1)

    Returns:
        Dict: Results for all bins, all four variants
    """
    results = {
        "event_number": event.get("event_number", 0),
        "test_name": event.get("test_name", ""),
        "config_key": event.get("config_key", ""),
        "water_temp": event.get("water_temp", ""),
        "shower_head": event.get("shower_head", "Standard"),
        "spray_pattern": event.get("spray_pattern"),
        "mannequin": event.get("mannequin", False),
        "door_position": event.get("door_position", ""),
        "planned_fan": event.get("planned_fan", ""),
        "time_of_day": event.get("time_of_day", ""),
        "fan_during_test": event.get("fan_during_test", False),
        "replicate_num": event.get("replicate_num", 0),
        "shower_on": event["shower_on"],
        "shower_off": event["shower_off"],
        "shower_duration_min": event.get("shower_duration_min", event.get("duration_min", 0)),
        "lambda_outside": lambda_outside,
        "lambda_entry": lambda_entry,
        "co2_event_idx": event.get("co2_event_idx", None),
        "flow_rate": event.get("flow_rate"),
    }

    time_of_day = event.get("time_of_day", "")
    lambda_values = {"outside": lambda_outside, "entry": lambda_entry}
    beta_window_start = event["shower_off"] + timedelta(hours=1)
    beta_window_end = event["deposition_end"]

    shower_on = event["shower_on"]
    pre_cutover = pd.notna(shower_on) and shower_on < ROOM_CUTOVER

    for bin_num in PARTICLE_BINS.keys():
        # ---- Primary series: p and peak_time, shared by outside/entry/bathroom ----
        p_result = calculate_penetration_factor(particle_data, shower_on, time_of_day, bin_num)
        results[f"bin{bin_num}_p_mean"] = p_result.get("p_mean", np.nan)
        results[f"bin{bin_num}_p_std"] = p_result.get("p_std", np.nan)

        p_mean = p_result.get("p_mean", np.nan)
        if np.isnan(p_mean):
            results[f"bin{bin_num}_peak_time"] = None
            for source in ("outside", "entry", "bathroom"):
                _fill_source_nan(results, bin_num, source, p_result.get("skip_reason", "Unknown"))
        else:
            peak_result = find_peak_time(
                particle_data, event["deposition_start"], event["deposition_end"], bin_num
            )
            peak_time = peak_result.get("peak_time")
            results[f"bin{bin_num}_peak_time"] = peak_time
            if peak_time is None:
                peak_time = event["shower_off"]

            beta_outside_val = np.nan
            for source in _LAMBDA_SOURCES:
                lambda_val = lambda_values[source]
                prefix = f"bin{bin_num}_{source}"

                if lambda_val is None or np.isnan(lambda_val):
                    _fill_source_nan(
                        results, bin_num, source, f"No lambda_{source} for this event"
                    )
                    continue

                beta_val = _compute_variant(
                    results,
                    particle_data,
                    event,
                    bin_num,
                    prefix,
                    p_mean,
                    peak_time,
                    lambda_val,
                    beta_window_start,
                    beta_window_end,
                    BEDROOM_VOLUME_M3,
                )
                if source == "outside":
                    beta_outside_val = beta_val

            _compute_bathroom_variant(
                results, particle_data, event, bin_num, p_mean, peak_time,
                lambda_outside, beta_outside_val,
            )

        # ---- Adjusted series (E3): independent fit, pre-cutover events only ----
        if not pre_cutover:
            _fill_adjusted_nan(
                results, bin_num,
                "Post-cutover: C_adjusted_room not computed (single-monitor variant only)",
            )
        elif lambda_outside is None or np.isnan(lambda_outside):
            _fill_adjusted_nan(results, bin_num, "No lambda_outside for this event")
        else:
            p_result_adj = calculate_penetration_factor(
                particle_data_adjusted, shower_on, time_of_day, bin_num
            )
            p_mean_adj = p_result_adj.get("p_mean", np.nan)
            results[f"bin{bin_num}_adjusted_p_mean"] = p_mean_adj
            results[f"bin{bin_num}_adjusted_p_std"] = p_result_adj.get("p_std", np.nan)

            if np.isnan(p_mean_adj):
                results[f"bin{bin_num}_adjusted_peak_time"] = None
                _fill_source_nan(
                    results, bin_num, "adjusted", p_result_adj.get("skip_reason", "Unknown")
                )
            else:
                peak_result_adj = find_peak_time(
                    particle_data_adjusted,
                    event["deposition_start"],
                    event["deposition_end"],
                    bin_num,
                )
                peak_time_adj = peak_result_adj.get("peak_time")
                results[f"bin{bin_num}_adjusted_peak_time"] = peak_time_adj
                if peak_time_adj is None:
                    peak_time_adj = event["shower_off"]

                _compute_variant(
                    results,
                    particle_data_adjusted,
                    event,
                    bin_num,
                    f"bin{bin_num}_adjusted",
                    p_mean_adj,
                    peak_time_adj,
                    lambda_outside,
                    beta_window_start,
                    beta_window_end,
                    BEDROOM_VOLUME_M3,
                )

    return results
