#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle Decay & Emission Analysis
===================================

This script analyzes particle concentration decay data from QuantAQ MODULAIR-PM
sensors to calculate particle penetration factors, other process rates, and shower
emission rates for the EPA Legionella study. The analysis uses a numerical
approach to solve the mass balance equation for twelve particle size bins,
producing four emission-rate variants per bin (E(t)1-E(t)4; see
src/particle_emission_variants.py).

Results characterize how shower-generated aerosols of different sizes penetrate,
deposit, and are emitted under controlled experimental conditions at varying water
temperatures and shower head configurations. Rather than commit to a single set of
inputs, every lambda-dependent metric (beta, E, Ct, R²) is computed across four
variants, each fixing a full combination of air-change-rate source, concentration
series, and control volume, per the report's Table 7 / Equation 10:

    E(t)1 (bin{n}_outside_*):  lambda_outside, primary concentration series,
        36.1086 m3 (bedroom). The primary reported value.
    E(t)2 (bin{n}_entry_*):    lambda_entry, primary series, 36.1086 m3.
        Brackets E(t)1 against the entry-zone CO2 air-change-rate estimate.
    E(t)3 (bin{n}_adjusted_*): lambda_outside, ratio-corrected concentration
        series, 36.1086 m3. Pre-cutover (single-monitor) events only.
    E(t)4 (bin{n}_bathroom_*): lambda_outside, primary series, 52.9903 m3
        (bedroom + bathroom). Reuses E(t)1's beta_other (volume-independent).

E(t)1 and E(t)2 (outside/entry) are the values consumed by the ~45 downstream
water-temp/spray-pattern/door/etc. boxplots and comparison figures
(src/plot_particle.py, plot_particle_boxplots.py, plot_comparison.py). E(t)3
and E(t)4 are written to all_results and compared against E(t)1 in
scripts/particle_emission_variant_figures.py.

Particle size bins analyzed (um):
    - Bin 0:  0.35-0.46
    - Bin 1:  0.46-0.66
    - Bin 2:  0.66-1.0
    - Bin 3:  1.0-1.3
    - Bin 4:  1.3-1.7
    - Bin 5:  1.7-2.3
    - Bin 6:  2.3-3.0
    - Bin 7:  3.0-4.0
    - Bin 8:  4.0-5.2
    - Bin 9:  5.2-6.5
    - Bin 10: 6.5-8.0
    - Bin 11: 8.0-10.0

Key Metrics Calculated:
    - p: Particle penetration factor (dimensionless, 0-1 range); one value for
      the primary concentration series (shared by E1/E2/E4), one for the
      adjusted series (E3 only)
    - beta_other: Effective other process loss rate (h-1); one value per bin
      per variant (E4 reuses E1's, since beta has no volume term)
    - E: Shower emission rate (particles/minute); one value per bin per variant
    - lambda_outside, lambda_entry: Air change rates from the CO2 decay analysis
      (co2_decay_analysis.py), used independently to bound E1 vs. E2 (h-1)

Analysis Features:
    - Numerical solution of time-dependent mass balance equation
    - Integration with CO2-derived air change rates
    - Per-bin analysis for size-dependent behavior (12 bins, 0.35-10.0 µm)
    - Statistical summaries across all shower events
    - Per-event decay-curve figures (this script); summary/boxplot/comparison
      figures (scripts/particle_emission_variant_figures.py, run after this script)
    - OUTDOOR_PM_EVENTS list: events for which outdoor PM concentrations are overlaid
      on the concentration panel (dotted lines)

Methodology:
    The mass balance equation for indoor particle concentration:
        V dC/dt = pQC_out - QC - beta_deposition CV + E
        dC/dt = p*lambda*C_out - lambda*C - beta_deposition*C + E/V

    1. Calculate penetration factor (p):
       - Use two averaging windows around each shower event (before and after)
       - For Night events:
           Before: 9pm (day before) to 2am (day of)
           After:  9am (day of) to 2pm (day of)
       - For Day events:
           Before: 9am (day of) to 2pm (day of)
           After:  9pm (day of) to 2am (next day)
       - p = C_inside / C_outside (averaged over each window, zeros excluded)
       - Final p = average of before and after window p values
       - Allowable range: 0-1 (values > 1 are capped at 1)
       - Computed independently on the primary series (E1/E2/E4) and the
         adjusted series (E3)

    2. Obtain air change rate (lambda):
       - Load lambda_outside and lambda_entry from CO2 decay analysis results
       - Steps 3-6 below run once per variant; every lambda-dependent result
         is stored per variant rather than averaged
       - Units: h-1

    3. Calculate other process rate (beta_other) when E=0:
       - Use 2-hour window after shower ends (DEPOSITION_WINDOW_HOURS)
       - Start time from peak concentration within the window to end of window
       - Solve numerically for each time step:
           beta = 1/dt - lambda - C_{t+1}/(C_t*dt) + p*lambda*(C_out,t/C_t)
       - Collect all estimates <= MAX_OTHER_PROCESS_RATE (no lower bound to
         avoid upward bias from excluding negative/noisy steps)
       - Apply 5th-95th percentile trim to remove extreme outliers symmetrically
       - Beta selected via R²-based single-step procedure (threshold 0.75):
         unclamped trimmed mean → keep if R² ≥ 0.75; otherwise beta = NaN
         (bin invalid, no Ct prediction)
       - Report selected beta as mean beta for the event/bin
       - Has no volume term: E4's beta_other equals E1's

    4. Calculate emission rate (E) from shower start to peak concentration:
       - Use shower ON to peak concentration time within analysis window
       - Solve numerically for E_t at each time step by rearranging the mass balance:
           (C_{t+1} - C_t)/dt = p*lambda*C_out,t - lambda*C_t - beta*C_t + E_t/V
           E_t/V = (C_{t+1} - C_t)/dt - p*lambda*C_out,t + lambda*C_t + beta*C_t
           E_t = V*(C_{t+1} - C_t)/dt - p*lambda*V*C_out,t + lambda*V*C_t + beta*V*C_t
       - V is 36.1086 m3 for E1/E2/E3, 52.9903 m3 for E4
       - Report E_mean and E_std from positive E_t values over the window
       - E_times and E_per_step (all steps including negative) stored for plotting

    5. Predict concentration Ct using forward Euler simulation:
       - Window: shower ON to 2 hours after shower OFF
       - Single continuous simulation using time-varying outdoor concentration:
           C_t(i+1) = C_t + dt*[p*lambda*C_out,t - C_t*(lambda + beta) + E_t/V]
       - E_t = E_mean from shower ON to peak_time, then E_t = 0
       - When E_mean is unavailable (emission calc failed), E_t = 0 throughout
         so a decay-only prediction is still generated for valid-beta bins
       - C_0 = measured bin concentration at shower ON
       - Returned as two continuous segments: emission phase (shower ON to peak)
         and decay phase (peak to deposition end); decay starts from predicted
         concentration at peak, not from the measured peak value
       - Plot both phases as a single continuous predicted Ct curve on figures
       - E_r_squared: R² of emission-phase forward Euler vs. measured concentration

    6. Calculate total emission (E_total) for each bin:
       - Area under the E_t vs. time curve using the trapezoidal rule:
           E_total = dt * sum[(E_t + E_t(i+1)) / 2]
       - Summed over all time steps from shower ON to peak concentration
       - Negative per-step contributions clipped to 0 before integration

Output Files:
    - particle_analysis_summary.xlsx: Multi-sheet workbook with:
        * all_results: Full results table (all metrics, all four variants,
          per event and bin)
        * p_penetration: Penetration factors per event and bin (includes test_name);
          primary series, one column per bin (E1/E2/E4 share this; E3's
          independent p is in all_results only, as bin{n}_adjusted_p_mean)
        * beta_deposition: Other process rates per event and bin (includes test_name);
          outside/entry only, columns doubled as bin{n}_outside_beta_other / bin{n}_entry_beta_other
        * beta_r_squared: R² of forward Euler decay simulation (includes test_name);
          doubled per source like beta_deposition
        * E_emission: Emission rates per event and bin (includes test_name); doubled per source
        * E_total_particles: Total emitted particle counts (E_total) per bin (includes
          test_name); doubled per source
        * E_r_squared: R² of forward Euler emission-phase simulation (includes test_name);
          doubled per source
        * peak_comparison: Measured vs. predicted concentration at peak_time and
          deposition_end for each bin, with percent difference (wide format,
          one row per event); predicted columns doubled per source
    - plots/event_figures/pm_decay/event_NN-YYYYYY_pm_decay.png: Individual event decay
      curves (four-panel): top panel shows measured concentrations and continuous
      predicted Ct (emission + decay phases) with decay R² in text box; three emission
      panels (bins 0–2, 3–6, 7–11) show per-step E_t, E_mean dashed lines, and emission
      R² annotation (first emission panel); shower markers use dotted lines; emission
      x-axes match concentration panel; emission y-axes clipped to 2nd-98th percentile;
      events in OUTDOOR_PM_EVENTS also overlay outdoor PM concentrations (dotted lines).
      Uses the E1 (outside) variant's result only.
    - plots/event_figures/excluded_events/: Figures for excluded events (duration-excluded
      showers saved here for reference rather than in type-specific subdirs)

    Summary/boxplot/comparison figures (penetration/deposition/emission bar charts,
    the fixed- and metric-axis boxplots, the shower-head and condition comparison
    families, and the new E1-vs-E3/E1-vs-E4 comparison figures) moved to
    scripts/particle_emission_variant_figures.py -- run it after this script.

Applications:
    - Characterizing size-resolved particle penetration and deposition in residential
      settings for Legionella exposure and aerosol transmission risk assessment
    - Quantifying shower-generated aerosol emission rates under varying water
      temperature and shower head configurations
    - Supporting development of building ventilation and filtration guidance for
      pathogen-carrying aerosol mitigation
    - Providing experimental data for validation of indoor air quality and
      airborne transmission models

Module Structure:
    - src/particle_calculations.py: Pure computation functions (p, beta, E, Ct)
    - src/particle_data_loader.py: Data loading and event identification
    - src/particle_room_correction.py: Builds the two "inside" concentration
      series (primary: raw C_bed1(t)/fleet C_room(t); adjusted: ratio-corrected
      C_adjusted_room(t)/fleet C_room(t))
    - src/particle_emission_variants.py: Per-event, per-bin analysis producing
      the four E(t)1-E(t)4 variants
    - scripts/particle_decay_analysis.py: Orchestration, per-event figures, and
      Excel save (this file)
    - scripts/particle_emission_variant_figures.py: Summary/boxplot/comparison
      figures, read from particle_analysis_summary.xlsx (run after this script)

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Date: 2026
Update log:
    2026-09-18  Room-concentration correction (src/particle_room_correction.py):
        "inside" data is now the position-weighted MODULAIR-PM fleet average
        (C_room) from 2026-06-03 to 2026-07-16, and a per-event, water-temp-
        bucketed ratio correction of MOD-PM-00195 before that date. Events are
        now loaded before particle data so the correction can use each event's
        window and water_temp.
    2026-09-22  Split into this pipeline script and
        scripts/particle_emission_variant_figures.py (summary/boxplot/comparison
        figures moved out). Added the E(t)3 (adjusted) and E(t)4 (bathroom)
        emission-rate variants (src/particle_emission_variants.py); the
        primary concentration series (E1/E2/E4) changed from the ratio-corrected
        series to the raw C_bed1(t)/fleet C_room(t) series -- the ratio-corrected
        series is now used for E3 only, per the report's Table 7. BEDROOM_VOLUME_M3
        precision updated to 36.1086 m3; BEDROOM_BATHROOM_VOLUME_M3 (52.9903 m3)
        added for E4.
    2026-09-22  Simplified beta selection (src/particle_calculations.py) from
        the R²-based four-step procedure (unclamped -> clamped >=0 -> 0 ->
        invalid; threshold 0.80) to a single step: unclamped trimmed mean,
        kept if R² >= 0.75, else beta = NaN (bin invalid). No clamping or
        forced beta=0 fallback.
"""

import sys
import warnings
from datetime import timedelta
from pathlib import Path
from typing import Optional

# Ensure stdout/stderr use UTF-8 on Windows (log files default to cp1252)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import src.sig_figs as sf  # noqa: E402
from src.data_paths import (  # noqa: E402
    get_data_root,
    get_event_figures_dir,
    get_event_figures_subdir,
)
from src.event_manager import (  # noqa: E402
    is_event_excluded,
    process_events_with_management,
)
from src.particle_calculations import (  # noqa: E402
    BEDROOM_BATHROOM_VOLUME_M3,
    BEDROOM_VOLUME_M3,
    DEPOSITION_WINDOW_HOURS,
    MAX_OTHER_PROCESS_RATE,
    MIN_POINTS_EMISSION,
    MIN_POINTS_OTHER_PROCESS,
    MIN_POINTS_PENETRATION,
    PARTICLE_BINS,
    TIME_STEP_MINUTES,
)
from src.particle_data_loader import (  # noqa: E402
    get_events_from_registry,
    identify_shower_events,
    load_and_merge_quantaq_data,
    load_co2_lambda_results,
    load_shower_log,
)
from src.particle_emission_variants import analyze_event_all_bins  # noqa: E402
from src.particle_room_correction import build_corrected_inside_data  # noqa: E402

# =============================================================================
# User-Configurable: Outdoor PM Overlay
# =============================================================================
# List event numbers for which the outdoor PM concentration should be overlaid
# on the pm_decay figure alongside the indoor concentration.  Useful for events
# where outdoor aerosol infiltration context is important (e.g. event 77,
# which occurred during an open-window night test).
# Set to an empty list [] to disable for all events.
OUTDOOR_PM_EVENTS: list = [75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85]

# =============================================================================
# Event Analysis Orchestration
# =============================================================================

# Sheets in particle_analysis_summary.xlsx that stay outside/entry-only
# (E1/E2), matching the ~45 downstream water-temp/spray-pattern/door/etc.
# boxplot and comparison figures in src/plot_particle.py, plot_particle_boxplots.py,
# and plot_comparison.py. The adjusted (E3) and bathroom (E4) variants
# computed by src.particle_emission_variants are added to all_results only;
# see scripts/particle_emission_variant_figures.py for their comparison figures.
_LAMBDA_SOURCES = ("outside", "entry")


# =============================================================================
# Main Analysis Pipeline
# =============================================================================


def run_particle_analysis(
    output_dir: Optional[Path] = None,
    generate_plots: bool = True,
    apply_sig_figs: bool = True,
) -> pd.DataFrame:
    """
    Run the complete particle decay and emission analysis.

    Parameters:
        output_dir (Path): Optional output directory (defaults to data_root/output)
        generate_plots (bool): If True, generate the per-event pm_decay figure for
            each event (summary/boxplot/comparison figures are generated separately
            by scripts/particle_emission_variant_figures.py)
        apply_sig_figs (bool): If True (default), round calculated float columns to
            SIG_FIGS_DATA significant figures before writing Excel output files and
            apply SIG_FIGS_FIGURE significant figures to figure annotations.
            Pass False (via --no-sig-figs) to preserve full floating-point precision.

    Returns:
        pd.DataFrame: DataFrame with analysis results for all events and bins
    """
    sf.set_enabled(apply_sig_figs)
    print("=" * 80)
    print("Particle Decay & Emission Analysis")
    print("Numerical Approach - Twelve Particle Size Bins, Four Emission Variants")
    print("=" * 80)
    print(f"Bedroom volume (E1/E2/E3): {BEDROOM_VOLUME_M3} m^3")
    print(f"Bedroom+bathroom volume (E4): {BEDROOM_BATHROOM_VOLUME_M3} m^3")
    print(f"Time step: {TIME_STEP_MINUTES} minute(s)")
    print("Penetration factor: averaged before/after windows (p capped at 1)")
    print(f"Deposition window: {DEPOSITION_WINDOW_HOURS} hour(s) after shower")
    print("Beta selection: R²-based single-step (unclamped trimmed mean → invalid; threshold 0.75)")
    print("\nValidation thresholds:")
    print(f"  Max other process rate (beta_other): {MAX_OTHER_PROCESS_RATE} h^-1")
    print(
        f"  Min data points: p={MIN_POINTS_PENETRATION}, beta_other={MIN_POINTS_OTHER_PROCESS}, E={MIN_POINTS_EMISSION}"
    )

    # Set output directory
    if output_dir is None:
        output_dir = get_data_root() / "output"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load events first: the room-concentration correction applied while loading
    # particle data (src.particle_room_correction) needs each pre-2026-06-03
    # event's shower_on/deposition_end window and water_temp bucket.
    events, co2_results, used_registry = get_events_from_registry(output_dir)

    if used_registry:
        print("  Using unified event registry for consistent event numbering")
    else:
        # Fall back to existing event management system
        print("\nNote: Registry not found. Using process_events_with_management().")
        print("  Run 'python scripts/event_registry.py' for unified numbering.\n")

        # Load shower log and identify events
        print("Loading shower log...")
        shower_log = load_shower_log()
        raw_events = identify_shower_events(shower_log)
        print(f"Found {len(raw_events)} raw shower events")

        # Load CO2 lambda results
        co2_results = load_co2_lambda_results()

        # Process events using the enhanced event management system
        print("\nProcessing events with event management system...")
        events, co2_events_processed, event_log = process_events_with_management(
            raw_events,
            [],  # CO2 events (will be loaded from co2_results)
            shower_log,
            co2_results,
            output_dir,
            create_synthetic=False,
        )

    # Load particle data: two parallel "inside" series (both merged with the
    # same "outside" data). particle_data is the primary series (raw C_bed1(t)
    # pre-cutover, fleet C_room(t)) feeding the outside/entry/bathroom
    # variants (E1/E2/E4); particle_data_adjusted is the ratio-corrected
    # series (C_adjusted_room(t) pre-cutover) feeding the adjusted variant
    # (E3) only. See src/particle_emission_variants.py.
    particle_data = load_and_merge_quantaq_data(events)
    particle_data_adjusted = load_and_merge_quantaq_data(
        events, inside_builder=build_corrected_inside_data
    )

    # Print event matching summary
    print("\nEvent Matching Summary:")
    matched_count = 0
    excluded_count = 0
    missing_lambda_count = 0

    for event in events:
        shower_time = event["shower_on"]

        # Check if excluded (time-based or duration-based)
        is_excluded_flag, exclusion_reason = is_event_excluded(shower_time)
        if not is_excluded_flag:
            is_excluded_flag = event.get("is_excluded", False)
            exclusion_reason = event.get("exclusion_reason", "")
        if is_excluded_flag:
            excluded_count += 1
            print(
                f"  Event {event.get('event_number', '?')} "
                f"({shower_time.strftime('%Y-%m-%d %H:%M')}): "
                f"EXCLUDED - {exclusion_reason}"
            )
            continue

        # Check if has at least one lambda source (outside or entry)
        lambda_outside_val = event.get("lambda_outside_mean", np.nan)
        lambda_entry_val = event.get("lambda_entry_mean", np.nan)
        if not (np.isnan(lambda_outside_val) and np.isnan(lambda_entry_val)):
            matched_count += 1
            co2_idx = event.get("co2_event_idx")
            if co2_idx is not None and co2_idx < len(co2_results):
                co2_time = co2_results.iloc[co2_idx]["injection_start"]
                print(
                    f"  {event.get('test_name', 'Event ' + str(event.get('event_number', '?')))} "
                    f"({shower_time.strftime('%m/%d %H:%M')}) "
                    f"-> CO2 {co2_idx + 1} ({co2_time.strftime('%H:%M')}), "
                    f"lambda_outside={lambda_outside_val:.4f}, "
                    f"lambda_entry={lambda_entry_val:.4f} h^-1"
                )
        else:
            missing_lambda_count += 1
            print(
                f"  {event.get('test_name', 'Event ' + str(event.get('event_number', '?')))} "
                f"({shower_time.strftime('%m/%d %H:%M')}): "
                f"No lambda value available"
            )

    print(
        f"\nTotal: {len(events)} events | Matched: {matched_count} | "
        f"Excluded: {excluded_count} | Missing lambda: {missing_lambda_count}"
    )

    # Analyze each event
    print("\nAnalyzing shower events...")
    results = []

    # Setup plot directory
    plot_dir = output_dir / "plots"
    pm_decay_dir = get_event_figures_subdir(output_dir, "pm_decay")
    if generate_plots:
        plot_dir.mkdir(exist_ok=True)
        pm_decay_dir.mkdir(parents=True, exist_ok=True)

    for event in events:
        event_num = event.get("event_number", 0)
        test_name = event.get("test_name", f"Event_{event_num}")
        shower_time = event["shower_on"]
        lambda_outside = event.get("lambda_outside_mean", np.nan)
        lambda_entry = event.get("lambda_entry_mean", np.nan)

        # Skip excluded events (time-based or duration-based)
        is_excluded_flag, exclusion_reason = is_event_excluded(shower_time)
        if not is_excluded_flag:
            is_excluded_flag = event.get("is_excluded", False)
            exclusion_reason = event.get("exclusion_reason", "")
        if is_excluded_flag:
            print(f"  {test_name}: Skipped (excluded: {exclusion_reason})")
            # Generate raw PM plot for excluded events that have a real event_number
            if generate_plots and event_num:
                try:
                    from src.plot_particle import plot_particle_decay_event
                    from src.plot_style import format_test_name_for_filename

                    excluded_dir = get_event_figures_subdir(output_dir, "excluded_events")
                    excluded_dir.mkdir(parents=True, exist_ok=True)

                    # Ensure deposition_end is set (fall back to shower_off + 2h)
                    # Use pd.isna() to catch both None and pd.NaT from the registry
                    if pd.isna(event.get("deposition_end")):
                        event = dict(event)
                        event["deposition_end"] = event["shower_off"] + timedelta(hours=2)

                    empty_result = {}
                    formatted_name = format_test_name_for_filename(test_name)
                    plot_path = (
                        excluded_dir / f"event_{event_num:02d}-{formatted_name}_pm_decay.png"
                    )
                    plot_particle_decay_event(
                        particle_data=particle_data,
                        event=event,
                        particle_bins=PARTICLE_BINS,
                        result=empty_result,
                        output_path=plot_path,
                        event_number=event_num,
                        test_name=test_name,
                    )
                except Exception as e:
                    print(f"    Warning: Failed to generate excluded plot for {test_name}: {e}")
            continue

        # Skip events with no lambda from either source
        if np.isnan(lambda_outside) and np.isnan(lambda_entry):
            print(f"  {test_name}: Skipped (no lambda from CO2 analysis)")
            continue

        print(
            f"  {test_name} ({shower_time.strftime('%m/%d %H:%M')}): "
            f"lambda_outside={lambda_outside:.4f}, lambda_entry={lambda_entry:.4f} h^-1"
        )

        result = analyze_event_all_bins(
            particle_data, particle_data_adjusted, event, lambda_outside, lambda_entry
        )
        results.append(result)

        # Print summary for this event with detailed skip reasons. A bin
        # counts as valid if either lambda source produced an emission rate.
        valid_bins = 0
        skipped_bins = []
        for bin_num in PARTICLE_BINS.keys():
            source_e_means = [
                result.get(f"bin{bin_num}_{source}_E_mean", np.nan) for source in _LAMBDA_SOURCES
            ]
            if any(not np.isnan(v) for v in source_e_means):
                valid_bins += 1
            else:
                skip_reason = result.get(
                    f"bin{bin_num}_outside_skip_reason",
                    result.get(f"bin{bin_num}_entry_skip_reason", "Unknown"),
                )
                skipped_bins.append((bin_num, skip_reason))

        print(f"    Successfully analyzed {valid_bins}/{len(PARTICLE_BINS)} bins")

        # Print skip reasons for failed bins (up to 3 for brevity)
        if skipped_bins and valid_bins < len(PARTICLE_BINS):
            for bin_num, reason in skipped_bins[:3]:
                bin_name = PARTICLE_BINS[bin_num]["name"]
                # Truncate long reasons
                if len(reason) > 80:
                    reason = reason[:77] + "..."
                print(f"      Bin {bin_num} ({bin_name} um): {reason}")
            if len(skipped_bins) > 3:
                print(f"      ... and {len(skipped_bins) - 3} more bins skipped")

        # Generate individual event plot if enabled (all bins on one plot)
        if generate_plots and valid_bins > 0:
            try:
                from src.plot_particle import plot_particle_decay_event
                from src.plot_style import format_test_name_for_filename

                # Format filename: event_01-0114_hw_morning_pm_decay.png
                formatted_name = format_test_name_for_filename(test_name)
                plot_path = pm_decay_dir / f"event_{event_num:02d}-{formatted_name}_pm_decay.png"
                plot_particle_decay_event(
                    particle_data=particle_data,
                    event=event,
                    particle_bins=PARTICLE_BINS,
                    result=result,
                    output_path=plot_path,
                    event_number=event_num,
                    test_name=test_name,
                    show_outdoor=event_num in OUTDOOR_PM_EVENTS,
                )
            except ImportError:
                pass  # Already warned about missing plot module
            except Exception as e:
                print(f"    Warning: Failed to generate plot for {test_name}: {e}")

    # Create results DataFrame
    results_df = pd.DataFrame(results)

    # Print overall statistics
    _print_overall_summary(results_df, results)

    # Save results
    _save_results(results_df, output_dir)

    print(
        "\nSummary/boxplot/comparison figures moved to "
        "scripts/particle_emission_variant_figures.py -- run it next."
    )

    return results_df


def _print_overall_summary(results_df: pd.DataFrame, results: list) -> None:
    """Print overall statistics summary to console."""
    print("\n" + "=" * 80)
    print("Overall Results Summary")
    print("=" * 80)

    if results_df.empty:
        print("\nNo events were analyzed (all skipped due to missing lambda or exclusions).")
        return

    for bin_num, bin_info in PARTICLE_BINS.items():
        bin_name = bin_info["name"]
        p_col = f"bin{bin_num}_p_mean"
        valid_p = (
            results_df[p_col].dropna() if p_col in results_df.columns else pd.Series(dtype=float)
        )

        print(f"\nBin {bin_num} ({bin_name} um):")
        if len(valid_p) > 0:
            print(f"  p (penetration):     {valid_p.mean():.3f} +/- {valid_p.std():.3f}")

        for source in _LAMBDA_SOURCES:
            beta_col = f"bin{bin_num}_{source}_beta_other"
            E_col = f"bin{bin_num}_{source}_E_mean"
            if beta_col not in results_df.columns or E_col not in results_df.columns:
                continue
            valid_beta = results_df[beta_col].dropna()
            valid_E = results_df[E_col].dropna()
            if len(valid_beta) > 0:
                print(f"  beta ({source}): {valid_beta.mean():.3f} +/- {valid_beta.std():.3f} h^-1")
            if len(valid_E) > 0:
                print(f"  E ({source}): {valid_E.mean():.2e} +/- {valid_E.std():.2e} #/min")
            print(f"  Valid events ({source}): {len(valid_E)}/{len(results)}")

        # E3 (adjusted) and E4 (bathroom): brief valid-event counts only, since
        # the detailed per-bin stats above are for the reported E1/E2 values.
        for variant in ("adjusted", "bathroom"):
            E_col = f"bin{bin_num}_{variant}_E_mean"
            if E_col not in results_df.columns:
                continue
            valid_E = results_df[E_col].dropna()
            print(f"  Valid events ({variant}): {len(valid_E)}/{len(results)}")


def _save_results(results_df: pd.DataFrame, output_dir: Path) -> None:
    """Save analysis results to Excel workbook with one sheet per metric.

    Every lambda-dependent metric is written twice per bin, once per source
    (bin{n}_outside_..., bin{n}_entry_...), bounding the result instead of
    blending outside/entry into one average. Penetration factor (bin{n}_p_mean)
    and peak_time are lambda-independent and stay single columns.

    all_results additionally carries the adjusted (E3, bin{n}_adjusted_...,
    including its own independent bin{n}_adjusted_p_mean/p_std/peak_time) and
    bathroom (E4, bin{n}_bathroom_...) emission-rate variants from
    src.particle_emission_variants; see scripts/particle_emission_variant_figures.py
    for their comparison figures. The per-metric sheets below (p_penetration
    through peak_comparison) stay outside/entry-only to match the ~45
    downstream water-temp/spray-pattern/etc. boxplots, which are unaffected.

    Sheets written:
        all_results       - Full results table (all metrics, all four
                            variants, per event and bin)
        p_penetration     - Penetration factors (outside/entry primary series
                            only)
        beta_deposition   - Other process rates (h⁻¹) per source: beta_other
                            (R²-selected, NaN if invalid) and
                            beta_other_raw_mean (unclamped trimmed mean,
                            identical to beta_other when valid) per bin
        beta_r_squared    - R² of forward Euler decay simulation, per source
        E_emission        - Mean emission rates (#/min), per source
        E_total_particles - Total emitted particle counts per bin (E_total, #),
                            per source
        E_r_squared       - R² of forward Euler emission-phase simulation, per
                            source
        peak_comparison   - Measured vs. predicted at peak_time and
                            deposition_end, per source
    """
    output_file = output_dir / "particle_analysis_summary.xlsx"

    if results_df.empty:
        print(f"\nNo results to save - skipping {output_file}")
        return

    # Create column rename mapping for units
    column_rename = {
        "shower_duration_min": "shower_duration (min)",
        "lambda_outside": "lambda_outside (h-1)",
        "lambda_entry": "lambda_entry (h-1)",
    }
    for bin_num in PARTICLE_BINS.keys():
        column_rename[f"bin{bin_num}_p_mean"] = f"bin{bin_num}_p_mean (-)"
        column_rename[f"bin{bin_num}_p_std"] = f"bin{bin_num}_p_std (-)"
        # adjusted (E3) uses its own independent p (different concentration
        # series than outside/entry/bathroom's shared primary-series p).
        column_rename[f"bin{bin_num}_adjusted_p_mean"] = f"bin{bin_num}_adjusted_p_mean (-)"
        column_rename[f"bin{bin_num}_adjusted_p_std"] = f"bin{bin_num}_adjusted_p_std (-)"
        for source in _LAMBDA_SOURCES + ("adjusted", "bathroom"):
            prefix = f"bin{bin_num}_{source}"
            column_rename[f"{prefix}_beta_other"] = f"{prefix}_beta_other (h-1)"
            column_rename[f"{prefix}_beta_other_raw_mean"] = f"{prefix}_beta_other_raw_mean (h-1)"
            column_rename[f"{prefix}_beta_other_std"] = f"{prefix}_beta_other_std (h-1)"
            column_rename[f"{prefix}_E_mean"] = f"{prefix}_E_mean (#/min)"
            column_rename[f"{prefix}_E_std"] = f"{prefix}_E_std (#/min)"
            column_rename[f"{prefix}_E_total"] = f"{prefix}_E_total (#)"

    results_df_export = results_df.rename(columns=column_rename)
    results_df_export = sf.apply_sig_figs_to_df(results_df_export)

    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
        # Main results
        results_df_export.to_excel(writer, sheet_name="all_results", index=False)

        # Shared ID columns present on every sheet
        id_cols = ["event_number", "test_name", "shower_on"]

        # Penetration factor is lambda-independent: one column per bin, no
        # source suffix.
        p_cols = id_cols + [f"bin{i}_p_mean (-)" for i in PARTICLE_BINS.keys()]

        # Every other metric is doubled: one column per bin per source.
        beta_cols = id_cols + [
            col
            for i in PARTICLE_BINS.keys()
            for source in _LAMBDA_SOURCES
            for col in (
                f"bin{i}_{source}_beta_other (h-1)",
                f"bin{i}_{source}_beta_other_raw_mean (h-1)",
                f"bin{i}_{source}_beta_step",
            )
        ]
        beta_r2_cols = id_cols + [
            f"bin{i}_{source}_beta_other_r_squared"
            for i in PARTICLE_BINS.keys()
            for source in _LAMBDA_SOURCES
        ]
        E_cols = id_cols + [
            f"bin{i}_{source}_E_mean (#/min)"
            for i in PARTICLE_BINS.keys()
            for source in _LAMBDA_SOURCES
        ]
        E_total_cols = id_cols + [
            f"bin{i}_{source}_E_total (#)"
            for i in PARTICLE_BINS.keys()
            for source in _LAMBDA_SOURCES
        ]
        E_total_cols = [c for c in E_total_cols if c in results_df_export.columns]

        E_r2_cols = id_cols + [
            f"bin{i}_{source}_E_r_squared"
            for i in PARTICLE_BINS.keys()
            for source in _LAMBDA_SOURCES
        ]
        E_r2_cols = [c for c in E_r2_cols if c in results_df_export.columns]

        results_df_export[p_cols].to_excel(writer, sheet_name="p_penetration", index=False)
        results_df_export[beta_cols].to_excel(writer, sheet_name="beta_deposition", index=False)
        results_df_export[beta_r2_cols].to_excel(writer, sheet_name="beta_r_squared", index=False)
        results_df_export[E_cols].to_excel(writer, sheet_name="E_emission", index=False)
        results_df_export[E_total_cols].to_excel(
            writer, sheet_name="E_total_particles", index=False
        )
        if E_r2_cols:
            results_df_export[E_r2_cols].to_excel(writer, sheet_name="E_r_squared", index=False)

        # Peak comparison sheet: measured vs. predicted at peak_time and deposition_end.
        # Wide format: one row per event, bins x source as column groups.
        peak_df = results_df[["event_number", "test_name"]].copy()
        for bin_num in PARTICLE_BINS.keys():
            for source in _LAMBDA_SOURCES:
                prefix = f"bin{bin_num}_{source}"
                meas_pk = f"{prefix}_peak_measured"
                pred_pk = f"{prefix}_peak_predicted"
                meas_de = f"{prefix}_deposition_end_measured"
                pred_de = f"{prefix}_deposition_end_predicted"

                if meas_pk not in results_df.columns:
                    continue

                peak_df[f"{prefix}_peak_measured (#/cm3)"] = results_df[meas_pk].values
                peak_df[f"{prefix}_peak_predicted (#/cm3)"] = results_df[pred_pk].values
                with np.errstate(invalid="ignore", divide="ignore"):
                    pct_pk = (
                        (results_df[pred_pk] - results_df[meas_pk]) / results_df[meas_pk] * 100.0
                    )
                peak_df[f"{prefix}_peak_pct_diff (%)"] = pct_pk.values

                peak_df[f"{prefix}_deposition_end_measured (#/cm3)"] = results_df[meas_de].values
                peak_df[f"{prefix}_deposition_end_predicted (#/cm3)"] = results_df[pred_de].values
                with np.errstate(invalid="ignore", divide="ignore"):
                    pct_de = (
                        (results_df[pred_de] - results_df[meas_de]) / results_df[meas_de] * 100.0
                    )
                peak_df[f"{prefix}_deposition_end_pct_diff (%)"] = pct_de.values

        peak_df = sf.apply_sig_figs_to_df(peak_df)
        peak_df.to_excel(writer, sheet_name="peak_comparison", index=False)

    print(f"\nResults saved to: {output_file}")


def main():
    """Main entry point for command-line usage."""
    import argparse

    parser = argparse.ArgumentParser(description="Particle Decay & Emission Analysis")
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for results (default: data_root/output)",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Disable per-event pm_decay figure generation (summary/boxplot/"
        "comparison figures are in scripts/particle_emission_variant_figures.py)",
    )
    parser.add_argument(
        "--no-sig-figs",
        action="store_true",
        help="Disable significant figure rounding on output data and figure annotations "
        "(default: 3 sig figs for files, 2 sig figs for figures)",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else None

    run_particle_analysis(
        output_dir=output_dir,
        generate_plots=not args.no_plot,
        apply_sig_figs=not args.no_sig_figs,
    )


if __name__ == "__main__":
    main()
