# Data Analysis

> **Author note:** Equation numbers and section cross-references require PI review before
> submission. Equation numbering uses the format (3-N) for Section 3; renumber when
> assembled into the full report. Items marked `[CHECK: ...]` are unresolved conflicts
> between this repository and the draft report; each names both values so the correct one
> can be confirmed rather than guessed.
>
> The analysis pipeline changed substantially between June and September 2026. Section 3.4
> describes the current aerosol analysis, built on the MODULAIR-PM fleet. Section 3.9
> records what the draft report describes that no script in this repository implements.

---

## 3.1 Overview

Raw sensor data, shower event logs, and CO₂ injection logs are processed through a Python pipeline (Python 3.14.2; NumPy 2.4.1, SciPy 1.17.0, pandas 2.3.3, Bokeh for interactive figures). All scripts and source modules are archived in the companion data and code repository (DOI: 10.18434/mds2-4153). Four analysis domains are treated:

1. **Air change rate (CO₂ decay).** The air exchange rate between Bedroom #1 and its surroundings is determined from the exponential decay of injected CO₂ after each event (Section 3.3).
2. **Room uniformity and the room-average concentration.** The MODULAIR-PM fleet installed 2026-06-04 establishes whether the bedroom can be treated as a single well-mixed zone, and produces the position-weighted room average C_room and the correction that transfers it to the single-monitor record (Sections 3.4.1 to 3.4.3).
3. **Aerosol mass balance.** Penetration factor, aerosol loss rate, and shower emission rate are estimated per event and per size bin from the room-average concentration (Sections 3.4.4 to 3.4.7).
4. **Environmental conditions.** Relative humidity, temperature, and wind are summarized over pre- and post-shower windows and over the onset and decay windows (Section 3.5).

All four are driven by a unified event registry (Section 3.2) that assigns a consistent event number, test name, and metadata record to every shower event.

### 3.1.1 Which Scripts Are Current

The aerosol analysis was rebuilt around the MODULAIR-PM fleet between June and September 2026. Not every script in the repository is part of the current pipeline.

**Current, and the basis for the draft report:**

| Script | Role |
|---|---|
| `scripts/event_registry.py` | Event registry, naming, exclusions |
| `scripts/co2_decay_analysis.py` | Air change rate λ |
| `src/moduair_loader.py` | MODULAIR-PM fleet chunk loader |
| `scripts/moduair_event_peak_times.py` | Per-monitor delta peak times (room uniformity) |
| `scripts/moduair_bin2_timeseries.py` | Per-bin fleet time series with C_room overlay |
| `scripts/moduair_correction_factor.py` | Inter-sensor ratios and the 195/813 Deming fit |
| `scripts/moduair_cave_ratio.py` | C_room, and the C_bed1 to C_room comparison |
| `scripts/export_config_timeseries_fleet.py` | Per-sensor and C_room config time-series exports |
| `scripts/rh_temp_other_analysis.py` | RH, temperature, and wind summaries |
| `scripts/hobo_onset_decay_figures.py` | HOBO onset/decay and pre/post figures |

**Superseded.** `scripts/particle_decay_analysis.py` and the calculation functions in `src/particle_calculations.py` implement the earlier single-sensor mass balance: one air change rate (λ_average), one control volume (`BEDROOM_VOLUME_M3 = 36.1`), a four-step R²-based acceptance hierarchy for β, and one emission rate. That path does not use C_room and is not what informs the draft report. It is documented in Section 3.8 because its outputs are still referenced by older figures.

[CHECK: `run_analysis_workflow.py` still lists `scripts/particle_decay_analysis.py` as the final pipeline step, so a full workflow run regenerates superseded outputs alongside current ones. Confirm whether that step should be removed, or the script updated to consume C_room.]

---

## 3.2 Event Identification and Registry

### 3.2.1 Shower Event Detection

Shower start and stop times are recorded by the automated control system in `shower_log_file.csv`. Each row is a discrete state change: shower on, shower off, bathroom fan on, or bathroom fan off, with a high-resolution timestamp. Consecutive on/off pairs per device are parsed into individual events characterized by start time, stop time, and duration. Only events beginning on or after 2026-01-15 15:00 are retained.

Analysis events run at 03:00 and 15:00 (Section 2.3.2). `scripts/moduair_event_peak_times.py` applies this directly, keeping only ON transitions within ±5 minutes of those two times so that off-hour water temperature verification runs are dropped before any fleet figure is built.

### 3.2.2 CO₂ Event Matching

CO₂ injection events are logged independently in `CO2_log_file.csv`, which records valve state changes and mixing fan activations. Each injection event is matched to the nearest shower event within a bidirectional tolerance of ±10 minutes, one to one. Unmatched shower events receive a synthetic CO₂ placeholder record with a duration inferred from neighboring matched events, so the registry stays complete and consistently numbered regardless of CO₂ data availability.

### 3.2.3 Event Numbering and Naming

Matched shower and CO₂ pairs receive sequential integer event numbers from 1. Each event also receives a structured test name (Section 2.2.6) based on the configuration active at its start time, resolved from the time-stamped transition tables in `src/event_manager.py`. The air change rate λ from the CO₂ decay analysis is stored in the registry for every event with a valid CO₂ match.

[CHECK: events relying on a synthetic CO₂ record are flagged, but whether they carry a λ forward from the nearest valid event of the same configuration, or are dropped from aerosol analysis entirely, is not stated in the code comments. Confirm the intended behavior and make the code match.]

### 3.2.4 Exclusion Criteria

Excluded events stay in the registry with `is_excluded = True`, a null event number, and an `exclusion_reason`. Exclusions are applied before any analysis step and excluded events do not appear in summary outputs.

**Duration-based exclusion.** Any shower event measuring outside 9 min 55 s to 10 min 5 s is treated as a manual water temperature verification run or a control-system test and excluded, with `exclusion_reason = "Water temperature testing (duration: X.X min)"`.

**Predefined individual event exclusions.** Four events are excluded for documented confounding activity:

| Date and time | Reason |
|---|---|
| 2026-01-22 15:00 | Tour in house during test |
| 2026-01-29 15:00 | People in house |
| 2026-05-13 15:00 | LVP flooring installation |
| 2026-05-21 15:00 | Bathroom flooring removal |

The last two are also permanent changes to the interior surfaces of the test rooms, not one-day disturbances. See the note at the end of Section 2.8.

**Predefined date range exclusions.** Four ranges are excluded:

| Excluded range | Reason |
|---|---|
| 2026-02-11 08:00 to 2026-02-13 11:00 | Conflicting log entries during a water temperature configuration change |
| 2026-03-08 00:00 to 2026-03-10 00:00 | Daylight saving clock change, instrument data misalignment, and elevated bedroom RH of uncertain origin |
| 2026-03-14 00:00 to 2026-03-15 12:00 | CO₂ injection system failure |
| 2026-05-11 00:00 to 2026-05-15 10:00 | CO₂ injection system failure |

**Aerosol analysis-specific exclusion.** Events whose CO₂ decay regression R² falls below 0.65 are excluded from aerosol analysis; the CO₂ result itself is retained. An unreliable λ propagates into the penetration factor, loss rate, and emission rate.

[CHECK: the threshold is 0.65 in code (`scripts/event_registry.py`, line 992). Three docstrings in that same file still state 0.75, including the printed exclusion message on line 997, and the draft report states 0.75. Confirm the intended threshold, then correct whichever of the three does not match.]

---

## 3.3 Air Change Rate Determination (CO₂ Decay Analysis)

### 3.3.1 Tracer Gas Decay Model

Bedroom #1 is treated as a single well-mixed zone. After the mixing fan stops and before the shower begins there are no internal CO₂ sources, and the CO₂ mass balance on the bedroom reduces to:

$$V \frac{dC}{dt} = Q\,C_{bg} - Q\,C \quad\Longrightarrow\quad \frac{dC}{dt} = \lambda \left( C_{bg} - C \right) \tag{3-1}$$

where $C(t)$ [ppm] is the bedroom CO₂ concentration, $\lambda = Q/V$ [h⁻¹] is the air change rate, and $C_{bg}$ [ppm] is a time-constant background concentration representing the mixed composition of air entering the bedroom. Three assumptions underlie Eq. (3-1): CO₂ outside and in the entryway does not change appreciably during the decay, airflow into and out of the room is balanced, and there are no interior CO₂ sources or sinks.

Integrating with $C(0) = C_0$ and rearranging into a form suitable for linear regression:

$$y(t) \equiv -\ln\!\left[\frac{C(t) - C_{bg}}{C_0 - C_{bg}}\right] = \lambda \cdot t \tag{3-2}$$

A linear regression of $y$ against $t$, with the intercept forced through the origin, gives $\lambda$ as the slope:

$$\hat{\lambda} = \frac{\displaystyle\sum_i t_i \cdot y_i}{\displaystyle\sum_i t_i^2} \tag{3-3}$$

### 3.3.2 Source Concentration Methods

$C_{bg}$ is not directly measured. Air enters the bedroom through the gap under the bedroom door, which draws from the entryway, and through envelope cracks, which draw from outdoors. Both contribute in unknown proportion, so three values are computed for every event to bracket the result:

| Method | $C_{bg}$ | Symbol |
|---|---|---|
| Outside only | $C_{out}$ | $\lambda_{out}$ |
| Entry only | $C_{ent}$ | $\lambda_{ent}$ |
| Average | $0.5\,C_{out} + 0.5\,C_{ent}$ | $\lambda_{avg}$ |

Each is averaged over the full decay window to give a single scalar for the regression. All three values, with their R², are written to `co2_lambda_summary.csv` and stored in the event registry.

The draft report uses $\lambda_{out}$ and $\lambda_{ent}$ as a bracket on the true air change rate and carries both separately through the mass balance. The superseded single-sensor pipeline used $\lambda_{avg}$ alone.

### 3.3.3 Data Preprocessing

Aranet4 data from Bedroom, Entry, and Outside are loaded from manufacturer-exported Excel files, parsed, resampled to a regular 1-minute grid by linear interpolation, then smoothed with a 6-minute centered rolling average before the regression. Bathroom Aranet4 data appear on event concentration plots for spatial context and are used in no λ calculation.

### 3.3.4 Regression Window

The decay window runs from shower_on + 10 min to shower_on + 2 h 10 min. Because every shower is exactly 10 minutes, this is identical to the window the draft report describes as shower-off to 2 hours after shower-off. The two descriptions agree; only the anchor differs.

The 10-minute offset from shower onset lets the airflow transient from shower start settle before the fit begins. $C_0$ is the first data point at window start, and $C_{bg}$ is computed once over the full 2-hour window.

A minimum initial concentration excess of 50 ppm is required for a valid regression.

[CHECK: the code applies the 50 ppm criterion at the start of the window, as $C_0 - C_{bg} \geq 50$ ppm. The draft report states the criterion at the end, as the bedroom concentration at $t$ = 2 h being less than 50 ppm above the outside or entryway average. These reject different events. Confirm which was applied to the reported λ values.]

An optional mode (`--entry-stop`) truncates the window once the bedroom concentration decays to within 100 ppm of the mean of entry and outside, preventing low-signal points from influencing the fit. It was not enabled for the primary analysis.

### 3.3.5 Uncertainty Considerations

1. **Source concentration assumption.** The spread between $\lambda_{out}$ and $\lambda_{ent}$ is the dominant inter-method uncertainty and is reported as a bracket rather than collapsed to a single value.
2. **Statistical fitting uncertainty.** Standard error of the regression slope, from sensor noise and short-term concentration fluctuation.
3. **Well-mixed zone assumption.** The single-zone model assumes spatial uniformity of bedroom CO₂ after the mixing fan stops. The 10-minute offset at window start mitigates residual gradients but does not eliminate them. Note that Section 3.4.1 shows the bedroom is *not* uniform for aerosol during the first hour after a warm shower; CO₂ is injected and mixed 15 minutes before the shower and decays over a 2-hour window, which is a different and more favorable mixing problem, but the assumption is not automatically safe.
4. **CO₂ sensor accuracy.** ±50 ppm ±3 % of reading. At the 1500 ppm to 2000 ppm injection peak, approximately ±95 ppm to ±110 ppm.
5. **Ambient variability.** Outdoor wind speed and direction drive infiltration. The reported λ is a 2-hour time average, not an instantaneous value.

---

## 3.4 Aerosol Analysis

All aerosol analyses run independently on each of the 12 OPC-N3 size bins in Table 2-3 (0.35 µm to 10.0 µm, indexed 0 to 11). Concentrations are number per cm³. Bin 2 (0.66 µm to 1.0 µm) is used for illustrative figures.

The mass balance in Section 3.4.4 requires that the bedroom be representable as a single uniform concentration. Sections 3.4.1 to 3.4.3 test that requirement, and the answer determines what concentration series the mass balance is given.

### 3.4.1 Room Uniformity Assessment

From 2026-06-04 the bedroom carried 10 monitors, and from 2026-07-08 it carried 14 (Section 2.5.1). Two analyses use that period.

**Per-bin concentration time series.** `scripts/moduair_bin2_timeseries.py` plots each bin's concentration for the co-located bedroom monitors over 2026-06-04 to 2026-07-16 as one interactive Bokeh figure per bin, one trace per monitor on a shared datetime axis, with the position-weighted C_room overlaid as a black trace. MOD-PM-00785 (outdoor) is excluded, and MOD-PM-00401 is dropped because it reads exactly zero for a large share of the window. Traces break at genuine data gaps so that missing periods are visible rather than interpolated across.

**Delta peak time.** `scripts/moduair_event_peak_times.py` computes, for every event and every monitor, the elapsed time from shower on to the maximum of the summed bins 0 to 11 within [shower_on, shower_off + 2 h]. Monitors installed late are gated to their install time, so events before a monitor went live are blank for that monitor only rather than being filled with a spurious value. The output is a CSV of event by monitor in minutes and an interactive figure.

The uniformity result is stratified by water temperature and is the reason a single-monitor concentration cannot be used directly:

- At 38 °C to 41 °C, the two High monitors peak approximately 20 min after shower on, the Mid monitors at 30 min to 40 min, and the Low monitors at 40 min to 120 min. Peak concentration at High 1, under the ceiling pitch, reaches roughly 1.5 times the peak at the Middle Bathroom monitor less than 1 m from the shower head, and reaches it sooner. A thermal plume carries aerosol through the doorway above 1.1 m, bypassing the bathroom monitor, then spills down the walls. Thermal stratification then prevents vertical mixing for the rest of the analysis window.
- At 23 °C the pattern inverts. The bathroom monitor peaks first and highest, the second monitor to respond is a Low monitor rather than a High one, bedroom peaks are similar across monitors, and decay is uniform with no second rise. There is little evidence of a thermal plume when water temperature is near air temperature.
- At 51 °C the plume is stronger and longer lived. Concentrations peak and fall at one monitor before peaking at another, some monitors peak twice, and the larger bins above 2.3 µm carry much higher concentrations than at the other two temperatures.

No single monitor represents the room average at any water temperature during the first hour.

### 3.4.2 MOD-PM-00195 Campaign Correction

Installing the fleet in June 2026 revealed that the historical MOD-PM-00195 record was biased. MOD-PM-00813 was placed alongside it at the Bed position, and `scripts/moduair_correction_factor.py` builds the co-located correlation over 2026-06-04 to 2026-07-16.

Three diagnostic ratios are computed per bin on a shared time base:

| Ratio | Definition |
|---|---|
| 813 | 195 / 813 |
| quad | 195 / mean(515, 465, 943, 516) |
| others | 195 / mean(remaining bedroom monitors, excluding 195, 785, 813, and the quad set) |

A value near 1.0 means 195 agrees with the reference; a sustained offset is a candidate multiplicative correction.

The correction itself is a per-bin orthogonal-distance (Deming) regression with 195 on the x-axis and 813 on the y-axis:

$$C_{195,\text{corrected}} = m \cdot C_{195,\text{measured}} + b \tag{3-4}$$

Deming rather than ordinary least squares, because OLS assumes the x variable is error free and would bias the slope. Both instruments are the same MODULAIR-PM model, so equal x and y error variance is assumed (delta = 1). Per-bin slope, intercept, and their standard errors are written to `moduair_correction_195_813_fit_jun04_jul16.csv`. The script fits and reports only; it does not modify the historical record.

[CHECK: the fit is produced but no script applies it. Confirm whether the reported single-monitor results have the 195 correction applied, and if so, where that application happens.]

### 3.4.3 Position-Weighted Room Average and Transfer to the Single-Monitor Record

`scripts/moduair_cave_ratio.py` implements the room average. The bedroom is divided into four elevation zones, with each zone's upper and lower bound set at half the distance between the mean monitor heights of adjacent zones:

| Zone | Monitors | High (m) | Low (m) | Weight |
|---|---|---|---|---|
| High | 402, 816 | 2.34 | 1.59 | 0.32 |
| Mid | 814, 554, 942, 401, 815 | 1.59 | 0.94 | 0.28 |
| Bed | 195, 813 | 0.94 | 0.47 | 0.20 |
| Low | 515, 465, 943, 516, 467 | 0.47 | 0.00 | 0.20 |

The High zone's upper bound of 2.34 m is the average height of the sloping ceiling. MOD-PM-00555 is deliberately excluded from the Mid zone because it is in the bathroom, not the bedroom. Weights are the fraction of total room height each zone spans. Each zone average is taken over its working monitors, and:

$$C_{room} = 0.32\,C_{high} + 0.28\,C_{mid} + 0.20\,C_{bed} + 0.20\,C_{low} \tag{3-5}$$

**Reporting-monitor filter.** C_room is computed only at 1-minute timestamps where more than eight fleet monitors report a finite, positive reading (`MIN_QUANTS = 8`, applied as strictly greater than). Minutes below that threshold are dropped rather than averaged over a partial fleet. For the aligned per-event curves the filter is applied per timestamp before averaging across events.

**Comparison of C_bed1 to C_room.** Concentrations are the raw `opc_bin{N}` values on a shared 1-minute grid with no smoothing. Three comparisons are produced per bin:

1. A scatter of C_bed1 against C_room over all qualifying 1-minute samples in the window, with a Deming fit and a 1:1 line, and a per-bin fit table of slope, intercept, r², and n.
2. The ratio of the event-averaged C_bed1 to the event-averaged C_room against time, aligned on minute index with 0 at shower on, from 15 min before shower on through the 2-hour deposition window. Events are grouped by registry water temperature code into W38-W41, W49, and W24. A ±1 standard deviation band is drawn only at minutes where at least three events contribute.
3. The same scatter restricted to W38-W41 events and split into two windows, each with its own Deming fit: **onset**, shower_on to shower_on + 60 min inclusive, and **decay**, shower_off + 60 min to shower_off + 120 min inclusive.

The onset and decay split is the substantive result. Correlation between C_bed1 and C_room is weak during onset and strong during decay, and C_bed1 reads low against C_room in both. The draft report gives average r² = 0.60 for onset and 0.91 for decay across bins up to 3.0 µm, with average slopes of 1.29 and 1.14 respectively. The room is therefore not uniform during the first hour, and is close to uniform during the second.

[CHECK: the draft report's Fig. 9 caption gives the onset slope as 1.25 while its body text gives 1.29, and its body text describes 73 experiments at 38 °C to 41 °C while its Fig. 10 caption describes 77. Confirm all four numbers against `moduair_room_ratio_fit_w38_w41.csv` and the event count the script reports.]

**Transfer to the single-monitor record.** The bin-specific ratio of the averaged curves is applied to every event that had only Bed 1, to estimate a room average:

$$C_{adjusted\,room}(t) = C_{bed\,1}(t) \cdot \frac{C_{room,average}(t)}{C_{bed\,1,average}(t)} \tag{3-6}$$

Three ratio curves are derived and applied by water temperature: the W38-W41 curve to events between 38 °C and 41 °C, the 51 °C curve to events above 41 °C, and the 23 °C curve to events below 38 °C.

[CHECK: Eq. (3-6) and the temperature-banded application appear in the draft report but in no script in this repository. `moduair_cave_ratio.py` produces the ratio curves; nothing consumes them to build C_adjusted room. Confirm where the adjustment is applied and whether that code should be added here.]

Two limits on the transfer are worth stating explicitly. The ratio curves for temperatures other than 38 °C to 41 °C rest on two events each, so no uncertainty band is computed for them. And the ratio is derived entirely from the June to July fleet period, then applied backwards across a five-month single-monitor record during which the room's floor surfaces changed twice (Section 2.8).

### 3.4.4 Mass Balance Model

Indoor particle concentration in Bedroom #1 for size bin $k$ follows a first-order single-zone mass balance:

$$V \frac{dC}{dt} = p\,Q\,C_{out}(t) - Q\,C(t) - \beta_{loss} V C(t) + E(t) \tag{3-7}$$

Dividing by $V$ and substituting $\lambda = Q/V$:

$$\frac{dC}{dt} = p\,\lambda\,C_{out}(t) - \lambda\,C(t) - \beta_{loss}\,C(t) + \frac{E(t)}{V} \tag{3-8}$$

| Symbol | Description | Units |
|---|---|---|
| $C(t)$ | Bedroom particle concentration, bin $k$ | # cm⁻³ |
| $C_{out}(t)$ | Outdoor particle concentration, bin $k$ | # cm⁻³ |
| $p$ | Penetration factor, fraction of outdoor particles crossing the envelope | dimensionless |
| $\lambda$ | Air change rate from CO₂ decay (Section 3.3) | h⁻¹ |
| $\beta_{loss}$ | First-order aerosol loss rate, surface deposition | h⁻¹ |
| $E(t)$ | Shower aerosol emission rate into the bedroom | # h⁻¹ |
| $V$ | Control volume | m³ |
| $Q$ | Balanced airflow to and from outside the room | m³ h⁻¹ |

Eq. (3-8) assumes all particle loss is first order and does not separate surface deposition from growth of particles into a larger bin. That second point matters: a negative $\beta_{loss}$ is physically interpretable here as net growth of smaller particles into the bin of interest, not as an error.

$p$, $\lambda$, and $\beta_{loss}$ are determined first, then $E(t)$ is solved for.

### 3.4.5 Penetration Factor

One MODULAIR-PM (MOD-PM-00785) sat outdoors for the whole campaign, approximately 5 m from the north wall at 1.5 m above grade. Its concentrations define $C_{out}$.

The penetration factor is the average ratio of indoor to outdoor concentration over a quiet window with no shower emission, taken from 6 hours before shower start to 1 hour before shower start:

$$p = \left. \frac{\overline{C(t)}}{\overline{C_{out}(t)}} \right|_{-6\,\text{h}}^{-1\,\text{h}} \tag{3-9}$$

By definition $0 \leq p \leq 1$. Where the measured indoor concentration exceeded the outdoor concentration, $p$ is set to 1.

Bins above 4.0 µm carry too few valid data for a penetration factor across most events.

[CHECK: two different penetration windows are documented. Eq. (3-9) is the draft report's definition, a single window from 6 h to 1 h before the shower. The superseded `src/particle_calculations.py` uses paired 6-hour windows before and after the event, chosen by time of day (before 20:00 to 02:00 and after 08:00 to 14:00 for night events, the reverse for day events), averaged and then capped at 1.0. Confirm which window produced the reported penetration factors, and note that the reported table stops at 3.0 µm to 4.0 µm while this document's bin table runs to 10.0 µm.]

### 3.4.6 Aerosol Loss Rate

With the shower off there is no emission term, and Eq. (3-8) reduces to:

$$\frac{dC}{dt} = p\,\lambda\,C_{out}(t) - \lambda\,C(t) - \beta_{loss}\,C(t) \tag{3-10}$$

$\beta_{loss}$ is assumed constant through an event. Discretizing with a forward difference:

$$\frac{C_{t_{i+1}} - C_t}{\Delta t} = p\,\lambda\,C_{out,t} - \lambda\,C_t - \beta_{loss}\,C_t \tag{3-11}$$

and rearranging for the loss rate at each time step:

$$\beta_{loss} = \frac{1}{\Delta t} - \lambda - \frac{C_{t_{i+1}}}{C_t\,\Delta t} + \frac{p\,\lambda\,C_{out,t}}{C_t} \tag{3-12}$$

Per-step values are averaged over 1 hour after shower off to 2 hours after shower off. That window is chosen deliberately: Section 3.4.3 shows the room reaches a near-uniform concentration only in the second hour, so the well-mixed assumption behind Eq. (3-10) holds there and not in the first hour.

Two loss rates are computed per event, one from $\lambda_{out}$ giving $\beta_{loss,out}$ and one from $\lambda_{ent}$ giving $\beta_{loss,entry}$, so the λ bracket carries through.

The loss rate is not computed when there is no air change rate for the event, or when the bin lacks continuous measurement for the full 2 hours after the shower. The second condition mostly affects bins above 3.0 µm.

Across events the loss rate runs roughly 5 to 10 times smaller than the air change rate. Some values are negative, which reflects growth of smaller particles into the bin of interest rather than a fitting failure.

[CHECK: the draft report gives the mean $\beta_{loss,out}$ for the 0.66 µm to 1.0 µm bin as −0.13 h⁻¹ ± 0.4 h⁻¹, then gives $\beta_{loss,entry}$ as 0.93 h⁻¹ ± 0.89 h⁻¹, which is character for character the $\lambda_{ent}$ value stated two paragraphs earlier. One of the two is a copy-paste error. Confirm the real $\beta_{loss,entry}$.]

### 3.4.7 Emission Rate

Solving Eq. (3-8) for the emission term at each time step:

$$E_t = p\,\lambda\,V\,C_{out,t} + \frac{V\left(C_t - C_{t_{i+1}}\right)}{\Delta t} - \lambda\,V\,C_t - \beta_{loss}\,V\,C_t \tag{3-13}$$

$E_t$ is computed from shower on until the concentration in the bin first falls below the previous time step, that is, over the rising limb only. The average over those $n$ rising steps is:

$$E_{average} = \Delta t \cdot \frac{E_{t_i} + E_{t_{i+1}} + \dots + E_{t_n}}{n} \tag{3-14}$$

[CHECK: Eq. (3-14) as written in the draft report multiplies a mean rate by a single time step, so with $E_t$ in # h⁻¹ the result carries units of particles, not particles per hour. Either the $\Delta t$ factor does not belong and the quantity is a mean rate, or the $1/n$ does not belong and the quantity is a time-integrated total. Confirm which is intended and correct the equation and the axis labels together.]

**Eight emission variants.** Rather than commit to a single set of inputs, eight emission rates are computed per bin per event from the combinations in Table 3-1. The maximum and minimum $E_{average}$ across the eight are the values presented in the results figures, so the reported emission rate is an interval rather than a point estimate.

**Table 3-1. Input combinations for the eight emission calculations.**

| Parameter | | E₁ | E₂ | E₃ | E₄ | E₅ | E₆ | E₇ | E₈ |
|---|---|---|---|---|---|---|---|---|---|
| Air change rate | λ_entry | X | X | X | X | | | | |
| | λ_out | | | | | X | X | X | X |
| Aerosol loss rate | β_loss,entry | X | X | X | X | | | | |
| | β_loss,out | | | | | X | X | X | X |
| Concentration | C_room(t) | X | X | | | X | X | | |
| | C_adjusted room(t) | | | X | X | | | X | X |
| Volume | 35.9 m³ | X | | X | | X | | X | |
| | 54.6 m³ | | X | | X | | X | | X |

The two concentration series are alternatives only for events with a single monitor. For events with nine or more monitors, $C_{room}$ is used and the variant count halves. The larger volume includes the bathroom inside the mass balance control volume; the smaller is the bedroom alone.

[CHECK: neither volume appears in this repository. `src/particle_calculations.py` carries `BEDROOM_VOLUME_M3 = 36.1`, annotated as 36.10859771 m³ from CAD, against the report's 35.9 m³. The 54.6 m³ combined volume has no stated derivation anywhere. Confirm both, and state whether 54.6 m³ accounts for the shower enclosure and the bathroom's own ceiling height.]

[CHECK: no script in this repository computes Eq. (3-13), Eq. (3-14), or Table 3-1. Confirm where the eight emission variants are calculated so the code can be brought into the repository or the external tool can be cited.]

---

## 3.5 Environmental Conditions Analysis

### 3.5.1 Pre- and Post-Shower Windows

`scripts/rh_temp_other_analysis.py` summarizes RH, temperature, and wind over two windows per event:

- **Pre-shower baseline:** 30 minutes before shower onset (`PRE_SHOWER_MINUTES = 30`).
- **Post-shower response:** 2 hours after shower offset (`POST_SHOWER_HOURS = 2`).

Mean and standard deviation are computed per sensor per window. The post-shower window matches the aerosol deposition window, so environmental conditions and aerosol behavior can be compared directly.

### 3.5.2 Onset and Decay Windows

`scripts/hobo_onset_decay_figures.py` uses the same two windows as the C_bed1 to C_room scatter in Section 3.4.3, so that HOBO temperature and RH can be read against the aerosol uniformity result:

- **onset:** shower_on to shower_on + 60 min, inclusive
- **decay:** shower_off + 60 min to shower_off + 120 min, inclusive

Only HOBO measurements falling inside an onset or decay window of some event are retained, over 2026-01-15 00:00 to 2026-07-16 23:59:59. Onset points are drawn as circles and decay points as squares, one color per sensor.

Sensor display labels follow the report section and intentionally differ from the loader room configuration for MB_E; see the CHECK in Section 2.5.3.

### 3.5.3 Paired Pre and Post Change

A second figure pair reduces each event to two points: the pooled mean across all five HOBO sensors' raw points in the 30-minute pre window and in the 2-hour post window, with pooled standard deviation whiskers, plotted at the window midpoints and colored by window.

Below each figure, a summary reports the paired change from pre to post: the mean difference and standard deviation across all events with data in both windows, and a paired t-test (`scipy.stats.ttest_rel`) of post against pre. The t statistic and p value are reported as NaN when fewer than two events qualify.

The pairing matters. Each event is its own control, so the test measures the within-event shift caused by the shower rather than the between-event variation caused by the room warming from January to July.

### 3.5.4 Bedroom Reference Conditions

For boxplot annotations and environmental stratification, bedroom RH and temperature per event are computed from five channels: Vaisala HMP155 Bed1, HOBO MB_Bed, HOBO MB_F, HOBO MB_C, and Aranet4 Bedroom.

The MODULAIR-PM `met_rh` and `met_temp` channels are excluded from all bedroom characterization. They read the instrument's internal flow cell, not room air.

Mean and standard deviation across the five sensors are computed over the 30-minute pre-shower window, with combined uncertainty estimated as:

$$u_{RH} = \frac{1.96\,\sigma_{RH}}{\sqrt{n}} \tag{3-15}$$

where $\sigma_{RH}$ is the sample standard deviation across the $n$ available sensors and 1.96 is the 97.5th percentile of the standard normal distribution. Temperature is treated identically. These values are written to the `Bedroom_Conditions` sheet of `rh_temp_wind_summary.xlsx` and are exempt from the significant-figure rounding in Section 3.7.1.

### 3.5.5 Sensors Excluded from RH Time-Series Figures

For event-level RH time-series figures, sensors with limited relevance to the bathroom-to-bedroom pathway are omitted to reduce clutter: Vaisala MBa RH, Vaisala Liv RH, Aranet4 Entry RH, Aranet4 Outside RH, and AIO2 outdoor RH. All remain in the full summary tables.

---

## 3.6 Quality Assurance and Data Exclusions

### 3.6.1 Consolidated Exclusion Criteria

**Table 3-2. Exclusion criteria across all analysis domains.**

| Criterion | Applied at | Disposition |
|---|---|---|
| Shower duration outside 10.0 min ± 5 s | All analysis | `is_excluded = True`; `event_number = NaN`; retained in registry |
| Four predefined individual events | All analysis | `is_excluded = True`; retained in registry |
| Four predefined date ranges | All analysis | `is_excluded = True`; retained in registry |
| Initial CO₂ concentration excess below 50 ppm | CO₂ decay (λ) | Event excluded from λ only |
| CO₂ decay regression R² below 0.65 | Aerosol analysis | CO₂ result retained; event excluded from aerosol analysis |
| Fewer than nine fleet monitors reporting at a minute | C_room | Minute dropped; no partial-fleet average |
| Monitor not yet installed at event time | Fleet per-sensor exports and delta peak times | That monitor blank for that event; other monitors unaffected |
| No continuous measurement over the 2 h after shower off | Aerosol loss rate | Bin excluded; primarily affects bins above 3.0 µm |
| No air change rate for the event | Aerosol loss rate and emission rate | Both excluded |

Two monitor-level exclusions are permanent rather than per-event. MOD-PM-00401 reads exactly zero for roughly 90 % of the fleet window and is dropped from the per-bin time series and flagged in the delta-peak output. MOD-PM-00555 is in the bathroom and is excluded from the C_mid zone average.

### 3.6.2 Flow Rate Filter for Summary Figures

For summary boxplots and categorical comparisons, events with measured flow outside 4.1 L min⁻¹ to 5.6 L min⁻¹ are excluded, so that head type or water temperature effects are not confounded with flow rate. The filter does not affect per-event time-series figures or the event registry.

`scripts/export_config_timeseries_fleet.py` treats flow differently on purpose: standard flow and 4.1 L min⁻¹ to 5.6 L min⁻¹ tagged events are pooled under the base configuration key, while restricted-flow events are reported as separate groups.

### 3.6.3 Registry Exclusion Flag and the Fleet Exports

The registry `is_excluded` flag is no longer applied to the MODULAIR-PM configuration exports. Excluded events still appear in those workbooks. The exports are a data product for inspection rather than a summary statistic, and dropping events silently made gaps hard to distinguish from missing data.

---

## 3.7 Reporting Conventions

### 3.7.1 Significant Figures

Numerical results are reported to three significant figures in data output files and to two in figure annotations, applied programmatically through `src/sig_figs.py`. Full precision is available by passing `--no-sig-figs` to any analysis script. The `Bedroom_Conditions` sheet of `rh_temp_wind_summary.xlsx` is exempt and is always written at full precision.

### 3.7.2 Summary Statistics

Unless a figure caption or table header says otherwise:

- **Mean:** arithmetic mean across replicate events within a configuration, or across configurations within a water temperature group.
- **Standard deviation:** sample standard deviation, divisor $n - 1$.
- **Uncertainty bars:** ±1.96 × standard error of the mean, an approximate 95 % confidence interval assuming approximate normality.
- **Boxplots:** center line is the median; box edges are the 25th and 75th percentiles; whiskers reach the most extreme observation within 1.5 × IQR of the nearest box edge; observations beyond the whiskers are plotted as open circles.
- **Standard deviation bands on aligned time series:** drawn only at minutes where at least three events contribute.

Water temperature groupings in boxplots include only baseline configurations, meaning standard nominal flow rate, no mannequin, and standard door and fan positions, unless stated otherwise. Variant runs conducted at W40 appear only in the categorical comparison figures specific to those variables.

### 3.7.3 Output Files

**Table 3-3. Principal analysis outputs, current pipeline.**

| Script | Output | Contents |
|---|---|---|
| `event_registry.py` | `event_log.csv` | Event registry: numbers, test names, λ, exclusion flags |
| `co2_decay_analysis.py` | `co2_lambda_summary.csv` | Per-event λ for all three source methods, R², window details |
| `co2_decay_analysis.py` | `co2_lambda_overall_summary.csv` | λ statistics aggregated by configuration |
| `co2_decay_analysis.py` | `plots/event_figures/co2_decay/event_NN-*.png` | Per-event CO₂ concentration and decay fit |
| `co2_decay_analysis.py` | `plots/air_change_rate_boxplot.png` | λ by water temperature, baseline configurations |
| `moduair_correction_factor.py` | `moduair_correction_195_813_fit_jun04_jul16.csv` | Per-bin Deming fit of 813 against 195: slope, intercept, standard errors |
| `moduair_correction_factor.py` | `moduair_correction_factor_ratios_jun04_jul16.csv` | The three diagnostic ratios per bin |
| `moduair_correction_factor.py` | `moduair_correction_factor_summary_jun04_jul16.csv` | Ratio summary statistics |
| `moduair_correction_factor.py` | `plots/moduair_correction/jun04_jul16/correction_factor_bin{N}.html` | Three ratios per bin on shared axes |
| `moduair_cave_ratio.py` | `moduair_room_ratio_fit.csv` | Per-bin C_bed1 against C_room Deming fit |
| `moduair_cave_ratio.py` | `moduair_room_ratio_fit_w38_w41.csv` | Per-bin fit split by onset and decay window |
| `moduair_cave_ratio.py` | `plots/moduair_room/c_bed1_vs_c_room_bin{N}[_onset\|_decay].html` | Scatter with Deming fit and 1:1 line |
| `moduair_cave_ratio.py` | `plots/moduair_room/c_bed1_c_room_ratio_time_bin{N}.html` | Ratio against aligned time by temperature group |
| `moduair_bin2_timeseries.py` | `plots/moduair_correction/bin{0-10}_timeseries.html` | Fleet concentration per bin with C_room overlay |
| `moduair_event_peak_times.py` | `moduair_event_peak_times.csv` | Event by monitor delta peak time, minutes |
| `moduair_event_peak_times.py` | `plots/moduair_correction/event_delta_peak_times.html` | Delta peak time figure |
| `export_config_timeseries_fleet.py` | `output/event_config_timeseries_fleet/MOD-PM-<sn>/*.xlsx` | Per-monitor aggregated and raw config time series |
| `export_config_timeseries_fleet.py` | `output/event_config_timeseries_fleet/C_room/*.xlsx` | C_room aggregated and raw config time series |
| `export_config_timeseries.py` | `output/event_config_timeseries.xlsx` | Single-sensor 1-min config time series |
| `export_event_timeseries.py` | `output/event_N_timeseries.xlsx` | Single-event predicted and measured concentration, all bins |
| `rh_temp_other_analysis.py` | `rh_temp_wind_summary.xlsx` | RH, temperature, wind statistics; `Bedroom_Conditions`; `Event_Log` |
| `rh_temp_other_analysis.py` | `plots/event_figures/{rh,temperature,wind}_timeseries/event_NN-*.png` | Per-event environmental figures |
| `hobo_onset_decay_figures.py` | `output/plots/hobo/hobo_{temp,rh}_onset_decay.html` | HOBO onset and decay scatter |
| `hobo_onset_decay_figures.py` | `output/plots/hobo/hobo_{temp,rh}_pre_post.html` | Paired pre and post per event with t-test summary |

The fleet exports cover events whose shower_on falls in 2026-06-03 through 2026-07-16, with per-monitor install gating. The C_room pass uses the full window event set with no gating, since a not-yet-installed monitor is simply absent from the reporting count at those minutes.

---

## 3.8 Superseded Single-Sensor Pipeline

`scripts/particle_decay_analysis.py`, with the calculation functions in `src/particle_calculations.py`, implements the aerosol mass balance as it stood before the fleet was installed. It is retained because older figures reference its outputs, and is recorded here so those figures can be dated correctly. It is not the basis of the draft report.

It differs from Sections 3.4.4 to 3.4.7 in five ways:

1. One air change rate, $\lambda_{avg}$, rather than the $\lambda_{out}$ and $\lambda_{ent}$ bracket.
2. One control volume, `BEDROOM_VOLUME_M3 = 36.1`, rather than 35.9 m³ and 54.6 m³.
3. Measured $C_{in}$ from MOD-PM-00195 alone, with no C_room and no adjustment.
4. Penetration from paired 6-hour windows before and after the event, chosen by time of day, averaged and capped at 1.0, rather than a single 6 h to 1 h pre-shower window.
5. $\beta$ from a trimmed mean over the full 2-hour post-shower window with a 5.0 h⁻¹ upper cap and a 5th to 95th percentile trim, then a four-step acceptance hierarchy testing forward Euler R² ≥ 0.80 against the unclamped mean, then a non-negative clamp, then $\beta = 0$, then NaN. The current approach instead restricts the average to the second hour, where Section 3.4.3 shows the room is near-uniform, and admits negative values as particle growth.

Its principal output is `particle_analysis_summary.xlsx`.

---

## 3.9 Divergence Between This Repository and the Draft Report

Three parts of the draft report's analysis are not implemented by any script here. They are listed together so they can be resolved as one decision rather than found one at a time.

| Report element | Where it appears | Status in this repository |
|---|---|---|
| $C_{adjusted\,room}(t)$, Eq. (3-6), and its temperature-banded application | Report Section 3.3.2, Equation 2 | Ratio curves are produced by `moduair_cave_ratio.py`; nothing consumes them to build the adjusted series |
| Eight emission variants, Table 3-1, and Eq. (3-13) and (3-14) | Report Section 3.3.2, Table 6 and Equations 10 and 11 | Not implemented anywhere |
| Control volumes 35.9 m³ and 54.6 m³ | Report Section 3.3.2, Table 6 | Only `BEDROOM_VOLUME_M3 = 36.1` exists |

Applying the 195 correction from Section 3.4.2 is a fourth case: the per-bin Deming fit is produced, but no script applies it to the historical record.

Resolving these means either bringing the calculation into this repository, or citing the external tool that performs it and archiving its inputs and outputs alongside the code.
