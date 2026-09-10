# Materials and Methods

> **Author note:** Section numbers, table numbers, and figure references require PI review
> before submission. Items marked `[CHECK: ...]` are unresolved conflicts between this
> repository, the draft report, and the instrument records; each names both values so the
> correct one can be confirmed rather than guessed. Items marked `[PLACEHOLDER: ...]` are
> not documented anywhere yet.
>
> Campaign period: 2026-01-14 (first W48 test day) through 2026-07-16 (last shower event).
> The experiment start cutoff enforced in software is 2026-01-15 15:00
> (`EXPERIMENT_START_DATE` in `src/event_manager.py`).
>
> The test configuration changed repeatedly over the six-month campaign. Section 2.8
> consolidates every change with its date so that any single-state description elsewhere
> in this document can be read against the period it applies to.

---

## 2.1 Test Facility

### 2.1.1 Building Description

The study was conducted in a single-story, double-wide manufactured home on the campus of the National Institute of Standards and Technology (NIST) in Gaithersburg, Maryland. The home has an interior floor area of approximately 140 m² (1,510 ft²) and contains three bedrooms, two full bathrooms, a family room, living room, kitchen, dining area, morning room, utility room, and an attached garage (Fig. 2-1). The home is supplied with NIST campus municipal water and heated by a residential electric hot water heater [PLACEHOLDER: confirm heater type, capacity, and recovery rate]. The building envelope is light-wood-frame construction typical of HUD-code manufactured housing.

The home is served by a forced-air HVAC system controlled by an Ecobee smart thermostat. The HVAC supply vents in Bathroom #1 and Bedroom #1 were sealed with tape for the full duration of the experiments. The two test rooms were therefore heated and cooled passively, by conduction and infiltration from the adjacent conditioned areas of the home and from outdoors. Room air temperature ranged from approximately 15 °C to 25 °C across the January to July campaign.

### 2.1.2 Test Rooms

The test area comprised two adjacent rooms: Bathroom #1 and Bedroom #1 (Fig. 2-1).

**Bathroom #1** has a floor area of 7.09 m² (76.33 ft²) and is in the north-central section of the home. [CHECK: the draft report gives the bathroom as 7 m²; this document gives 7.09 m² from the CAD floor plan. Confirm which is reported.] It contains a [PLACEHOLDER: describe the shower fixture, e.g. fiberglass tub/shower combination or dedicated shower stall, and interior shower dimensions]. A shower curtain was drawn closed for every shower event. The bathroom is equipped with a ceiling-mounted exhaust fan with a measured flow rate of 35.6 m³ h⁻¹ [PLACEHOLDER: confirm fan manufacturer, model, and rated airflow, and the method used to measure 35.6 m³ h⁻¹]. A linen closet (0.72 m², 7.73 ft²) and a utility room (3.66 m², 39.39 ft²) adjoin the bathroom on its west side.

**Bedroom #1** has a floor area of 14.71 m² (158.31 ft²) and a volume of 36.1 m³ (36.10859771 m³ from architectural CAD drawings; `BEDROOM_VOLUME_M3` in `src/particle_calculations.py`). [CHECK: the draft report gives the bedroom floor area as 16.5 m² and uses control volumes of 35.9 m³ (bedroom) and 54.6 m³ (bedroom plus bathroom) in the emission calculation. Neither 35.9 m³ nor 54.6 m³ appears in any script in this repository. Confirm the reported floor area, and the basis and provenance of both volumes.] The ceiling is not flat: it rises from 2.1 m at the windows to 2.7 m at the closet wall. This pitch matters for the aerosol analysis, because the highest-mounted monitors sit under the peak (Section 2.5.1).

Bedroom #1 shares a common wall with Bathroom #1 and the two rooms are connected by an interior [PLACEHOLDER: solid-core / hollow-core] door nominally [PLACEHOLDER: door width × height, mm]. A separate interior door on the east wall opens to the main hallway. Each door has a ventilation gap of 4 cm to 7 cm at the bottom. With the HVAC vents sealed and the bedroom door closed, most of the air exchange between the bedroom and the rest of the home occurred through the gap under the bedroom door; the remainder occurred through cracks and other openings in the exterior and interior walls.

A walk-in closet (2.37 m², 25.5 ft²) is attached to the south side of Bedroom #1. The closet doors were kept closed during all tests.

### 2.1.3 Room Contents and Coordinate System

All monitor positions in this document use a right-handed coordinate system with the origin (0 m, 0 m, 0 m) on the floor at the corner of the bedroom entrance door. The three coordinates are north, west, and vertical, all in meters.

Bedroom #1 contained the following during testing:

- A 0.8 m by 2.4 m table along the wall under the windows, holding the aerosol, temperature, and relative humidity monitors. Two desk chairs and a filing cabinet were stored under this table.
- A 0.6 m by 0.9 m table on the wall opposite the bathtub, holding one desk chair and the tracer gas injection system.
- Two sitting chairs (0.55 m by 0.74 m by 0.74 m) on the wall with the bedroom exit door, present for most but not all of the campaign [CHECK: identify the dates the sitting chairs were removed; the report's Appendix A.1 photograph caption notes their removal but gives no date].
- A 0.4 m by 0.4 m by 0.4 m chemical cabinet on the closet wall between the exit and closet doors.

---

## 2.2 Experimental Design

Seven independent experimental variables were varied across the test period (Table 2-1). Variables were changed sequentially rather than in a randomized or full-factorial design. Each unique combination held constant for one or more complete shower events is termed a *configuration*, and all shower events within a configuration are treated as replicates. The authoritative transition tables are the module-level lists in `src/event_manager.py`; Table 2-2 is generated from them and Section 2.8 summarizes the non-configuration changes.

**Table 2-1. Independent experimental variables and levels investigated.**

| Variable | Levels | Values |
|---|---|---|
| Water temperature code | 17 | W11, W14, W22, W23, W24, W25, W30, W37, W38, W40, W43, W48, W49, W52, W53 °C (W23 is retained for documentation only; its events fall inside an excluded range) |
| Shower head type | 4 | Standard, Pepco, FilterWand, Used |
| Spray pattern | Up to 4 per head | Wide, Narrow, Mid (Pepco); rainfall, 12Nozzle, SingleWide, SingleNarrow (Used); none (Standard, FilterWand) |
| Mannequin presence | 2 | Present, Absent |
| Bathroom door position | 3 | Open, Closed, Ajar |
| Bedroom door position | 2 | Closed, Ajar |
| Bathroom exhaust fan | 2 | Off; On for 12 min from shower onset |
| Measured flow rate | continuous | 1.4 L min⁻¹ to 5.8 L min⁻¹ (Section 2.2.3) |

[CHECK: the draft report's Table 1 identifies the four heads as Shower Head 1 through Shower Head 4 with wand, setting count, and approximate age. This repository identifies them as Standard, Pepco, FilterWand, and Used. No mapping between the two schemes exists in either place. A mapping table is needed before the report and the event registry can be read together.]

### 2.2.1 Water Temperature

Water temperature was set manually at the shower mixing valve and verified by measurement in the collection bucket used for the flow rate determination (Section 2.2.3), in triplicate at the beginning and end of each shower setting. The first of the three buckets read 1 °C to 2 °C cooler than the following two; that bucket was typically collected within the first two minutes of the shower turning on. Water temperature measurements had a mean relative standard deviation of 2.5 %. For experiments in which temperature was not the variable under test, the mean water temperature was 39.0 °C with a relative standard deviation of 2.3 %.

Water temperature codes span W11 (unheated municipal supply, 11 °C) to W53 (53 °C). The code recorded in a configuration key is the nominal setpoint, not the measured bucket temperature.

### 2.2.2 Shower Head Types and Spray Patterns

Four shower heads were used. Two have a wand; two do not. Several have multiple selectable spray patterns:

- **Standard:** installed at campaign start, single fixed spray pattern, no spray pattern recorded in the configuration key.
- **Pepco:** installed 2026-02-24, three selectable patterns recorded as Wide, Narrow, and Mid.
- **FilterWand:** filtered shower wand, installed 2026-04-10 and again 2026-05-22, no spray pattern recorded.
- **Used:** installed 2026-04-13, 2026-05-26, four patterns recorded as rainfall, 12Nozzle, SingleWide, and SingleNarrow.

[PLACEHOLDER: brand, model, rated flow rate, and any flow-restrictor specification for each of the four heads. The draft report gives approximate age and prior usage (Shower Head 1 about 20 years, rarely used; Shower Head 3 new, never used; Shower Head 4 about 14 years, daily use; Shower Head 2 unknown) but not brand or model.]

### 2.2.3 Water Flow Rate

The shower was turned on and off by a solenoid-controlled valve. Flow rate was set manually and measured gravimetrically: a tared 19 L bucket was placed in the water stream, the collected water was massed after exactly 1 minute, and the volume was computed assuming a water density of 1 kg L⁻¹. Measurements were taken in triplicate at the beginning and end of each shower setting. The shower was adjusted to the maximum flow that still held the desired temperature.

The mean flow rate for non-restricted experiments was 4.9 L min⁻¹ with a mean relative standard deviation of 1.1 %. The maximum across all heads was 5.8 L min⁻¹.

A restricted-flow series ran from 2026-05-07 to 2026-05-22 at four settings recorded in `FLOW_RATE_TRANSITIONS`: 1.4 L min⁻¹ (from 2026-05-07 08:00), 2.1 L min⁻¹ (from 2026-05-11 08:35), 2.9 L min⁻¹ (from 2026-05-17 09:00), and 4.2 L min⁻¹ (from 2026-05-20 09:45). [CHECK: three different value sets are in circulation for this series. `src/event_manager.py` records 1.4, 2.1, 2.9, 4.2 L min⁻¹. The draft report states 1.4, 2.0, 2.9, 4.2 L min⁻¹. The docstring of `scripts/export_config_timeseries_fleet.py` groups "1.4 LPM and 2.2 LPM" as separate config groups. Confirm the measured values and reconcile the fleet export grouping.]

The standard analysis range is 4.1 L min⁻¹ to 5.6 L min⁻¹. A configuration key carries a `_FlowRateX.XLPM` suffix only when the measured rate falls outside that range, so that atypical flow events remain distinguishable within an otherwise identical configuration.

### 2.2.4 Mannequin

A human-analog mannequin [PLACEHOLDER: material, approximate size and mass, posture, and placement within the shower enclosure] was present for selected configurations to evaluate the influence of a body analog on aerosol generation and transport. The mannequin was added and removed eight times across the campaign (Table 2-2).

### 2.2.5 Door and Fan Configurations

The bathroom door, between Bathroom #1 and Bedroom #1, was tested Open, Closed, and Ajar [PLACEHOLDER: define "ajar" as an approximate angle or latch-side gap width]. The bedroom door, between Bedroom #1 and the main hallway, was tested Closed and Ajar. The single Ajar bedroom-door period (2026-05-01 17:10 to 2026-05-07 08:00) was unintended and is recorded as a distinct configuration rather than discarded.

The bathroom exhaust fan was tested Off and On. When On, it ran for 12 minutes beginning at shower onset under automated control. The fan was active from 2026-03-31 12:45 to 2026-04-10 08:55; the 12-minute planned duration was recorded from 2026-03-31 12:45 and cleared 2026-04-07 17:00.

### 2.2.6 Configuration Naming Convention

Each configuration is identified by a structured key:

```
W##[_HeadType[_SprayPattern]][_Mannequin]_BathDoorXxx_BdrmDoorXxx_FanXxx[_FlowRateX.XLPM]
```

`W##` is the water temperature code in °C, `HeadType` and `SprayPattern` identify the head and pattern, `_Mannequin` is appended when the mannequin was present, `BathDoorXxx` and `BdrmDoorXxx` give door positions (Open, Closed, Ajar), `FanXxx` gives fan state (On, Off), and the `_FlowRateX.XLPM` suffix appears only for measured flow outside 4.1 L min⁻¹ to 5.6 L min⁻¹.

Individual events within a configuration receive a sequential event number and a structured test name:

```
MMDD_W##[_Head[_Pattern]][_Mannequin]_DoorPos[_Fan]_RNN
```

`MMDD` is the calendar date, `DoorPos` is an abbreviated door position code, and `RNN` is the replicate index within that configuration.

### 2.2.7 Test Timeline and Configuration Transitions

Table 2-2 lists every configuration transition recorded in `src/event_manager.py`, in chronological order. Each row states only what changed; every other variable carried forward from the preceding row. The 2025-12-08 flow rate entry predates the campaign and is a pre-experiment baseline. The W23 entry on 2026-02-11 is retained for documentation only, because its events fall inside an excluded range (Section 3.2.4).

**Table 2-2. Chronological configuration transitions, 2025-12-08 through 2026-07-15.**

| Date and time | Change |
|---|---|
| 2025-12-08 08:00 | Flow rate = 4.6 L/min |
| 2026-01-14 00:00 | Bath door = Open; Bedroom door = Closed; Exhaust fan = Off; Fan duration = none planned; Mannequin = Absent; Shower head = Standard; Spray pattern = none; Water temp = W48 |
| 2026-01-22 14:00 | Flow rate = 4.9 L/min; Water temp = W11 |
| 2026-02-02 17:00 | Flow rate = 5.0 L/min; Water temp = W25 |
| 2026-02-05 10:00 | Water temp = W30 |
| 2026-02-09 10:00 | Flow rate = 4.9 L/min; Water temp = W37 |
| 2026-02-11 08:00 | Flow rate = 5.1 L/min; Water temp = W23 |
| 2026-02-13 11:00 | Water temp = W22 |
| 2026-02-16 08:00 | Flow rate = 4.8 L/min |
| 2026-02-16 11:00 | Water temp = W43 |
| 2026-02-18 10:23 | Flow rate = 5.0 L/min; Water temp = W14 |
| 2026-02-20 08:00 | Flow rate = 4.6 L/min; Water temp = W53 |
| 2026-02-24 08:00 | Flow rate = 4.4 L/min; Shower head = Pepco; Spray pattern = Wide; Water temp = W52 |
| 2026-02-26 10:23 | Flow rate = 4.5 L/min; Water temp = W49 |
| 2026-03-02 09:23 | Flow rate = 4.3 L/min; Spray pattern = Narrow; Water temp = W40 |
| 2026-03-04 10:23 | Flow rate = 4.5 L/min; Spray pattern = Wide |
| 2026-03-06 08:47 | Flow rate = 4.8 L/min; Spray pattern = Mid |
| 2026-03-09 08:47 | Flow rate = 4.2 L/min; Spray pattern = Narrow |
| 2026-03-11 09:47 | Flow rate = 4.3 L/min; Mannequin = Present |
| 2026-03-12 10:47 | Flow rate = 4.2 L/min; Mannequin = Absent; Water temp = W38 |
| 2026-03-13 08:45 | Mannequin = Present |
| 2026-03-17 08:25 | Mannequin = Absent |
| 2026-03-18 09:45 | Flow rate = 4.6 L/min; Mannequin = Present |
| 2026-03-19 08:25 | Mannequin = Absent; Spray pattern = Wide |
| 2026-03-22 09:25 | Mannequin = Present |
| 2026-03-27 08:25 | Mannequin = Absent |
| 2026-03-31 12:45 | Exhaust fan = On; Fan duration = 12 min |
| 2026-04-07 17:00 | Fan duration = none planned |
| 2026-04-08 09:15 | Bath door = Closed; Flow rate = 4.7 L/min |
| 2026-04-10 08:55 | Bath door = Open; Exhaust fan = Off; Flow rate = 5.0 L/min; Shower head = FilterWand; Spray pattern = none |
| 2026-04-13 11:35 | Flow rate = 5.6 L/min; Shower head = Used; Spray pattern = rainfall |
| 2026-04-15 08:25 | Spray pattern = 12Nozzle |
| 2026-04-17 08:35 | Flow rate = 5.5 L/min; Spray pattern = SingleWide |
| 2026-04-20 08:35 | Flow rate = 5.4 L/min; Spray pattern = SingleNarrow |
| 2026-04-23 14:15 | Flow rate = 4.3 L/min; Mannequin = Present; Shower head = Pepco; Spray pattern = Narrow |
| 2026-04-29 08:35 | Flow rate = 4.2 L/min; Mannequin = Absent |
| 2026-05-01 17:10 | Bedroom door = Ajar; Flow rate = 4.3 L/min |
| 2026-05-04 08:30 | Bath door = Ajar |
| 2026-05-07 08:00 | Bath door = Open; Bedroom door = Closed; Flow rate = 1.4 L/min |
| 2026-05-11 08:35 | Flow rate = 2.1 L/min |
| 2026-05-17 09:00 | Flow rate = 2.9 L/min |
| 2026-05-20 09:45 | Flow rate = 4.2 L/min |
| 2026-05-22 09:30 | Flow rate = 4.7 L/min; Shower head = FilterWand |
| 2026-05-22 10:30 | Spray pattern = none |
| 2026-05-26 10:30 | Flow rate = 5.6 L/min; Shower head = Used; Spray pattern = rainfall |
| 2026-06-01 08:00 | Flow rate = 5.7 L/min; Spray pattern = 12Nozzle |
| 2026-06-03 08:30 | Flow rate = 5.6 L/min; Spray pattern = SingleNarrow |
| 2026-06-08 13:15 | Flow rate = 5.5 L/min; Spray pattern = SingleWide |
| 2026-06-11 10:30 | Flow rate = 5.8 L/min; Spray pattern = rainfall |
| 2026-06-22 11:45 | Flow rate = 5.8 L/min; Mannequin = Present |
| 2026-06-26 07:45 | Flow rate = 5.5 L/min; Spray pattern = SingleNarrow |
| 2026-07-02 09:45 | Flow rate = 5.5 L/min; Mannequin = Absent |
| 2026-07-06 10:15 | Flow rate = 5.6 L/min |
| 2026-07-14 10:15 | Flow rate = 5.2 L/min; Water temp = W49 |
| 2026-07-15 10:15 | Flow rate = 5.5 L/min; Water temp = W24 |

---

## 2.3 Shower Test Protocol

### 2.3.1 Automated Control System

The shower solenoid valve, bathroom exhaust fan, CO₂ tracer gas injection valve, and CO₂ mixing fan were all actuated by a computer-based control system [PLACEHOLDER: describe the control hardware and software, e.g. NI cDAQ relay outputs, LabVIEW state machine, or Python script with USB relay board]. Automation held event timing consistent across replicates and removed operator variability from shower start and stop times.

### 2.3.2 Twice-Daily Test Cycle

Shower events ran twice per day, at 03:00 and 15:00 local time. Event selection in `scripts/moduair_event_peak_times.py` gates shower ON transitions to within ±5 minutes of 03:00 or 15:00 for exactly this reason: off-hour ON transitions are water temperature verification runs, not analysis events. The twelve-hour spacing let the room return to near-ambient conditions before the next event.

[CHECK: the draft report states that the protocol ran at 02:00 and 14:00. The repository uses 03:00 and 15:00, and the predefined event exclusions in `src/event_manager.py` are timestamped at 15:00. The campaign crossed the 2026-03-08 daylight saving transition, and 2026-03-08 00:00 to 2026-03-10 00:00 is an excluded range for instrument data misalignment. Confirm whether the clock time changed at the DST transition or whether one of the two documents is simply wrong.]

The cycle below is given as minutes relative to shower onset:

| Time relative to shower on | Action |
|---|---|
| −60 min | Bathroom exhaust fan on |
| −30 min | CO₂ mixing fan on |
| −22 min | CO₂ injection begins; bathroom exhaust fan off |
| −16 min | CO₂ injection ends |
| −15 min | CO₂ mixing fan off; quiet equilibration period begins |
| 0 | Shower on; bathroom exhaust fan on where the fan configuration is On |
| +10 min | Shower off |
| +12 min | Bathroom exhaust fan off where applicable |
| +2 h 10 min | End of analysis window |

The injection raised the bedroom CO₂ concentration to roughly 1500 ppm to 2000 ppm.

**Injection duration changed during the campaign.** The initial protocol injected for 4 minutes. From 2026-01-22 the duration was extended to 6 minutes to raise the initial bedroom concentration and improve the signal-to-noise ratio of the decay fit.

[CHECK: two readings of the change are in circulation and they disagree about which end of the injection moved. The docstring of `scripts/co2_decay_analysis.py` describes injection from :40 to :44 past the hour with the mixing fan running to :45, which is the 4-minute protocol. The draft report describes a 6-minute injection starting 22 minutes before the shower, which ends at :44 and is consistent with the same :45 fan stop, implying the start moved earlier from :40 to :38. An earlier revision of this document instead recorded the end moving later, from :44 to :46. Confirm against the CO₂ injection log (`CO2_log_file.csv`) which end moved, then correct `scripts/co2_decay_analysis.py`, which still documents only the 4-minute protocol.]

The 10-minute shower duration is longer than the 5 min to 8 min typical of reported shower events but falls inside the reported distribution. Duration was enforced by the control system; any event measuring outside 9 min 55 s to 10 min 5 s is treated as a water temperature verification or control-system test and excluded (Section 3.2.4).

### 2.3.3 Water Temperature and Flow Rate Verification

Water temperature and volumetric flow rate were measured in triplicate at the beginning and end of each shower setting, by the gravimetric method in Section 2.2.3. The measured rates are recorded in `FLOW_RATE_TRANSITIONS` and are used to identify events outside the 4.1 L min⁻¹ to 5.6 L min⁻¹ standard analysis range.

---

## 2.4 Carbon Dioxide Tracer Gas System

### 2.4.1 Equipment

Air change rate between Bedroom #1 and the rest of the home and outdoors was determined using pure CO₂ as an inert tracer gas. The system consisted of a compressed CO₂ cylinder [PLACEHOLDER: cylinder size and gas purity] connected through a pressure regulator and a computer-controlled solenoid valve [PLACEHOLDER: manufacturer, model, and whether a mass flow controller set the injection rate]. The injection point and the mixing fan were both on the smaller 0.6 m by 0.9 m table in Bedroom #1, just outside the bathroom door [PLACEHOLDER: injection height above floor; mixing fan manufacturer, model, and blade diameter].

### 2.4.2 Injection and Mixing Protocol

Timing is given in Section 2.3.2. The mixing fan ran through the injection and for a further period afterwards, stopping 15 minutes before shower onset. That quiet interval let the tracer approach a spatially uniform concentration in the bedroom before shower-induced airflow began. The single-zone well-mixed assumption underlying the decay fit (Section 3.3.1) rests on this interval.

---

## 2.5 Instrumentation

### 2.5.1 Particle Counters (QuantAQ MODULAIR-PM)

Particle number concentrations were measured with QuantAQ MODULAIR-PM optical particle counters (QuantAQ, Inc., Somerville, MA). Each unit contains an Alphasense OPC-N3 that reports particle number concentration in 24 size bins spanning 0.35 µm to 40 µm. This study analyzes the 12 lowest bins, 0.35 µm to 10.0 µm (Table 2-3), the range most relevant to respiratory deposition and to Legionella transmission risk. Concentrations are reported as number per cm³.

**The monitor count changed twice during the campaign, and this is the single most consequential change in the dataset.**

| Period | Monitors in Bedroom #1 | Notes |
|---|---|---|
| 2026-01-14 to 2026-06-03 | 1 (Bed 1, MOD-PM-00195) | Plus MOD-PM-00785 outdoors throughout |
| 2026-06-04 to 2026-07-07 | 10 | Fleet installed; MOD-PM-00555 live from 2026-06-26 01:00 |
| 2026-07-08 to 2026-07-16 | 14 | MOD-PM-00465, 00515, 00516 live from 2026-07-08 13:00 |

Before 2026-06-04 a single bedroom monitor cannot establish whether the room was uniformly mixed. The fleet period exists to answer that question and to derive a correction applicable to the single-monitor record (Sections 3.4.1 through 3.4.3).

Table 2-4 gives the fleet layout. Positions use the coordinate system of Section 2.1.3. Labels are those defined in `MODUAIR_SENSOR_LABELS` in `src/plot_style.py`, and the elevation zone assignment is that used by `scripts/moduair_cave_ratio.py`.

**Table 2-4. MODULAIR-PM fleet layout, positions in meters.**

| Label | Serial | North | West | Vertical | Zone | Live from |
|---|---|---|---|---|---|---|
| Low 1 | MOD-PM-00515 | 1.98 | 3.43 | 0.24 | Low | 2026-07-08 13:00 |
| Low 2 | MOD-PM-00465 | 0.61 | 0.58 | 0.25 | Low | 2026-07-08 13:00 |
| Low 3 | MOD-PM-00943 | 1.98 | 1.91 | 0.28 | Low | fleet install |
| Low 4 | MOD-PM-00516 | 1.91 | 0.64 | 0.30 | Low | 2026-07-08 13:00 |
| Low 5 | MOD-PM-00467 | 0.71 | 3.23 | 0.39 | Low | fleet install |
| Bed 1 | MOD-PM-00195 | 0.25 | 1.55 | 0.64 | Bed | 2026-01-14 |
| Bed 2 | MOD-PM-00813 | 0.48 | 1.80 | 0.64 | Bed | fleet install |
| Middle Bathroom | MOD-PM-00555 | 2.24 | 5.59 | 1.14 | excluded from C_mid | 2026-06-26 01:00 |
| Middle 1 | MOD-PM-00814 | 2.24 | 3.45 | 1.17 | Mid | fleet install |
| Middle 2 | MOD-PM-00554 | 1.88 | 0.58 | 1.22 | Mid | fleet install |
| Middle 3 | MOD-PM-00942 | 1.98 | 1.93 | 1.22 | Mid | fleet install |
| Middle 4 | MOD-PM-00401 | 3.48 | 0.46 | 1.32 | Mid | fleet install |
| Middle 5 | MOD-PM-00815 | 3.56 | 1.75 | 1.32 | Mid | fleet install |
| High 1 | MOD-PM-00402 | 0.33 | 1.93 | 1.93 | High | fleet install |
| High 2 | MOD-PM-00816 | 2.03 | 1.91 | 1.98 | High | fleet install |
| Outdoor | MOD-PM-00785 | n/a | n/a | 1.5 | outdoor reference | 2026-01-14 |

Notes on individual monitors:

- **MOD-PM-00195 (Bed 1)** is the only bedroom monitor present for the whole campaign. It sits 0.25 m from the south wall at 0.64 m, chosen to represent exposure for a person in a bed in a room adjoining an active shower. Its historical record was found to be biased when the fleet was installed; the correction is described in Section 3.4.2.
- **MOD-PM-00401 (Middle 4)** reads exactly zero for roughly 90 % of the 2026-06-04 to 2026-07-16 window and is treated as a near-dead sensor. It is dropped from `scripts/moduair_bin2_timeseries.py` and flagged with a printed warning in `scripts/moduair_event_peak_times.py`.
- **MOD-PM-00555 (Middle Bathroom)** is in the bathroom, not the bedroom, and is excluded from the C_mid zone average and from the delta-peak figure. It is retained on the per-bin time series because it shows when aerosol first appears at the source.
- **MOD-PM-00467 and MOD-PM-00785** are excluded from the per-sensor fleet config exports in `scripts/export_config_timeseries_fleet.py`.

The outdoor monitor was approximately 5 m from the north exterior wall of the home at 1.5 m above grade. Its concentrations define C_out in the mass balance.

The MODULAIR-PM also reports temperature and relative humidity from inside the instrument flow cell as `met_temp` and `met_rh`. These do not represent room air and are excluded from all environmental analyses.

Data were acquired at native cadence, retrieved from the QuantAQ cloud API, and archived as weekly raw chunk files by the separate NIST_moduair-pm repository. This repository reads those chunks through `src/moduair_loader.py` and does not download them. Two loading modes exist: `load_sensor_bins()` resamples to a regular 1-minute grid and applies a 10-minute centered rolling average, and `load_sensor_bins_raw()` returns native-cadence records with no resampling and no smoothing.

**Table 2-3. Alphasense OPC-N3 size bins analyzed in this study.**

| Bin | Lower bound (µm) | Upper bound (µm) | Reporting group |
|---|---|---|---|
| 0 | 0.35 | 0.46 | Fine (bins 0 to 2) |
| 1 | 0.46 | 0.66 | Fine (bins 0 to 2) |
| 2 | 0.66 | 1.0 | Fine (bins 0 to 2) |
| 3 | 1.0 | 1.3 | Accumulation and coarse (bins 3 to 6) |
| 4 | 1.3 | 1.7 | Accumulation and coarse (bins 3 to 6) |
| 5 | 1.7 | 2.3 | Accumulation and coarse (bins 3 to 6) |
| 6 | 2.3 | 3.0 | Accumulation and coarse (bins 3 to 6) |
| 7 | 3.0 | 4.0 | Coarse (bins 7 to 11) |
| 8 | 4.0 | 5.2 | Coarse (bins 7 to 11) |
| 9 | 5.2 | 6.5 | Coarse (bins 7 to 11) |
| 10 | 6.5 | 8.0 | Coarse (bins 7 to 11) |
| 11 | 8.0 | 10.0 | Coarse (bins 7 to 11) |

Bin 2 (0.66 µm to 1.0 µm) is the bin used for illustrative figures in the draft report. Bins above 4.0 µm carry too few valid data for a penetration factor across most events.

### 2.5.2 Carbon Dioxide, Temperature, and Relative Humidity (Aranet4 PRO)

CO₂ concentration, air temperature, and relative humidity were measured at four locations with Aranet4 PRO wireless sensors (SAF Tehnika, Riga, Latvia):

| Designation | Location (north, west, vertical, m) | Purpose |
|---|---|---|
| Aranet4 Bedroom | 1.6, 0.5, 0.6 | Primary CO₂ decay tracer zone |
| Aranet4 Bathroom | 5.6, 3.0, 1.5 | Spatial context; not used in λ |
| Aranet4 Entry | Building entry, just outside the bedroom entrance door | Background CO₂ reference for the decay model |
| Aranet4 Outside | Approximately 5 m from the north exterior wall, 1.5 m above grade | Outdoor CO₂ and ambient reference |

[CHECK: the bathroom Aranet4 coordinates above are read from the draft report as (5.6 north, 3.0 west, 1.5 vertical), but the bathroom is west of the bedroom and the report's own monitor tables put bathroom instruments at a west coordinate near 5.6 m. The north and west values may be transposed. Confirm against the floor plan.]

Each sensor uses a non-dispersive infrared detector for CO₂, with a manufacturer-stated accuracy of ±50 ppm ±3 % of reading, and capacitive sensors for temperature and RH. Data were recorded at 1-minute intervals and exported from the Aranet PRO base station as Excel files.

The Bathroom Aranet4 was installed in March 2026. For events before that, Bathroom #1 conditions are characterized from the co-located HOBO UX100 and Vaisala sensors (Sections 2.5.3 and 2.5.4). [PLACEHOLDER: confirm the exact installation date of the Bathroom Aranet4 so that events can be split cleanly on it.]

### 2.5.3 Temperature and Relative Humidity (HOBO UX100 Data Loggers)

Six HOBO UX100 data loggers (Onset Computer Corporation, Bourne, MA) recorded temperature and relative humidity throughout the campaign. Two models were deployed: UX100-011A, a two-channel logger with external probe, and UX100-011, a single integrated probe. Raw temperatures are recorded in °F and converted to °C during processing, after the factory calibration offsets in Table 2-5 are applied additively.

**Table 2-5. HOBO UX100 deployment, calibration offsets, and naming.**

| Loader key | Report label | Model | Serial | Temp offset (°F) | RH offset (%) | Position (north, west, vertical, m) |
|---|---|---|---|---|---|---|
| MB_C | Bedroom East | UX100-011 | not recorded | +0.03 | +0.05 | 1.9, 0.6, 1.4 |
| MB_Bed | Bedroom Central | UX100-011A | 21904272 | −0.23 | −0.47 | 2.0, 1.9, 1.4 |
| MB_E | Bedroom West | UX100-011 | 20244192 | −0.12 | +0.50 | 2.2, 3.5, 1.4 |
| MB_F | Bedroom North | UX100-011 | not recorded | +0.24 | −0.11 | 3.6, 1.8, 1.4 |
| MB_Bath | Bathroom | UX100-011 | 10355906 | +0.13 | +0.57 | 2.2, 5.6, 1.4 |
| MB_D (also MB_G) | not in the report table | UX100-011A | 21904275 | −0.04 | −0.55 | Bathroom #1 |

Two naming conflicts sit in this table and both need resolution before the report is assembled.

[CHECK: `src/env_data_loader.py` describes MB_E as "HOBO Bath/Bed, doorway between bath and bedroom", and `data_config.json` carries the same description. The draft report and `scripts/hobo_onset_decay_figures.py` both call it Bedroom West and place it at 2.2 m north, 3.5 m west, 1.4 m vertical, which is inside the bedroom rather than in the doorway. The docstring of `hobo_onset_decay_figures.py` states the difference is intentional and follows the report. Confirm the physical location, then make the loader description and the report agree.]

[CHECK: MB_D is loaded as "HOBO Bathroom1" and is the primary bathroom logger in this repository, with MB_G as its early-period filename alias for the same serial number 21904275. The draft report's Table 2 lists only five monitors and has no entry for it. Confirm whether MB_D was removed, whether the report's single "Bathroom" row is MB_D or MB_Bath, and which of the two the reported bathroom conditions come from.]

All five monitors in the report's table are at 1.4 m above the floor, with four of them on a line running through the middle of the bedroom, roughly centered on the bathroom door. [PLACEHOLDER: confirm the HOBO logging interval and the mounting method.]

### 2.5.4 Temperature and Relative Humidity (Vaisala HMP Series, Continuous DAQ)

Two Vaisala probe series provided continuous temperature and RH through the NI cDAQ-9178 chassis (Section 2.5.7).

**Vaisala HMP155** probes (Vaisala Oyj, Helsinki, Finland) were deployed in Bedroom #1 (Bed1) and the living room (Liv), each providing one RH channel and one temperature channel. They output 0 V to 10 V DC, require 24 V DC excitation, and were factory calibrated on 2023-05-11. The floor plan instrument legend places them at approximately 1.2 m above the floor.

**Vaisala HMP45A** probes were deployed at several locations including Bathroom #1 (MBa). They output 0 V to 1 V DC and require 7 V to 35 V DC excitation. [PLACEHOLDER: HMP45A calibration dates.]

[PLACEHOLDER: confirm the DAQ scan rate for both probe series and any software averaging applied before writing to file.]

Analog voltages were converted to engineering units with the manufacturer transfer functions and the per-sensor coefficients in the DAQ channel map (`docs/MH cDAQ channels MAP and wiring notes 123125.xlsx`).

### 2.5.5 Outdoor Meteorology (Met One AIO2)

Outdoor wind speed, wind direction, air temperature, relative humidity, and barometric pressure were measured by a Met One AIO2 all-in-one weather sensor (Met One Instruments, Inc., Grants Pass, OR) on a 10 m meteorological tower [PLACEHOLDER: tower position relative to the home, cardinal direction and distance]. Manufacturer specifications recorded in `data_config.json` give a wind speed range of 0 m s⁻¹ to 75 m s⁻¹ (calibrated 0 m s⁻¹ to 60 m s⁻¹), accuracy of ±0.5 m s⁻¹ or 5 % of reading, whichever is greater, and resolution of 0.1 m s⁻¹. The AIO2 communicates over RS-232 to a dedicated logger independent of the NI cDAQ system. [PLACEHOLDER: confirm the AIO2 logging interval.]

### 2.5.6 Differential Pressure (Setra 264)

A Setra Model 264 differential pressure transducer array (Setra Systems, Inc., Boxborough, MA) was installed to measure pressure differences across the envelope and between zones. Eleven channels (0 V to 5 V DC output) were wired to the NI cDAQ-9178 through NI 9201 modules in slots 1 and 2. The array was later disconnected. Differential pressure is not reported in this study. [PLACEHOLDER: confirm the disconnection date and whether any usable DP data were collected beforehand.]

### 2.5.7 Data Acquisition System

Analog sensors (Vaisala HMP155, Vaisala HMP45A, Setra 264) were sampled by a National Instruments CompactDAQ-9178 chassis (National Instruments, Austin, TX), an 8-slot USB chassis rated for 9 V to 30 V DC input at up to 15 W. Five NI 9201 analog input modules, 8 channels each with a ±10 V input range, occupied slots 1 through 5; slots 6 through 8 were unoccupied. The NI 9201 is specified for ±0.04 % full-scale accuracy at 25 °C.

---

## 2.6 Sensor Placement Summary

Positions for the MODULAIR-PM fleet are in Table 2-4, for the HOBO loggers in Table 2-5, and for the Aranet4 sensors in Section 2.5.2. All use the coordinate origin defined in Section 2.1.3. Table 2-6 summarizes the remaining instruments.

**Table 2-6. Remaining instrument placement.**

| Sensor | Location | Approx. height (m) | Notes |
|---|---|---|---|
| Vaisala HMP155 Bed1 | Bedroom #1 | 1.2 | RH and temperature, continuous DAQ |
| Vaisala HMP155 Liv | Living room | 1.2 | RH and temperature, continuous DAQ |
| Vaisala HMP45A MBa | Bathroom #1 | [PLACEHOLDER] | RH and temperature, continuous DAQ |
| Met One AIO2 | Outdoors, 10 m tower | 10.0 | Wind, temperature, RH, pressure |
| Setra 264 | Multiple zones | Various | Disconnected; not reported |

---

## 2.7 Instrument Calibration and Applied Corrections

HOBO UX100 loggers were factory calibrated by Onset before deployment. The offsets in Table 2-5 are all within 0.24 °F for temperature and 0.57 % for RH. They are applied additively to raw data during processing, before the °F to °C conversion.

Vaisala HMP155 probes were factory calibrated on 2023-05-11, with records on file. HMP45A calibration dates are outstanding. Transfer functions for both series are applied from the manufacturer documentation using the per-sensor coefficients in the DAQ channel map.

Aranet4 PRO sensors carry a manufacturer CO₂ accuracy of ±50 ppm ±3 % of reading. At the 1500 ppm to 2000 ppm peak bedroom concentrations produced by the injection protocol, that is approximately ±95 ppm to ±110 ppm. No field zero or span check was performed. [PLACEHOLDER: confirm whether any co-location check of the four Aranet4 sensors was performed, and whether the CO₂ injection mass flow was ever verified.]

**MODULAIR-PM MOD-PM-00195 campaign correction.** When the fleet was installed in June 2026, the historical MOD-PM-00195 record was found to be biased. MOD-PM-00813 was placed alongside it at the Bed position so a co-located correlation could be built, and an orthogonal-distance (Deming) regression is fit per bin with 195 on the x-axis and 813 on the y-axis. The fit maps a measured 195 value to its 813-equivalent. Deming rather than ordinary least squares is used because both instruments are the same model and carry comparable measurement error, so equal x and y error variance is assumed (delta = 1). Per-bin slope, intercept, and standard errors are written by `scripts/moduair_correction_factor.py`; the script fits and reports only and does not modify the historical record. The analysis procedure is in Section 3.4.2.

---

## 2.8 Campaign Change Log

The test configuration, the instrument suite, and the room itself all changed during the six-month campaign. Any statement elsewhere in this document that describes a single state should be read against this table. Configuration-variable changes are in Table 2-2 and are not repeated here.

**Table 2-7. Non-configuration changes over the campaign.**

| Date | Change | Consequence for analysis |
|---|---|---|
| 2026-01-14 | Campaign begins. One bedroom PM monitor (MOD-PM-00195) and one outdoor monitor (MOD-PM-00785) | Room uniformity cannot be assessed from measurement before 2026-06-04 |
| 2026-01-15 15:00 | Software experiment start cutoff (`EXPERIMENT_START_DATE`) | Shower events before this are not registered |
| 2026-01-22 | CO₂ injection duration extended from 4 min to 6 min | Higher initial bedroom CO₂; better decay signal-to-noise after this date. See the CHECK in Section 2.3.2 on which end of the injection moved |
| 2026-03 (day outstanding) | Bathroom Aranet4 installed | Bathroom CO₂, temperature, and RH available only after this date |
| 2026-03-08 to 2026-03-10 | Daylight saving transition; instrument data misalignment and elevated bedroom RH of uncertain origin | Excluded date range |
| 2026-03-14 to 2026-03-15 12:00 | CO₂ injection system failure | Excluded date range |
| 2026-05-11 to 2026-05-15 10:00 | CO₂ injection system failure | Excluded date range |
| 2026-05-13 | Luxury vinyl plank flooring installed in the bedroom | Single event excluded. The bedroom floor surface changed for the remainder of the campaign |
| 2026-05-21 | Bathroom flooring removed | Single event excluded. The bathroom floor surface changed for the remainder of the campaign |
| 2026-06-04 | MODULAIR-PM fleet installed: nine additional bedroom monitors plus MOD-PM-00813 co-located with MOD-PM-00195 at the Bed position | Room uniformity assessment and the C_room average become possible. Start of the fixed 2026-06-04 to 2026-07-16 fleet analysis window |
| 2026-06-26 01:00 | MOD-PM-00555 (Middle Bathroom) live | Bathroom aerosol source timing available from this date |
| 2026-07-08 13:00 | MOD-PM-00465, 00515, 00516 live, bringing the bedroom fleet to 14 | The Low zone is fully populated only after this date. Events before it are gated per sensor |
| 2026-07-16 | Last shower event | End of campaign |
| Date outstanding | Two sitting chairs removed from the bedroom | Bedroom furniture and surface area changed. [CHECK: identify the date; the report's Appendix A.1 photographs show both states] |

Two of these deserve emphasis because they cut across the whole dataset.

**The flooring work on 2026-05-13 and 2026-05-21 changed the interior surfaces of both test rooms partway through the campaign.** Only the two shower events on those days are excluded. The aerosol loss rate βloss is a surface deposition term, so a change in floor material and area is a plausible step change in βloss between the pre- and post-flooring periods. [CHECK: no analysis in this repository splits events on the flooring dates. Confirm whether βloss should be tested for a step change across 2026-05-13 and 2026-05-21 before results are pooled across the full campaign.]

**The single-monitor period (2026-01-14 to 2026-06-03) covers most of the events, and the fleet period (2026-06-04 to 2026-07-16) covers the room uniformity evidence.** Every conclusion about room-average concentration for the single-monitor period rests on a correction derived from the fleet period and transferred backwards. Section 3.4.3 describes the transfer and its limits.
