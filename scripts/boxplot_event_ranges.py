#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Event-Range Boxplot Figures
============================

Generate emission_etotal boxplot figures for specific event number ranges.

Mannequin impacts:
  - events 317-327 versus 329-348
  - events 295-307 versus 309-315

Shower head settings:
  - events 277-281 versus 283-286 versus 295-307

Figures use emission_etotal metric, temperature filter 36-42 C, and outside/entry
air-change-rate sources. Output is saved to results/figures/.

Author: Nathan Lima
Institution: National Institute of Standards and Technology (NIST)
Created: 2026-09-24
Update log:
  2026-09-24  Initial version.
"""

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_paths import get_data_root
from src.particle_calculations import PARTICLE_BINS
from src.plot_comparison import _draw_categorical_comparison_boxplot
from src.plot_particle_boxplots import _TEMP_BOXPLOT_CONFIG, _extract_config_temp
from src.plot_style import apply_style

_UNIT_SUFFIX_RE = re.compile(r" \([^)]*\)$")


def _strip_unit_suffixes(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns={c: _UNIT_SUFFIX_RE.sub("", c) for c in df.columns})


def _load_results(output_dir: Path) -> pd.DataFrame:
    path = output_dir / "particle_analysis_summary.xlsx"
    if not path.exists():
        raise FileNotFoundError(f"Particle analysis summary not found: {path}")
    df = pd.read_excel(path, sheet_name="all_results", engine="openpyxl")
    df = _strip_unit_suffixes(df)
    return df


def _temp_in_range(row, lo=36.0, hi=42.0) -> bool:
    t = _extract_config_temp(row.get("config_key", ""))
    return t is not None and lo <= t <= hi


def _make_group_defs(df, ranges):
    """
    Build group definitions for _draw_categorical_comparison_boxplot.

    Parameters
    ----------
    df : pd.DataFrame
        Temperature-filtered results dataframe.
    ranges : list of tuple
        Each entry is (label, start_event, end_event).

    Returns
    -------
    list
        List of (group_key, tick_label, filter_fn) tuples.
    """
    group_defs = []
    for label, start_ev, end_ev in ranges:
        sub = df[df["event_number"].between(start_ev, end_ev)]
        config_set = set(sub["config_key"].astype(str).unique())
        # filter_fn receives config_key string
        def _fn(ck, s=config_set):
            return str(ck) in s
        group_defs.append((label, label.replace("_", " "), _fn))
    return group_defs


def _plot_event_range_comparison(
    df_temp: pd.DataFrame,
    ranges,
    title_base: str,
    x_label: str,
    stem_prefix: str,
    plot_dir: Path,
):
    """
    Plot emission_etotal boxplots for a set of event ranges.

    Parameters
    ----------
    df_temp : pd.DataFrame
        Temperature-filtered results.
    ranges : list
        List of (label, start, end) tuples.
    title_base : str
        Figure title prefix.
    x_label : str
        X-axis label.
    stem_prefix : str
        Output filename stem.
    plot_dir : Path
        Directory for output figures.
    """
    group_defs = _make_group_defs(df_temp, ranges)
    cfg = _TEMP_BOXPLOT_CONFIG["emission_etotal"]
    base_path = plot_dir / f"{stem_prefix}_emission_etotal_boxplot.png"
    _draw_categorical_comparison_boxplot(
        results_df=df_temp,
        particle_bins=PARTICLE_BINS,
        output_path=base_path,
        cfg=cfg,
        group_defs=group_defs,
        title_base=title_base,
        x_label=x_label,
        rh_data=None,
        temp_filter=None,
    )


def main():
    apply_style()
    data_root = get_data_root()
    output_dir = data_root / "output"
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    print("Loading particle analysis summary...")
    df = _load_results(output_dir)
    print(f"  Loaded {len(df)} rows")

    print("Applying temperature filter 36-42 C...")
    df_temp = df[df.apply(lambda r: _temp_in_range(r, 36.0, 42.0), axis=1)].copy()
    print(f"  {len(df_temp)} rows after temperature filter")

    # Mannequin impacts comparison 1
    print("\nMannequin impacts: events 317-327 vs 329-348")
    ranges1 = [("317-327", 317, 327), ("329-348", 329, 348)]
    _plot_event_range_comparison(
        df_temp,
        ranges1,
        title_base="Mannequin Impacts — Emission",
        x_label="Event Range",
        stem_prefix="mannequin_impacts_317-327_vs_329-348",
        plot_dir=plot_dir,
    )

    # Mannequin impacts comparison 2
    print("\nMannequin impacts: events 295-307 vs 309-315")
    ranges2 = [("295-307", 295, 307), ("309-315", 309, 315)]
    _plot_event_range_comparison(
        df_temp,
        ranges2,
        title_base="Mannequin Impacts — Emission",
        x_label="Event Range",
        stem_prefix="mannequin_impacts_295-307_vs_309-315",
        plot_dir=plot_dir,
    )

    # Shower head settings
    print("\nShower head settings: 277-281 vs 283-286 vs 295-307")
    ranges3 = [("277-281", 277, 281), ("283-286", 283, 286), ("295-307", 295, 307)]
    _plot_event_range_comparison(
        df_temp,
        ranges3,
        title_base="Shower Head Settings — Emission",
        x_label="Event Range",
        stem_prefix="shower_head_settings_277-281_vs_283-286_vs_295-307",
        plot_dir=plot_dir,
    )

    print("\nDone. Figures saved to:", plot_dir)


if __name__ == "__main__":
    main()
