from __future__ import annotations

from typing import Tuple

import numpy as np
import pandas as pd

from .math_utils import circular_day_distance, doy_noleap
from .year_config import resolve_reference_years


def compute_daily_thresholds(
    df: pd.DataFrame,
    station_col: str,
    variable_col: str,
    ref_mask: pd.Series,
    doy_col: str,
    cfg: dict,
) -> pd.DataFrame:
    window = int(cfg["index_construction"]["percentile_window_days"])
    low_p = float(cfg["index_construction"]["lower_percentile"])
    high_p = float(cfg["index_construction"]["upper_percentile"])
    min_samples = int(cfg["index_construction"]["min_reference_samples_per_doy"])

    ref_df = df.loc[ref_mask & df[variable_col].notna(), [station_col, doy_col, variable_col]].copy()
    stations = sorted(ref_df[station_col].dropna().unique().tolist())
    all_days = np.arange(1, 366)
    rows = []

    for st in stations:
        sdf = ref_df.loc[ref_df[station_col] == st]
        day_vals = sdf[doy_col].to_numpy(dtype=int)
        vals = sdf[variable_col].to_numpy(dtype=float)
        for day in all_days:
            dist = circular_day_distance(day_vals, day, 365)
            sample = vals[dist <= window]
            if len(sample) < min_samples:
                sample = vals
            rows.append(
                {
                    station_col: st,
                    doy_col: day,
                    f"{variable_col}_p10": np.nanpercentile(sample, low_p) if len(sample) else np.nan,
                    f"{variable_col}_p90": np.nanpercentile(sample, high_p) if len(sample) else np.nan,
                }
            )
    return pd.DataFrame(rows)




def create_extreme_indices(df: pd.DataFrame, cfg: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    dcfg = cfg["data"]
    station_col = dcfg["station_id_col"]
    station_name_col = dcfg["station_name_col"]
    year_col = dcfg["year_col"]
    month_col = dcfg["month_col"]
    day_col = dcfg["day_col"]
    tmin_col = dcfg["tmin_col"]
    tmax_col = dcfg["tmax_col"]

    out = df.copy()
    out["date"] = pd.to_datetime(out[[year_col, month_col, day_col]])
    if bool(cfg["index_construction"]["drop_feb29"]):
        out = out.loc[~((out["date"].dt.month == 2) & (out["date"].dt.day == 29))].copy()
    out["doy"] = doy_noleap(out["date"])

    ref_years, ref_label = resolve_reference_years(cfg, out)
    if ref_years is None:
        ref_mask = pd.Series(True, index=out.index)
    else:
        y0, y1 = ref_years
        ref_mask = out[year_col].between(y0, y1)

    thresholds_tmax = compute_daily_thresholds(out, station_col, tmax_col, ref_mask, "doy", cfg)
    thresholds_tmin = compute_daily_thresholds(out, station_col, tmin_col, ref_mask, "doy", cfg)
    thresholds = thresholds_tmax.merge(thresholds_tmin, on=[station_col, "doy"], how="outer")

    out = out.merge(thresholds, on=[station_col, "doy"], how="left")
    annual_min_coverage_pct = float(cfg["index_construction"].get("annual_min_valid_coverage_pct", 80.0))

    def _event_indicator(value_col: str, threshold_col: str, comparator) -> pd.Series:
        valid = out[value_col].notna() & out[threshold_col].notna()
        indicator = pd.Series(np.nan, index=out.index, dtype=float)
        indicator.loc[valid] = comparator(out.loc[valid, value_col], out.loc[valid, threshold_col]).astype(float)
        return indicator

    out["warm_days"] = _event_indicator(tmax_col, f"{tmax_col}_p90", lambda x, y: x > y)
    out["cool_days"] = _event_indicator(tmax_col, f"{tmax_col}_p10", lambda x, y: x < y)
    out["warm_nights"] = _event_indicator(tmin_col, f"{tmin_col}_p90", lambda x, y: x > y)
    out["cool_nights"] = _event_indicator(tmin_col, f"{tmin_col}_p10", lambda x, y: x < y)

    annual = out.groupby([station_col, station_name_col, year_col], as_index=False).agg(
        warm_days=("warm_days", lambda s: s.sum(min_count=1)),
        warm_nights=("warm_nights", lambda s: s.sum(min_count=1)),
        cool_days=("cool_days", lambda s: s.sum(min_count=1)),
        cool_nights=("cool_nights", lambda s: s.sum(min_count=1)),
        valid_days_tmax=("warm_days", lambda s: int(s.notna().sum())),
        valid_days_tmin=("warm_nights", lambda s: int(s.notna().sum())),
    )
    annual["coverage_pct_tmax"] = annual["valid_days_tmax"] / 365.0 * 100.0
    annual["coverage_pct_tmin"] = annual["valid_days_tmin"] / 365.0 * 100.0
    annual.loc[annual["coverage_pct_tmax"] < annual_min_coverage_pct, ["warm_days", "cool_days"]] = np.nan
    annual.loc[annual["coverage_pct_tmin"] < annual_min_coverage_pct, ["warm_nights", "cool_nights"]] = np.nan
    annual["reference_period_label"] = ref_label
    return out, annual
