from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import QuantileRegressor

from .compound_dry_hot import (
    _add_empirical_return_periods,
    _build_definition_aggregates,
    _make_yearly_extent,
    _periods_from_cfg,
    _trend_summary,
)
from .config_utils import get_focus_quantiles, get_time_scale_years
from .indices import create_extreme_indices


def _fit_ols_slope(years: np.ndarray, values: np.ndarray, time_scale_years: float) -> float:
    x = (years - years.min()) / float(time_scale_years)
    return float(np.polyfit(x, values, 1)[0])


def _fit_quantile_slope(years: np.ndarray, values: np.ndarray, tau: float, time_scale_years: float) -> float:
    x = ((years - years.min()) / float(time_scale_years)).reshape(-1, 1)
    model = QuantileRegressor(quantile=float(tau), alpha=0.0, fit_intercept=True, solver="highs")
    model.fit(x, values)
    return float(model.coef_[0])


def mask_internal_temperature_inconsistencies(data: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mask logically inconsistent temperature observations for sensitivity analysis."""
    dcfg = cfg["data"]
    tmin_col = dcfg["tmin_col"]
    tmax_col = dcfg["tmax_col"]
    tmean_col = dcfg["tmean_col"]
    station_col = dcfg["station_id_col"]
    cleaned = data.copy()
    for col in (tmin_col, tmax_col, tmean_col):
        cleaned[col] = pd.to_numeric(cleaned[col], errors="coerce")

    range_conflict = cleaned[tmin_col].notna() & cleaned[tmax_col].notna() & (cleaned[tmin_col] > cleaned[tmax_col])
    mean_outside = (
        ~range_conflict
        & cleaned[tmean_col].notna()
        & cleaned[tmin_col].notna()
        & cleaned[tmax_col].notna()
        & ((cleaned[tmean_col] < cleaned[tmin_col]) | (cleaned[tmean_col] > cleaned[tmax_col]))
    )
    cleaned.loc[range_conflict, [tmin_col, tmax_col, tmean_col]] = np.nan
    cleaned.loc[mean_outside, tmean_col] = np.nan

    summary = pd.DataFrame(
        [
            {
                "screen": "tmin_greater_than_tmax",
                "rows_flagged": int(range_conflict.sum()),
                "stations_affected": int(data.loc[range_conflict, station_col].nunique()),
                "action_in_sensitivity": "set tmin, tmax, and tmean to missing",
            },
            {
                "screen": "tmean_outside_valid_min_max_range",
                "rows_flagged": int(mean_outside.sum()),
                "stations_affected": int(data.loc[mean_outside, station_col].nunique()),
                "action_in_sensitivity": "set tmean to missing",
            },
        ]
    )
    summary["total_input_rows"] = int(len(data))
    summary["rows_flagged_pct"] = 100.0 * summary["rows_flagged"] / len(data)
    return cleaned, summary


def _network_quantile_summary(annual: pd.DataFrame, cfg: dict, variant: str) -> pd.DataFrame:
    year_col = cfg["data"]["year_col"]
    time_scale = get_time_scale_years(cfg)
    quantiles = get_focus_quantiles(cfg)
    rows = []
    for index_cfg in cfg["indices"]:
        idx = index_cfg["name"]
        regional = annual.groupby(year_col, as_index=False)[idx].mean().dropna()
        years = regional[year_col].to_numpy(dtype=float)
        values = regional[idx].to_numpy(dtype=float)
        row = {
            "variant": variant,
            "index_name": idx,
            "n_years": int(len(values)),
            "ols_slope": _fit_ols_slope(years, values, time_scale),
        }
        for tau in quantiles:
            row[f"slope_{tau:0.2f}"] = _fit_quantile_slope(years, values, tau, time_scale)
        row["Delta1"] = row[f"slope_{max(quantiles):0.2f}"] - row[f"slope_{min(quantiles):0.2f}"]
        rows.append(row)
    return pd.DataFrame(rows)


def _compound_qc_comparison(raw: pd.DataFrame, cleaned: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    variants = []
    for label, frame in (("as_observed", raw), ("consistency_screened", cleaned)):
        aggregates = pd.concat(
            [_build_definition_aggregates(frame, cfg, definition) for definition in cfg["advanced_analyses"]["compound_dry_hot"]["definitions"]],
            ignore_index=True,
        )
        station_year = _add_empirical_return_periods(aggregates, cfg)
        yearly_extent = _make_yearly_extent(station_year, cfg)
        baseline, comparison = _periods_from_cfg(cfg, yearly_extent["year"])
        trend = _trend_summary(yearly_extent, cfg, baseline, comparison)
        trend.insert(0, "variant", label)
        variants.append(trend)
    combined = pd.concat(variants, ignore_index=True)
    keys = ["definition", "definition_title", "return_period_threshold_years"]
    metrics = [
        "kendall_tau",
        "kendall_p_value",
        "linear_slope_stations_per_decade",
        "baseline_mean_affected_stations",
        "comparison_mean_affected_stations",
    ]
    left = combined.loc[combined["variant"] == "as_observed", keys + metrics]
    right = combined.loc[combined["variant"] == "consistency_screened", keys + metrics]
    out = left.merge(right, on=keys, suffixes=("_as_observed", "_screened"))
    for metric in metrics:
        out[f"{metric}_difference"] = out[f"{metric}_screened"] - out[f"{metric}_as_observed"]
    return out


def run_internal_consistency_sensitivity(
    data: pd.DataFrame,
    annual_as_observed: pd.DataFrame,
    cfg: dict,
    outdir: Path,
    progress_callback=None,
) -> dict[str, pd.DataFrame]:
    module_cfg = cfg.get("advanced_analyses", {}).get("internal_consistency_sensitivity", {})
    if not module_cfg.get("enabled", True):
        return {}
    if progress_callback:
        progress_callback("Running internal-temperature-consistency sensitivity analysis...")
    cleaned, screening = mask_internal_temperature_inconsistencies(data, cfg)
    _, annual_screened = create_extreme_indices(cleaned, cfg, progress_callback=None)
    raw_summary = _network_quantile_summary(annual_as_observed, cfg, "as_observed")
    screened_summary = _network_quantile_summary(annual_screened, cfg, "consistency_screened")
    quantile = raw_summary.merge(screened_summary, on=["index_name", "n_years"], suffixes=("_as_observed", "_screened"))
    for metric in ["ols_slope", *[f"slope_{q:0.2f}" for q in get_focus_quantiles(cfg)], "Delta1"]:
        quantile[f"{metric}_difference"] = quantile[f"{metric}_screened"] - quantile[f"{metric}_as_observed"]
    compound = _compound_qc_comparison(data, cleaned, cfg)

    tables_dir = outdir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    screening.to_csv(tables_dir / "temperature_internal_consistency_screening.csv", index=False)
    quantile.to_csv(tables_dir / "temperature_internal_consistency_quantile_sensitivity.csv", index=False)
    compound.to_csv(tables_dir / "temperature_internal_consistency_compound_sensitivity.csv", index=False)
    if progress_callback:
        progress_callback("Saved internal-temperature-consistency sensitivity tables.")
    return {
        "temperature_internal_consistency_screening": screening,
        "temperature_internal_consistency_quantile_sensitivity": quantile,
        "temperature_internal_consistency_compound_sensitivity": compound,
    }
