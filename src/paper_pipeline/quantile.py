from __future__ import annotations

from typing import Dict, List, Tuple
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import IterationLimitWarning

from .math_utils import (
    iid_bootstrap,
    make_quantile_grid,
    maximum_entropy_bootstrap,
    moving_block_bootstrap,
    residual_bootstrap,
)


def quantile_loss(u: np.ndarray, tau: float) -> float:
    u = np.asarray(u, dtype=float)
    return float(np.sum(np.where(u >= 0, tau * u, (tau - 1.0) * u)))


def exact_small_sample_quantile_slope(years: np.ndarray, values: np.ndarray, tau: float) -> float:
    x = ((np.asarray(years, dtype=float) - np.min(years)) / 10.0).astype(float)
    y = np.asarray(values, dtype=float)
    n = len(y)
    if n < 3:
        return np.nan
    candidates = set()
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[j] - x[i]
            if dx != 0:
                candidates.add(float((y[j] - y[i]) / dx))
    if not candidates:
        return 0.0
    best_b, best_obj = None, None
    for b in candidates:
        resid = y - b * x
        a = float(np.quantile(resid, tau))
        obj = quantile_loss(y - (a + b * x), tau)
        if (best_obj is None) or (obj < best_obj - 1e-12) or (abs(obj - best_obj) <= 1e-12 and abs(b) < abs(best_b)):
            best_obj = obj
            best_b = b
    return float(best_b)


def fit_quantile_slope(years: np.ndarray, values: np.ndarray, tau: float, max_iter: int = 5000) -> float:
    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(years) & np.isfinite(values)
    years = years[mask]
    values = values[mask]
    if len(values) < 3:
        return np.nan
    if np.nanstd(values) == 0:
        return 0.0
    if len(values) <= 12:
        return exact_small_sample_quantile_slope(years, values, tau)
    x = (years - years.min()) / 10.0
    X = sm.add_constant(x)
    try:
        result = _fit_quantreg_with_retry(values, X, tau=tau, max_iter=max_iter)
        return float(result.params[1])
    except Exception:
        return np.nan


def fit_ols_slope(years: np.ndarray, values: np.ndarray) -> float:
    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(years) & np.isfinite(values)
    years = years[mask]
    values = values[mask]
    if len(values) < 3:
        return np.nan
    x = (years - years.min()) / 10.0
    X = sm.add_constant(x)
    try:
        return float(sm.OLS(values, X).fit().params[1])
    except Exception:
        return np.nan


def fit_quantile_line(years: np.ndarray, values: np.ndarray, tau: float, max_iter: int = 5000) -> Dict[str, float | np.ndarray]:
    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(years) & np.isfinite(values)
    years = years[mask]
    values = values[mask]
    if len(values) < 3:
        return {
            "years": years,
            "values": values,
            "x": np.array([], dtype=float),
            "fitted": np.array([], dtype=float),
            "slope": np.nan,
            "intercept": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
        }

    order = np.argsort(years)
    years = years[order]
    values = values[order]
    x = (years - years.min()) / 10.0

    try:
        X = sm.add_constant(x)
        res = _fit_quantreg_with_retry(values, X, tau=tau, max_iter=max_iter)
        ci = res.conf_int()
        ci_low = float(ci[1, 0]) if np.ndim(ci) == 2 else np.nan
        ci_high = float(ci[1, 1]) if np.ndim(ci) == 2 else np.nan
        return {
            "years": years,
            "values": values,
            "x": x,
            "fitted": np.asarray(res.predict(X), dtype=float),
            "slope": float(res.params[1]),
            "intercept": float(res.params[0]),
            "ci_low": ci_low,
            "ci_high": ci_high,
        }
    except Exception:
        slope = fit_quantile_slope(years, values, tau=tau, max_iter=max_iter)
        intercept = float(np.quantile(values - slope * x, tau)) if np.isfinite(slope) else np.nan
        fitted = intercept + slope * x if np.isfinite(intercept) and np.isfinite(slope) else np.array([], dtype=float)
        return {
            "years": years,
            "values": values,
            "x": x,
            "fitted": np.asarray(fitted, dtype=float),
            "slope": float(slope),
            "intercept": intercept,
            "ci_low": np.nan,
            "ci_high": np.nan,
        }


def _fit_quantreg_with_retry(values: np.ndarray, X: np.ndarray, tau: float, max_iter: int):
    model = sm.QuantReg(values, X)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", IterationLimitWarning)
        result = model.fit(q=tau, max_iter=max_iter)
    reached_limit = any(isinstance(w.message, IterationLimitWarning) for w in caught)
    if reached_limit:
        retry_iter = max(int(max_iter) * 10, 5000)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IterationLimitWarning)
            result = model.fit(q=tau, max_iter=retry_iter)
    return result


def fit_ols_line(years: np.ndarray, values: np.ndarray) -> Dict[str, float | np.ndarray]:
    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(years) & np.isfinite(values)
    years = years[mask]
    values = values[mask]
    if len(values) < 3:
        return {
            "years": years,
            "values": values,
            "x": np.array([], dtype=float),
            "fitted": np.array([], dtype=float),
            "slope": np.nan,
            "intercept": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
        }

    order = np.argsort(years)
    years = years[order]
    values = values[order]
    x = (years - years.min()) / 10.0
    X = sm.add_constant(x)
    try:
        res = sm.OLS(values, X).fit()
        ci = res.conf_int()
        ci_low = float(ci[1, 0]) if np.ndim(ci) == 2 else np.nan
        ci_high = float(ci[1, 1]) if np.ndim(ci) == 2 else np.nan
        return {
            "years": years,
            "values": values,
            "x": x,
            "fitted": np.asarray(res.predict(X), dtype=float),
            "slope": float(res.params[1]),
            "intercept": float(res.params[0]),
            "ci_low": ci_low,
            "ci_high": ci_high,
        }
    except Exception:
        slope = fit_ols_slope(years, values)
        intercept = float(np.nanmean(values) - slope * np.nanmean(x)) if np.isfinite(slope) else np.nan
        fitted = intercept + slope * x if np.isfinite(intercept) and np.isfinite(slope) else np.array([], dtype=float)
        return {
            "years": years,
            "values": values,
            "x": x,
            "fitted": np.asarray(fitted, dtype=float),
            "slope": float(slope),
            "intercept": intercept,
            "ci_low": np.nan,
            "ci_high": np.nan,
        }


def bootstrap_qr(years: np.ndarray, values: np.ndarray, focus_quantiles: List[float], cfg: dict, rng: np.random.Generator) -> pd.DataFrame:
    n_reps = int(cfg["bootstrap"]["n_reps"])
    method = str(cfg["bootstrap"]["method"]).lower()
    block_length = int(cfg["bootstrap"]["block_length"])
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    records = []
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    years_valid = np.asarray(years, dtype=float)[valid]
    values_valid = values[valid]
    if len(values_valid) < 3:
        return pd.DataFrame()
    x_decades = (years_valid - years_valid.min()) / 10.0

    valid_methods = {"meboot", "residual", "moving_block", "iid"}
    if method not in valid_methods:
        raise ValueError(f"Unknown bootstrap method '{method}'. Valid options: {sorted(valid_methods)}")

    for rep in range(n_reps):
        if method == "meboot":
            yb = maximum_entropy_bootstrap(values_valid, rng)
        elif method == "residual":
            yb = residual_bootstrap(x_decades, values_valid, rng)
        elif method == "moving_block":
            yb = moving_block_bootstrap(values_valid, block_length, rng)
        else:
            yb = iid_bootstrap(values_valid, rng)
        row = {"replicate": rep}
        for tau in focus_quantiles:
            row[f"slope_{tau:0.2f}"] = fit_quantile_slope(years_valid, yb, tau=tau, max_iter=max_iter)
        records.append(row)

    boot = pd.DataFrame(records)
    if boot.empty:
        return boot
    boot["Delta1"] = boot["slope_0.95"] - boot["slope_0.05"]
    boot["Delta2"] = boot["slope_0.95"] - boot["slope_0.50"]
    boot["Delta3"] = boot["slope_0.50"] - boot["slope_0.05"]
    return boot


def summarize_bootstrap(boot: pd.DataFrame, alpha: float) -> Dict[str, float]:
    out = {}
    if boot.empty:
        return out
    lo, hi = alpha / 2, 1 - alpha / 2
    for col in [c for c in boot.columns if c != "replicate"]:
        series = pd.to_numeric(boot[col], errors="coerce").dropna()
        suffix = col.replace("slope_", "")
        out[f"boot_mean_{suffix}"] = float(series.mean()) if len(series) else np.nan
        out[f"boot_sd_{suffix}"] = float(series.std(ddof=1)) if len(series) > 1 else (0.0 if len(series) == 1 else np.nan)
        out[f"boot_median_{suffix}"] = float(series.median()) if len(series) else np.nan
        out[f"boot_ci_low_{suffix}"] = float(series.quantile(lo)) if len(series) else np.nan
        out[f"boot_ci_high_{suffix}"] = float(series.quantile(hi)) if len(series) else np.nan
    return out


def run_station_qr(annual: pd.DataFrame, cfg: dict) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    dcfg = cfg["data"]
    station_col = dcfg["station_id_col"]
    station_name_col = dcfg["station_name_col"]
    year_col = dcfg["year_col"]
    full_q = make_quantile_grid(**cfg["quantile_regression"]["full_quantiles"])
    focus_q = [float(x) for x in cfg["quantile_regression"]["focus_quantiles"]]
    min_years = int(cfg["quantile_regression"]["min_years_required_to_run"])
    recommended_years = int(cfg["quantile_regression"]["min_years_recommended_for_publication"])
    alpha = float(cfg["bootstrap"]["alpha"])
    rng = np.random.default_rng(int(cfg["project"]["random_seed"]))

    all_quant_records, summary_records, boot_records = [], [], []

    for index_cfg in cfg["indices"]:
        idx_name = index_cfg["name"]
        for (station_id, station_name), sdf in annual.groupby([station_col, station_name_col]):
            sdf = sdf[[year_col, idx_name]].sort_values(year_col).copy()
            years, values = sdf[year_col].to_numpy(dtype=float), sdf[idx_name].to_numpy(dtype=float)
            n_years = int(np.isfinite(values).sum())

            all_slopes = {}
            for tau in full_q:
                slope = fit_quantile_slope(years, values, tau=tau, max_iter=int(cfg["quantile_regression"]["max_iter"]))
                all_quant_records.append({"index_name": idx_name, "station_id": station_id, "station_name": station_name, "tau": tau, "slope": slope, "n_years": n_years})
                all_slopes[tau] = slope

            row = {
                "index_name": idx_name,
                "station_id": station_id,
                "station_name": station_name,
                "n_years": n_years,
                "publication_warning_short_record": bool(n_years < recommended_years),
                "ols_slope": fit_ols_slope(years, values),
                "slope_0.05": all_slopes.get(0.05, np.nan),
                "slope_0.50": all_slopes.get(0.50, np.nan),
                "slope_0.95": all_slopes.get(0.95, np.nan),
            }
            row["Delta1"] = row["slope_0.95"] - row["slope_0.05"]
            row["Delta2"] = row["slope_0.95"] - row["slope_0.50"]
            row["Delta3"] = row["slope_0.50"] - row["slope_0.05"]

            if cfg["bootstrap"]["enabled"] and n_years >= min_years:
                boot = bootstrap_qr(years, values, focus_q, cfg, rng)
                if not boot.empty:
                    boot.insert(0, "index_name", idx_name)
                    boot.insert(1, "station_id", station_id)
                    boot.insert(2, "station_name", station_name)
                    boot_records.append(boot)
                    row.update(summarize_bootstrap(boot.drop(columns=["index_name", "station_id", "station_name"]), alpha))

            summary_records.append(row)

    return pd.DataFrame(all_quant_records), pd.DataFrame(summary_records), (pd.concat(boot_records, ignore_index=True) if boot_records else pd.DataFrame())
