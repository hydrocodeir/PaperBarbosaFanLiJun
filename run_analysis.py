from __future__ import annotations

import argparse
import json
import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm


def make_quantile_grid(start: float, stop: float, step: float) -> List[float]:
    vals = []
    x = start
    while x <= stop + 1e-9:
        vals.append(round(float(x), 2))
        x += step
    return vals




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
    best_b = None
    best_obj = None
    for b in candidates:
        resid = y - b * x
        a = float(np.quantile(resid, tau))
        obj = quantile_loss(y - (a + b * x), tau)
        if (best_obj is None) or (obj < best_obj - 1e-12) or (abs(obj - best_obj) <= 1e-12 and abs(b) < abs(best_b)):
            best_obj = obj
            best_b = b
    return float(best_b)

def doy_noleap(dt: pd.Series) -> pd.Series:
    doy = dt.dt.dayofyear.astype(int)
    leap = dt.dt.is_leap_year
    after_feb28 = (dt.dt.month > 2) | ((dt.dt.month == 2) & (dt.dt.day == 29))
    doy = doy - ((leap & after_feb28).astype(int))
    return doy


def circular_day_distance(days: np.ndarray, center: int, max_day: int = 365) -> np.ndarray:
    raw = np.abs(days - center)
    return np.minimum(raw, max_day - raw)


def moving_block_bootstrap(arr: np.ndarray, block_length: int, rng: np.random.Generator) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    n = len(arr)
    if n == 0:
        return arr.copy()
    block_length = max(1, min(int(block_length), n))
    starts = rng.integers(0, n, size=int(math.ceil(n / block_length)))
    out = []
    for s in starts:
        idx = [(s + i) % n for i in range(block_length)]
        out.extend(arr[idx].tolist())
        if len(out) >= n:
            break
    return np.asarray(out[:n], dtype=float)


def iid_bootstrap(arr: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    n = len(arr)
    if n == 0:
        return arr.copy()
    idx = rng.integers(0, n, size=n)
    return arr[idx]


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
    x = (years - years.min()) / 10.0  # per decade
    X = sm.add_constant(x)
    try:
        model = sm.QuantReg(values, X)
        res = model.fit(q=tau, max_iter=max_iter)
        return float(res.params[1])
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
        res = sm.OLS(values, X).fit()
        return float(res.params[1])
    except Exception:
        return np.nan


def bootstrap_qr(years: np.ndarray, values: np.ndarray, focus_quantiles: List[float], cfg: dict, rng: np.random.Generator):
    n_reps = int(cfg["bootstrap"]["n_reps"])
    method = cfg["bootstrap"]["method"]
    block_length = int(cfg["bootstrap"]["block_length"])
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    records = []
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    years_valid = np.asarray(years, dtype=float)[valid]
    values_valid = values[valid]
    if len(values_valid) < 3:
        return pd.DataFrame()

    for rep in range(n_reps):
        if method == "moving_block":
            yb = moving_block_bootstrap(values_valid, block_length, rng)
        elif method == "iid":
            yb = iid_bootstrap(values_valid, rng)
        else:
            raise ValueError(f"Unsupported bootstrap method: {method}")
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
    lo = alpha / 2
    hi = 1 - alpha / 2
    for col in [c for c in boot.columns if c != "replicate"]:
        series = pd.to_numeric(boot[col], errors="coerce").dropna()
        if len(series) == 0:
            out[f"boot_mean_{col.replace('slope_', '').replace('Delta', 'Delta')}"] = np.nan
            out[f"boot_sd_{col.replace('slope_', '').replace('Delta', 'Delta')}"] = np.nan
            out[f"boot_median_{col.replace('slope_', '').replace('Delta', 'Delta')}"] = np.nan
            out[f"boot_ci_low_{col.replace('slope_', '').replace('Delta', 'Delta')}"] = np.nan
            out[f"boot_ci_high_{col.replace('slope_', '').replace('Delta', 'Delta')}"] = np.nan
        else:
            suffix = col.replace("slope_", "")
            out[f"boot_mean_{suffix}"] = float(series.mean())
            out[f"boot_sd_{suffix}"] = float(series.std(ddof=1)) if len(series) > 1 else 0.0
            out[f"boot_median_{suffix}"] = float(series.median())
            out[f"boot_ci_low_{suffix}"] = float(series.quantile(lo))
            out[f"boot_ci_high_{suffix}"] = float(series.quantile(hi))
    return out


def compute_daily_thresholds(df: pd.DataFrame, station_col: str, date_col: str, variable_col: str, ref_mask: pd.Series, doy_col: str, cfg: dict) -> pd.DataFrame:
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
                    f"{variable_col}_nref": int(len(sample)),
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

    ref_years = cfg["index_construction"]["reference_years"]
    if ref_years is None:
        ref_mask = pd.Series(True, index=out.index)
        ref_label = "all_available_years"
    else:
        y0, y1 = int(ref_years[0]), int(ref_years[1])
        ref_mask = out[year_col].between(y0, y1)
        ref_label = f"{y0}_{y1}"

    thresholds_tmax = compute_daily_thresholds(out, station_col, "date", tmax_col, ref_mask, "doy", cfg)
    thresholds_tmin = compute_daily_thresholds(out, station_col, "date", tmin_col, ref_mask, "doy", cfg)
    thresholds = thresholds_tmax.merge(thresholds_tmin, on=[station_col, "doy"], how="outer")

    out = out.merge(thresholds, on=[station_col, "doy"], how="left")
    out["warm_days"] = (out[tmax_col] > out[f"{tmax_col}_p90"]).astype(float)
    out["cool_days"] = (out[tmax_col] < out[f"{tmax_col}_p10"]).astype(float)
    out["warm_nights"] = (out[tmin_col] > out[f"{tmin_col}_p90"]).astype(float)
    out["cool_nights"] = (out[tmin_col] < out[f"{tmin_col}_p10"]).astype(float)

    annual = (
        out.groupby([station_col, station_name_col, year_col], as_index=False)[["warm_days", "warm_nights", "cool_days", "cool_nights"]]
        .sum(min_count=1)
    )
    annual["reference_period_label"] = ref_label
    return out, annual


def run_station_qr(annual: pd.DataFrame, cfg: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
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

    all_quant_records = []
    summary_records = []
    boot_records = []

    for index_cfg in cfg["indices"]:
        idx_name = index_cfg["name"]
        for (station_id, station_name), sdf in annual.groupby([station_col, station_name_col]):
            sdf = sdf[[year_col, idx_name]].sort_values(year_col).copy()
            years = sdf[year_col].to_numpy(dtype=float)
            values = sdf[idx_name].to_numpy(dtype=float)
            n_years = int(np.isfinite(values).sum())

            all_slopes = {}
            for tau in full_q:
                slope = fit_quantile_slope(years, values, tau=tau, max_iter=int(cfg["quantile_regression"]["max_iter"]))
                all_quant_records.append(
                    {
                        "index_name": idx_name,
                        "station_id": station_id,
                        "station_name": station_name,
                        "tau": tau,
                        "slope": slope,
                        "n_years": n_years,
                    }
                )
                all_slopes[tau] = slope

            row = {
                "index_name": idx_name,
                "station_id": station_id,
                "station_name": station_name,
                "n_years": n_years,
                "publication_warning_short_record": bool(n_years < recommended_years),
                "ols_slope": fit_ols_slope(years, values),
            }
            for tau in focus_q:
                row[f"slope_{tau:0.2f}"] = all_slopes.get(tau, np.nan)
            row["Delta1"] = row.get("slope_0.95", np.nan) - row.get("slope_0.05", np.nan)
            row["Delta2"] = row.get("slope_0.95", np.nan) - row.get("slope_0.50", np.nan)
            row["Delta3"] = row.get("slope_0.50", np.nan) - row.get("slope_0.05", np.nan)

            if cfg["bootstrap"]["enabled"] and n_years >= min_years:
                boot = bootstrap_qr(years, values, focus_q, cfg, rng)
                if not boot.empty:
                    boot.insert(0, "index_name", idx_name)
                    boot.insert(1, "station_id", station_id)
                    boot.insert(2, "station_name", station_name)
                    boot_records.append(boot)
                    row.update(summarize_bootstrap(boot.drop(columns=["index_name", "station_id", "station_name"]), alpha))
            summary_records.append(row)

    summary = pd.DataFrame(summary_records)
    all_quant = pd.DataFrame(all_quant_records)
    boot_long = pd.concat(boot_records, ignore_index=True) if boot_records else pd.DataFrame()
    return all_quant, summary, boot_long


def build_feature_table(summary: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    feature_rows = []
    for _, row in summary.iterrows():
        out = row.to_dict()
        if "slope_0.05" in out and "slope_0.95" in out:
            out["Delta1"] = out["slope_0.95"] - out["slope_0.05"]
        if "slope_0.50" in out and "slope_0.95" in out:
            out["Delta2"] = out["slope_0.95"] - out["slope_0.50"]
        if "slope_0.05" in out and "slope_0.50" in out:
            out["Delta3"] = out["slope_0.50"] - out["slope_0.05"]
        feature_rows.append(out)
    return pd.DataFrame(feature_rows)


def run_clustering(features: pd.DataFrame, cfg: dict) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    if features.empty or not cfg["clustering"]["enabled"]:
        return pd.DataFrame(), {}
    mode = cfg["clustering"]["feature_mode"]
    if mode == "simple":
        feature_cols = cfg["clustering"]["simple_features"]
    elif mode == "uncertainty":
        feature_cols = cfg["clustering"]["uncertainty_features"]
    else:
        raise ValueError(f"Unsupported feature mode: {mode}")

    cluster_rows = []
    artifacts = {}
    for idx_name, sdf in features.groupby("index_name"):
        X = sdf[feature_cols].apply(pd.to_numeric, errors="coerce")
        X = X.replace([np.inf, -np.inf], np.nan)
        # fill missing values column-wise with median so clustering can still run
        X = X.apply(lambda col: col.fillna(col.median()), axis=0)
        scaler = StandardScaler() if cfg["clustering"]["standardize"] else None
        Xv = scaler.fit_transform(X) if scaler is not None else X.to_numpy()

        if cfg["clustering"]["algorithm"] == "hierarchical":
            Z = linkage(Xv, method=cfg["clustering"]["linkage"], metric=cfg["clustering"]["metric"])
            labels = fcluster(Z, t=int(cfg["clustering"]["n_clusters"]), criterion="maxclust")
            artifacts[idx_name] = Z
        elif cfg["clustering"]["algorithm"] == "kmeans":
            km = KMeans(n_clusters=int(cfg["clustering"]["n_clusters"]), random_state=int(cfg["project"]["random_seed"]), n_init=20)
            labels = km.fit_predict(Xv) + 1
            artifacts[idx_name] = km.cluster_centers_
        else:
            raise ValueError(f"Unsupported clustering algorithm: {cfg['clustering']['algorithm']}")

        cdf = sdf[["index_name", "station_id", "station_name"]].copy()
        cdf["cluster"] = labels
        cluster_rows.append(cdf)

    cluster_df = pd.concat(cluster_rows, ignore_index=True) if cluster_rows else pd.DataFrame()
    return cluster_df, artifacts


def _sort_station_names_for_heatmap(plot_df: pd.DataFrame, mode: str) -> List[str]:
    if mode == "alphabetic":
        return sorted(plot_df["station_name"].unique().tolist())
    else:
        # default: sort by Delta1 descending
        order_df = plot_df.groupby("station_name", as_index=False)["Delta1"].mean().sort_values("Delta1", ascending=False)
        return order_df["station_name"].tolist()


def plot_data_coverage(annual: pd.DataFrame, outdir: Path, cfg: dict):
    counts = annual.groupby("station_name")["year"].agg(["min", "max", "count"]).sort_values("count", ascending=False)
    fig, ax = plt.subplots(figsize=(10, 8))
    y = np.arange(len(counts))
    ax.barh(y, counts["count"], color="steelblue")
    ax.set_yticks(y)
    ax.set_yticklabels(counts.index)
    ax.invert_yaxis()
    ax.set_xlabel("Number of annual index values")
    ax.set_title("Data coverage by station")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / f"data_coverage_by_station.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def plot_region_quantile_slopes(annual: pd.DataFrame, outdir: Path, cfg: dict):
    years = np.sort(annual["year"].unique())
    mean_df = annual.groupby("year", as_index=False)[[idx["name"] for idx in cfg["indices"]]].mean()
    qgrid = make_quantile_grid(**cfg["quantile_regression"]["full_quantiles"])
    focus = [float(x) for x in cfg["quantile_regression"]["focus_quantiles"]]
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        vals = mean_df[idx_name].to_numpy(dtype=float)
        slopes = [fit_quantile_slope(years, vals, tau=q, max_iter=int(cfg["quantile_regression"]["max_iter"])) for q in qgrid]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(qgrid, slopes, marker="o", markersize=2.5, linewidth=1.2)
        for fq in focus:
            fq_slope = fit_quantile_slope(years, vals, tau=fq, max_iter=int(cfg["quantile_regression"]["max_iter"]))
            ax.scatter([fq], [fq_slope], s=40, zorder=3)
            ax.annotate(f"q={fq:.2f}\n{fq_slope:.2f}", (fq, fq_slope), textcoords="offset points", xytext=(5, 5), fontsize=8)
        ax.axhline(0, color="black", linewidth=0.8, alpha=0.7)
        ax.set_xlabel("Quantile")
        ax.set_ylabel(f"Slope ({idx_cfg['unit']} per decade)")
        ax.set_title(f"Regional mean quantile-regression slopes - {idx_cfg['title']}")
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(outdir / f"region_quantile_slopes_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def plot_station_heatmap(features: pd.DataFrame, outdir: Path, cfg: dict):
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        sdf = features.loc[features["index_name"] == idx_name, ["station_name", "slope_0.05", "slope_0.50", "slope_0.95", "Delta1"]].copy()
        if sdf.empty:
            continue
        order = _sort_station_names_for_heatmap(sdf, cfg["plots"]["heatmap_station_order"])
        sdf["station_name"] = pd.Categorical(sdf["station_name"], categories=order, ordered=True)
        sdf = sdf.sort_values("station_name")
        mat = sdf.set_index("station_name")[["slope_0.05", "slope_0.50", "slope_0.95", "Delta1"]]
        fig, ax = plt.subplots(figsize=(8, max(6, len(mat) * 0.28)))
        im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="coolwarm")
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels(["q0.05", "q0.50", "q0.95", "Delta1"], rotation=0)
        ax.set_yticks(range(mat.shape[0]))
        ax.set_yticklabels(mat.index.tolist(), fontsize=8)
        ax.set_title(f"Station-level QR slopes and asymmetry - {idx_cfg['title']}")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                val = mat.iloc[i, j]
                if np.isfinite(val):
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label(f"Slope / delta ({idx_cfg['unit']} per decade)")
        fig.tight_layout()
        fig.savefig(outdir / f"station_focus_heatmap_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def plot_delta_uncertainty(features: pd.DataFrame, outdir: Path, cfg: dict):
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        cols = ["station_name", "boot_mean_Delta1", "boot_ci_low_Delta1", "boot_ci_high_Delta1"]
        if not set(cols).issubset(features.columns):
            continue
        sdf = features.loc[features["index_name"] == idx_name, cols].copy().dropna()
        if sdf.empty:
            continue
        sdf = sdf.sort_values("boot_mean_Delta1", ascending=False)
        fig, ax = plt.subplots(figsize=(9, max(6, len(sdf) * 0.28)))
        y = np.arange(len(sdf))
        x = sdf["boot_mean_Delta1"].to_numpy()
        xerr = np.vstack([x - sdf["boot_ci_low_Delta1"].to_numpy(), sdf["boot_ci_high_Delta1"].to_numpy() - x])
        ax.errorbar(x, y, xerr=xerr, fmt="o", capsize=3)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(sdf["station_name"].tolist(), fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Bootstrap mean Delta1 with 95% CI")
        ax.set_title(f"Asymmetric trend uncertainty by station - {idx_cfg['title']}")
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(outdir / f"delta_uncertainty_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def _scatter_station_map(df: pd.DataFrame, value_col: str, title: str, outpath: Path, cfg: dict, discrete: bool = False):
    fig, ax = plt.subplots(figsize=(7.5, 6))
    x = df["longitude"].to_numpy()
    y = df["latitude"].to_numpy()
    c = df[value_col].to_numpy()
    scatter = ax.scatter(x, y, c=c, s=70, cmap="viridis" if discrete else "coolwarm", edgecolor="black", linewidth=0.4)
    if cfg["plots"]["annotate_stations_on_map"]:
        for _, r in df.iterrows():
            ax.annotate(r["station_name"], (r["longitude"], r["latitude"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)
    ax.grid(alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.85)
    cbar.set_label(value_col)
    fig.tight_layout()
    fig.savefig(outpath, dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def plot_maps(features: pd.DataFrame, stations: pd.DataFrame, cluster_df: pd.DataFrame, outdir: Path, cfg: dict):
    merged = features.merge(stations, on=["station_id", "station_name"], how="left")
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        sdf = merged.loc[merged["index_name"] == idx_name].copy()
        if sdf.empty:
            continue
        for col, lab in [("slope_0.05", "q0.05 slope"), ("slope_0.50", "q0.50 slope"), ("slope_0.95", "q0.95 slope"), ("Delta1", "Delta1 (q0.95-q0.05)")]:
            if col in sdf.columns:
                _scatter_station_map(
                    sdf,
                    value_col=col,
                    title=f"Spatial distribution of {lab} - {idx_cfg['title']}",
                    outpath=outdir / f"map_{idx_name}_{col}.{cfg['plots']['save_format']}",
                    cfg=cfg,
                )
        if not cluster_df.empty:
            if "cluster" in sdf.columns:
                csdf = sdf.copy()
            else:
                csdf = sdf.merge(cluster_df.loc[cluster_df["index_name"] == idx_name, ["station_id", "cluster"]], on="station_id", how="left")
            if "cluster" in csdf.columns:
                _scatter_station_map(
                    csdf,
                    value_col="cluster",
                    title=f"Station clusters - {idx_cfg['title']}",
                    outpath=outdir / f"map_{idx_name}_cluster.{cfg['plots']['save_format']}",
                    cfg=cfg,
                    discrete=True,
                )


def plot_dendrograms(features: pd.DataFrame, artifacts: Dict[str, np.ndarray], outdir: Path, cfg: dict):
    if not artifacts:
        return
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        if idx_name not in artifacts or cfg["clustering"]["algorithm"] != "hierarchical":
            continue
        sdf = features.loc[features["index_name"] == idx_name].copy()
        labels = sdf["station_name"].tolist()
        Z = artifacts[idx_name]
        fig, ax = plt.subplots(figsize=(10, 6))
        dendrogram(Z, labels=labels, orientation="right", leaf_font_size=8, ax=ax)
        ax.set_title(f"Hierarchical clustering dendrogram - {idx_cfg['title']}")
        ax.set_xlabel("Distance")
        fig.tight_layout()
        fig.savefig(outdir / f"dendrogram_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def generate_report(annual: pd.DataFrame, features: pd.DataFrame, cluster_df: pd.DataFrame, cfg: dict, outdir: Path):
    n_years_total = annual["year"].nunique()
    station_count = annual["station_id"].nunique()
    lines = []
    lines.append("# Quantile Regression + Bootstrap + Clustering Report")
    lines.append("")
    lines.append("## Data audit")
    lines.append(f"- Number of stations: **{station_count}**")
    lines.append(f"- Number of years in uploaded data: **{n_years_total}**")
    lines.append(f"- Year range: **{annual['year'].min()}-{annual['year'].max()}**")
    rec = int(cfg["quantile_regression"]["min_years_recommended_for_publication"])
    if n_years_total < rec:
        lines.append(f"- WARNING: Uploaded data are shorter than the recommended length for publication-grade trend inference (**{rec}+ years**).")
        lines.append("- Current outputs are therefore *computationally correct for the uploaded sample* but should be treated as a **pilot / code validation run**, not a final Q1 result.")
    lines.append("")
    lines.append("## Highest Delta1 stations by index")
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        sdf = features.loc[features["index_name"] == idx_name, ["station_name", "Delta1"]].dropna().sort_values("Delta1", ascending=False).head(5)
        lines.append(f"### {idx_cfg['title']}")
        if sdf.empty:
            lines.append("- No results")
        else:
            for _, r in sdf.iterrows():
                lines.append(f"- {r['station_name']}: Delta1 = {r['Delta1']:.3f}")
        lines.append("")
    if not cluster_df.empty:
        lines.append("## Cluster sizes")
        for idx_cfg in cfg["indices"]:
            idx_name = idx_cfg["name"]
            counts = cluster_df.loc[cluster_df["index_name"] == idx_name, "cluster"].value_counts().sort_index()
            lines.append(f"### {idx_cfg['title']}")
            for k, v in counts.items():
                lines.append(f"- Cluster {int(k)}: {int(v)} stations")
            lines.append("")
    (outdir / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Quantile regression + bootstrap + clustering pipeline for temperature extremes.")
    parser.add_argument("--config", type=str, default="config.yaml")
    args = parser.parse_args()

    cfg_path = Path(args.config)
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))

    outdir = Path(cfg["paths"]["output_dir"])
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(cfg["paths"]["data_csv"])
    stations = pd.read_csv(cfg["paths"]["station_csv"])

    # Build annual extreme indices
    daily_with_indices, annual = create_extreme_indices(data, cfg)

    # Run QR and bootstrap
    qr_all, qr_summary, boot_long = run_station_qr(annual, cfg)

    # Feature engineering + clustering
    feature_table = build_feature_table(qr_summary, cfg)
    cluster_df, artifacts = run_clustering(feature_table, cfg)

    # Join cluster assignments back to feature table
    if not cluster_df.empty:
        feature_table = feature_table.merge(cluster_df, on=["index_name", "station_id", "station_name"], how="left")

    # Export tables
    annual.to_csv(tables_dir / "annual_extreme_indices.csv", index=False)
    qr_all.to_csv(tables_dir / "qr_all_quantiles_long.csv", index=False)
    qr_summary.to_csv(tables_dir / "qr_focus_slopes_and_bootstrap_summary.csv", index=False)
    feature_table.to_csv(tables_dir / "clustering_feature_table.csv", index=False)
    cluster_df.to_csv(tables_dir / "cluster_assignments.csv", index=False)
    if cfg["bootstrap"]["save_long_table"] and not boot_long.empty:
        boot_long.to_csv(tables_dir / "bootstrap_distributions_long.csv", index=False)

    # Compact publication-ready summary tables
    focus_cols = [
        "index_name", "station_id", "station_name", "n_years",
        "slope_0.05", "slope_0.50", "slope_0.95",
        "Delta1", "Delta2", "Delta3",
        "boot_mean_0.05", "boot_sd_0.05", "boot_ci_low_0.05", "boot_ci_high_0.05",
        "boot_mean_0.50", "boot_sd_0.50", "boot_ci_low_0.50", "boot_ci_high_0.50",
        "boot_mean_0.95", "boot_sd_0.95", "boot_ci_low_0.95", "boot_ci_high_0.95",
        "boot_mean_Delta1", "boot_sd_Delta1", "boot_ci_low_Delta1", "boot_ci_high_Delta1",
        "cluster",
    ]
    focus_cols = [c for c in focus_cols if c in feature_table.columns]
    feature_table[focus_cols].to_csv(tables_dir / "publication_summary_table.csv", index=False)

    # Figures
    plot_data_coverage(annual, figs_dir, cfg)
    plot_region_quantile_slopes(annual, figs_dir, cfg)
    plot_station_heatmap(feature_table, figs_dir, cfg)
    plot_delta_uncertainty(feature_table, figs_dir, cfg)
    plot_dendrograms(feature_table, artifacts, figs_dir, cfg)
    plot_maps(feature_table, stations, cluster_df, figs_dir, cfg)

    # Metadata + report
    metadata = {
        "project_name": cfg["project"]["name"],
        "year_range": [int(annual["year"].min()), int(annual["year"].max())],
        "n_years": int(annual["year"].nunique()),
        "n_stations": int(annual["station_id"].nunique()),
        "bootstrap_reps": int(cfg["bootstrap"]["n_reps"]),
        "focus_quantiles": cfg["quantile_regression"]["focus_quantiles"],
        "publication_warning_short_record": bool(annual["year"].nunique() < int(cfg["quantile_regression"]["min_years_recommended_for_publication"])),
    }
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    generate_report(annual, feature_table, cluster_df, cfg, outdir)

    print(f"Done. Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
