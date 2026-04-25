from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Callable, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

from .config_utils import get_focus_quantiles, get_primary_delta, get_time_scale_years, slope_col, tau_label
from .indices import create_extreme_indices
from .plotting import (
    _draw_boundary,
    _format_geo_axis,
    _get_boundary_path_from_cfg,
    _get_plot_extent,
    _load_boundary_geometry,
    _panel_label,
    apply_publication_theme,
)
from .progress_utils import ProgressTracker
from .quantile import fit_ols_line, fit_quantile_slope, run_station_qr
from .year_config import build_split_periods, format_year_range_label, get_effective_year_range


def _expected_sign(index_name: str) -> int:
    return -1 if str(index_name).startswith("cool") else 1


def _is_expected_direction(index_name: str, value: float) -> bool:
    if not np.isfinite(value):
        return False
    return _expected_sign(index_name) * float(value) > 0


def _fast_no_bootstrap_cfg(cfg: dict) -> dict:
    cfg_alt = deepcopy(cfg)
    focus_quantiles = get_focus_quantiles(cfg_alt)
    cfg_alt["quantile_regression"]["full_quantiles"] = {
        "start": min(focus_quantiles),
        "stop": max(focus_quantiles),
        "step": max(0.01, round((max(focus_quantiles) - min(focus_quantiles)) / max(1, len(focus_quantiles) - 1), 2)),
    }
    cfg_alt["bootstrap"]["enabled"] = False
    return cfg_alt


def _default_signal_cfg(cfg: dict) -> dict:
    year_range = get_effective_year_range(cfg)
    periods = build_split_periods(year_range)
    baseline = [periods[0][1], periods[0][2]] if periods else list(year_range or (1991, 2007))
    comparison = [periods[1][1], periods[1][2]] if len(periods) > 1 else baseline
    user_cfg = cfg.get("advanced_analyses", {}).get("climate_change_signal", {})
    return {
        "enabled": bool(user_cfg.get("enabled", False)),
        "fixed_baseline_years": list(user_cfg.get("fixed_baseline_years", baseline)),
        "comparison_years": list(user_cfg.get("comparison_years", comparison)),
        "emergence_snr_threshold": float(user_cfg.get("emergence_snr_threshold", 2.0)),
        "fingerprint_permutations": int(user_cfg.get("fingerprint_permutations", 499)),
    }


def _annual_temperature_anomaly(data: pd.DataFrame, cfg: dict, baseline_years: list[int]) -> pd.DataFrame:
    dcfg = cfg["data"]
    station_col = dcfg["station_id_col"]
    station_name_col = dcfg["station_name_col"]
    year_col = dcfg["year_col"]
    tmean_col = dcfg.get("tmean_col", "tmean")
    tmin_col = dcfg.get("tmin_col", "tmin")
    tmax_col = dcfg.get("tmax_col", "tmax")

    out = data[[station_col, station_name_col, year_col, tmin_col, tmax_col, tmean_col]].copy()
    out[tmean_col] = pd.to_numeric(out[tmean_col], errors="coerce")
    fallback = (pd.to_numeric(out[tmin_col], errors="coerce") + pd.to_numeric(out[tmax_col], errors="coerce")) / 2.0
    out["temperature_mean"] = out[tmean_col].where(out[tmean_col].notna(), fallback)
    annual = (
        out.groupby([station_col, station_name_col, year_col], as_index=False)
        .agg(temperature_mean=("temperature_mean", "mean"), valid_temperature_days=("temperature_mean", lambda s: int(s.notna().sum())))
    )
    y0, y1 = int(baseline_years[0]), int(baseline_years[1])
    baseline = (
        annual.loc[annual[year_col].between(y0, y1)]
        .groupby([station_col, station_name_col], as_index=False)["temperature_mean"]
        .mean()
        .rename(columns={"temperature_mean": "station_baseline_temperature_c"})
    )
    annual = annual.merge(baseline, on=[station_col, station_name_col], how="left")
    annual["station_temperature_anomaly_c"] = annual["temperature_mean"] - annual["station_baseline_temperature_c"]
    regional = (
        annual.groupby(year_col, as_index=False)
        .agg(
            regional_temperature_anomaly_c=("station_temperature_anomaly_c", "mean"),
            regional_temperature_mean_c=("temperature_mean", "mean"),
            n_stations=("station_temperature_anomaly_c", lambda s: int(s.notna().sum())),
        )
        .rename(columns={year_col: "year"})
    )
    return regional


def _fit_response_to_warming(x: np.ndarray, y: np.ndarray, tau: float | None, max_iter: int) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(y) < 3 or np.nanstd(x) <= 0:
        return {"slope_per_c": np.nan, "ci_low_per_c": np.nan, "ci_high_per_c": np.nan}
    X = sm.add_constant(x)
    try:
        if tau is None:
            res = sm.OLS(y, X).fit()
            ci = res.conf_int()
        else:
            res = sm.QuantReg(y, X).fit(q=float(tau), max_iter=max_iter)
            ci = res.conf_int()
        return {
            "slope_per_c": float(res.params[1]),
            "ci_low_per_c": float(ci[1, 0]) if np.ndim(ci) == 2 else np.nan,
            "ci_high_per_c": float(ci[1, 1]) if np.ndim(ci) == 2 else np.nan,
        }
    except Exception:
        return {"slope_per_c": np.nan, "ci_low_per_c": np.nan, "ci_high_per_c": np.nan}


def _build_fixed_baseline_outputs(
    data: pd.DataFrame,
    cfg: dict,
    signal_cfg: dict,
    outdir: Path,
    progress_callback: Callable[[str], None] | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures" / "advanced_climate_change_signal"
    cfg_fixed = _fast_no_bootstrap_cfg(cfg)
    cfg_fixed["index_construction"]["reference_years"] = signal_cfg["fixed_baseline_years"]
    nested_callback = (
        (lambda message: progress_callback(f"Climate signal [fixed baseline] -> {message}"))
        if progress_callback is not None
        else None
    )
    _, fixed_annual = create_extreme_indices(data, cfg_fixed, progress_callback=nested_callback)
    fixed_annual.to_csv(tables_dir / "fixed_baseline_annual_extreme_indices.csv", index=False)

    _, fixed_qr, _ = run_station_qr(fixed_annual, cfg_fixed, progress_callback=nested_callback)
    fixed_qr.to_csv(tables_dir / "fixed_baseline_qr_summary.csv", index=False)

    station_col = cfg["data"]["station_id_col"]
    station_name_col = cfg["data"]["station_name_col"]
    year_col = cfg["data"]["year_col"]
    b0, b1 = [int(x) for x in signal_cfg["fixed_baseline_years"]]
    c0, c1 = [int(x) for x in signal_cfg["comparison_years"]]
    rows = []
    for index_cfg in cfg["indices"]:
        idx = index_cfg["name"]
        for (station_id, station_name), sdf in fixed_annual.groupby([station_col, station_name_col]):
            early = pd.to_numeric(sdf.loc[sdf[year_col].between(b0, b1), idx], errors="coerce").dropna()
            late = pd.to_numeric(sdf.loc[sdf[year_col].between(c0, c1), idx], errors="coerce").dropna()
            if early.empty or late.empty:
                continue
            delta = float(late.mean() - early.mean())
            rows.append(
                {
                    "index_name": idx,
                    "station_id": station_id,
                    "station_name": station_name,
                    "baseline_period": f"{b0}-{b1}",
                    "comparison_period": f"{c0}-{c1}",
                    "baseline_mean_days": float(early.mean()),
                    "comparison_mean_days": float(late.mean()),
                    "late_minus_baseline_days": delta,
                    "expected_direction": "increase" if _expected_sign(idx) > 0 else "decrease",
                    "direction_consistent": bool(_is_expected_direction(idx, delta)),
                }
            )
    station_change = pd.DataFrame(rows)
    summary_rows = []
    if not station_change.empty:
        for idx, sdf in station_change.groupby("index_name"):
            vals = pd.to_numeric(sdf["late_minus_baseline_days"], errors="coerce").dropna()
            summary_rows.append(
                {
                    "index_name": idx,
                    "n_stations": int(len(vals)),
                    "mean_late_minus_baseline_days": float(vals.mean()),
                    "median_late_minus_baseline_days": float(vals.median()),
                    "direction_consistent_stations": int(sdf["direction_consistent"].sum()),
                    "direction_consistency_fraction": float(sdf["direction_consistent"].mean()),
                }
            )
        station_change.to_csv(tables_dir / "fixed_baseline_period_change_station_level.csv", index=False)
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary.to_csv(tables_dir / "fixed_baseline_period_change_summary.csv", index=False)
        _plot_fixed_baseline_summary(summary, figs_dir / f"fixed_baseline_period_change_summary.{cfg['plots']['save_format']}", cfg)
    return fixed_annual, fixed_qr, summary


def _plot_fixed_baseline_summary(summary: pd.DataFrame, outpath: Path, cfg: dict) -> None:
    apply_publication_theme(cfg)
    order = [idx["name"] for idx in cfg["indices"]]
    plot_df = summary.set_index("index_name").reindex(order).dropna(subset=["mean_late_minus_baseline_days"])
    if plot_df.empty:
        return
    colors = ["#b2182b" if str(idx).startswith("warm") else "#2166ac" for idx in plot_df.index]
    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    bars = ax.bar(range(len(plot_df)), plot_df["mean_late_minus_baseline_days"], color=colors, edgecolor="black", linewidth=0.35)
    ax.axhline(0, color="#333333", linewidth=0.8)
    ax.set_xticks(range(len(plot_df)))
    ax.set_xticklabels([x.replace("_", " ").title() for x in plot_df.index], rotation=15, ha="right")
    ax.set_ylabel("Late minus baseline mean (days/year)")
    ax.set_title("Fixed-Baseline Extreme-Frequency Shift", loc="left")
    for bar, value, frac in zip(bars, plot_df["mean_late_minus_baseline_days"], plot_df["direction_consistency_fraction"]):
        va = "bottom" if value >= 0 else "top"
        y = value + (0.8 if value >= 0 else -0.8)
        ax.text(bar.get_x() + bar.get_width() / 2, y, f"{value:+.1f}\n{frac:.0%}", ha="center", va=va, fontsize=8.5)
    fig.savefig(outpath, dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def _build_warming_link_outputs(
    annual: pd.DataFrame,
    temp_anomaly: pd.DataFrame,
    cfg: dict,
    outdir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures" / "advanced_climate_change_signal"
    focus_quantiles = get_focus_quantiles(cfg)
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    station_col = cfg["data"]["station_id_col"]
    station_name_col = cfg["data"]["station_name_col"]
    year_col = cfg["data"]["year_col"]

    merged = annual.merge(temp_anomaly[["year", "regional_temperature_anomaly_c"]], left_on=year_col, right_on="year", how="left")
    station_rows = []
    for idx_cfg in cfg["indices"]:
        idx = idx_cfg["name"]
        for (station_id, station_name), sdf in merged.groupby([station_col, station_name_col]):
            x = pd.to_numeric(sdf["regional_temperature_anomaly_c"], errors="coerce").to_numpy(dtype=float)
            y = pd.to_numeric(sdf[idx], errors="coerce").to_numpy(dtype=float)
            ols = _fit_response_to_warming(x, y, tau=None, max_iter=max_iter)
            row = {
                "index_name": idx,
                "station_id": station_id,
                "station_name": station_name,
                "ols_slope_per_c": ols["slope_per_c"],
                "ols_expected_direction": bool(_is_expected_direction(idx, ols["slope_per_c"])),
            }
            for tau in focus_quantiles:
                fit = _fit_response_to_warming(x, y, tau=tau, max_iter=max_iter)
                suffix = f"{tau:0.2f}"
                row[f"slope_per_c_{suffix}"] = fit["slope_per_c"]
                row[f"ci_low_per_c_{suffix}"] = fit["ci_low_per_c"]
                row[f"ci_high_per_c_{suffix}"] = fit["ci_high_per_c"]
                row[f"expected_direction_{suffix}"] = bool(_is_expected_direction(idx, fit["slope_per_c"]))
            station_rows.append(row)
    station_df = pd.DataFrame(station_rows)
    if not station_df.empty:
        station_df.to_csv(tables_dir / "warming_link_station_quantile_response.csv", index=False)

    network_rows = []
    network = annual.groupby(year_col, as_index=False)[[idx["name"] for idx in cfg["indices"]]].mean()
    network = network.merge(temp_anomaly[["year", "regional_temperature_anomaly_c"]], left_on=year_col, right_on="year", how="left")
    for idx_cfg in cfg["indices"]:
        idx = idx_cfg["name"]
        x = pd.to_numeric(network["regional_temperature_anomaly_c"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(network[idx], errors="coerce").to_numpy(dtype=float)
        ols = _fit_response_to_warming(x, y, tau=None, max_iter=max_iter)
        network_rows.append(
            {
                "index_name": idx,
                "model": "OLS",
                "tau": np.nan,
                "slope_per_c": ols["slope_per_c"],
                "ci_low_per_c": ols["ci_low_per_c"],
                "ci_high_per_c": ols["ci_high_per_c"],
                "expected_direction": bool(_is_expected_direction(idx, ols["slope_per_c"])),
            }
        )
        for tau in focus_quantiles:
            fit = _fit_response_to_warming(x, y, tau=tau, max_iter=max_iter)
            network_rows.append(
                {
                    "index_name": idx,
                    "model": "QR",
                    "tau": tau,
                    "slope_per_c": fit["slope_per_c"],
                    "ci_low_per_c": fit["ci_low_per_c"],
                    "ci_high_per_c": fit["ci_high_per_c"],
                    "expected_direction": bool(_is_expected_direction(idx, fit["slope_per_c"])),
                }
            )
    network_df = pd.DataFrame(network_rows)
    if not network_df.empty:
        network_df.to_csv(tables_dir / "warming_link_network_quantile_response.csv", index=False)
        _plot_warming_link_response(network_df, figs_dir / f"warming_link_network_quantile_response.{cfg['plots']['save_format']}", cfg)
    return station_df, network_df


def _plot_warming_link_response(network_df: pd.DataFrame, outpath: Path, cfg: dict) -> None:
    apply_publication_theme(cfg)
    qr = network_df.loc[network_df["model"] == "QR"].copy()
    if qr.empty:
        return
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.8), constrained_layout=True)
    for i, idx_cfg in enumerate(cfg["indices"]):
        ax = axes.flat[i]
        sdf = qr.loc[qr["index_name"] == idx_cfg["name"]].sort_values("tau")
        if sdf.empty:
            continue
        ax.plot(sdf["tau"], sdf["slope_per_c"], color="#202020", linewidth=1.8)
        ax.scatter(sdf["tau"], sdf["slope_per_c"], s=46, color="#202020", zorder=3)
        ols = network_df.loc[(network_df["index_name"] == idx_cfg["name"]) & (network_df["model"] == "OLS")]
        if not ols.empty and np.isfinite(float(ols.iloc[0]["slope_per_c"])):
            ax.axhline(float(ols.iloc[0]["slope_per_c"]), color="#c62828", linewidth=1.5, linestyle="--", label="OLS")
        ax.axhline(0, color="#555555", linewidth=0.8)
        ax.set_xticks(sdf["tau"].tolist())
        ax.set_xticklabels([tau_label(t) for t in sdf["tau"]])
        ax.set_title(f"{_panel_label(i)} {idx_cfg['title']}", loc="left")
        ax.set_ylabel("Days/year per deg C")
        ax.set_xlabel("Quantile")
        if i == 0:
            ax.legend(frameon=False)
    fig.savefig(outpath, dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def _build_emergence_outputs(
    qr_summary: pd.DataFrame,
    stations: pd.DataFrame,
    cfg: dict,
    signal_cfg: dict,
    outdir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures" / "advanced_climate_change_signal"
    focus_quantiles = get_focus_quantiles(cfg)
    primary_delta = get_primary_delta(cfg)
    threshold = float(signal_cfg["emergence_snr_threshold"])
    rows = []
    metrics = [f"{tau:0.2f}" for tau in focus_quantiles] + [primary_delta]
    for _, row in qr_summary.iterrows():
        idx = row["index_name"]
        for metric in metrics:
            mean_col = f"boot_mean_{metric}"
            sd_col = f"boot_sd_{metric}"
            if mean_col not in qr_summary.columns or sd_col not in qr_summary.columns:
                continue
            signal = pd.to_numeric(pd.Series([row.get(mean_col)]), errors="coerce").iloc[0]
            sd = pd.to_numeric(pd.Series([row.get(sd_col)]), errors="coerce").iloc[0]
            snr = abs(signal / sd) if np.isfinite(signal) and np.isfinite(sd) and sd > 0 else np.nan
            rows.append(
                {
                    "index_name": idx,
                    "station_id": row["station_id"],
                    "station_name": row["station_name"],
                    "metric": metric,
                    "signal": signal,
                    "noise_sd": sd,
                    "signal_to_noise": snr,
                    "expected_direction": bool(_is_expected_direction(idx, signal)),
                    "emerged": bool(_is_expected_direction(idx, signal) and np.isfinite(snr) and snr >= threshold),
                    "snr_threshold": threshold,
                }
            )
    station_df = pd.DataFrame(rows)
    summary_rows = []
    if not station_df.empty:
        station_df.to_csv(tables_dir / "climate_signal_emergence_station_level.csv", index=False)
        for (idx, metric), sdf in station_df.groupby(["index_name", "metric"]):
            valid = sdf.dropna(subset=["signal_to_noise"])
            summary_rows.append(
                {
                    "index_name": idx,
                    "metric": metric,
                    "n_stations": int(len(valid)),
                    "expected_direction_stations": int(valid["expected_direction"].sum()),
                    "emerged_stations": int(valid["emerged"].sum()),
                    "expected_direction_fraction": float(valid["expected_direction"].mean()) if len(valid) else np.nan,
                    "emergence_fraction": float(valid["emerged"].mean()) if len(valid) else np.nan,
                    "median_signal_to_noise": float(valid["signal_to_noise"].median()) if len(valid) else np.nan,
                }
            )
    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        summary_df.to_csv(tables_dir / "climate_signal_emergence_summary.csv", index=False)
        _plot_emergence_map(station_df, stations, cfg, figs_dir / f"climate_signal_emergence_q50_maps.{cfg['plots']['save_format']}")
    return station_df, summary_df


def _plot_emergence_map(station_df: pd.DataFrame, stations: pd.DataFrame, cfg: dict, outpath: Path) -> None:
    apply_publication_theme(cfg)
    metric = "0.50"
    plot_df = station_df.loc[station_df["metric"] == metric].merge(stations, on=["station_id", "station_name"], how="left")
    if plot_df.empty:
        return
    boundary_geom = _load_boundary_geometry(_get_boundary_path_from_cfg(cfg))
    vmax = max(2.0, float(np.nanmax(plot_df["signal_to_noise"]))) if np.isfinite(plot_df["signal_to_noise"]).any() else 2.0
    fig, axes = plt.subplots(2, 2, figsize=(12.3, 9.8), constrained_layout=True)
    for i, idx_cfg in enumerate(cfg["indices"]):
        ax = axes.flat[i]
        sdf = plot_df.loc[plot_df["index_name"] == idx_cfg["name"]].copy()
        _draw_boundary(ax, boundary_geom, linewidth=0.95, zorder=2)
        scatter = ax.scatter(
            sdf["longitude"],
            sdf["latitude"],
            c=sdf["signal_to_noise"],
            s=np.where(sdf["emerged"], 86, 42),
            cmap="viridis",
            vmin=0,
            vmax=vmax,
            edgecolor=np.where(sdf["emerged"], "black", "#777777"),
            linewidth=0.45,
            zorder=3,
        )
        xmin, xmax, ymin, ymax = _get_plot_extent(sdf, boundary_geom=boundary_geom, pad_deg=0.15)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")
        _format_geo_axis(ax, show_x=(i >= 2), show_y=(i % 2 == 0))
        ax.set_title(f"{_panel_label(i)} {idx_cfg['title']}", loc="left")
        ax.grid(False)
    cbar = fig.colorbar(scatter, ax=axes.ravel().tolist(), shrink=0.86, pad=0.02)
    cbar.set_label("Signal-to-noise ratio at q0.50")
    fig.savefig(outpath, dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def _network_fingerprint_score(annual: pd.DataFrame, cfg: dict) -> tuple[float, list[dict[str, object]]]:
    focus_quantiles = get_focus_quantiles(cfg)
    time_scale_years = get_time_scale_years(cfg)
    year_col = cfg["data"]["year_col"]
    mean_df = annual.groupby(year_col, as_index=False)[[idx["name"] for idx in cfg["indices"]]].mean()
    rows = []
    for idx_cfg in cfg["indices"]:
        idx = idx_cfg["name"]
        years = mean_df[year_col].to_numpy(dtype=float)
        values = mean_df[idx].to_numpy(dtype=float)
        slopes = {}
        for tau in focus_quantiles:
            slopes[f"{tau:0.2f}"] = fit_quantile_slope(
                years,
                values,
                tau=tau,
                time_scale_years=time_scale_years,
                max_iter=int(cfg["quantile_regression"]["max_iter"]),
            )
        delta = slopes[f"{max(focus_quantiles):0.2f}"] - slopes[f"{min(focus_quantiles):0.2f}"]
        for tau, slope in slopes.items():
            rows.append({"index_name": idx, "metric": tau, "value": slope, "expected": _is_expected_direction(idx, slope)})
        rows.append({"index_name": idx, "metric": get_primary_delta(cfg), "value": delta, "expected": _is_expected_direction(idx, delta)})
    score = float(np.mean([row["expected"] for row in rows])) if rows else np.nan
    return score, rows


def _fingerprint_permutation_test(
    annual: pd.DataFrame,
    cfg: dict,
    n_permutations: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed_score, observed_components = _network_fingerprint_score(annual, cfg)
    year_col = cfg["data"]["year_col"]
    rng = np.random.default_rng(int(cfg["project"]["random_seed"]) + 2024)
    mean_df = annual.groupby(year_col, as_index=False)[[idx["name"] for idx in cfg["indices"]]].mean().sort_values(year_col)
    perm_rows = []
    for rep in range(int(n_permutations)):
        shifted = mean_df.copy()
        for idx_cfg in cfg["indices"]:
            vals = shifted[idx_cfg["name"]].to_numpy(dtype=float)
            if len(vals) > 1:
                offset = int(rng.integers(1, len(vals)))
                shifted[idx_cfg["name"]] = np.roll(vals, offset)
        score, _ = _network_fingerprint_score(shifted, cfg)
        perm_rows.append({"replicate": rep, "fingerprint_score": score})
    perm_df = pd.DataFrame(perm_rows)
    p_value = float((1 + (perm_df["fingerprint_score"] >= observed_score).sum()) / (len(perm_df) + 1)) if not perm_df.empty else np.nan
    observed_df = pd.DataFrame(observed_components)
    observed_df["observed_fingerprint_score"] = observed_score
    observed_df["permutation_p_value"] = p_value
    return observed_df, perm_df


def _safe_fraction(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator else np.nan


def _build_fingerprint_outputs(
    qr_summary: pd.DataFrame,
    fixed_summary: pd.DataFrame,
    warming_station: pd.DataFrame,
    emergence_summary: pd.DataFrame,
    fixed_annual: pd.DataFrame,
    cfg: dict,
    signal_cfg: dict,
    outdir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures" / "advanced_climate_change_signal"
    primary_delta = get_primary_delta(cfg)
    components = []

    for idx_cfg in cfg["indices"]:
        idx = idx_cfg["name"]
        sdf = qr_summary.loc[qr_summary["index_name"] == idx]
        if not sdf.empty:
            vals = pd.to_numeric(sdf[slope_col(0.50)], errors="coerce").dropna()
            components.append(
                {
                    "component": "station_q50_trend_direction",
                    "index_name": idx,
                    "component_value": _safe_fraction(sum(_is_expected_direction(idx, v) for v in vals), len(vals)),
                    "numerator": int(sum(_is_expected_direction(idx, v) for v in vals)),
                    "denominator": int(len(vals)),
                }
            )
            delta_vals = pd.to_numeric(sdf[primary_delta], errors="coerce").dropna()
            components.append(
                {
                    "component": "tail_asymmetry_direction",
                    "index_name": idx,
                    "component_value": _safe_fraction(sum(_is_expected_direction(idx, v) for v in delta_vals), len(delta_vals)),
                    "numerator": int(sum(_is_expected_direction(idx, v) for v in delta_vals)),
                    "denominator": int(len(delta_vals)),
                }
            )

    if not fixed_summary.empty:
        for _, row in fixed_summary.iterrows():
            components.append(
                {
                    "component": "fixed_baseline_period_shift",
                    "index_name": row["index_name"],
                    "component_value": float(row["direction_consistency_fraction"]),
                    "numerator": int(row["direction_consistent_stations"]),
                    "denominator": int(row["n_stations"]),
                }
            )

    if not warming_station.empty:
        q50_col = "expected_direction_0.50"
        for idx, sdf in warming_station.groupby("index_name"):
            components.append(
                {
                    "component": "warming_link_q50_direction",
                    "index_name": idx,
                    "component_value": float(sdf[q50_col].mean()) if q50_col in sdf.columns else np.nan,
                    "numerator": int(sdf[q50_col].sum()) if q50_col in sdf.columns else 0,
                    "denominator": int(len(sdf)),
                }
            )

    if not emergence_summary.empty:
        es = emergence_summary.loc[emergence_summary["metric"] == "0.50"]
        for _, row in es.iterrows():
            components.append(
                {
                    "component": "q50_signal_emergence",
                    "index_name": row["index_name"],
                    "component_value": float(row["emergence_fraction"]),
                    "numerator": int(row["emerged_stations"]),
                    "denominator": int(row["n_stations"]),
                }
            )

    fdr_path = tables_dir / "station_significance_fdr.csv"
    if fdr_path.exists():
        fdr = pd.read_csv(fdr_path)
        fdr = fdr.loc[np.isclose(pd.to_numeric(fdr["tau"], errors="coerce"), 0.50)].copy()
        for idx, sdf in fdr.groupby("index_name"):
            valid = sdf.loc[pd.to_numeric(sdf["fdr_reject"], errors="coerce").fillna(0) > 0].copy()
            numerator = sum(_is_expected_direction(idx, v) for v in pd.to_numeric(valid["slope"], errors="coerce").dropna())
            components.append(
                {
                    "component": "q50_fdr_field_significance",
                    "index_name": idx,
                    "component_value": _safe_fraction(numerator, sdf["station_id"].nunique()),
                    "numerator": int(numerator),
                    "denominator": int(sdf["station_id"].nunique()),
                }
            )

    hom_path = tables_dir / "homogeneity_flag_exclusion_sensitivity.csv"
    if hom_path.exists():
        hom = pd.read_csv(hom_path)
        all_col = f"{primary_delta}_all_stations"
        ex_col = f"{primary_delta}_exclude_flagged"
        if {all_col, ex_col, "index_name"}.issubset(hom.columns):
            numer = 0
            denom = 0
            for _, row in hom.iterrows():
                idx = row["index_name"]
                all_ok = _is_expected_direction(idx, float(row[all_col]))
                ex_ok = _is_expected_direction(idx, float(row[ex_col]))
                numer += int(all_ok and ex_ok)
                denom += 1
            components.append(
                {
                    "component": "homogeneity_exclusion_sign_stability",
                    "index_name": "all_indices",
                    "component_value": _safe_fraction(numer, denom),
                    "numerator": int(numer),
                    "denominator": int(denom),
                }
            )

    boot_depth_path = tables_dir / "bootstrap_depth_sensitivity_summary.csv"
    if boot_depth_path.exists():
        boot_depth = pd.read_csv(boot_depth_path)
        corr = pd.to_numeric(
            boot_depth.loc[boot_depth["metric"] == f"boot_mean_{primary_delta}", "station_correlation"],
            errors="coerce",
        ).dropna()
        if not corr.empty:
            components.append(
                {
                    "component": "bootstrap_depth_rank_stability",
                    "index_name": "all_indices",
                    "component_value": float(np.clip(corr.mean(), 0, 1)),
                    "numerator": float(corr.mean()),
                    "denominator": 1,
                }
            )

    component_df = pd.DataFrame(components)
    if not component_df.empty:
        overall = float(pd.to_numeric(component_df["component_value"], errors="coerce").dropna().mean())
        component_df["overall_fingerprint_score"] = overall
        component_df.to_csv(tables_dir / "climate_fingerprint_component_scores.csv", index=False)
    observed_df, perm_df = _fingerprint_permutation_test(
        fixed_annual if not fixed_annual.empty else qr_summary,
        cfg,
        int(signal_cfg["fingerprint_permutations"]),
    )
    observed_df.to_csv(tables_dir / "climate_fingerprint_network_components.csv", index=False)
    perm_df.to_csv(tables_dir / "climate_fingerprint_permutation_null.csv", index=False)
    if not component_df.empty:
        _plot_fingerprint_score(component_df, observed_df, perm_df, figs_dir / f"climate_fingerprint_score.{cfg['plots']['save_format']}", cfg)
    return component_df, observed_df, perm_df


def _plot_fingerprint_score(component_df: pd.DataFrame, observed_df: pd.DataFrame, perm_df: pd.DataFrame, outpath: Path, cfg: dict) -> None:
    apply_publication_theme(cfg)
    agg = component_df.groupby("component", as_index=False)["component_value"].mean().sort_values("component_value")
    overall = float(component_df["overall_fingerprint_score"].iloc[0])
    observed = float(observed_df["observed_fingerprint_score"].iloc[0]) if not observed_df.empty else np.nan
    p_value = float(observed_df["permutation_p_value"].iloc[0]) if not observed_df.empty else np.nan
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2), constrained_layout=True)
    ax = axes[0]
    bars = ax.barh(range(len(agg)), agg["component_value"], color="#3b7ea1", edgecolor="black", linewidth=0.3)
    ax.axvline(overall, color="#b2182b", linewidth=1.8, linestyle="--", label=f"Overall={overall:.2f}")
    ax.set_yticks(range(len(agg)))
    ax.set_yticklabels(agg["component"].str.replace("_", " ").str.title().tolist(), fontsize=8.5)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Component score")
    ax.set_title("(a) Climate-Fingerprint Component Scores", loc="left")
    ax.legend(frameon=False)
    for bar, value in zip(bars, agg["component_value"]):
        ax.text(value + 0.015, bar.get_y() + bar.get_height() / 2, f"{value:.2f}", va="center", fontsize=8.5)

    ax = axes[1]
    if not perm_df.empty:
        ax.hist(perm_df["fingerprint_score"], bins=np.linspace(0, 1, 16), color="#bdbdbd", edgecolor="white")
    if np.isfinite(observed):
        ax.axvline(observed, color="#b2182b", linewidth=2.0)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Network fingerprint score under circular-shift null")
    ax.set_ylabel("Frequency")
    ax.set_title(f"(b) Permutation Test (observed={observed:.2f}, p={p_value:.3f})", loc="left")
    fig.savefig(outpath, dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def run_climate_change_signal_analysis(
    data: pd.DataFrame,
    annual: pd.DataFrame,
    qr_summary: pd.DataFrame,
    feature_table: pd.DataFrame,
    stations: pd.DataFrame,
    cfg: dict,
    outdir: Path,
    progress_callback: Callable[[str], None] | None = None,
) -> Dict[str, pd.DataFrame]:
    signal_cfg = _default_signal_cfg(cfg)
    if not signal_cfg["enabled"]:
        return {}

    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures" / "advanced_climate_change_signal"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)
    results: Dict[str, pd.DataFrame] = {}

    tasks = [
        "regional temperature anomaly",
        "fixed-baseline indices",
        "warming-linked response",
        "signal emergence",
        "fingerprint score",
    ]
    tracker = ProgressTracker("Climate-change signal analysis", len(tasks), progress_callback)

    tracker.emit(1, detail="regional temperature anomaly")
    temp_anomaly = _annual_temperature_anomaly(data, cfg, signal_cfg["fixed_baseline_years"])
    temp_anomaly.to_csv(tables_dir / "regional_temperature_anomaly.csv", index=False)
    temp_trend = fit_ols_line(
        temp_anomaly["year"].to_numpy(dtype=float),
        temp_anomaly["regional_temperature_anomaly_c"].to_numpy(dtype=float),
        time_scale_years=get_time_scale_years(cfg),
    )
    temp_summary = pd.DataFrame(
        [
            {
                "baseline_period": format_year_range_label(tuple(signal_cfg["fixed_baseline_years"])),
                "temperature_trend_c_per_decade": float(temp_trend["slope"]),
                "temperature_trend_ci_low": float(temp_trend["ci_low"]),
                "temperature_trend_ci_high": float(temp_trend["ci_high"]),
                "first_year_anomaly_c": float(temp_anomaly.sort_values("year").iloc[0]["regional_temperature_anomaly_c"]),
                "last_year_anomaly_c": float(temp_anomaly.sort_values("year").iloc[-1]["regional_temperature_anomaly_c"]),
                "last_minus_first_anomaly_c": float(
                    temp_anomaly.sort_values("year").iloc[-1]["regional_temperature_anomaly_c"]
                    - temp_anomaly.sort_values("year").iloc[0]["regional_temperature_anomaly_c"]
                ),
            }
        ]
    )
    temp_summary.to_csv(tables_dir / "regional_temperature_anomaly_summary.csv", index=False)
    results["regional_temperature_anomaly"] = temp_anomaly
    results["regional_temperature_anomaly_summary"] = temp_summary

    tracker.emit(2, detail="fixed-baseline indices")
    fixed_annual, fixed_qr, fixed_summary = _build_fixed_baseline_outputs(data, cfg, signal_cfg, outdir, progress_callback)
    results["fixed_baseline_annual_extreme_indices"] = fixed_annual
    results["fixed_baseline_qr_summary"] = fixed_qr
    results["fixed_baseline_period_change_summary"] = fixed_summary

    tracker.emit(3, detail="warming-linked response")
    warming_station, warming_network = _build_warming_link_outputs(fixed_annual, temp_anomaly, cfg, outdir)
    results["warming_link_station_quantile_response"] = warming_station
    results["warming_link_network_quantile_response"] = warming_network

    tracker.emit(4, detail="signal emergence")
    emergence_station, emergence_summary = _build_emergence_outputs(qr_summary, stations, cfg, signal_cfg, outdir)
    results["climate_signal_emergence_station_level"] = emergence_station
    results["climate_signal_emergence_summary"] = emergence_summary

    tracker.emit(5, detail="fingerprint score")
    component_df, network_components, perm_df = _build_fingerprint_outputs(
        qr_summary,
        fixed_summary,
        warming_station,
        emergence_summary,
        fixed_annual,
        cfg,
        signal_cfg,
        outdir,
    )
    results["climate_fingerprint_component_scores"] = component_df
    results["climate_fingerprint_network_components"] = network_components
    results["climate_fingerprint_permutation_null"] = perm_df
    return results
