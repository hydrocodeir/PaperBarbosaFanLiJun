from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm, spearmanr
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

from .indices import create_extreme_indices
from .plotting import (
    FIG3_CMAP,
    _build_interpolation_grid,
    _draw_boundary,
    _fig3_levels_from_values,
    _format_geo_axis,
    _get_boundary_path_from_cfg,
    _get_plot_extent,
    _interpolate_station_surface,
    _load_boundary_geometry,
    _mask_surface_to_boundary,
    _panel_label,
    _short_tau_label,
    apply_publication_theme,
)
from .quantile import add_sensitivity_check_columns, run_station_qr


def _analytic_pvalue_from_ci(slope: float, ci_low: float, ci_high: float) -> float:
    if not (np.isfinite(slope) and np.isfinite(ci_low) and np.isfinite(ci_high)):
        return np.nan
    width = ci_high - ci_low
    if not np.isfinite(width) or width <= 0:
        return np.nan
    se = width / (2 * 1.96)
    if not np.isfinite(se) or se <= 0:
        return np.nan
    z_score = abs(float(slope) / se)
    return float(2 * (1 - norm.cdf(z_score)))


def _knn_weights(coords: np.ndarray, k: int) -> np.ndarray:
    n = len(coords)
    if n < 2:
        return np.zeros((n, n), dtype=float)
    dmat = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(dmat, np.inf)
    k = max(1, min(k, n - 1))
    weights = np.zeros((n, n), dtype=float)
    for i in range(n):
        nn_idx = np.argsort(dmat[i])[:k]
        inv_d = 1.0 / np.maximum(dmat[i, nn_idx], 1e-12)
        inv_sum = inv_d.sum()
        if inv_sum > 0:
            weights[i, nn_idx] = inv_d / inv_sum
    return weights


def _moran_i(values: np.ndarray, coords: np.ndarray, k_neighbors: int, permutations: int, rng: np.random.Generator) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    coords = np.asarray(coords, dtype=float)
    mask = np.isfinite(values) & np.all(np.isfinite(coords), axis=1)
    values = values[mask]
    coords = coords[mask]
    if len(values) < 4:
        return np.nan, np.nan

    weights = _knn_weights(coords, k_neighbors)
    w_sum = weights.sum()
    if w_sum <= 0:
        return np.nan, np.nan
    centered = values - values.mean()
    denom = np.sum(centered ** 2)
    if denom <= 0:
        return np.nan, np.nan
    obs = float((len(values) / w_sum) * ((weights * np.outer(centered, centered)).sum() / denom))

    permuted = []
    for _ in range(permutations):
        shuffled = rng.permutation(centered)
        permuted.append(float((len(values) / w_sum) * ((weights * np.outer(shuffled, shuffled)).sum() / denom)))
    permuted_arr = np.asarray(permuted, dtype=float)
    p_value = float((1 + np.sum(np.abs(permuted_arr) >= abs(obs))) / (permutations + 1))
    return obs, p_value


def _plot_two_heatmaps(left: pd.DataFrame, right: pd.DataFrame, left_title: str, right_title: str, outpath: Path, cmap: str = "viridis") -> None:
    apply_publication_theme()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    for ax, data, title in zip(axes, [left, right], [left_title, right_title]):
        im = ax.imshow(data.to_numpy(dtype=float), aspect="auto", cmap=cmap)
        ax.set_xticks(range(data.shape[1]))
        ax.set_xticklabels([f"{float(c):.2f}" for c in data.columns], rotation=0)
        ax.set_yticks(range(data.shape[0]))
        ax.set_yticklabels(data.index.tolist())
        ax.set_xlabel("Quantile τ")
        ax.set_title(title)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                val = data.iloc[i, j]
                if pd.notna(val):
                    ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=8.5, color="white")
        fig.colorbar(im, ax=ax, shrink=0.82)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def _plot_single_heatmap(data: pd.DataFrame, title: str, outpath: Path, cmap: str = "coolwarm", fmt: str = "{:.2f}") -> None:
    apply_publication_theme()
    fig, ax = plt.subplots(figsize=(7.8, 5.8), constrained_layout=True)
    im = ax.imshow(data.to_numpy(dtype=float), aspect="auto", cmap=cmap)
    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels([f"{float(c):.2f}" if isinstance(c, (float, int)) else str(c) for c in data.columns], rotation=0)
    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels(data.index.tolist())
    ax.set_title(title)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data.iloc[i, j]
            if pd.notna(val):
                ax.text(j, i, fmt.format(float(val)), ha="center", va="center", fontsize=8.5, color="black")
    fig.colorbar(im, ax=ax, shrink=0.86)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def run_spatial_inference(qr_summary: pd.DataFrame, stations: pd.DataFrame, cfg: dict, outdir: Path) -> Dict[str, pd.DataFrame]:
    analysis_cfg = cfg.get("advanced_analyses", {}).get("spatial_inference", {})
    if not analysis_cfg.get("enabled", False) or qr_summary.empty:
        return {}

    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures" / "advanced_spatial_inference"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)
    quantiles = [float(x) for x in analysis_cfg.get("quantiles", cfg["quantile_regression"]["focus_quantiles"])]
    alpha = float(analysis_cfg.get("fdr_alpha", 0.05))
    moran_perms = int(analysis_cfg.get("moran_permutations", 499))
    k_neighbors = int(analysis_cfg.get("moran_k_neighbors", 5))
    rng = np.random.default_rng(int(cfg["project"]["random_seed"]))

    merged = qr_summary.merge(stations, on=["station_id", "station_name"], how="left")
    station_rows = []
    moran_rows = []

    for tau in quantiles:
        suffix = f"{tau:0.2f}"
        slope_col = f"slope_{suffix}"
        ci_low_col = f"ci_low_{suffix}"
        ci_high_col = f"ci_high_{suffix}"
        if slope_col not in merged.columns:
            continue

        for idx_name, sdf in merged.groupby("index_name"):
            local = sdf.copy()
            if ci_low_col in local.columns and ci_high_col in local.columns:
                local["analytic_p"] = [
                    _analytic_pvalue_from_ci(slope, lo, hi)
                    for slope, lo, hi in zip(local[slope_col], local[ci_low_col], local[ci_high_col])
                ]
            else:
                local["analytic_p"] = np.nan

            valid_p = pd.to_numeric(local["analytic_p"], errors="coerce")
            qvals = np.full(len(local), np.nan, dtype=float)
            reject = np.full(len(local), np.nan, dtype=float)
            mask = np.isfinite(valid_p.to_numpy(dtype=float))
            if mask.any():
                rej_mask, q_mask, _, _ = multipletests(valid_p.to_numpy(dtype=float)[mask], alpha=alpha, method="fdr_bh")
                qvals[np.where(mask)[0]] = q_mask
                reject[np.where(mask)[0]] = rej_mask.astype(float)
            local["fdr_q"] = qvals
            local["fdr_reject"] = reject
            local["tau"] = tau
            station_rows.append(
                local[[
                    "index_name", "station_id", "station_name", "latitude", "longitude", "tau",
                    slope_col, ci_low_col, ci_high_col, "analytic_p", "fdr_q", "fdr_reject",
                ]].rename(columns={slope_col: "slope", ci_low_col: "ci_low", ci_high_col: "ci_high"})
            )

            moran_local = local[["latitude", "longitude", slope_col]].dropna().copy()
            moran_i, moran_p = _moran_i(
                moran_local[slope_col].to_numpy(dtype=float),
                moran_local[["longitude", "latitude"]].to_numpy(dtype=float),
                k_neighbors=k_neighbors,
                permutations=moran_perms,
                rng=rng,
            )
            moran_rows.append(
                {
                    "index_name": idx_name,
                    "tau": tau,
                    "n_stations": int(len(moran_local)),
                    "moran_i": moran_i,
                    "moran_p_perm": moran_p,
                }
            )

    station_df = pd.concat(station_rows, ignore_index=True) if station_rows else pd.DataFrame()
    moran_df = pd.DataFrame(moran_rows)
    if not station_df.empty:
        station_df.to_csv(tables_dir / "station_significance_fdr.csv", index=False)
    if not moran_df.empty:
        moran_df.to_csv(tables_dir / "spatial_autocorrelation_moran.csv", index=False)

    if not station_df.empty:
        raw_counts = (
            station_df.assign(raw_sig=lambda d: (pd.to_numeric(d["analytic_p"], errors="coerce") < alpha).astype(float))
            .groupby(["index_name", "tau"], as_index=False)[["raw_sig", "fdr_reject"]]
            .sum()
        )
        raw_pivot = raw_counts.pivot(index="index_name", columns="tau", values="raw_sig").sort_index()
        fdr_pivot = raw_counts.pivot(index="index_name", columns="tau", values="fdr_reject").sort_index()
        _plot_two_heatmaps(
            raw_pivot,
            fdr_pivot,
            "Raw Significant Counts",
            "FDR-Rejected Counts",
            figs_dir / f"raw_vs_fdr_counts.{cfg['plots']['save_format']}",
            cmap="YlGnBu",
        )
    if not moran_df.empty:
        moran_pivot = moran_df.pivot(index="index_name", columns="tau", values="moran_i").sort_index()
        _plot_single_heatmap(
            moran_pivot,
            "Moran's I of Station Slopes",
            figs_dir / f"moran_i_heatmap.{cfg['plots']['save_format']}",
            cmap="RdBu_r",
            fmt="{:.2f}",
        )

    return {"station_significance_fdr": station_df, "spatial_autocorrelation_moran": moran_df}


def _prepare_fast_sensitivity_cfg(cfg: dict, bootstrap_enabled: bool, bootstrap_method: str | None = None) -> dict:
    cfg_alt = deepcopy(cfg)
    cfg_alt["quantile_regression"]["full_quantiles"] = {"start": 0.05, "stop": 0.95, "step": 0.45}
    cfg_alt["bootstrap"]["enabled"] = bool(bootstrap_enabled)
    if bootstrap_method is not None:
        cfg_alt["bootstrap"]["method"] = bootstrap_method
    return cfg_alt


def _summarize_metric_comparison(base_df: pd.DataFrame, alt_df: pd.DataFrame, metrics: list[str], label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    keys = ["index_name", "station_id", "station_name"]
    merged = base_df[keys + metrics].merge(
        alt_df[keys + metrics],
        on=keys,
        how="inner",
        suffixes=("_base", "_alt"),
    )
    summary_rows = []
    for idx_name, sdf in merged.groupby("index_name"):
        for metric in metrics:
            base_col = f"{metric}_base"
            alt_col = f"{metric}_alt"
            valid = sdf[[base_col, alt_col]].dropna()
            if valid.empty:
                continue
            summary_rows.append(
                {
                    "comparison": label,
                    "index_name": idx_name,
                    "metric": metric,
                    "n_pairs": int(len(valid)),
                    "mean_base": float(valid[base_col].mean()),
                    "mean_alt": float(valid[alt_col].mean()),
                    "mean_abs_diff": float(np.mean(np.abs(valid[alt_col] - valid[base_col]))),
                    "correlation": float(np.corrcoef(valid[base_col], valid[alt_col])[0, 1]) if len(valid) >= 2 else np.nan,
                }
            )
            sdf[f"{metric}_diff"] = sdf[alt_col] - sdf[base_col]
    merged["comparison"] = label
    return merged, pd.DataFrame(summary_rows)


def _plot_metric_sensitivity(summary_df: pd.DataFrame, title: str, outpath: Path) -> None:
    if summary_df.empty:
        return
    apply_publication_theme()
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    plot_df = summary_df.copy()
    plot_df["label"] = plot_df["index_name"] + " | " + plot_df["metric"]
    ax.barh(np.arange(len(plot_df)), plot_df["mean_abs_diff"], color="#31688e")
    ax.set_yticks(np.arange(len(plot_df)))
    ax.set_yticklabels(plot_df["label"].tolist(), fontsize=8.5)
    ax.set_xlabel("Mean absolute difference")
    ax.set_title(title)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def _plot_interpolation_comparison(qr_summary: pd.DataFrame, stations: pd.DataFrame, cfg: dict, methods: list[str], outdir: Path) -> pd.DataFrame:
    figs_dir = outdir / "figures" / "advanced_method_sensitivity"
    figs_dir.mkdir(parents=True, exist_ok=True)
    boundary_geom = _load_boundary_geometry(_get_boundary_path_from_cfg(cfg))
    taus = [float(x) for x in cfg["quantile_regression"].get("sensitivity_check_quantiles", [0.05, 0.95])]
    rows = []
    for tau in taus:
        suffix = f"{tau:0.2f}"
        slope_col = f"slope_{suffix}"
        fig, axes = plt.subplots(len(cfg["indices"]), len(methods), figsize=(5.0 * len(methods), 3.8 * len(cfg["indices"])), constrained_layout=True)
        if len(cfg["indices"]) == 1:
            axes = np.asarray([axes])
        if len(methods) == 1:
            axes = axes[:, None]

        for row_idx, idx_cfg in enumerate(cfg["indices"]):
            idx_name = idx_cfg["name"]
            sdf = qr_summary.loc[qr_summary["index_name"] == idx_name, ["station_id", "station_name", slope_col]].merge(
                stations,
                on=["station_id", "station_name"],
                how="left",
            )
            interp_df = sdf[["longitude", "latitude", slope_col]].rename(columns={slope_col: "slope_value"}).dropna().copy()
            if interp_df.empty:
                continue
            surfaces = {}
            grid_x, grid_y = _build_interpolation_grid(interp_df, boundary_geom=boundary_geom, nx=220, ny=220)
            levels = _fig3_levels_from_values(sdf[slope_col])
            for col_idx, method in enumerate(methods):
                ax = axes[row_idx, col_idx]
                surface = _interpolate_station_surface(
                    interp_df,
                    grid_x,
                    grid_y,
                    method=method,
                    smooth=float(cfg.get("spatial_visualization", {}).get("interpolation_smooth", 0.35)),
                )
                surface = _mask_surface_to_boundary(grid_x, grid_y, surface, boundary_geom=boundary_geom)
                surfaces[method] = surface
                ax.contourf(grid_x, grid_y, surface, levels=levels, cmap=FIG3_CMAP, extend="both")
                _draw_boundary(ax, boundary_geom, linewidth=0.8)
                ax.scatter(sdf["longitude"], sdf["latitude"], c=sdf[slope_col], cmap=FIG3_CMAP, s=14, edgecolor="none", zorder=4)
                xmin, xmax, ymin, ymax = _get_plot_extent(sdf, boundary_geom=boundary_geom, pad_deg=0.1)
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.set_aspect("equal", adjustable="box")
                _format_geo_axis(ax, show_x=(row_idx == len(cfg["indices"]) - 1), show_y=(col_idx == 0))
                ax.set_title(f"{_panel_label(row_idx * len(methods) + col_idx)} {idx_cfg['title']} | {method}", loc="left", fontsize=10.5)

            method_pairs = [(methods[i], methods[j]) for i in range(len(methods)) for j in range(i + 1, len(methods))]
            for left, right in method_pairs:
                a = surfaces[left].ravel()
                b = surfaces[right].ravel()
                mask = np.isfinite(a) & np.isfinite(b)
                rows.append(
                    {
                        "index_name": idx_name,
                        "tau": tau,
                        "method_left": left,
                        "method_right": right,
                        "n_cells": int(mask.sum()),
                        "surface_correlation": float(np.corrcoef(a[mask], b[mask])[0, 1]) if mask.sum() >= 2 else np.nan,
                        "surface_rmse": float(np.sqrt(np.mean((a[mask] - b[mask]) ** 2))) if mask.sum() else np.nan,
                    }
                )
        fig.suptitle(f"Interpolation Sensitivity Maps at τ = {_short_tau_label(tau)}", fontsize=14, fontweight="bold")
        fig.savefig(figs_dir / f"interpolation_sensitivity_tau_{tau:0.2f}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)
    return pd.DataFrame(rows)


def run_method_sensitivity(data: pd.DataFrame, annual: pd.DataFrame, qr_summary: pd.DataFrame, stations: pd.DataFrame, cfg: dict, outdir: Path) -> Dict[str, pd.DataFrame]:
    analysis_cfg = cfg.get("advanced_analyses", {}).get("method_sensitivity", {})
    if not analysis_cfg.get("enabled", False):
        return {}

    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures" / "advanced_method_sensitivity"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)
    results: Dict[str, pd.DataFrame] = {}

    metric_list = ["slope_0.05", "slope_0.50", "slope_0.95", "Delta1"]
    ref_periods = analysis_cfg.get("reference_periods", {})
    current_label = "current"
    current_ref = ref_periods.get(current_label, cfg["index_construction"]["reference_years"])
    base_qr = qr_summary.copy()

    reference_station_parts = []
    reference_summary_parts = []
    for label, ref_years in ref_periods.items():
        if label == current_label or ref_years == current_ref:
            continue
        cfg_alt = _prepare_fast_sensitivity_cfg(cfg, bootstrap_enabled=False)
        cfg_alt["index_construction"]["reference_years"] = ref_years
        _, annual_alt = create_extreme_indices(data, cfg_alt)
        _, qr_alt, _ = run_station_qr(annual_alt, cfg_alt)
        station_cmp, summary_cmp = _summarize_metric_comparison(base_qr, qr_alt, metric_list, f"reference_period:{label}")
        station_cmp["alternative"] = label
        summary_cmp["alternative"] = label
        reference_station_parts.append(station_cmp)
        reference_summary_parts.append(summary_cmp)

    reference_station_df = pd.concat(reference_station_parts, ignore_index=True) if reference_station_parts else pd.DataFrame()
    reference_summary_df = pd.concat(reference_summary_parts, ignore_index=True) if reference_summary_parts else pd.DataFrame()
    if not reference_station_df.empty:
        reference_station_df.to_csv(tables_dir / "reference_period_sensitivity_station_level.csv", index=False)
        results["reference_period_sensitivity_station_level"] = reference_station_df
    if not reference_summary_df.empty:
        reference_summary_df.to_csv(tables_dir / "reference_period_sensitivity_summary.csv", index=False)
        _plot_metric_sensitivity(reference_summary_df, "Reference-Period Sensitivity", figs_dir / f"reference_period_sensitivity.{cfg['plots']['save_format']}")
        results["reference_period_sensitivity_summary"] = reference_summary_df

    bootstrap_station_parts = []
    bootstrap_summary_parts = []
    bootstrap_methods = [str(x) for x in analysis_cfg.get("bootstrap_methods", [cfg["bootstrap"]["method"]])]
    current_method = str(cfg["bootstrap"]["method"]).lower()
    boot_metrics = [
        "boot_mean_0.05", "boot_ci_low_0.05", "boot_ci_high_0.05",
        "boot_mean_0.95", "boot_ci_low_0.95", "boot_ci_high_0.95",
    ]
    for method in bootstrap_methods:
        if method.lower() == current_method:
            continue
        cfg_alt = _prepare_fast_sensitivity_cfg(cfg, bootstrap_enabled=True, bootstrap_method=method)
        _, qr_alt, _ = run_station_qr(annual, cfg_alt)
        qr_alt = add_sensitivity_check_columns(qr_alt, cfg_alt)
        station_cmp, summary_cmp = _summarize_metric_comparison(base_qr, qr_alt, boot_metrics, f"bootstrap_method:{method}")
        station_cmp["alternative"] = method
        summary_cmp["alternative"] = method
        bootstrap_station_parts.append(station_cmp)
        bootstrap_summary_parts.append(summary_cmp)

    bootstrap_station_df = pd.concat(bootstrap_station_parts, ignore_index=True) if bootstrap_station_parts else pd.DataFrame()
    bootstrap_summary_df = pd.concat(bootstrap_summary_parts, ignore_index=True) if bootstrap_summary_parts else pd.DataFrame()
    if not bootstrap_station_df.empty:
        bootstrap_station_df.to_csv(tables_dir / "bootstrap_method_sensitivity_station_level.csv", index=False)
        results["bootstrap_method_sensitivity_station_level"] = bootstrap_station_df
    if not bootstrap_summary_df.empty:
        bootstrap_summary_df.to_csv(tables_dir / "bootstrap_method_sensitivity_summary.csv", index=False)
        _plot_metric_sensitivity(bootstrap_summary_df, "Bootstrap-Method Sensitivity", figs_dir / f"bootstrap_method_sensitivity.{cfg['plots']['save_format']}")
        results["bootstrap_method_sensitivity_summary"] = bootstrap_summary_df

    interpolation_methods = [str(x) for x in analysis_cfg.get("interpolation_methods", [cfg["spatial_visualization"]["interpolation_method"]])]
    interpolation_df = _plot_interpolation_comparison(qr_summary, stations, cfg, interpolation_methods, outdir)
    if not interpolation_df.empty:
        interpolation_df.to_csv(tables_dir / "interpolation_method_sensitivity_summary.csv", index=False)
        results["interpolation_method_sensitivity_summary"] = interpolation_df

    return results


def run_driver_analysis(feature_table: pd.DataFrame, stations: pd.DataFrame, cfg: dict, outdir: Path) -> Dict[str, pd.DataFrame]:
    analysis_cfg = cfg.get("advanced_analyses", {}).get("driver_analysis", {})
    if not analysis_cfg.get("enabled", False) or feature_table.empty:
        return {}

    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures" / "advanced_driver_analysis"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)
    predictors = list(analysis_cfg.get("predictors", ["latitude", "longitude", "elevation"]))
    metrics = list(analysis_cfg.get("metrics", ["slope_0.05", "slope_0.50", "slope_0.95", "Delta1"]))
    merged = feature_table.merge(stations, on=["station_id", "station_name"], how="left")

    rows = []
    for idx_name, sdf in merged.groupby("index_name"):
        for metric in metrics:
            use_cols = predictors + [metric]
            local = sdf[use_cols].apply(pd.to_numeric, errors="coerce").dropna()
            if len(local) < len(predictors) + 2:
                continue
            X_raw = local[predictors].to_numpy(dtype=float)
            y_raw = local[metric].to_numpy(dtype=float)
            x_scaler = StandardScaler()
            y_scaler = StandardScaler()
            X_std = x_scaler.fit_transform(X_raw)
            y_std = y_scaler.fit_transform(y_raw.reshape(-1, 1)).ravel()
            res = sm.OLS(y_std, sm.add_constant(X_std)).fit()
            for pred_idx, predictor in enumerate(predictors):
                rho, rho_p = spearmanr(local[predictor], local[metric])
                rows.append(
                    {
                        "index_name": idx_name,
                        "metric": metric,
                        "predictor": predictor,
                        "std_beta": float(res.params[pred_idx + 1]),
                        "std_beta_pvalue": float(res.pvalues[pred_idx + 1]),
                        "model_r2": float(res.rsquared),
                        "spearman_rho": float(rho),
                        "spearman_pvalue": float(rho_p),
                        "n_stations": int(len(local)),
                    }
                )

    driver_df = pd.DataFrame(rows)
    if not driver_df.empty:
        driver_df.to_csv(tables_dir / "driver_analysis_summary.csv", index=False)

        delta_df = driver_df.loc[driver_df["metric"] == "Delta1"].pivot(index="index_name", columns="predictor", values="std_beta").sort_index()
        if not delta_df.empty:
            _plot_single_heatmap(
                delta_df,
                "Standardized Driver Effects on Delta1",
                figs_dir / f"driver_effects_delta1_heatmap.{cfg['plots']['save_format']}",
                cmap="RdBu_r",
                fmt="{:.2f}",
            )

        apply_publication_theme()
        for idx_cfg in cfg["indices"]:
            idx_name = idx_cfg["name"]
            sdf = merged.loc[merged["index_name"] == idx_name].copy()
            fig, axes = plt.subplots(1, len(predictors), figsize=(5.0 * len(predictors), 4.3), constrained_layout=True)
            if len(predictors) == 1:
                axes = [axes]
            for i, predictor in enumerate(predictors):
                ax = axes[i]
                local = sdf[[predictor, "Delta1"]].apply(pd.to_numeric, errors="coerce").dropna()
                ax.scatter(local[predictor], local["Delta1"], s=42, color="#2a6f97", alpha=0.85)
                if len(local) >= 2:
                    coeff = np.polyfit(local[predictor], local["Delta1"], deg=1)
                    xs = np.linspace(local[predictor].min(), local[predictor].max(), 100)
                    ax.plot(xs, coeff[0] * xs + coeff[1], color="#c1121f", linewidth=1.8)
                driver_row = driver_df.loc[
                    (driver_df["index_name"] == idx_name) & (driver_df["metric"] == "Delta1") & (driver_df["predictor"] == predictor)
                ]
                rho_text = ""
                if not driver_row.empty:
                    rho_text = f"rho={float(driver_row.iloc[0]['spearman_rho']):+.2f}\np={float(driver_row.iloc[0]['spearman_pvalue']):.3f}"
                ax.set_title(f"{_panel_label(i)} {predictor}")
                ax.set_xlabel(predictor)
                ax.set_ylabel("Delta1")
                ax.text(0.02, 0.98, rho_text, transform=ax.transAxes, va="top", ha="left", fontsize=8.5)
            fig.suptitle(f"Driver Relationships for {idx_cfg['title']}", fontsize=14, fontweight="bold")
            fig.savefig(figs_dir / f"driver_scatter_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
            plt.close(fig)
    return {"driver_analysis_summary": driver_df}


def run_regionalization_analysis(feature_table: pd.DataFrame, cfg: dict, outdir: Path) -> Dict[str, pd.DataFrame]:
    analysis_cfg = cfg.get("advanced_analyses", {}).get("regionalization", {})
    if not analysis_cfg.get("enabled", False) or feature_table.empty or "cluster" not in feature_table.columns:
        return {}

    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures" / "advanced_regionalization"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)
    metrics = list(analysis_cfg.get("summary_metrics", ["slope_0.05", "slope_0.50", "slope_0.95", "Delta1"]))
    cluster_df = feature_table.dropna(subset=["cluster"]).copy()
    if cluster_df.empty:
        return {}

    rows = []
    for (idx_name, cluster), sdf in cluster_df.groupby(["index_name", "cluster"]):
        for metric in metrics:
            vals = pd.to_numeric(sdf[metric], errors="coerce").dropna()
            if vals.empty:
                continue
            rows.append(
                {
                    "index_name": idx_name,
                    "cluster": int(cluster),
                    "metric": metric,
                    "n_stations": int(len(vals)),
                    "mean": float(vals.mean()),
                    "median": float(vals.median()),
                    "std": float(vals.std(ddof=1)) if len(vals) > 1 else 0.0,
                    "min": float(vals.min()),
                    "max": float(vals.max()),
                }
            )
    regional_df = pd.DataFrame(rows)
    if not regional_df.empty:
        regional_df.to_csv(tables_dir / "regional_cluster_composites.csv", index=False)

        apply_publication_theme()
        for idx_cfg in cfg["indices"]:
            idx_name = idx_cfg["name"]
            sdf = cluster_df.loc[cluster_df["index_name"] == idx_name].copy()
            if sdf.empty:
                continue

            heat_df = regional_df.loc[regional_df["index_name"] == idx_name].pivot(index="cluster", columns="metric", values="median").sort_index()
            if not heat_df.empty:
                _plot_single_heatmap(
                    heat_df,
                    f"Cluster Median Composite - {idx_cfg['title']}",
                    figs_dir / f"cluster_composite_heatmap_{idx_name}.{cfg['plots']['save_format']}",
                    cmap="RdBu_r",
                    fmt="{:.2f}",
                )

            fig, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
            groups = []
            labels = []
            for cluster_id, cdf in sdf.groupby("cluster"):
                vals = pd.to_numeric(cdf["Delta1"], errors="coerce").dropna()
                if vals.empty:
                    continue
                groups.append(vals.to_numpy(dtype=float))
                labels.append(f"C{int(cluster_id)}\n(n={len(vals)})")
            if groups:
                ax.boxplot(groups, labels=labels, patch_artist=True, boxprops={"facecolor": "#a8dadc"}, medianprops={"color": "#1d3557"})
                ax.axhline(0, color="black", linewidth=0.8)
                ax.set_ylabel("Delta1")
                ax.set_title(f"Cluster-wise Delta1 Distribution - {idx_cfg['title']}")
                fig.savefig(figs_dir / f"cluster_delta1_boxplot_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
            plt.close(fig)

    return {"regional_cluster_composites": regional_df}
