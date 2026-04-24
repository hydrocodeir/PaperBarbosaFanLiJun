from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.interpolate import Rbf, griddata
from scipy.stats import gaussian_kde
from shapely import contains_xy
from shapely.geometry import Point
from shapely.ops import polygonize, unary_union as shp_unary_union

from .clustering import screen_clustering_features
from .config_utils import (
    get_delta_definitions,
    get_default_figure_dpi,
    get_focus_quantiles,
    get_full_quantiles,
    get_plot_quantiles,
    get_primary_delta,
    get_sensitivity_quantiles,
    get_time_scale_years,
    get_time_unit_label,
    metric_label,
    slope_col,
    tau_label,
)
from .quantile import fit_ols_line, fit_quantile_line
from .quantile import bootstrap_qr, build_station_seed, summarize_bootstrap
from .year_config import build_split_periods, format_year_range_label, get_effective_year_range


def apply_publication_theme(cfg: dict | None = None) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    figure_dpi = 120 if cfg is None else get_default_figure_dpi(cfg)
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif", "Times"],
        "axes.titlesize": 13,
        "axes.titleweight": "semibold",
        "axes.labelsize": 11,
        "axes.linewidth": 0.9,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.fontsize": 9,
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "#d0d0d0",
        "figure.dpi": figure_dpi,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })


def _panel_label(i: int) -> str:
    return f"({chr(97 + i)})"


def _station_slug(station_name: str) -> str:
    safe = "".join(ch.lower() if ch.isalnum() else "_" for ch in str(station_name))
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe.strip("_")


def _format_slope(value: float) -> str:
    return "NA" if not np.isfinite(value) else f"{value:+.2f}"


def _format_metric_title(metric: str) -> str:
    return metric_label(metric).replace("_", " ")


def _resolve_quantile_axis_limits(quantiles: list[float]) -> tuple[float, float]:
    return (max(0.01, min(quantiles) - 0.01), min(0.99, max(quantiles) + 0.01))


def _quantile_axis_ticks(quantiles: list[float], step: float = 0.05) -> list[float]:
    if not quantiles:
        return []
    qmin = min(quantiles)
    qmax = max(quantiles)
    start = np.ceil(qmin / step) * step
    stop = np.floor(qmax / step) * step
    ticks = np.arange(start, stop + step * 0.5, step, dtype=float)
    return [round(float(t), 2) for t in ticks]


def _time_unit_ylabel(idx_cfg: dict, cfg: dict) -> str:
    return f"Trend ({idx_cfg['unit']} per {get_time_unit_label(cfg)})"


def _get_plot_year_range(cfg: dict, df: pd.DataFrame | None = None) -> tuple[int, int] | None:
    return get_effective_year_range(cfg, df)


FIG3_CMAP = LinearSegmentedColormap.from_list(
    "paper2_fig3",
    [
        "#10207a",
        "#1f4aa8",
        "#3b7d2a",
        "#8cbf26",
        "#f2f0d0",
        "#f0dd3d",
        "#ef9b1a",
        "#d6281f",
    ],
    N=256,
)


def _load_boundary_geometry(boundary_path: Path | None):
    if boundary_path is None:
        return None
    path = Path(boundary_path)
    if not path.exists():
        return None
    gdf = gpd.read_file(path)
    if gdf.empty:
        return None
    gdf = gdf.to_crs(4326)
    geom = gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all") else gdf.unary_union
    if geom.geom_type in {"LineString", "MultiLineString"}:
        polys = list(polygonize(geom))
        if polys:
            geom = shp_unary_union(polys)
    return geom


def _get_boundary_path_from_cfg(cfg: dict) -> Path:
    spatial_cfg = cfg.get("spatial_visualization", {})
    return Path(spatial_cfg.get("iran_boundary_geojson", "data/Iran.geojson"))


def _get_plot_extent(meta: pd.DataFrame, boundary_geom=None, pad_deg: float = 0.4):
    if boundary_geom is not None:
        minx, miny, maxx, maxy = boundary_geom.bounds
    else:
        minx = float(meta["longitude"].min())
        maxx = float(meta["longitude"].max())
        miny = float(meta["latitude"].min())
        maxy = float(meta["latitude"].max())
    return (minx - pad_deg, maxx + pad_deg, miny - pad_deg, maxy + pad_deg)


def _build_interpolation_grid(merged: pd.DataFrame, boundary_geom=None, nx: int = 320, ny: int = 320, pad_deg: float = 0.0):
    xmin, xmax, ymin, ymax = _get_plot_extent(merged, boundary_geom=boundary_geom, pad_deg=pad_deg)
    grid_x, grid_y = np.meshgrid(np.linspace(xmin, xmax, nx), np.linspace(ymin, ymax, ny))
    return grid_x, grid_y


def _interpolate_station_surface(merged: pd.DataFrame, grid_x, grid_y, method: str = "thin_plate_spline", smooth: float = 0.35):
    x = merged["longitude"].to_numpy(dtype=float)
    y = merged["latitude"].to_numpy(dtype=float)
    z = merged["slope_value"].to_numpy(dtype=float)
    if len(merged) < 4:
        return np.full_like(grid_x, np.nan, dtype=float)

    method = (method or "thin_plate_spline").lower()
    rbf_aliases = {
        "thin_plate_spline": "thin_plate",
        "rbf": "thin_plate",
        "spline": "thin_plate",
        "multiquadric": "multiquadric",
        "inverse": "inverse",
        "gaussian": "gaussian",
        "linear_rbf": "linear",
        "quintic": "quintic",
    }
    if method in rbf_aliases:
        rbf = Rbf(x, y, z, function=rbf_aliases[method], smooth=smooth)
        return rbf(grid_x, grid_y)
    if method in {"linear", "cubic", "nearest"}:
        surface = griddata(np.column_stack([x, y]), z, (grid_x, grid_y), method=method)
        if method != "nearest" and np.isnan(surface).any():
            nearest = griddata(np.column_stack([x, y]), z, (grid_x, grid_y), method="nearest")
            surface = np.where(np.isnan(surface), nearest, surface)
        return surface
    rbf = Rbf(x, y, z, function="thin_plate", smooth=smooth)
    return rbf(grid_x, grid_y)


def _mask_surface_to_boundary(grid_x, grid_y, surface, boundary_geom=None):
    if boundary_geom is None:
        return surface
    try:
        mask = contains_xy(boundary_geom, grid_x.ravel(), grid_y.ravel()).reshape(grid_x.shape)
    except Exception:
        pts = [Point(float(x), float(y)) for x, y in zip(grid_x.ravel(), grid_y.ravel())]
        mask = np.array([boundary_geom.covers(pt) for pt in pts], dtype=bool).reshape(grid_x.shape)
    return np.where(mask, surface, np.nan)


def _draw_boundary(ax, boundary_geom, linewidth: float = 1.0, zorder: float = 4):
    if boundary_geom is None:
        return
    geos = getattr(boundary_geom, "geoms", [boundary_geom])
    for geom in geos:
        try:
            x, y = geom.exterior.xy
            ax.plot(x, y, color="black", linewidth=linewidth, zorder=zorder)
            for ring in geom.interiors:
                xi, yi = ring.xy
                ax.plot(xi, yi, color="black", linewidth=max(0.35, linewidth * 0.4), zorder=zorder)
        except Exception:
            continue


def _format_geo_axis(ax, show_x: bool, show_y: bool) -> None:
    xticks = np.arange(44, 64.1, 2.5)
    yticks = np.arange(25, 40.1, 2.0)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    if show_x:
        ax.set_xticklabels([f"{x:.1f}°E" for x in xticks], fontsize=8)
    else:
        ax.set_xticklabels([])
    if show_y:
        ax.set_yticklabels([f"{y:.0f}°N" for y in yticks], fontsize=8)
    else:
        ax.set_yticklabels([])
    ax.tick_params(direction="out", length=3.5, width=0.8)


def _fig3_levels_from_values(values: pd.Series) -> np.ndarray:
    finite = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return np.arange(-10, 12, 2, dtype=float)
    vmax = float(np.nanmax(np.abs(finite)))
    vmax = max(vmax, 2.0)
    step = 1.0 if vmax <= 6 else 2.0
    bound = step * np.ceil(vmax / step)
    return np.arange(-bound, bound + step, step, dtype=float)


def _style_timeseries_axis(ax, years: np.ndarray, values: np.ndarray, idx_cfg: dict, panel_idx: int) -> None:
    ax.set_title(f"{_panel_label(panel_idx)} {idx_cfg['title']}", loc="left")
    ax.set_xlabel("Year")
    ax.set_ylabel(f"{idx_cfg['title']} ({idx_cfg['unit']})")
    if len(years):
        ax.set_xlim(years.min() - 1, years.max() + 1)
    if len(values):
        ymin = np.nanmin(values)
        ymax = np.nanmax(values)
        pad = max(3.0, 0.08 * max(1.0, ymax - ymin))
        ax.set_ylim(ymin - pad, ymax + pad)
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.grid(False, axis="x")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _style_coeff_axis(ax, idx_cfg: dict, panel_idx: int, cfg: dict, quantiles: list[float]) -> None:
    ax.set_title(f"{_panel_label(panel_idx)} {idx_cfg['title']}", loc="left")
    ax.set_xlabel("Quantile, τ")
    ax.set_ylabel(_time_unit_ylabel(idx_cfg, cfg))
    ax.set_xlim(*_resolve_quantile_axis_limits(quantiles))
    ax.set_xticks(_quantile_axis_ticks(quantiles))
    ax.axhline(0, color="#3a3a3a", linewidth=0.8, linestyle="--", alpha=0.85, zorder=1)
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.grid(False, axis="x")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _draw_quantile_coefficient_panel(
    ax,
    coeff_df: pd.DataFrame,
    ols_fit: dict,
    idx_cfg: dict,
    cfg: dict,
    quantiles: list[float],
    panel_idx: int | None = None,
) -> None:
    if coeff_df.empty:
        return

    shade_df = coeff_df.dropna(subset=["ci_low", "ci_high"])
    if not shade_df.empty:
        ax.fill_between(
            shade_df["tau"].to_numpy(),
            shade_df["ci_low"].to_numpy(),
            shade_df["ci_high"].to_numpy(),
            color="#bdbdbd",
            alpha=0.45,
            linewidth=0,
            zorder=1,
            label="QR 95% CI",
        )

    ax.plot(coeff_df["tau"], coeff_df["slope"], color="#1d1d1d", linewidth=0.9, alpha=0.7, zorder=2)
    ax.scatter(coeff_df["tau"], coeff_df["slope"], s=13, color="#111111", zorder=3, label="QR slope")

    if np.isfinite(float(ols_fit["slope"])):
        ax.axhline(float(ols_fit["slope"]), color="#c62828", linewidth=2.0, zorder=4, label="LSM mean trend")
    if np.isfinite(float(ols_fit["ci_low"])):
        ax.axhline(float(ols_fit["ci_low"]), color="#c62828", linewidth=1.2, linestyle="--", zorder=4, label="LSM 95% CI")
    if np.isfinite(float(ols_fit["ci_high"])):
        ax.axhline(float(ols_fit["ci_high"]), color="#c62828", linewidth=1.2, linestyle="--", zorder=4)

    if panel_idx is None:
        ax.set_xlabel("Quantile, τ")
        ax.set_ylabel(_time_unit_ylabel(idx_cfg, cfg))
        ax.set_xlim(*_resolve_quantile_axis_limits(quantiles))
        ax.set_xticks(_quantile_axis_ticks(quantiles))
        ax.axhline(0, color="#3a3a3a", linewidth=0.8, linestyle="--", alpha=0.85, zorder=1)
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
        ax.grid(False, axis="x")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    else:
        _style_coeff_axis(ax, idx_cfg, panel_idx, cfg, quantiles)


def _annotate_textbox(ax, text: str) -> None:
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.8,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#c9c9c9", "alpha": 0.95},
    )


def _annotate_textbox_bottom_left(ax, text: str) -> None:
    ax.text(
        0.02,
        0.03,
        text,
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=8.8,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#c9c9c9", "alpha": 0.95},
    )


def _build_quantile_profile_df(
    years: np.ndarray,
    values: np.ndarray,
    quantiles: list[float],
    cfg: dict,
    *,
    series_key: str,
) -> pd.DataFrame:
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    time_scale_years = get_time_scale_years(cfg)
    coeff_rows = []
    missing_ci = []

    for tau in quantiles:
        qr_fit = fit_quantile_line(years, values, tau=float(tau), time_scale_years=time_scale_years, max_iter=max_iter)
        row = {
            "tau": float(tau),
            "slope": float(qr_fit["slope"]),
            "ci_low": float(qr_fit["ci_low"]),
            "ci_high": float(qr_fit["ci_high"]),
        }
        coeff_rows.append(row)
        if np.isfinite(row["slope"]) and not (np.isfinite(row["ci_low"]) and np.isfinite(row["ci_high"])):
            missing_ci.append(float(tau))

    if missing_ci and cfg.get("bootstrap", {}).get("enabled", False):
        rng = np.random.default_rng(build_station_seed(int(cfg["project"]["random_seed"]), series_key, "profile_ci"))
        boot = bootstrap_qr(years, values, missing_ci, cfg, rng)
        if not boot.empty:
            boot_summary = summarize_bootstrap(boot, float(cfg["bootstrap"]["alpha"]))
            for row in coeff_rows:
                tau_suffix = f"{float(row['tau']):0.2f}"
                if float(row["tau"]) in missing_ci:
                    boot_low = boot_summary.get(f"boot_ci_low_{tau_suffix}", np.nan)
                    boot_high = boot_summary.get(f"boot_ci_high_{tau_suffix}", np.nan)
                    if np.isfinite(boot_low) and np.isfinite(boot_high):
                        row["ci_low"] = float(boot_low)
                        row["ci_high"] = float(boot_high)

    return pd.DataFrame(coeff_rows).dropna(subset=["tau", "slope"])


def plot_station_paper2_figures(annual: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    qgrid = get_full_quantiles(cfg)
    display_quantiles = get_plot_quantiles(cfg, "station_timeseries")
    time_scale_years = get_time_scale_years(cfg)
    fig_root = outdir / "paper2_station_figures"
    fig1_dir = fig_root / "figure1_timeseries"
    fig2_dir = fig_root / "figure2_quantile_coefficients"
    fig1_dir.mkdir(parents=True, exist_ok=True)
    fig2_dir.mkdir(parents=True, exist_ok=True)

    palette = ["#c62828", "#1f1f1f", "#2e7d32", "#1565c0", "#6a1b9a"]
    qr_colors = {tau: palette[i % len(palette)] for i, tau in enumerate(display_quantiles)}

    for (station_id, station_name), sdf in annual.groupby(["station_id", "station_name"]):
        sdf = sdf.sort_values("year").copy()
        if sdf.empty:
            continue

        years = sdf["year"].to_numpy(dtype=float)
        fig1, axes1 = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
        fig2, axes2 = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
        fig1.suptitle(f"Figure 1. Station-Level Quantile and Mean Trends - {station_name} ({int(station_id)})", fontsize=15, fontweight="bold", y=1.02)
        fig2.suptitle(f"Figure 2. Quantile Regression Coefficients - {station_name} ({int(station_id)})", fontsize=15, fontweight="bold", y=1.02)

        for i, idx_cfg in enumerate(cfg["indices"]):
            idx_name = idx_cfg["name"]
            values = sdf[idx_name].to_numpy(dtype=float)
            ax1 = axes1.flat[i]
            ax2 = axes2.flat[i]

            ax1.plot(years, values, color="#c7c7c7", linewidth=1.0, alpha=0.95, zorder=1)
            ax1.scatter(years, values, s=18, color="#d9d9d9", edgecolor="white", linewidth=0.45, alpha=0.95, zorder=2)

            annotation_lines = []
            for tau in display_quantiles:
                qr_fit = fit_quantile_line(years, values, tau=tau, time_scale_years=time_scale_years, max_iter=max_iter)
                if len(qr_fit["fitted"]):
                    ax1.plot(
                        qr_fit["years"],
                        qr_fit["fitted"],
                        color=qr_colors[tau],
                        linewidth=2.4 if tau == display_quantiles[len(display_quantiles) // 2] else 2.1,
                        zorder=3,
                        label=f"QR τ={tau:.2f}",
                    )
                annotation_lines.append(f"τ={tau:.2f}: {_format_slope(float(qr_fit['slope']))}")

            ols_fit = fit_ols_line(years, values, time_scale_years=time_scale_years)
            if len(ols_fit["fitted"]):
                ax1.plot(
                    ols_fit["years"],
                    ols_fit["fitted"],
                    color="#111111",
                    linewidth=2.0,
                    linestyle=(0, (1.2, 1.8)),
                    zorder=4,
                    label="LSM mean",
                )
            annotation_lines.append(f"Mean: {_format_slope(float(ols_fit['slope']))}")
            _style_timeseries_axis(ax1, years, values, idx_cfg, i)
            _annotate_textbox(ax1, f"Slope (days/{get_time_unit_label(cfg)})\n" + "\n".join(annotation_lines))
            if i == 0:
                ax1.legend(loc="upper right", ncol=2)

            coeff_df = _build_quantile_profile_df(
                years,
                values,
                qgrid,
                cfg,
                series_key=f"station:{station_id}:{station_name}:{idx_name}",
            )

            if not coeff_df.empty:
                _draw_quantile_coefficient_panel(ax2, coeff_df, ols_fit, idx_cfg, cfg, qgrid, panel_idx=i)
            _annotate_textbox(
                ax2,
                f"LSM slope: {_format_slope(float(ols_fit['slope']))}\n95% CI: [{_format_slope(float(ols_fit['ci_low']))}, {_format_slope(float(ols_fit['ci_high']))}]",
            )
            if i == 0:
                ax2.legend(loc="upper right")

        slug = _station_slug(station_name)
        fig1.savefig(fig1_dir / f"{int(station_id)}_{slug}_figure1.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        fig2.savefig(fig2_dir / f"{int(station_id)}_{slug}_figure2.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig1)
        plt.close(fig2)


def plot_ijoc_station_comparisons(annual: pd.DataFrame, feature_table: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    fig_dir = outdir / "ijoc_station_comparisons"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = outdir.parent / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    qgrid = get_full_quantiles(cfg)
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    focus_quantiles = get_plot_quantiles(cfg, "station_comparison")
    time_scale_years = get_time_scale_years(cfg)
    primary_delta = get_primary_delta(cfg)
    cluster_colors = {1: "#1b5e20", 2: "#1565c0", 3: "#ef6c00", 4: "#8e24aa"}
    station_short = {
        "Tehran (Mehrabad Airport)": "Tehran",
        "Bushehr (Airport)": "Bushehr",
        "Torbat-E Heydariyeh": "Torbat-e H.",
    }

    def _configured_feature_cols() -> list[str]:
        mode = str(cfg.get("clustering", {}).get("feature_mode", "simple"))
        if mode == "uncertainty":
            return list(cfg["clustering"].get("uncertainty_features", []))
        return list(cfg["clustering"].get("simple_features", []))

    screening_df = screen_clustering_features(
        feature_table,
        cfg,
        feature_cols=_configured_feature_cols(),
        feature_set_label="representative_station_selection",
    )

    def _kept_feature_cols(idx_name: str) -> list[str]:
        if screening_df.empty:
            return _configured_feature_cols()
        kept = screening_df.loc[
            (screening_df["index_name"] == idx_name)
            & (screening_df["feature_set"] == "representative_station_selection")
            & (screening_df["status"].astype(str).str.startswith("kept")),
            "feature",
        ].tolist()
        return kept or _configured_feature_cols()

    def _select_representative_rows(idx_name: str) -> list[dict]:
        sdf = feature_table.loc[feature_table["index_name"] == idx_name].copy()
        if sdf.empty or "cluster" not in sdf.columns:
            return []
        sdf = sdf.dropna(subset=["cluster"]).copy()
        if sdf.empty:
            return []

        feature_cols = [col for col in _kept_feature_cols(idx_name) if col in sdf.columns]
        if not feature_cols:
            feature_cols = [slope_col(tau) for tau in focus_quantiles if slope_col(tau) in sdf.columns]
        if not feature_cols:
            return []

        X = sdf[feature_cols].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
        X = X.apply(lambda col: col.fillna(col.median()), axis=0)
        means = X.mean(axis=0)
        stds = X.std(axis=0, ddof=0).replace(0, 1.0)
        Xz = (X - means) / stds
        cluster_series = pd.to_numeric(sdf["cluster"], errors="coerce").astype("Int64")

        rows = []
        unique_clusters = sorted(cluster_series.dropna().astype(int).unique().tolist())
        for cluster_id in unique_clusters:
            cluster_mask = cluster_series == cluster_id
            cluster_sdf = sdf.loc[cluster_mask].copy()
            cluster_Xz = Xz.loc[cluster_mask].copy()
            if cluster_sdf.empty or cluster_Xz.empty:
                continue

            centroid = cluster_Xz.mean(axis=0).to_numpy(dtype=float)
            distances = np.sqrt(((cluster_Xz.to_numpy(dtype=float) - centroid) ** 2).sum(axis=1))
            candidate = (
                cluster_sdf.assign(distance_to_cluster_centroid=distances)
                .sort_values(["distance_to_cluster_centroid", "station_name"], ascending=[True, True])
                .iloc[0]
            )
            rows.append(
                {
                    "index_name": idx_name,
                    "station_id": candidate["station_id"],
                    "station_name": candidate["station_name"],
                    "cluster": int(cluster_id),
                    "selection_method": "closest_to_cluster_centroid",
                    "feature_space": ",".join(feature_cols),
                    "distance_to_cluster_centroid": float(candidate["distance_to_cluster_centroid"]),
                    **{slope_col(tau): float(candidate.get(slope_col(tau), np.nan)) for tau in focus_quantiles},
                    primary_delta: float(candidate.get(primary_delta, np.nan)),
                }
            )
        return rows

    representative_rows = []
    for idx_cfg in cfg["indices"]:
        representative_rows.extend(_select_representative_rows(idx_cfg["name"]))
    representative_df = pd.DataFrame(representative_rows)
    if not representative_df.empty:
        representative_df = representative_df.sort_values(["index_name", "cluster", "distance_to_cluster_centroid", "station_name"]).reset_index(drop=True)
        representative_df.to_csv(tables_dir / "representative_station_selection.csv", index=False)

    def build_curve_rows(idx_name: str):
        rows = []
        local_reps = representative_df.loc[representative_df["index_name"] == idx_name].copy()
        if local_reps.empty:
            return rows
        local_reps = local_reps.sort_values(["cluster", "distance_to_cluster_centroid", "station_name"])
        for _, rep in local_reps.iterrows():
            station_name = str(rep["station_name"])
            sdf = annual.loc[annual["station_name"] == station_name].sort_values("year").copy()
            if sdf.empty:
                continue
            values = sdf[idx_name].to_numpy(dtype=float)
            years = sdf["year"].to_numpy(dtype=float)
            coeff_rows = []
            for tau in qgrid:
                qr_fit = fit_quantile_line(years, values, tau=float(tau), time_scale_years=time_scale_years, max_iter=max_iter)
                coeff_rows.append({"tau": float(tau), "slope": float(qr_fit["slope"])})
            coeff_df = pd.DataFrame(coeff_rows).dropna(subset=["tau", "slope"])
            if coeff_df.empty:
                continue
            meta = feature_table.loc[
                (feature_table["index_name"] == idx_name) & (feature_table["station_name"] == station_name)
            ].iloc[0]
            cluster = int(meta["cluster"]) if pd.notna(meta["cluster"]) else -1
            rows.append({
                "station_name": station_name,
                "label": station_short.get(station_name, station_name),
                "cluster": cluster,
                "color": cluster_colors.get(cluster, "#333333"),
                "coeff_df": coeff_df,
                "meta": meta,
                "distance_to_cluster_centroid": float(rep["distance_to_cluster_centroid"]),
            })
        return rows

    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        curve_rows = build_curve_rows(idx_name)
        if not curve_rows:
            continue

        fig, ax = plt.subplots(figsize=(9.5, 6.2))
        legend_handles = []
        legend_labels = []

        for row in curve_rows:
            coeff_df = row["coeff_df"]
            meta = row["meta"]
            cluster = row["cluster"]
            color = row["color"]
            line, = ax.plot(
                coeff_df["tau"],
                coeff_df["slope"],
                color=color,
                linewidth=2.3,
                alpha=0.92,
            )
            for fq in focus_quantiles:
                point = coeff_df.loc[np.isclose(coeff_df["tau"], fq)]
                if not point.empty:
                    ax.scatter(
                        point["tau"],
                        point["slope"],
                        s=34,
                        color=color,
                        edgecolor="white",
                        linewidth=0.4,
                        zorder=4,
                    )

            legend_handles.append(line)
            legend_labels.append(
                f"{row['label']} (C{cluster}; " + ", ".join(
                    f"{tau_label(tau)}={meta[slope_col(tau)]:+.2f}" for tau in focus_quantiles if slope_col(tau) in meta.index
                ) + ")"
            )

        ax.axhline(0, color="#444444", linewidth=0.8, linestyle="--", alpha=0.9)
        ax.set_xlim(0.01, 0.99)
        ax.set_xticks(_quantile_axis_ticks(qgrid))
        ax.set_xlabel("Quantile, τ")
        ax.set_ylabel(_time_unit_ylabel(idx_cfg, cfg))
        ax.set_title(f"Representative station comparison - {idx_cfg['title']}")
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
        ax.grid(False, axis="x")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(legend_handles, legend_labels, loc="upper left", fontsize=8.5, framealpha=0.95)
        _annotate_textbox_bottom_left(
            ax,
            "Representative stations were chosen algorithmically as the closest\n"
            "observed member to each cluster centroid in the screened feature space.",
        )
        fig.tight_layout()
        fig.savefig(fig_dir / f"station_comparison_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)

    # Main multi-panel figure for the manuscript
    idx_order = [idx["name"] for idx in cfg["indices"]]
    idx_cfg_map = {idx["name"]: idx for idx in cfg["indices"]}
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 10.4), sharex=True, constrained_layout=True)
    fig.suptitle(
        "Representative station comparisons across cluster-defined quantile structures",
        fontsize=16,
        fontweight="semibold",
        y=1.01,
    )
    focus_quantiles = get_plot_quantiles(cfg, "station_comparison")

    for panel_idx, idx_name in enumerate(idx_order):
        ax = axes.flat[panel_idx]
        idx_cfg = idx_cfg_map[idx_name]
        curve_rows = build_curve_rows(idx_name)
        if not curve_rows:
            continue

        handles = []
        labels = []
        for row in curve_rows:
            coeff_df = row["coeff_df"]
            meta = row["meta"]
            color = row["color"]
            label = row["label"]
            line, = ax.plot(
                coeff_df["tau"],
                coeff_df["slope"],
                color=color,
                linewidth=2.15,
                alpha=0.95,
                solid_capstyle="round",
                zorder=2,
            )
            handles.append(line)
            labels.append(f"{label} (C{row['cluster']})")

            markers = ["o", "s", "D", "^", "v", "P"]
            for fq, marker in zip(focus_quantiles, markers):
                point = coeff_df.loc[np.isclose(coeff_df["tau"], fq)]
                if not point.empty:
                    ax.scatter(
                        point["tau"],
                        point["slope"],
                        s=44 if fq == focus_quantiles[len(focus_quantiles) // 2] else 38,
                        color=color,
                        marker=marker,
                        edgecolor="white",
                        linewidth=0.55,
                        zorder=4,
                    )

            end_point = coeff_df.iloc[-1]
            ax.text(
                min(0.985, float(end_point["tau"]) + 0.006),
                float(end_point["slope"]),
                label,
                color=color,
                fontsize=8.6,
                va="center",
                ha="left",
                clip_on=True,
            )

        ax.set_title(f"{_panel_label(panel_idx)} {idx_cfg['title']}", loc="left")
        ax.axhline(0, color="#5b5b5b", linewidth=0.85, linestyle=(0, (3, 2)), alpha=0.95, zorder=1)
        ax.set_xlim(0.01, 0.99)
        ax.set_xticks(_quantile_axis_ticks(qgrid))
        ax.grid(True, axis="y", alpha=0.18, linewidth=0.6)
        ax.grid(False, axis="x")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#888888")
        ax.spines["bottom"].set_color("#888888")
        if panel_idx % 2 == 0:
            ax.set_ylabel(_time_unit_ylabel(idx_cfg, cfg))
        if panel_idx >= 2:
            ax.set_xlabel("Quantile, τ")

        summary_lines = [
            f"{row['label']}: " + " / ".join(
                f"{row['meta'][slope_col(tau)]:+.1f}" for tau in focus_quantiles if slope_col(tau) in row["meta"].index
            )
            for row in curve_rows
        ]
        ax.text(
            0.02,
            0.03,
            " / ".join(tau_label(tau) for tau in focus_quantiles) + "\n" + "\n".join(summary_lines),
            transform=ax.transAxes,
            va="bottom",
            ha="left",
            fontsize=7.7,
            bbox={"boxstyle": "round,pad=0.32", "facecolor": "white", "edgecolor": "#d0d0d0", "alpha": 0.94},
        )

    fig.text(
        0.5,
        -0.01,
        f"Markers denote configured focal quantiles: {', '.join(tau_label(tau) for tau in focus_quantiles)}. "
        "Stations are the closest observed member to each cluster centroid in the screened clustering feature space.",
        ha="center",
        va="top",
        fontsize=9,
    )
    fig.savefig(fig_dir / f"main_figure_representative_stations.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def plot_paper2_figure3_maps(annual: pd.DataFrame, stations: pd.DataFrame, outdir: Path, cfg: dict, taus: List[float] | None = None):
    apply_publication_theme()
    fig_dir = outdir / "paper2_figure3_quantile_maps"
    fig_dir.mkdir(parents=True, exist_ok=True)
    taus = taus or get_plot_quantiles(cfg, "paper2_maps")
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    time_scale_years = get_time_scale_years(cfg)
    spatial_cfg = cfg.get("spatial_visualization", {})
    boundary_path = _get_boundary_path_from_cfg(cfg)
    interpolation_method = spatial_cfg.get("interpolation_method", "thin_plate_spline")
    interpolation_smooth = float(spatial_cfg.get("interpolation_smooth", 0.35))
    boundary_geom = _load_boundary_geometry(boundary_path)
    for tau in taus:
        rows = []
        slope_col = f"slope_{tau:0.2f}"

        for idx_cfg in cfg["indices"]:
            idx_name = idx_cfg["name"]
            sdf = annual[["station_id", "station_name", "year", idx_name]].copy().sort_values(["station_id", "year"])
            for (station_id, station_name), gdf in sdf.groupby(["station_id", "station_name"]):
                years = gdf["year"].to_numpy(dtype=float)
                values = gdf[idx_name].to_numpy(dtype=float)
                fit = fit_quantile_line(years, values, tau=float(tau), time_scale_years=time_scale_years, max_iter=max_iter)
                rows.append(
                    {
                        "index_name": idx_name,
                        "station_id": station_id,
                        "station_name": station_name,
                        slope_col: float(fit["slope"]),
                    }
                )

        map_df = pd.DataFrame(rows).merge(stations, on=["station_id", "station_name"], how="left")
        levels = _fig3_levels_from_values(map_df[slope_col])
        norm = BoundaryNorm(levels, FIG3_CMAP.N, clip=False)

        fig, axes = plt.subplots(2, 2, figsize=(15, 11), constrained_layout=True)
        fig.suptitle(
            f"Figure 3. Quantile-Regression Slope Maps at τ = {_short_tau_label(float(tau))}",
            fontsize=15,
            fontweight="bold",
            y=1.02,
        )

        mappable = None
        for i, idx_cfg in enumerate(cfg["indices"]):
            idx_name = idx_cfg["name"]
            ax = axes.flat[i]
            sdf = map_df.loc[map_df["index_name"] == idx_name].copy()
            ax.set_title(f"{_panel_label(i)} {idx_cfg['title']}", loc="left")
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            ax.grid(False)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_facecolor("#f4f4f4")

            if sdf.empty:
                ax.text(0.5, 0.5, "No station data", transform=ax.transAxes, ha="center", va="center")
                continue

            interp_df = sdf[["longitude", "latitude", slope_col]].rename(columns={slope_col: "slope_value"}).dropna().copy()
            grid_x, grid_y = _build_interpolation_grid(interp_df, boundary_geom=boundary_geom, nx=320, ny=320)
            surface = _interpolate_station_surface(
                interp_df,
                grid_x,
                grid_y,
                method=interpolation_method,
                smooth=interpolation_smooth,
            )
            surface = _mask_surface_to_boundary(grid_x, grid_y, surface, boundary_geom=boundary_geom)
            mappable = ax.contourf(
                grid_x,
                grid_y,
                surface,
                levels=levels,
                cmap=FIG3_CMAP,
                norm=norm,
                extend="both",
                antialiased=True,
                zorder=1,
            )
            _draw_boundary(ax, boundary_geom, linewidth=1.0)

            ax.scatter(
                sdf["longitude"],
                sdf["latitude"],
                c=sdf[slope_col],
                s=38,
                cmap=FIG3_CMAP,
                norm=norm,
                edgecolor="white",
                linewidth=0.45,
                alpha=0.96,
                zorder=5,
            )

            xmin, xmax, ymin, ymax = _get_plot_extent(sdf, boundary_geom=boundary_geom, pad_deg=0.1)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_aspect("equal", adjustable="box")
            row_idx, col_idx = divmod(i, 2)
            _format_geo_axis(ax, show_x=(row_idx == 1), show_y=(col_idx == 0))

            _annotate_textbox_bottom_left(
                ax,
                f"Range: {_format_slope(float(pd.to_numeric(sdf[slope_col], errors='coerce').min()))} to {_format_slope(float(pd.to_numeric(sdf[slope_col], errors='coerce').max()))}\n"
                f"Station count: {len(sdf)} | Surface: visual interpolation only",
            )

        if mappable is not None:
            cbar = fig.colorbar(mappable, ax=axes.ravel().tolist(), shrink=0.88, pad=0.02, ticks=levels)
            cbar.set_label(f"Slope (days per {get_time_unit_label(cfg)})")

        fig.text(
            0.5,
            -0.01,
            "Colored station points are the primary evidence. Background shading shows an interpolated slope surface for visual orientation only and is not used as a standalone inferential layer.",
            ha="center",
            fontsize=10,
        )
        fig.savefig(fig_dir / f"figure3_tau_{tau:0.2f}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def _bootstrap_panel_facecolor(tau: float) -> str:
    palette = {
        0.05: "#fbe9e7",
        0.50: "#f4f4f4",
        0.95: "#e8f5e9",
    }
    return palette.get(float(tau), "#fafafa")


def _plot_bootstrap_distribution(ax, samples: np.ndarray, punctual_estimate: float, tau: float, idx_cfg: dict, panel_idx: int, show_legend: bool) -> None:
    samples = np.asarray(samples, dtype=float)
    samples = samples[np.isfinite(samples)]
    ax.set_facecolor(_bootstrap_panel_facecolor(tau))

    if len(samples) == 0:
        ax.text(0.5, 0.5, "No bootstrap samples", transform=ax.transAxes, ha="center", va="center", fontsize=10)
        ax.set_title(f"{_panel_label(panel_idx)} {idx_cfg['title']} | τ = {tau:.2f}", loc="left")
        return

    n_bins = min(24, max(12, int(np.sqrt(len(samples)) * 1.8)))
    hist_color = "#90caf9"
    edge_color = "#356a95"
    ax.hist(samples, bins=n_bins, density=True, color=hist_color, edgecolor=edge_color, linewidth=0.7, alpha=0.72, zorder=1, label="Bootstrap density" if show_legend else None)

    unique_count = len(np.unique(np.round(samples, 6)))
    if unique_count > 1:
        xs = np.linspace(samples.min(), samples.max(), 300)
        kde = gaussian_kde(samples)
        ax.plot(xs, kde(xs), color="#0d47a1", linewidth=2.0, zorder=2, label="KDE" if show_legend else None)

    mean_val = float(np.mean(samples))
    median_val = float(np.median(samples))
    ci_low, ci_high = np.quantile(samples, [0.025, 0.975])

    ax.axvline(punctual_estimate, color="#1a1a1a", linewidth=2.0, linestyle=(0, (4, 3)), zorder=3, label="Punctual estimate" if show_legend else None)
    ax.axvspan(ci_low, ci_high, color="#b0bec5", alpha=0.28, zorder=0, label="Bootstrap 95% interval" if show_legend else None)
    ax.axvline(mean_val, color="#c62828", linewidth=1.4, linestyle="-", alpha=0.9, zorder=2, label="Bootstrap mean" if show_legend else None)

    ax.set_title(f"{_panel_label(panel_idx)} {idx_cfg['title']} | τ = {tau:.2f}", loc="left")
    ax.set_xlabel(f"Slope ({idx_cfg['unit']} per decade)")
    ax.set_ylabel("Density")
    ax.grid(True, axis="y", alpha=0.22, linewidth=0.6)
    ax.grid(False, axis="x")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    textbox = (
        f"Estimate: {_format_slope(punctual_estimate)}\n"
        f"Mean: {_format_slope(mean_val)}\n"
        f"Median: {_format_slope(median_val)}\n"
        f"95% CI: [{_format_slope(float(ci_low))}, {_format_slope(float(ci_high))}]"
    )
    _annotate_textbox(ax, textbox)


def plot_station_paper1_figure4(boot_long: pd.DataFrame, qr_summary: pd.DataFrame, outdir: Path, cfg: dict):
    if boot_long.empty or qr_summary.empty:
        return

    apply_publication_theme()
    fig_dir = outdir / "paper1_station_figures" / "figure4_bootstrap_distributions"
    fig_dir.mkdir(parents=True, exist_ok=True)
    taus = get_plot_quantiles(cfg, "bootstrap_distributions")

    for (station_id, station_name), station_boot in boot_long.groupby(["station_id", "station_name"]):
        station_summary = qr_summary.loc[
            (qr_summary["station_id"] == station_id) & (qr_summary["station_name"] == station_name)
        ].copy()
        if station_summary.empty:
            continue

        fig, axes = plt.subplots(len(taus), len(cfg["indices"]), figsize=(18, 12.5), constrained_layout=True)
        fig.suptitle(
            f"Figure 4. Bootstrap Distributions of Quantile Slopes - {station_name} ({int(station_id)})",
            fontsize=15,
            fontweight="bold",
            y=1.01,
        )

        for row_idx, tau in enumerate(taus):
            slope_col = f"slope_{tau:0.2f}"
            for col_idx, idx_cfg in enumerate(cfg["indices"]):
                ax = axes[row_idx, col_idx]
                idx_name = idx_cfg["name"]
                idx_boot = station_boot.loc[station_boot["index_name"] == idx_name, slope_col].to_numpy(dtype=float)
                summary_row = station_summary.loc[station_summary["index_name"] == idx_name]
                estimate = float(summary_row.iloc[0][slope_col]) if (not summary_row.empty and slope_col in summary_row.columns) else np.nan
                panel_idx = row_idx * len(cfg["indices"]) + col_idx
                _plot_bootstrap_distribution(ax, idx_boot, estimate, tau, idx_cfg, panel_idx, show_legend=(row_idx == 0 and col_idx == 0))

        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper center", ncol=min(4, len(handles)), bbox_to_anchor=(0.5, 1.02))

        fig.text(
            0.5,
            -0.01,
            "Dashed vertical lines indicate the punctual quantile-regression estimates. Shaded bands show the bootstrap 95% interval.",
            ha="center",
            fontsize=10,
        )
        slug = _station_slug(station_name)
        fig.savefig(fig_dir / f"{int(station_id)}_{slug}_figure4.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def _short_tau_label(tau: float) -> str:
    return f"{tau:.2f}".rstrip("0").rstrip(".")


def _plot_single_quantile_dendrogram(ax, sdf: pd.DataFrame, idx_cfg: dict, tau: float, panel_idx: int) -> None:
    slope_col = f"slope_{tau:0.2f}"
    plot_df = sdf[["station_name", slope_col]].copy()
    plot_df[slope_col] = pd.to_numeric(plot_df[slope_col], errors="coerce")
    plot_df = plot_df.dropna().sort_values(slope_col, ascending=True)

    ax.set_title(f"{_panel_label(panel_idx)} {idx_cfg['title']}", loc="left")
    if len(plot_df) < 2:
        ax.text(0.5, 0.5, "Insufficient stations for clustering", transform=ax.transAxes, ha="center", va="center", fontsize=10)
        ax.set_axis_off()
        return

    X = plot_df[[slope_col]].to_numpy(dtype=float)
    Z = linkage(X, method="average", metric="euclidean")
    dendrogram(
        Z,
        labels=plot_df["station_name"].tolist(),
        orientation="right",
        leaf_font_size=8,
        color_threshold=None,
        above_threshold_color="#375a7f",
        ax=ax,
    )
    ax.set_xlabel(f"Slope distance at τ = {_short_tau_label(tau)}")
    ax.set_ylabel("Station")
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)

    stats_text = (
        f"n = {len(plot_df)}\n"
        f"Range: {_format_slope(float(plot_df[slope_col].min()))} to {_format_slope(float(plot_df[slope_col].max()))}\n"
        f"Median: {_format_slope(float(plot_df[slope_col].median()))}"
    )
    ax.text(
        0.98,
        0.02,
        stats_text,
        transform=ax.transAxes,
        va="bottom",
        ha="right",
        fontsize=8.3,
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "#d0d0d0", "alpha": 0.94},
    )


def plot_paper1_quantile_dendrograms(qr_summary: pd.DataFrame, outdir: Path, cfg: dict):
    if qr_summary.empty:
        return

    apply_publication_theme()
    fig_dir = outdir / "paper1_quantile_dendrograms"
    fig_dir.mkdir(parents=True, exist_ok=True)
    figure_map = list(enumerate(get_plot_quantiles(cfg, "paper1_dendrograms"), start=5))

    for fig_no, tau in figure_map:
        fig, axes = plt.subplots(2, 2, figsize=(16, 12), constrained_layout=True)
        fig.suptitle(
            f"Figure {fig_no}. Average-Linkage Dendrograms of Station Slopes at Quantile τ = {_short_tau_label(tau)}",
            fontsize=15,
            fontweight="bold",
            y=1.01,
        )

        for i, idx_cfg in enumerate(cfg["indices"]):
            idx_name = idx_cfg["name"]
            ax = axes.flat[i]
            sdf = qr_summary.loc[qr_summary["index_name"] == idx_name].copy()
            _plot_single_quantile_dendrogram(ax, sdf, idx_cfg, tau, i)

        fig.text(
            0.5,
            -0.01,
            f"Average linkage is applied to station slopes estimated at τ = {_short_tau_label(tau)}. Branch lengths represent Euclidean dissimilarity in slope magnitude.",
            ha="center",
            fontsize=10,
        )
        fig.savefig(fig_dir / f"figure_{fig_no}_tau_{tau:0.2f}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def _sort_station_names_for_heatmap(plot_df: pd.DataFrame, mode: str) -> List[str]:
    if mode == "alphabetic":
        return sorted(plot_df["station_name"].unique().tolist())
    delta_candidates = [col for col in plot_df.columns if col != "station_name"]
    primary_delta = plot_df.attrs.get("primary_delta", delta_candidates[-1] if delta_candidates else "")
    order_df = plot_df.groupby("station_name", as_index=False)[primary_delta].mean().sort_values(primary_delta, ascending=False)
    return order_df["station_name"].tolist()


def plot_station_heatmap(features: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    focus_quantiles = get_focus_quantiles(cfg)
    primary_delta = get_primary_delta(cfg)
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        cols = ["station_name"] + [slope_col(tau) for tau in focus_quantiles] + [primary_delta]
        sdf = features.loc[features["index_name"] == idx_name, cols].copy()
        if sdf.empty:
            continue
        sdf.attrs["primary_delta"] = primary_delta
        order = _sort_station_names_for_heatmap(sdf, cfg["plots"]["heatmap_station_order"])
        sdf["station_name"] = pd.Categorical(sdf["station_name"], categories=order, ordered=True)
        mat_cols = [slope_col(tau) for tau in focus_quantiles] + [primary_delta]
        mat = sdf.sort_values("station_name").set_index("station_name")[mat_cols]

        fig, ax = plt.subplots(figsize=(8, max(6, len(mat) * 0.3)))
        im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="RdBu_r")
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels([metric_label(col) for col in mat.columns])
        ax.set_yticks(range(mat.shape[0]))
        ax.set_yticklabels(mat.index.tolist())
        ax.set_title(f"Station quantile trend signatures - {idx_cfg['title']}")
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label(_time_unit_ylabel(idx_cfg, cfg))
        fig.tight_layout()
        fig.savefig(outdir / f"station_focus_heatmap_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def plot_region_quantile_slopes(annual: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    years = np.sort(annual["year"].unique())
    mean_df = annual.groupby("year", as_index=False)[[idx["name"] for idx in cfg["indices"]]].mean()
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    qgrid = get_full_quantiles(cfg)
    time_scale_years = get_time_scale_years(cfg)

    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        vals = mean_df[idx_name].to_numpy(dtype=float)

        coeff_df = _build_quantile_profile_df(
            years,
            vals,
            qgrid,
            cfg,
            series_key=f"regional_mean:{idx_name}",
        )
        ols_fit = fit_ols_line(years, vals, time_scale_years=time_scale_years)

        fig, ax = plt.subplots(figsize=(8, 5))
        _draw_quantile_coefficient_panel(ax, coeff_df, ols_fit, idx_cfg, cfg, qgrid)
        ax.set_title(f"Descriptive regional-average quantile profile - {idx_cfg['title']}")
        _annotate_textbox_bottom_left(
            ax,
            "Annual values were averaged across stations before fitting QR.\n"
            "Use as a descriptive pooled summary; station-level QR remains the main evidence.\n"
            f"OLS: {_format_slope(float(ols_fit['slope']))} | 95% CI: [{_format_slope(float(ols_fit['ci_low']))}, {_format_slope(float(ols_fit['ci_high']))}]",
        )
        ax.legend(loc="upper right")
        fig.tight_layout()
        fig.savefig(outdir / f"region_quantile_slopes_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def plot_ijoc_regional_quantile_panels(annual: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    years = np.sort(annual["year"].unique())
    mean_df = annual.groupby("year", as_index=False)[[idx["name"] for idx in cfg["indices"]]].mean()
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    qgrid = get_full_quantiles(cfg)
    time_scale_years = get_time_scale_years(cfg)
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 10.2), constrained_layout=True)

    for i, idx_cfg in enumerate(cfg["indices"]):
        idx_name = idx_cfg["name"]
        vals = mean_df[idx_name].to_numpy(dtype=float)
        coeff_df = _build_quantile_profile_df(
            years,
            vals,
            qgrid,
            cfg,
            series_key=f"ijoc_regional_mean:{idx_name}",
        )
        ols_fit = fit_ols_line(years, vals, time_scale_years=time_scale_years)
        ax = axes.flat[i]
        _draw_quantile_coefficient_panel(ax, coeff_df, ols_fit, idx_cfg, cfg, qgrid, panel_idx=i)
        ax.text(
            0.02,
            0.03,
            "Descriptive regional-average panel\n"
            "Station-level QR and cross-station summaries are primary.\n"
            f"OLS: {_format_slope(float(ols_fit['slope']))}",
            transform=ax.transAxes,
            fontsize=8.3,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#d0d0d0", "alpha": 0.94},
        )
        if i == 0:
            ax.legend(loc="upper left")
    fig.text(
        0.5,
        -0.01,
        "Panels summarize quantile-regression fits applied to annual station-mean series and are intended as descriptive pooled views only. "
        "Primary inference in the manuscript is based on station-level quantile regression, bootstrap summaries, and cross-station diagnostics.",
        ha="center",
        va="top",
        fontsize=9,
    )
    fig.savefig(outdir / f"ijoc_regional_quantile_panels.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def plot_dendrograms(features: pd.DataFrame, artifacts: Dict[str, np.ndarray], outdir: Path, cfg: dict):
    if not artifacts or cfg["clustering"]["algorithm"] != "hierarchical":
        return
    apply_publication_theme()
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        if idx_name not in artifacts:
            continue
        labels = features.loc[features["index_name"] == idx_name, "station_name"].tolist()
        fig, ax = plt.subplots(figsize=(10, 6))
        dendrogram(artifacts[idx_name], labels=labels, orientation="right", leaf_font_size=8, ax=ax, color_threshold=None)
        ax.set_title(f"Hierarchical clustering dendrogram - {idx_cfg['title']}")
        ax.set_xlabel("Distance")
        fig.tight_layout()
        fig.savefig(outdir / f"dendrogram_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def _scatter_station_map(df: pd.DataFrame, value_col: str, title: str, outpath: Path, cfg: dict, discrete: bool = False):
    boundary_geom = _load_boundary_geometry(_get_boundary_path_from_cfg(cfg))
    fig, ax = plt.subplots(figsize=(7.5, 6))
    _draw_boundary(ax, boundary_geom, linewidth=1.0, zorder=2)
    scatter = ax.scatter(
        df["longitude"],
        df["latitude"],
        c=df[value_col],
        s=75,
        cmap=("tab10" if discrete else "RdBu_r"),
        edgecolor="black",
        linewidth=0.35,
        zorder=3,
    )
    xmin, xmax, ymin, ymax = _get_plot_extent(df, boundary_geom=boundary_geom, pad_deg=0.1)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    _format_geo_axis(ax, show_x=True, show_y=True)
    ax.set_facecolor("#f4f4f4")
    ax.grid(False)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)
    fig.colorbar(scatter, ax=ax, shrink=0.82, label=value_col)
    fig.tight_layout()
    fig.savefig(outpath, dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def plot_ijoc_study_area(stations: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    boundary_geom = _load_boundary_geometry(_get_boundary_path_from_cfg(cfg))
    year_label = format_year_range_label(_get_plot_year_range(cfg))
    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    _draw_boundary(ax, boundary_geom, linewidth=1.0, zorder=2)
    scatter = ax.scatter(
        stations["longitude"],
        stations["latitude"],
        c=stations["elevation"],
        s=72,
        cmap="terrain",
        edgecolor="black",
        linewidth=0.35,
        zorder=3,
    )
    xmin, xmax, ymin, ymax = _get_plot_extent(stations, boundary_geom=boundary_geom, pad_deg=0.25)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    _format_geo_axis(ax, show_x=True, show_y=True)
    ax.set_facecolor("#f6f6f6")
    ax.grid(False)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Study area and meteorological station network", loc="left")
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.84)
    cbar.set_label("Elevation (m)")
    ax.text(
        0.02,
        0.03,
        f"{len(stations)} stations | {year_label}",
        transform=ax.transAxes,
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#d0d0d0", "alpha": 0.96},
    )
    fig.tight_layout()
    fig.savefig(outdir / f"ijoc_study_area.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def plot_ijoc_main_delta_maps(feature_table: pd.DataFrame, stations: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    boundary_geom = _load_boundary_geometry(_get_boundary_path_from_cfg(cfg))
    primary_delta = get_primary_delta(cfg)
    merged = feature_table.merge(stations, on=["station_id", "station_name"], how="left")
    vals = pd.to_numeric(merged[primary_delta], errors="coerce")
    vmax = max(4.0, float(np.nanmax(np.abs(vals)))) if np.isfinite(vals).any() else 4.0
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 10.2), constrained_layout=True)
    for i, idx_cfg in enumerate(cfg["indices"]):
        ax = axes.flat[i]
        sdf = merged.loc[merged["index_name"] == idx_cfg["name"]].copy()
        _draw_boundary(ax, boundary_geom, linewidth=0.95, zorder=2)
        scatter = ax.scatter(
            sdf["longitude"],
            sdf["latitude"],
            c=sdf[primary_delta],
            s=84,
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            edgecolor="black",
            linewidth=0.35,
            zorder=3,
        )
        xmin, xmax, ymin, ymax = _get_plot_extent(sdf, boundary_geom=boundary_geom, pad_deg=0.15)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")
        _format_geo_axis(ax, show_x=(i >= 2), show_y=(i % 2 == 0))
        ax.set_facecolor("#f7f7f7")
        ax.grid(False)
        ax.set_title(f"{_panel_label(i)} {idx_cfg['title']}", loc="left")
        if i >= 2:
            ax.set_xlabel("Longitude")
        if i % 2 == 0:
            ax.set_ylabel("Latitude")
    cbar = fig.colorbar(scatter, ax=axes.ravel().tolist(), shrink=0.86, pad=0.02)
    cbar.set_label(primary_delta)
    fig.savefig(outdir / f"ijoc_main_{primary_delta.lower()}_maps.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def plot_ijoc_robustness_synthesis(outdir: Path, cfg: dict):
    apply_publication_theme()
    tables_dir = outdir.parent / "tables"
    focus_quantiles = get_focus_quantiles(cfg)
    primary_delta = get_primary_delta(cfg)
    sensitivity_quantiles = get_sensitivity_quantiles(cfg)
    def _safe_read_csv(path: Path) -> pd.DataFrame:
        return pd.read_csv(path) if path.exists() else pd.DataFrame()

    ref = _safe_read_csv(tables_dir / "reference_period_sensitivity_summary.csv")
    boot = _safe_read_csv(tables_dir / "bootstrap_method_sensitivity_summary.csv")
    interp = _safe_read_csv(tables_dir / "interpolation_method_sensitivity_summary.csv")
    cluster = _safe_read_csv(tables_dir / "cluster_robustness_summary.csv")

    ref_mat = pd.DataFrame()
    if not ref.empty and {"alternative", "metric", "index_name", "mean_abs_diff"}.issubset(ref.columns):
        ref = ref[ref["alternative"] == "all_available"].copy()
        ref["metric"] = ref["metric"].map(metric_label).fillna(ref["metric"])
        ref_mat = ref.pivot(index="index_name", columns="metric", values="mean_abs_diff").reindex(
            index=["warm_days", "warm_nights", "cool_days", "cool_nights"],
            columns=[metric_label(slope_col(tau)) for tau in focus_quantiles] + [primary_delta],
        )

    boot_mat = pd.DataFrame()
    if not boot.empty and {"alternative", "metric", "index_name", "mean_abs_diff"}.issubset(boot.columns):
        current_method = str(cfg.get("bootstrap", {}).get("method", "")).lower()
        configured_methods = [
            str(x).lower()
            for x in cfg.get("advanced_analyses", {}).get("method_sensitivity", {}).get("bootstrap_methods", [])
        ]
        candidate_methods = [m for m in configured_methods if m != current_method]
        preferred_alt = candidate_methods[0] if candidate_methods else None
        if preferred_alt and preferred_alt in set(boot["alternative"].astype(str).str.lower()):
            boot = boot[boot["alternative"].astype(str).str.lower() == preferred_alt].copy()
        else:
            boot = boot.iloc[:0].copy() if boot.empty else boot[boot["alternative"] == boot["alternative"].iloc[0]].copy()
        metric_map = {}
        for tau in sensitivity_quantiles:
            qlab = tau_label(tau)
            suffix = f"{tau:0.2f}"
            metric_map[f"boot_mean_{suffix}"] = f"{qlab} mean"
            metric_map[f"boot_ci_low_{suffix}"] = f"{qlab} CI low"
            metric_map[f"boot_ci_high_{suffix}"] = f"{qlab} CI high"
        boot["metric"] = boot["metric"].replace(metric_map)
        boot_mat = boot.pivot(index="index_name", columns="metric", values="mean_abs_diff").reindex(
            index=["warm_days", "warm_nights", "cool_days", "cool_nights"],
            columns=[f"{tau_label(tau)} mean" for tau in sensitivity_quantiles] + [f"{tau_label(tau)} CI low" for tau in sensitivity_quantiles],
        )

    interp_mat = pd.DataFrame()
    if not interp.empty and {"method_left", "method_right", "tau", "index_name", "surface_correlation"}.issubset(interp.columns):
        interp = interp[(interp["method_left"] == "linear_rbf") & (interp["method_right"] == "linear")].copy()
        interp["tau_lab"] = interp["tau"].map(lambda x: tau_label(float(x)))
        interp_mat = interp.pivot(index="index_name", columns="tau_lab", values="surface_correlation").reindex(
            index=["warm_days", "warm_nights", "cool_days", "cool_nights"],
            columns=[tau_label(tau) for tau in sensitivity_quantiles],
        )

    show_ref_panel = not ref_mat.empty
    if show_ref_panel:
        fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.6), constrained_layout=True)
        ax_ref = axes[0, 0]
        ax_boot = axes[0, 1]
        ax_interp = axes[1, 0]
        ax_cluster = axes[1, 1]
        panel_prefix = {
            "ref": "(a)",
            "boot": "(b)",
            "interp": "(c)",
            "cluster": "(d)",
        }
    else:
        fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.8), constrained_layout=True)
        ax_ref = None
        ax_boot, ax_interp, ax_cluster = axes
        panel_prefix = {
            "boot": "(a)",
            "interp": "(b)",
            "cluster": "(c)",
        }

    def draw_heatmap(ax, data, title, cmap, vmin=None, vmax=None, fmt="{:.2f}", empty_message="Not available"):
        if data.empty:
            ax.axis("off")
            ax.text(0.5, 0.5, empty_message, ha="center", va="center", fontsize=10)
            ax.set_title(title, loc="left")
            return None
        im = ax.imshow(data.to_numpy(dtype=float), aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks(range(data.shape[1]))
        ax.set_xticklabels(data.columns.tolist(), rotation=0)
        ax.set_yticks(range(data.shape[0]))
        ax.set_yticklabels([x.replace("_", " ").title() for x in data.index.tolist()])
        ax.set_title(title, loc="left")
        for r in range(data.shape[0]):
            for c in range(data.shape[1]):
                value = data.iloc[r, c]
                if pd.notna(value):
                    ax.text(c, r, fmt.format(float(value)), ha="center", va="center", fontsize=8.5)
        return im

    im0 = None
    if show_ref_panel and ax_ref is not None:
        im0 = draw_heatmap(
            ax_ref,
            ref_mat,
            f"{panel_prefix['ref']} Reference-period sensitivity\nMean absolute difference",
            "YlOrBr",
        )
    bootstrap_panel_title = "(b) Bootstrap-method sensitivity\nMean absolute difference"
    current_method = str(cfg.get("bootstrap", {}).get("method", "")).lower()
    configured_methods = [
        str(x).lower()
        for x in cfg.get("advanced_analyses", {}).get("method_sensitivity", {}).get("bootstrap_methods", [])
    ]
    candidate_methods = [m for m in configured_methods if m != current_method]
    if candidate_methods:
        bootstrap_panel_title = f"{panel_prefix['boot']} Bootstrap-method sensitivity\n{current_method} vs {candidate_methods[0]}"
    else:
        bootstrap_panel_title = f"{panel_prefix['boot']} Bootstrap-method sensitivity\nMean absolute difference"
    im1 = draw_heatmap(ax_boot, boot_mat, bootstrap_panel_title, "YlOrRd")
    im2 = draw_heatmap(
        ax_interp,
        interp_mat,
        f"{panel_prefix['interp']} Interpolation agreement\nSurface correlation",
        "GnBu",
        vmin=0.5,
        vmax=1.0,
    )

    ax = ax_cluster
    if cluster.empty or "adjusted_rand_index" not in cluster.columns:
        ax.axis("off")
        ax.text(0.5, 0.5, "Not available", ha="center", va="center", fontsize=10)
        ax.set_title(f"{panel_prefix['cluster']} Clustering robustness", loc="left")
    else:
        cluster = cluster.set_index("index_name").reindex(["warm_days", "warm_nights", "cool_days", "cool_nights"])
        bars = ax.bar(
            range(len(cluster)),
            cluster["adjusted_rand_index"],
            color=["#b22222", "#2e7d32", "#1f77b4", "#9467bd"],
            edgecolor="black",
            linewidth=0.35,
        )
        ax.set_xticks(range(len(cluster)))
        ax.set_xticklabels([x.replace("_", " ").title() for x in cluster.index.tolist()], rotation=15, ha="right")
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Adjusted Rand index")
        ax.set_title(f"{panel_prefix['cluster']} Clustering robustness", loc="left")
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
        ax.grid(False, axis="x")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for bar, value in zip(bars, cluster["adjusted_rand_index"]):
            ax.text(bar.get_x() + bar.get_width()/2, value + 0.03, f"{value:.2f}", ha="center", va="bottom", fontsize=8.5)

    heatmap_items = []
    if show_ref_panel and im0 is not None and ax_ref is not None:
        heatmap_items.append((ax_ref, im0))
    if im1 is not None:
        heatmap_items.append((ax_boot, im1))
    if im2 is not None:
        heatmap_items.append((ax_interp, im2))
    for ax_, im in heatmap_items:
        shrink = 0.82 if show_ref_panel else 0.9
        cbar = fig.colorbar(im, ax=ax_, shrink=shrink)
        cbar.ax.tick_params(labelsize=8)

    fig.savefig(outdir / f"ijoc_robustness_synthesis.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def plot_ijoc_split_period_comparison(annual: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    periods = build_split_periods(_get_plot_year_range(cfg, annual))
    if len(periods) < 2:
        return
    time_scale_years = get_time_scale_years(cfg)
    split_quantiles = get_focus_quantiles(cfg)
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.6), constrained_layout=True)
    for i, idx_cfg in enumerate(cfg["indices"]):
        ax = axes.flat[i]
        vals_by_period = []
        for label, y0, y1 in periods:
            sub = annual.loc[annual["year"].between(y0, y1)]
            mean_df = sub.groupby("year", as_index=False)[idx_cfg["name"]].mean()
            years = mean_df["year"].to_numpy(dtype=float)
            values = mean_df[idx_cfg["name"]].to_numpy(dtype=float)
            ols = fit_ols_line(years, values, time_scale_years=time_scale_years)
            row_vals = [
                float(
                    fit_quantile_line(
                        years,
                        values,
                        tau=tau,
                        time_scale_years=time_scale_years,
                        max_iter=int(cfg["quantile_regression"]["max_iter"]),
                    )["slope"]
                )
                for tau in split_quantiles
            ]
            row_vals.append(float(ols["slope"]))
            vals_by_period.append(row_vals)
        vals = np.array(vals_by_period, dtype=float)
        x = np.arange(vals.shape[1])
        width = 0.34
        ax.bar(x - width/2, vals[0], width=width, color="#b0bec5", edgecolor="black", linewidth=0.3, label=periods[0][0])
        ax.bar(x + width/2, vals[1], width=width, color="#ef6c00", edgecolor="black", linewidth=0.3, label=periods[1][0])
        ax.axhline(0, color="#555555", linewidth=0.8, linestyle="--")
        ax.set_xticks(x)
        ax.set_xticklabels([tau_label(tau) for tau in split_quantiles] + ["OLS"])
        ax.set_title(f"{_panel_label(i)} {idx_cfg['title']}", loc="left")
        ax.set_ylabel(_time_unit_ylabel(idx_cfg, cfg))
        ax.grid(True, axis="y", alpha=0.22, linewidth=0.6)
        ax.grid(False, axis="x")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if i == 0:
            ax.legend(loc="upper left")
    fig.savefig(outdir / f"ijoc_split_period_comparison.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def plot_maps(features: pd.DataFrame, stations: pd.DataFrame, cluster_df: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    merged = features.merge(stations, on=["station_id", "station_name"], how="left")
    focus_quantiles = get_focus_quantiles(cfg)
    primary_delta = get_primary_delta(cfg)
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        sdf = merged.loc[merged["index_name"] == idx_name].copy()
        if sdf.empty:
            continue
        for col, lab in [(slope_col(tau), f"{tau_label(tau)} slope") for tau in focus_quantiles] + [(primary_delta, primary_delta)]:
            _scatter_station_map(sdf, col, f"Spatial distribution of {lab} - {idx_cfg['title']}", outdir / f"map_{idx_name}_{col}.{cfg['plots']['save_format']}", cfg)
        if not cluster_df.empty:
            # `features` may already include `cluster` from pipeline merge.
            # Merging again can create suffixes (cluster_x/cluster_y), so we resolve safely.
            csdf = sdf.merge(
                cluster_df.loc[cluster_df["index_name"] == idx_name, ["station_id", "cluster"]],
                on="station_id",
                how="left",
                suffixes=("", "_from_cluster_df"),
            )
            if "cluster" not in csdf.columns and "cluster_from_cluster_df" in csdf.columns:
                csdf["cluster"] = csdf["cluster_from_cluster_df"]
            if "cluster" in csdf.columns:
                csdf_cluster = csdf.dropna(subset=["cluster"]).copy()
                if csdf_cluster.empty:
                    continue
                _scatter_station_map(
                    csdf_cluster,
                    "cluster",
                    f"Station clusters - {idx_cfg['title']}",
                    outdir / f"map_{idx_name}_cluster.{cfg['plots']['save_format']}",
                    cfg,
                    True,
                )


def plot_data_coverage(annual: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    counts = annual.groupby("station_name")["year"].agg(["min", "max", "count"]).sort_values("count", ascending=False)
    fig, ax = plt.subplots(figsize=(10, 8))
    y = np.arange(len(counts))
    ax.barh(y, counts["count"], color="#2f6c8f")
    ax.set_yticks(y)
    ax.set_yticklabels(counts.index)
    ax.invert_yaxis()
    ax.set_xlabel("Number of annual index values")
    ax.set_title("Data coverage by station")
    fig.tight_layout()
    fig.savefig(outdir / f"data_coverage_by_station.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def plot_delta_uncertainty(features: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    primary_delta = get_primary_delta(cfg)
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        cols = ["station_name", f"boot_mean_{primary_delta}", f"boot_ci_low_{primary_delta}", f"boot_ci_high_{primary_delta}"]
        if not set(cols).issubset(features.columns):
            continue
        sdf = features.loc[features["index_name"] == idx_name, cols].copy().dropna().sort_values(f"boot_mean_{primary_delta}", ascending=False)
        if sdf.empty:
            continue
        fig, ax = plt.subplots(figsize=(9, max(6, len(sdf) * 0.28)))
        y = np.arange(len(sdf))
        x = sdf[f"boot_mean_{primary_delta}"].to_numpy()
        xerr = np.vstack([x - sdf[f"boot_ci_low_{primary_delta}"].to_numpy(), sdf[f"boot_ci_high_{primary_delta}"].to_numpy() - x])
        ax.errorbar(x, y, xerr=xerr, fmt="o", capsize=3, color="#1f4e79")
        ax.axvline(0, color="black", linewidth=0.9)
        ax.set_yticks(y)
        ax.set_yticklabels(sdf["station_name"].tolist())
        ax.invert_yaxis()
        ax.set_xlabel(f"Bootstrap mean {primary_delta} with 95% CI")
        ax.set_title(f"Asymmetric trend uncertainty - {idx_cfg['title']}")
        fig.tight_layout()
        fig.savefig(outdir / f"delta_uncertainty_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)
    split_quantiles = get_plot_quantiles(cfg, "split_period_bars")
    time_scale_years = get_time_scale_years(cfg)
