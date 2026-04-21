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

from .quantile import fit_ols_line, fit_quantile_line, fit_quantile_slope
from .math_utils import make_quantile_grid


def apply_publication_theme() -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
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
        "figure.dpi": 120,
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


def _style_coeff_axis(ax, idx_cfg: dict, panel_idx: int) -> None:
    ax.set_title(f"{_panel_label(panel_idx)} {idx_cfg['title']}", loc="left")
    ax.set_xlabel("Quantile, τ")
    ax.set_ylabel(f"Trend ({idx_cfg['unit']} per decade)")
    ax.set_xlim(0.01, 0.99)
    ax.set_xticks([0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99])
    ax.axhline(0, color="#3a3a3a", linewidth=0.8, linestyle="--", alpha=0.85, zorder=1)
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.grid(False, axis="x")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


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


def plot_station_paper2_figures(annual: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    qgrid = make_quantile_grid(0.01, 0.99, 0.01)
    fig_root = outdir / "paper2_station_figures"
    fig1_dir = fig_root / "figure1_timeseries"
    fig2_dir = fig_root / "figure2_quantile_coefficients"
    fig1_dir.mkdir(parents=True, exist_ok=True)
    fig2_dir.mkdir(parents=True, exist_ok=True)

    qr_colors = {
        0.10: "#c62828",
        0.50: "#1f1f1f",
        0.90: "#2e7d32",
    }

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
            for tau in (0.10, 0.50, 0.90):
                qr_fit = fit_quantile_line(years, values, tau=tau, max_iter=max_iter)
                if len(qr_fit["fitted"]):
                    ax1.plot(
                        qr_fit["years"],
                        qr_fit["fitted"],
                        color=qr_colors[tau],
                        linewidth=2.4 if tau == 0.50 else 2.1,
                        zorder=3,
                        label=f"QR τ={tau:.2f}",
                    )
                annotation_lines.append(f"τ={tau:.2f}: {_format_slope(float(qr_fit['slope']))}")

            ols_fit = fit_ols_line(years, values)
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
            _annotate_textbox(ax1, "Slope (days/decade)\n" + "\n".join(annotation_lines))
            if i == 0:
                ax1.legend(loc="upper right", ncol=2)

            coeff_rows = []
            for tau in qgrid:
                qr_fit = fit_quantile_line(years, values, tau=float(tau), max_iter=max_iter)
                coeff_rows.append({
                    "tau": float(tau),
                    "slope": float(qr_fit["slope"]),
                    "ci_low": float(qr_fit["ci_low"]),
                    "ci_high": float(qr_fit["ci_high"]),
                })
            coeff_df = pd.DataFrame(coeff_rows).dropna(subset=["tau", "slope"])

            if not coeff_df.empty:
                shade_df = coeff_df.dropna(subset=["ci_low", "ci_high"])
                if not shade_df.empty:
                    ax2.fill_between(
                        shade_df["tau"].to_numpy(),
                        shade_df["ci_low"].to_numpy(),
                        shade_df["ci_high"].to_numpy(),
                        color="#bdbdbd",
                        alpha=0.45,
                        linewidth=0,
                        zorder=1,
                        label="QR 95% CI" if i == 0 else None,
                    )
                ax2.plot(coeff_df["tau"], coeff_df["slope"], color="#1d1d1d", linewidth=1.2, alpha=0.9, zorder=2)
                ax2.scatter(coeff_df["tau"], coeff_df["slope"], s=12, color="#111111", zorder=3, label="QR slope" if i == 0 else None)

            if np.isfinite(float(ols_fit["slope"])):
                ax2.axhline(float(ols_fit["slope"]), color="#c62828", linewidth=2.0, zorder=4, label="LSM mean trend" if i == 0 else None)
            if np.isfinite(float(ols_fit["ci_low"])) and np.isfinite(float(ols_fit["ci_high"])):
                ax2.axhline(float(ols_fit["ci_low"]), color="#c62828", linewidth=1.2, linestyle="--", zorder=4, label="LSM 95% CI" if i == 0 else None)
                ax2.axhline(float(ols_fit["ci_high"]), color="#c62828", linewidth=1.2, linestyle="--", zorder=4)

            _style_coeff_axis(ax2, idx_cfg, i)
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


def plot_paper2_figure3_maps(annual: pd.DataFrame, stations: pd.DataFrame, outdir: Path, cfg: dict, taus: List[float] | None = None):
    apply_publication_theme()
    fig_dir = outdir / "paper2_figure3_quantile_maps"
    fig_dir.mkdir(parents=True, exist_ok=True)
    taus = taus or [0.05, 0.10, 0.50, 0.90, 0.95]
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    spatial_cfg = cfg.get("spatial_visualization", {})
    boundary_path = Path(spatial_cfg.get("iran_boundary_geojson", "data/Iran.geojson"))
    interpolation_method = spatial_cfg.get("interpolation_method", "thin_plate_spline")
    interpolation_smooth = float(spatial_cfg.get("interpolation_smooth", 0.35))
    boundary_geom = _load_boundary_geometry(boundary_path)
    for tau in taus:
        rows = []
        slope_col = f"slope_{tau:0.2f}"
        ci_low_col = f"ci_low_{tau:0.2f}"
        ci_high_col = f"ci_high_{tau:0.2f}"

        for idx_cfg in cfg["indices"]:
            idx_name = idx_cfg["name"]
            sdf = annual[["station_id", "station_name", "year", idx_name]].copy().sort_values(["station_id", "year"])
            for (station_id, station_name), gdf in sdf.groupby(["station_id", "station_name"]):
                years = gdf["year"].to_numpy(dtype=float)
                values = gdf[idx_name].to_numpy(dtype=float)
                fit = fit_quantile_line(years, values, tau=float(tau), max_iter=max_iter)
                rows.append(
                    {
                        "index_name": idx_name,
                        "station_id": station_id,
                        "station_name": station_name,
                        slope_col: float(fit["slope"]),
                        ci_low_col: float(fit["ci_low"]),
                        ci_high_col: float(fit["ci_high"]),
                    }
                )

        map_df = pd.DataFrame(rows).merge(stations, on=["station_id", "station_name"], how="left")
        sig_mask = (
            np.isfinite(map_df[ci_low_col])
            & np.isfinite(map_df[ci_high_col])
            & ((map_df[ci_low_col] > 0) | (map_df[ci_high_col] < 0))
        )
        map_df["significant_95"] = sig_mask
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

            if sdf["significant_95"].any():
                sig_df = sdf.loc[sdf["significant_95"]]
            else:
                sig_df = sdf.iloc[0:0]

            nonsig_df = sdf.loc[~sdf["significant_95"]]

            if not nonsig_df.empty:
                ax.scatter(
                    nonsig_df["longitude"],
                    nonsig_df["latitude"],
                    s=28,
                    facecolor="white",
                    edgecolor="#5f6368",
                    linewidth=0.55,
                    alpha=0.95,
                    zorder=5,
                    label="Not significant" if i == 0 else None,
                )

            if not sig_df.empty:
                ax.scatter(
                    sig_df["longitude"],
                    sig_df["latitude"],
                    s=68,
                    facecolor="#1fa3b8",
                    edgecolor="white",
                    linewidth=0.9,
                    alpha=0.98,
                    zorder=6,
                    label="Significant at 95% CI" if i == 0 else None,
                )

            ax.scatter(
                sdf["longitude"],
                sdf["latitude"],
                c=sdf[slope_col],
                s=18,
                cmap=FIG3_CMAP,
                norm=norm,
                edgecolor="none",
                alpha=0.9,
                zorder=4,
            )

            xmin, xmax, ymin, ymax = _get_plot_extent(sdf, boundary_geom=boundary_geom, pad_deg=0.1)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_aspect("equal", adjustable="box")
            row_idx, col_idx = divmod(i, 2)
            _format_geo_axis(ax, show_x=(row_idx == 1), show_y=(col_idx == 0))

            n_sig = int(sdf["significant_95"].sum())
            _annotate_textbox_bottom_left(
                ax,
                f"Range: {_format_slope(float(pd.to_numeric(sdf[slope_col], errors='coerce').min()))} to {_format_slope(float(pd.to_numeric(sdf[slope_col], errors='coerce').max()))}\n"
                f"95% significant: {n_sig}/{len(sdf)}",
            )

        if mappable is not None:
            cbar = fig.colorbar(mappable, ax=axes.ravel().tolist(), shrink=0.88, pad=0.02, ticks=levels)
            cbar.set_label("Slope (days per decade)")

        fig.text(
            0.5,
            -0.01,
            "Filled teal circles denote stations significant at the 95% level; open grey circles denote non-significant stations. Background shading shows the interpolated quantile-regression slope field.",
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
    taus = [0.05, 0.50, 0.95]

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
    figure_map = [
        (5, 0.05),
        (6, 0.50),
        (7, 0.95),
    ]

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
    order_df = plot_df.groupby("station_name", as_index=False)["Delta1"].mean().sort_values("Delta1", ascending=False)
    return order_df["station_name"].tolist()


def plot_station_heatmap(features: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        sdf = features.loc[features["index_name"] == idx_name, ["station_name", "slope_0.05", "slope_0.50", "slope_0.95", "Delta1"]].copy()
        if sdf.empty:
            continue
        order = _sort_station_names_for_heatmap(sdf, cfg["plots"]["heatmap_station_order"])
        sdf["station_name"] = pd.Categorical(sdf["station_name"], categories=order, ordered=True)
        mat = sdf.sort_values("station_name").set_index("station_name")[["slope_0.05", "slope_0.50", "slope_0.95", "Delta1"]]

        fig, ax = plt.subplots(figsize=(8, max(6, len(mat) * 0.3)))
        im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="RdBu_r")
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels(["q0.05", "q0.50", "q0.95", "Δ(q95-q05)"])
        ax.set_yticks(range(mat.shape[0]))
        ax.set_yticklabels(mat.index.tolist())
        ax.set_title(f"Station quantile trend signatures - {idx_cfg['title']}")
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label(f"Trend ({idx_cfg['unit']} per decade)")
        fig.tight_layout()
        fig.savefig(outdir / f"station_focus_heatmap_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)


def plot_region_quantile_slopes(annual: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    years = np.sort(annual["year"].unique())
    mean_df = annual.groupby("year", as_index=False)[[idx["name"] for idx in cfg["indices"]]].mean()
    qgrid = make_quantile_grid(**cfg["quantile_regression"]["full_quantiles"])
    focus = [float(x) for x in cfg["quantile_regression"]["focus_quantiles"]]

    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        vals = mean_df[idx_name].to_numpy(dtype=float)
        slopes = [fit_quantile_slope(years, vals, tau=q, max_iter=int(cfg["quantile_regression"]["max_iter"])) for q in qgrid]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(qgrid, slopes, color="#1f4e79", linewidth=1.6)
        for fq in focus:
            fq_slope = fit_quantile_slope(years, vals, tau=fq, max_iter=int(cfg["quantile_regression"]["max_iter"]))
            ax.scatter([fq], [fq_slope], s=35, color="#b22222", zorder=3)
        ax.axhline(0, color="black", linewidth=0.9)
        ax.set_xlabel("Quantile τ")
        ax.set_ylabel(f"Slope ({idx_cfg['unit']} / decade)")
        ax.set_title(f"Regional quantile-regression slope profile - {idx_cfg['title']}")
        fig.tight_layout()
        fig.savefig(outdir / f"region_quantile_slopes_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
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
    boundary_geom = _load_boundary_geometry(Path("data") / "Iran.geojson")
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


def plot_maps(features: pd.DataFrame, stations: pd.DataFrame, cluster_df: pd.DataFrame, outdir: Path, cfg: dict):
    apply_publication_theme()
    merged = features.merge(stations, on=["station_id", "station_name"], how="left")
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        sdf = merged.loc[merged["index_name"] == idx_name].copy()
        if sdf.empty:
            continue
        for col, lab in [("slope_0.05", "q0.05 slope"), ("slope_0.50", "q0.50 slope"), ("slope_0.95", "q0.95 slope"), ("Delta1", "Δ(q95-q05)")]:
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
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        cols = ["station_name", "boot_mean_Delta1", "boot_ci_low_Delta1", "boot_ci_high_Delta1"]
        if not set(cols).issubset(features.columns):
            continue
        sdf = features.loc[features["index_name"] == idx_name, cols].copy().dropna().sort_values("boot_mean_Delta1", ascending=False)
        if sdf.empty:
            continue
        fig, ax = plt.subplots(figsize=(9, max(6, len(sdf) * 0.28)))
        y = np.arange(len(sdf))
        x = sdf["boot_mean_Delta1"].to_numpy()
        xerr = np.vstack([x - sdf["boot_ci_low_Delta1"].to_numpy(), sdf["boot_ci_high_Delta1"].to_numpy() - x])
        ax.errorbar(x, y, xerr=xerr, fmt="o", capsize=3, color="#1f4e79")
        ax.axvline(0, color="black", linewidth=0.9)
        ax.set_yticks(y)
        ax.set_yticklabels(sdf["station_name"].tolist())
        ax.invert_yaxis()
        ax.set_xlabel("Bootstrap mean Δ(q95-q05) with 95% CI")
        ax.set_title(f"Asymmetric trend uncertainty - {idx_cfg['title']}")
        fig.tight_layout()
        fig.savefig(outdir / f"delta_uncertainty_{idx_name}.{cfg['plots']['save_format']}", dpi=int(cfg["plots"]["dpi"]))
        plt.close(fig)
