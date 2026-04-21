from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram

from .quantile import fit_quantile_slope
from .math_utils import make_quantile_grid


def apply_publication_theme() -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update({
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.dpi": 120,
    })


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
    fig, ax = plt.subplots(figsize=(7.5, 6))
    scatter = ax.scatter(df["longitude"], df["latitude"], c=df[value_col], s=75, cmap=("tab10" if discrete else "RdBu_r"), edgecolor="black", linewidth=0.35)
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
