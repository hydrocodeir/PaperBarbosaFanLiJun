from __future__ import annotations

from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config_utils import (
    boot_ci_high_col,
    boot_ci_low_col,
    boot_mean_col,
    get_focus_quantiles,
    get_plot_dpi,
    get_primary_delta,
    metric_label,
)
from .progress_utils import ProgressTracker
from .quantile import bootstrap_qr, build_station_seed, summarize_bootstrap


def _station_summary_from_boot(
    annual: pd.DataFrame,
    cfg: dict,
    index_name: str,
    n_reps: int,
    progress_callback: Callable[[str], None] | None = None,
) -> pd.DataFrame:
    station_col = cfg["data"]["station_id_col"]
    station_name_col = cfg["data"]["station_name_col"]
    year_col = cfg["data"]["year_col"]
    alpha = float(cfg["bootstrap"]["alpha"])
    min_years = int(cfg["quantile_regression"]["min_years_required_to_run"])
    base_seed = int(cfg["project"]["random_seed"])

    boot_cfg = {**cfg}
    boot_cfg["bootstrap"] = {**cfg["bootstrap"], "n_reps": int(n_reps)}
    focus_q = get_focus_quantiles(cfg)
    primary_delta = get_primary_delta(cfg)

    rows: list[dict] = []
    sub = annual[[station_col, station_name_col, year_col, index_name]].copy()

    station_groups = list(sub.groupby([station_col, station_name_col]))
    total_stations = len(station_groups)
    tracker = ProgressTracker(
        f"Bootstrap-depth sensitivity [{index_name}, {n_reps} reps]",
        total_stations,
        progress_callback,
    )
    for station_no, ((station_id, station_name), sdf) in enumerate(station_groups, start=1):
        tracker.emit(station_no, detail=f"station={station_name} ({station_id})")
        sdf = sdf.sort_values(year_col)
        years = sdf[year_col].to_numpy(dtype=float)
        values = sdf[index_name].to_numpy(dtype=float)
        n_years = int(np.isfinite(values).sum())
        if n_years < min_years:
            continue

        rng = np.random.default_rng(build_station_seed(base_seed, index_name, station_id))
        boot = bootstrap_qr(years, values, focus_q, boot_cfg, rng)
        if boot.empty:
            continue

        summary = summarize_bootstrap(boot, alpha)
        row = {
            "index_name": index_name,
            "station_id": station_id,
            "station_name": station_name,
            "n_years": n_years,
            "n_reps": int(n_reps),
        }
        row.update(summary)
        for tau in focus_q:
            suffix = f"{tau:0.2f}"
            row[f"boot_ci_width_{suffix}"] = row[boot_ci_high_col(suffix)] - row[boot_ci_low_col(suffix)]
        row[f"boot_ci_width_{primary_delta}"] = row[boot_ci_high_col(primary_delta)] - row[boot_ci_low_col(primary_delta)]
        rows.append(row)

    return pd.DataFrame(rows)


def _regional_comparison_table(
    baseline: pd.DataFrame,
    deeper: pd.DataFrame,
    baseline_reps: int,
    deeper_reps: int,
    target_indices: list[str],
    focus_quantiles: list[float],
    primary_delta: str,
) -> pd.DataFrame:
    records: list[dict] = []
    metrics = [boot_mean_col(f"{tau:0.2f}") for tau in focus_quantiles]
    metrics.append(boot_mean_col(primary_delta))
    metrics.extend([f"boot_ci_width_{tau:0.2f}" for tau in focus_quantiles])
    metrics.append(f"boot_ci_width_{primary_delta}")
    for index_name in target_indices:
        bdf = baseline[baseline["index_name"] == index_name].copy()
        ddf = deeper[deeper["index_name"] == index_name].copy()
        merged = bdf.merge(
            ddf,
            on=["index_name", "station_id", "station_name"],
            suffixes=(f"_{baseline_reps}", f"_{deeper_reps}"),
        )
        if merged.empty:
            continue

        for metric in metrics:
            old_col = f"{metric}_{baseline_reps}"
            new_col = f"{metric}_{deeper_reps}"
            old_vals = pd.to_numeric(merged[old_col], errors="coerce")
            new_vals = pd.to_numeric(merged[new_col], errors="coerce")
            diff = new_vals - old_vals
            corr = old_vals.corr(new_vals)
            records.append(
                {
                    "index_name": index_name,
                    "metric": metric,
                    f"mean_{baseline_reps}": float(old_vals.mean()),
                    f"mean_{deeper_reps}": float(new_vals.mean()),
                    "mean_difference": float(diff.mean()),
                    "median_abs_station_difference": float(diff.abs().median()),
                    "max_abs_station_difference": float(diff.abs().max()),
                    "station_correlation": float(corr) if pd.notna(corr) else np.nan,
                    "n_stations": int(len(merged)),
                }
            )
    return pd.DataFrame(records)


def _plot_bootstrap_depth_sensitivity(
    comparison: pd.DataFrame,
    outpath: Path,
    baseline_reps: int,
    deeper_reps: int,
    target_indices: list[str],
    focus_quantiles: list[float],
    primary_delta: str,
    cfg: dict,
) -> None:
    plot_metrics = [(boot_mean_col(f"{tau:0.2f}"), f"{metric_label(f'slope_{tau:0.2f}')} mean") for tau in focus_quantiles]
    plot_metrics.append((boot_mean_col(primary_delta), f"{primary_delta} mean"))
    n_panels = len(plot_metrics)
    ncols = 2
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(10, max(7, 3.2 * nrows)), constrained_layout=True)
    axes_arr = np.atleast_1d(axes).ravel()
    palette = ["#1f1f1f", "#b22222", "#1d3557", "#2a9d8f", "#6d597a", "#8d6e63"]
    color_map = {name: palette[i % len(palette)] for i, name in enumerate(target_indices)}

    for ax, (metric, label) in zip(axes_arr, plot_metrics):
        sub = comparison[comparison["metric"] == metric].copy()
        if sub.empty:
            ax.axis("off")
            continue
        x = np.arange(len(sub))
        old_col = f"mean_{baseline_reps}"
        new_col = f"mean_{deeper_reps}"
        width = 0.34
        for i, (_, row) in enumerate(sub.iterrows()):
            color = color_map.get(row["index_name"], "#555555")
            ax.bar(i - width / 2, row[old_col], width=width, color=color, alpha=0.55, edgecolor="black", linewidth=0.6)
            ax.bar(i + width / 2, row[new_col], width=width, color=color, alpha=0.9, edgecolor="black", linewidth=0.6)
            ax.text(
                i,
                max(row[old_col], row[new_col]) + (0.04 * max(1.0, abs(max(row[old_col], row[new_col])))),
                f"r={row['station_correlation']:.2f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )
        ax.set_xticks(x)
        ax.set_xticklabels([name.replace("_", " ").title() for name in sub["index_name"]], rotation=0)
        ax.set_title(label, fontsize=10)
        ax.grid(axis="y", alpha=0.25, linewidth=0.6)
        ax.set_axisbelow(True)
    for ax in axes_arr[len(plot_metrics):]:
        ax.axis("off")

    handles = [
        plt.Line2D([0], [0], color="black", lw=8, alpha=0.55, label=f"{baseline_reps} replicates"),
        plt.Line2D([0], [0], color="black", lw=8, alpha=0.9, label=f"{deeper_reps} replicates"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Bootstrap-depth sensitivity for configured regional summaries", fontsize=12, y=1.04)
    fig.savefig(outpath, dpi=get_plot_dpi(cfg), bbox_inches="tight")
    plt.close(fig)
def run_bootstrap_depth_sensitivity(
    annual: pd.DataFrame,
    cfg: dict,
    outdir: Path,
    progress_callback: Callable[[str], None] | None = None,
) -> dict[str, Path]:
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures"
    baseline_reps = int(cfg["bootstrap"]["n_reps"])
    deeper_reps = int(cfg.get("advanced_analyses", {}).get("bootstrap_depth_sensitivity", {}).get("n_reps", 800))
    target_indices = [idx_cfg["name"] for idx_cfg in cfg["indices"]]
    focus_quantiles = get_focus_quantiles(cfg)
    primary_delta = get_primary_delta(cfg)

    baseline_rows = []
    pub = pd.read_csv(tables_dir / "publication_summary_table.csv")
    if progress_callback is not None:
        progress_callback("Bootstrap-depth sensitivity: preparing baseline summary from publication table")
    for index_name in target_indices:
        sub = pub[pub["index_name"] == index_name].copy()
        if sub.empty:
            continue
        keep = ["index_name", "station_id", "station_name"]
        keep.extend([boot_mean_col(f"{tau:0.2f}") for tau in focus_quantiles])
        keep.append(boot_mean_col(primary_delta))
        for tau in focus_quantiles:
            suffix = f"{tau:0.2f}"
            keep.extend([boot_ci_low_col(suffix), boot_ci_high_col(suffix)])
        keep.extend([boot_ci_low_col(primary_delta), boot_ci_high_col(primary_delta)])
        tmp = sub[keep].copy()
        tmp["n_reps"] = baseline_reps
        for tau in focus_quantiles:
            suffix = f"{tau:0.2f}"
            tmp[f"boot_ci_width_{suffix}"] = tmp[boot_ci_high_col(suffix)] - tmp[boot_ci_low_col(suffix)]
        tmp[f"boot_ci_width_{primary_delta}"] = tmp[boot_ci_high_col(primary_delta)] - tmp[boot_ci_low_col(primary_delta)]
        baseline_rows.append(tmp)
    baseline = pd.concat(baseline_rows, ignore_index=True) if baseline_rows else pd.DataFrame()

    deeper_frames = []
    total_indices = len(target_indices)
    index_tracker = ProgressTracker(
        f"Bootstrap-depth sensitivity index loop [{deeper_reps} reps]",
        total_indices,
        progress_callback,
    )
    for index_no, index_name in enumerate(target_indices, start=1):
        index_tracker.emit(index_no, detail=f"index={index_name}")
        deeper_frames.append(
            _station_summary_from_boot(
                annual,
                cfg,
                index_name,
                deeper_reps,
                progress_callback=progress_callback,
            )
        )
    deeper = pd.concat(deeper_frames, ignore_index=True) if deeper_frames else pd.DataFrame()

    if progress_callback is not None:
        progress_callback("Bootstrap-depth sensitivity: summarizing baseline vs deeper-bootstrap comparison")
    comparison = _regional_comparison_table(
        baseline,
        deeper,
        baseline_reps,
        deeper_reps,
        target_indices,
        focus_quantiles,
        primary_delta,
    )
    station_compare = baseline.merge(
        deeper,
        on=["index_name", "station_id", "station_name"],
        suffixes=(f"_{baseline_reps}", f"_{deeper_reps}"),
    )
    comparison_path = tables_dir / "bootstrap_depth_sensitivity_summary.csv"
    station_path = tables_dir / "bootstrap_depth_sensitivity_station_comparison.csv"
    fig_path = figs_dir / "ijoc_bootstrap_depth_sensitivity.png"

    comparison.to_csv(comparison_path, index=False)
    station_compare.to_csv(station_path, index=False)
    if progress_callback is not None:
        progress_callback("Bootstrap-depth sensitivity: exporting summary figure and tables")
    _plot_bootstrap_depth_sensitivity(
        comparison,
        fig_path,
        baseline_reps,
        deeper_reps,
        target_indices,
        focus_quantiles,
        primary_delta,
        cfg,
    )

    return {
        "summary_table": comparison_path,
        "station_table": station_path,
        "figure": fig_path,
    }
