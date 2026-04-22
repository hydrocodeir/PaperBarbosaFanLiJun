from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .quantile import bootstrap_qr, build_station_seed, summarize_bootstrap


TARGET_INDICES = ["warm_days", "warm_nights"]
TARGET_QUANTILES = [0.05, 0.50, 0.95]
def _station_summary_from_boot(
    annual: pd.DataFrame,
    cfg: dict,
    index_name: str,
    n_reps: int,
) -> pd.DataFrame:
    station_col = cfg["data"]["station_id_col"]
    station_name_col = cfg["data"]["station_name_col"]
    year_col = cfg["data"]["year_col"]
    alpha = float(cfg["bootstrap"]["alpha"])
    min_years = int(cfg["quantile_regression"]["min_years_required_to_run"])
    base_seed = int(cfg["project"]["random_seed"])

    boot_cfg = {**cfg}
    boot_cfg["bootstrap"] = {**cfg["bootstrap"], "n_reps": int(n_reps)}
    focus_q = [float(x) for x in TARGET_QUANTILES]

    rows: list[dict] = []
    sub = annual[[station_col, station_name_col, year_col, index_name]].copy()

    for (station_id, station_name), sdf in sub.groupby([station_col, station_name_col]):
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
        for tau in TARGET_QUANTILES:
            suffix = f"{tau:0.2f}"
            row[f"boot_ci_width_{suffix}"] = row[f"boot_ci_high_{suffix}"] - row[f"boot_ci_low_{suffix}"]
        row["boot_ci_width_Delta1"] = row["boot_ci_high_Delta1"] - row["boot_ci_low_Delta1"]
        rows.append(row)

    return pd.DataFrame(rows)


def _regional_comparison_table(
    baseline: pd.DataFrame,
    deeper: pd.DataFrame,
    baseline_reps: int,
    deeper_reps: int,
) -> pd.DataFrame:
    records: list[dict] = []
    for index_name in TARGET_INDICES:
        bdf = baseline[baseline["index_name"] == index_name].copy()
        ddf = deeper[deeper["index_name"] == index_name].copy()
        merged = bdf.merge(
            ddf,
            on=["index_name", "station_id", "station_name"],
            suffixes=(f"_{baseline_reps}", f"_{deeper_reps}"),
        )
        if merged.empty:
            continue

        for metric in [
            "boot_mean_0.05",
            "boot_mean_0.50",
            "boot_mean_0.95",
            "boot_mean_Delta1",
            "boot_ci_width_0.05",
            "boot_ci_width_0.50",
            "boot_ci_width_0.95",
            "boot_ci_width_Delta1",
        ]:
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
) -> None:
    plot_metrics = [
        ("boot_mean_0.05", "q0.05 mean"),
        ("boot_mean_0.50", "q0.50 mean"),
        ("boot_mean_0.95", "q0.95 mean"),
        ("boot_mean_Delta1", "Delta1 mean"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    color_map = {"warm_days": "#1f1f1f", "warm_nights": "#b22222"}

    for ax, (metric, label) in zip(axes.ravel(), plot_metrics):
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
            ax.bar(
                i - width / 2,
                row[old_col],
                width=width,
                color=color,
                alpha=0.55,
                edgecolor="black",
                linewidth=0.6,
            )
            ax.bar(
                i + width / 2,
                row[new_col],
                width=width,
                color=color,
                alpha=0.9,
                edgecolor="black",
                linewidth=0.6,
            )
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

    handles = [
        plt.Line2D([0], [0], color="black", lw=8, alpha=0.55, label=f"{baseline_reps} replicates"),
        plt.Line2D([0], [0], color="black", lw=8, alpha=0.9, label=f"{deeper_reps} replicates"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Bootstrap-depth sensitivity for core regional summaries", fontsize=12, y=1.04)
    fig.savefig(outpath, dpi=400, bbox_inches="tight")
    plt.close(fig)


def run_bootstrap_depth_sensitivity(annual: pd.DataFrame, cfg: dict, outdir: Path) -> dict[str, Path]:
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures"
    baseline_reps = int(cfg["bootstrap"]["n_reps"])
    deeper_reps = int(cfg.get("advanced_analyses", {}).get("bootstrap_depth_sensitivity", {}).get("n_reps", 800))

    baseline_rows = []
    pub = pd.read_csv(tables_dir / "publication_summary_table.csv")
    for index_name in TARGET_INDICES:
        sub = pub[pub["index_name"] == index_name].copy()
        if sub.empty:
            continue
        keep = [
            "index_name",
            "station_id",
            "station_name",
            "boot_mean_0.05",
            "boot_mean_0.50",
            "boot_mean_0.95",
            "boot_mean_Delta1",
            "boot_ci_low_0.05",
            "boot_ci_high_0.05",
            "boot_ci_low_0.50",
            "boot_ci_high_0.50",
            "boot_ci_low_0.95",
            "boot_ci_high_0.95",
            "boot_ci_low_Delta1",
            "boot_ci_high_Delta1",
        ]
        tmp = sub[keep].copy()
        tmp["n_reps"] = baseline_reps
        for tau in ["0.05", "0.50", "0.95"]:
            tmp[f"boot_ci_width_{tau}"] = tmp[f"boot_ci_high_{tau}"] - tmp[f"boot_ci_low_{tau}"]
        tmp["boot_ci_width_Delta1"] = tmp["boot_ci_high_Delta1"] - tmp["boot_ci_low_Delta1"]
        baseline_rows.append(tmp)
    baseline = pd.concat(baseline_rows, ignore_index=True)

    deeper_frames = []
    for index_name in TARGET_INDICES:
        deeper_frames.append(_station_summary_from_boot(annual, cfg, index_name, deeper_reps))
    deeper = pd.concat(deeper_frames, ignore_index=True)

    comparison = _regional_comparison_table(baseline, deeper, baseline_reps, deeper_reps)
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
    _plot_bootstrap_depth_sensitivity(comparison, fig_path, baseline_reps, deeper_reps)

    return {
        "summary_table": comparison_path,
        "station_table": station_path,
        "figure": fig_path,
    }
