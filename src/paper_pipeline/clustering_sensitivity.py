from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .clustering import compare_clusterings, run_clustering
from .config_utils import get_plot_dpi


def _method_specs(cfg: dict) -> list[dict]:
    return list(cfg.get("advanced_analyses", {}).get("alternative_clustering_methods", {}).get("methods", []))


def _method_label(spec: dict) -> str:
    if spec["algorithm"] == "hierarchical":
        return f"hierarchical_{spec['linkage']}_{spec['metric']}"
    return str(spec["algorithm"])


def _plot_alternative_clustering_heatmap(summary: pd.DataFrame, outpath: Path, cfg: dict) -> None:
    if summary.empty:
        return

    pivot = summary.pivot(index="index_name", columns="method_label", values="adjusted_rand_index")
    pivot = pivot.reindex(sorted(pivot.index), axis=0)
    pivot = pivot.reindex(sorted(pivot.columns), axis=1)

    fig, ax = plt.subplots(figsize=(10, 4.8), constrained_layout=True)
    im = ax.imshow(pivot.to_numpy(dtype=float), cmap="Greys", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels([col.replace("_", "\n", 1).replace("_", "\n") for col in pivot.columns], fontsize=9)
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels([idx.replace("_", " ").title() for idx in pivot.index], fontsize=10)
    ax.set_title("Alternative-clustering sensitivity (Adjusted Rand index)", fontsize=12)

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.iat[i, j]
            ax.text(
                j,
                i,
                f"{value:.2f}",
                ha="center",
                va="center",
                color="white" if value < 0.45 else "black",
                fontsize=9,
            )

    cbar = fig.colorbar(im, ax=ax, shrink=0.9)
    cbar.set_label("Adjusted Rand index", fontsize=10)
    fig.savefig(outpath, dpi=get_plot_dpi(cfg), bbox_inches="tight")
    plt.close(fig)


def run_alternative_clustering_sensitivity(feature_table: pd.DataFrame, cfg: dict, outdir: Path) -> dict[str, Path]:
    specs = _method_specs(cfg)
    if feature_table.empty or not specs:
        return {}

    figs_dir = outdir / "figures"
    tables_dir = outdir / "tables"

    baseline = run_clustering(feature_table, cfg, label_col="cluster_baseline")[0]
    rows = []
    assignment_frames = [baseline]

    for spec in specs:
        label = _method_label(spec)
        alt_df, _ = run_clustering(
            feature_table,
            cfg,
            label_col=f"cluster_{label}",
            algorithm=spec["algorithm"],
            linkage_method=spec.get("linkage"),
            metric=spec.get("metric"),
        )
        if alt_df.empty:
            continue
        assignment_frames.append(alt_df)
        comp = compare_clusterings(
            baseline,
            alt_df,
            primary_label="cluster_baseline",
            secondary_label=f"cluster_{label}",
        )
        if comp.empty:
            continue
        comp["method_label"] = label
        comp["algorithm"] = spec["algorithm"]
        comp["linkage"] = spec.get("linkage", "")
        comp["metric"] = spec.get("metric", "")
        rows.append(comp)

    if not rows:
        return {}

    summary = pd.concat(rows, ignore_index=True)
    summary = summary[
        [
            "index_name",
            "method_label",
            "algorithm",
            "linkage",
            "metric",
            "n_stations_compared",
            "adjusted_rand_index",
            "same_label_fraction",
        ]
    ].sort_values(["index_name", "adjusted_rand_index"], ascending=[True, False])

    assignments = assignment_frames[0]
    for frame in assignment_frames[1:]:
        assignments = assignments.merge(frame, on=["index_name", "station_id", "station_name"], how="left")

    summary_path = tables_dir / "alternative_clustering_sensitivity_summary.csv"
    assignments_path = tables_dir / "alternative_clustering_assignments.csv"
    fig_path = figs_dir / "ijoc_alternative_clustering_sensitivity.png"

    summary.to_csv(summary_path, index=False)
    assignments.to_csv(assignments_path, index=False)
    _plot_alternative_clustering_heatmap(summary, fig_path, cfg)

    return {
        "summary_table": summary_path,
        "assignments_table": assignments_path,
        "figure": fig_path,
    }
