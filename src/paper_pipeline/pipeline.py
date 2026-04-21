from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from .clustering import build_feature_table, run_clustering
from .indices import create_extreme_indices
from .plotting import (
    plot_data_coverage,
    plot_delta_uncertainty,
    plot_dendrograms,
    plot_maps,
    plot_paper1_quantile_dendrograms,
    plot_paper2_figure3_maps,
    plot_station_paper1_figure4,
    plot_region_quantile_slopes,
    plot_station_paper2_figures,
    plot_station_heatmap,
)
from .quantile import run_station_qr
from .reporting import generate_report


def run_pipeline(config_path: str = "config.yaml") -> Path:
    cfg = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    outdir = Path(cfg["paths"]["output_dir"])
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures"
    for p in (outdir, tables_dir, figs_dir):
        p.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(cfg["paths"]["data_csv"])
    stations = pd.read_csv(cfg["paths"]["station_csv"])

    _, annual = create_extreme_indices(data, cfg)
    qr_all, qr_summary, boot_long = run_station_qr(annual, cfg)

    feature_table = build_feature_table(qr_summary)
    cluster_df, artifacts = run_clustering(feature_table, cfg)
    if not cluster_df.empty:
        feature_table = feature_table.merge(cluster_df, on=["index_name", "station_id", "station_name"], how="left")

    annual.to_csv(tables_dir / "annual_extreme_indices.csv", index=False)
    qr_all.to_csv(tables_dir / "qr_all_quantiles_long.csv", index=False)
    qr_summary.to_csv(tables_dir / "qr_focus_slopes_and_bootstrap_summary.csv", index=False)
    feature_table.to_csv(tables_dir / "clustering_feature_table.csv", index=False)
    cluster_df.to_csv(tables_dir / "cluster_assignments.csv", index=False)
    if cfg["bootstrap"]["save_long_table"] and not boot_long.empty:
        boot_long.to_csv(tables_dir / "bootstrap_distributions_long.csv", index=False)

    focus_cols = [
        "index_name", "station_id", "station_name", "n_years",
        "slope_0.05", "slope_0.50", "slope_0.95", "Delta1", "Delta2", "Delta3",
        "boot_mean_0.05", "boot_sd_0.05", "boot_ci_low_0.05", "boot_ci_high_0.05",
        "boot_mean_0.50", "boot_sd_0.50", "boot_ci_low_0.50", "boot_ci_high_0.50",
        "boot_mean_0.95", "boot_sd_0.95", "boot_ci_low_0.95", "boot_ci_high_0.95",
        "boot_mean_Delta1", "boot_sd_Delta1", "boot_ci_low_Delta1", "boot_ci_high_Delta1", "cluster",
    ]
    feature_table[[c for c in focus_cols if c in feature_table.columns]].to_csv(tables_dir / "publication_summary_table.csv", index=False)

    plot_data_coverage(annual, figs_dir, cfg)
    plot_region_quantile_slopes(annual, figs_dir, cfg)
    plot_station_heatmap(feature_table, figs_dir, cfg)
    plot_delta_uncertainty(feature_table, figs_dir, cfg)
    plot_dendrograms(feature_table, artifacts, figs_dir, cfg)
    plot_maps(feature_table, stations, cluster_df, figs_dir, cfg)
    plot_station_paper2_figures(annual, figs_dir, cfg)
    plot_paper2_figure3_maps(annual, stations, figs_dir, cfg)
    plot_station_paper1_figure4(boot_long, qr_summary, figs_dir, cfg)
    plot_paper1_quantile_dendrograms(qr_summary, figs_dir, cfg)

    metadata = {
        "project_name": cfg["project"]["name"],
        "year_range": [int(annual["year"].min()), int(annual["year"].max())],
        "n_years": int(annual["year"].nunique()),
        "n_stations": int(annual["station_id"].nunique()),
        "bootstrap_reps": int(cfg["bootstrap"]["n_reps"]),
        "focus_quantiles": cfg["quantile_regression"]["focus_quantiles"],
    }
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    generate_report(annual, feature_table, cluster_df, cfg, outdir)
    return outdir
