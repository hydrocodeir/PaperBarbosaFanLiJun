from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import yaml

from .advanced_analysis import (
    run_driver_analysis,
    run_method_sensitivity,
    run_regionalization_analysis,
    run_spatial_inference,
)
from .bootstrap_depth_sensitivity import run_bootstrap_depth_sensitivity
from .clustering import build_feature_table, compare_clusterings, run_clustering, screen_clustering_features
from .clustering_sensitivity import run_alternative_clustering_sensitivity
from .data_quality import run_data_quality_assessment
from .homogeneity_sensitivity import run_homogeneity_exclusion_sensitivity
from .indices import create_extreme_indices
from .plotting import (
    plot_data_coverage,
    plot_delta_uncertainty,
    plot_dendrograms,
    plot_ijoc_station_comparisons,
    plot_ijoc_main_delta_maps,
    plot_ijoc_regional_quantile_panels,
    plot_ijoc_robustness_synthesis,
    plot_ijoc_split_period_comparison,
    plot_ijoc_study_area,
    plot_maps,
    plot_paper1_quantile_dendrograms,
    plot_paper2_figure3_maps,
    plot_station_paper1_figure4,
    plot_region_quantile_slopes,
    plot_station_paper2_figures,
    plot_station_heatmap,
)
from .quantile import add_sensitivity_check_columns, run_station_qr
from .reporting import generate_report
from .year_config import filter_to_analysis_years, format_year_range_label, get_effective_year_range


def run_pipeline(config_path: str = "config.yaml") -> Path:
    cfg = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    outdir = Path(cfg["paths"]["output_dir"])
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures"
    status_path = outdir / "run_status.txt"
    for p in (outdir, tables_dir, figs_dir):
        p.mkdir(parents=True, exist_ok=True)

    def log_status(message: str):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{timestamp}] {message}"
        print(line, flush=True)
        with status_path.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")

    data = pd.read_csv(cfg["paths"]["data_csv"])
    data = filter_to_analysis_years(data, cfg)
    stations = pd.read_csv(cfg["paths"]["station_csv"])
    analysis_year_range = get_effective_year_range(cfg, data)
    log_status(
        f"Loaded input tables: data_rows={len(data)}, stations={stations.shape[0]}, "
        f"analysis_years={format_year_range_label(analysis_year_range)}"
    )

    log_status("Running data-quality and homogeneity diagnostics...")
    dq_results = run_data_quality_assessment(data, cfg, outdir)
    log_status("Saved data-quality and homogeneity diagnostics.")

    log_status("Building annual extreme indices...")
    _, annual = create_extreme_indices(data, cfg)
    annual.to_csv(tables_dir / "annual_extreme_indices.csv", index=False)
    log_status(f"Saved annual indices table with {len(annual)} rows.")

    log_status("Running sensitivity check excluding homogeneity-flagged stations...")
    run_homogeneity_exclusion_sensitivity(annual, dq_results["homogeneity"], cfg, outdir)
    log_status("Saved homogeneity-exclusion sensitivity diagnostics.")

    metadata = {
        "project_name": cfg["project"]["name"],
        "year_range": [int(annual["year"].min()), int(annual["year"].max())],
        "n_years": int(annual["year"].nunique()),
        "n_stations": int(annual["station_id"].nunique()),
        "bootstrap_reps": int(cfg["bootstrap"]["n_reps"]),
        "focus_quantiles": cfg["quantile_regression"]["focus_quantiles"],
    }
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    log_status("Wrote run metadata.")

    log_status("Starting quantile regression and bootstrap stage. This can take a while with current settings.")
    qr_all, qr_summary, boot_long = run_station_qr(annual, cfg, progress_callback=log_status)
    qr_summary = add_sensitivity_check_columns(qr_summary, cfg)
    log_status("Quantile regression stage completed.")

    feature_table = build_feature_table(qr_summary)
    baseline_screening_df = screen_clustering_features(
        feature_table,
        cfg,
        feature_set_label="baseline",
    )
    cluster_df, artifacts = run_clustering(feature_table, cfg)
    reduced_cluster_df = pd.DataFrame()
    cluster_robustness_df = pd.DataFrame()
    sensitivity_screening_df = pd.DataFrame()
    if not cluster_df.empty:
        feature_table = feature_table.merge(cluster_df, on=["index_name", "station_id", "station_name"], how="left")
    robustness_cfg = cfg["clustering"].get("robustness_check", {})
    if robustness_cfg.get("enabled", False):
        sensitivity_screening_df = screen_clustering_features(
            feature_table,
            cfg,
            feature_cols=robustness_cfg.get("reduced_features", cfg["clustering"]["simple_features"]),
            feature_set_label="sensitivity_rerun",
        )
        reduced_cluster_df, _ = run_clustering(
            feature_table,
            cfg,
            feature_cols=robustness_cfg.get("reduced_features", cfg["clustering"]["simple_features"]),
            label_col="cluster_reduced_features",
        )
        if not reduced_cluster_df.empty:
            feature_table = feature_table.merge(
                reduced_cluster_df,
                on=["index_name", "station_id", "station_name"],
                how="left",
            )
            cluster_robustness_df = compare_clusterings(cluster_df, reduced_cluster_df)

    qr_all.to_csv(tables_dir / "qr_all_quantiles_long.csv", index=False)
    qr_summary.to_csv(tables_dir / "qr_focus_slopes_and_bootstrap_summary.csv", index=False)
    feature_table.to_csv(tables_dir / "clustering_feature_table.csv", index=False)
    if not baseline_screening_df.empty or not sensitivity_screening_df.empty:
        pd.concat([df for df in [baseline_screening_df, sensitivity_screening_df] if not df.empty], ignore_index=True).to_csv(
            tables_dir / "clustering_feature_screening_summary.csv",
            index=False,
        )
    cluster_df.to_csv(tables_dir / "cluster_assignments.csv", index=False)
    if not reduced_cluster_df.empty:
        reduced_cluster_df.to_csv(tables_dir / "cluster_assignments_reduced_features.csv", index=False)
    if not cluster_robustness_df.empty:
        cluster_robustness_df.to_csv(tables_dir / "cluster_robustness_summary.csv", index=False)
    if cfg["bootstrap"]["save_long_table"] and not boot_long.empty:
        boot_long.to_csv(tables_dir / "bootstrap_distributions_long.csv", index=False)
    log_status("Saved quantile, bootstrap, and clustering tables.")

    log_status("Running alternative-clustering sensitivity analysis...")
    run_alternative_clustering_sensitivity(feature_table, cfg, outdir)
    log_status("Saved alternative-clustering sensitivity diagnostics.")

    focus_cols = [
        "index_name", "station_id", "station_name", "n_years",
        "insufficient_years_for_qr",
        "slope_0.05", "ci_low_0.05", "ci_high_0.05",
        "slope_0.50", "ci_low_0.50", "ci_high_0.50",
        "slope_0.95", "ci_low_0.95", "ci_high_0.95",
        "Delta1", "Delta2", "Delta3",
        "boot_mean_0.05", "boot_sd_0.05", "boot_ci_low_0.05", "boot_ci_high_0.05",
        "boot_mean_0.50", "boot_sd_0.50", "boot_ci_low_0.50", "boot_ci_high_0.50",
        "boot_mean_0.95", "boot_sd_0.95", "boot_ci_low_0.95", "boot_ci_high_0.95",
        "analytic_sig_0.05", "bootstrap_sig_0.05", "sig_agree_0.05", "sensitivity_status_0.05",
        "analytic_sig_0.95", "bootstrap_sig_0.95", "sig_agree_0.95", "sensitivity_status_0.95",
        "boot_mean_Delta1", "boot_sd_Delta1", "boot_ci_low_Delta1", "boot_ci_high_Delta1", "cluster",
        "cluster_reduced_features",
    ]
    feature_table[[c for c in focus_cols if c in feature_table.columns]].to_csv(tables_dir / "publication_summary_table.csv", index=False)
    log_status("Saved publication summary table.")

    log_status("Running limited higher-replication bootstrap-depth sensitivity...")
    run_bootstrap_depth_sensitivity(annual, cfg, outdir)
    log_status("Saved bootstrap-depth sensitivity diagnostics.")

    log_status("Rendering figures...")
    plot_data_coverage(annual, figs_dir, cfg)
    plot_region_quantile_slopes(annual, figs_dir, cfg)
    plot_ijoc_regional_quantile_panels(annual, figs_dir, cfg)
    plot_ijoc_study_area(stations, figs_dir, cfg)
    plot_station_heatmap(feature_table, figs_dir, cfg)
    plot_delta_uncertainty(feature_table, figs_dir, cfg)
    plot_dendrograms(feature_table, artifacts, figs_dir, cfg)
    plot_maps(feature_table, stations, cluster_df, figs_dir, cfg)
    plot_ijoc_main_delta_maps(feature_table, stations, figs_dir, cfg)
    plot_ijoc_split_period_comparison(annual, figs_dir, cfg)
    plot_ijoc_station_comparisons(annual, feature_table, figs_dir, cfg)
    plot_station_paper2_figures(annual, figs_dir, cfg)
    plot_paper2_figure3_maps(annual, stations, figs_dir, cfg)
    plot_station_paper1_figure4(boot_long, qr_summary, figs_dir, cfg)
    plot_paper1_quantile_dendrograms(qr_summary, figs_dir, cfg)
    log_status("Finished rendering figures.")

    log_status("Running advanced publication analyses...")
    advanced_results = {}
    advanced_results.update(run_spatial_inference(qr_summary, stations, cfg, outdir))
    advanced_results.update(run_method_sensitivity(data, annual, qr_summary, stations, cfg, outdir))
    advanced_results.update(run_driver_analysis(feature_table, stations, cfg, outdir))
    advanced_results.update(run_regionalization_analysis(feature_table, stations, cfg, outdir))
    log_status("Finished advanced publication analyses.")

    log_status("Rendering robustness synthesis figure from completed sensitivity outputs...")
    plot_ijoc_robustness_synthesis(figs_dir, cfg)
    log_status("Saved robustness synthesis figure.")

    generate_report(
        annual,
        feature_table,
        cluster_df,
        cfg,
        outdir,
        cluster_robustness_df=cluster_robustness_df,
        advanced_results=advanced_results,
    )
    log_status("Generated final report.")
    return outdir
