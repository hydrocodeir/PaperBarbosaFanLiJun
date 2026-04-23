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
from .config_utils import (
    boot_ci_high_col,
    boot_ci_low_col,
    boot_mean_col,
    boot_sd_col,
    ci_high_col,
    ci_low_col,
    get_delta_definitions,
    get_focus_quantiles,
    get_primary_delta,
    get_sensitivity_quantiles,
    slope_col,
    validate_analysis_config,
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
from .progress_utils import ProgressTracker
from .quantile import add_sensitivity_check_columns, run_station_qr
from .reporting import generate_report
from .year_config import filter_to_analysis_years, format_year_range_label, get_effective_year_range


def _read_cached_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def run_pipeline(config_path: str = "config.yaml", start_phase: int = 1) -> Path:
    cfg = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    validate_analysis_config(cfg)
    if start_phase not in {1, 8, 9, 10}:
        raise ValueError("start_phase must be one of: 1, 8, 9, 10")
    outdir = Path(cfg["paths"]["output_dir"])
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures"
    status_path = outdir / "run_status.txt"
    status_summary_path = outdir / "run_status_summary.txt"
    status_detail_path = outdir / "run_status_detail.txt"
    for p in (outdir, tables_dir, figs_dir):
        p.mkdir(parents=True, exist_ok=True)

    def _write_status_line(path: Path, line: str) -> None:
        with path.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")

    def log_status(message: str, *, level: str = "summary"):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{timestamp}] {message}"
        print(line, flush=True)
        _write_status_line(status_path, line)
        if level == "summary":
            _write_status_line(status_summary_path, line)
        elif level == "detail":
            _write_status_line(status_detail_path, line)
        else:
            _write_status_line(status_summary_path, line)
            _write_status_line(status_detail_path, line)

    def log_summary(message: str) -> None:
        log_status(message, level="summary")

    def log_detail(message: str) -> None:
        log_status(message, level="detail")

    phase_names = [
        "Load input tables",
        "Data-quality diagnostics",
        "Annual index construction",
        "Homogeneity sensitivity",
        "Quantile regression and bootstrap",
        "Clustering and feature tables",
        "Bootstrap-depth sensitivity",
        "Figure rendering",
        "Advanced publication analyses",
        "Final reporting",
    ]
    pipeline_tracker = ProgressTracker("Pipeline", len(phase_names), log_summary)

    def log_phase_start(phase_no: int, detail: str | None = None) -> None:
        phase_label = phase_names[phase_no - 1]
        extra = f"[PHASE {phase_no}/{len(phase_names)}] {phase_label}"
        if detail is not None:
            extra = f"{extra} | {detail}"
        pipeline_tracker.emit(phase_no, detail=extra)

    data = pd.read_csv(cfg["paths"]["data_csv"])
    data = filter_to_analysis_years(data, cfg)
    stations = pd.read_csv(cfg["paths"]["station_csv"])
    analysis_year_range = get_effective_year_range(cfg, data)
    log_phase_start(1, detail=f"data_rows={len(data)} | stations={stations.shape[0]}")
    log_summary(
        f"Loaded input tables: data_rows={len(data)}, stations={stations.shape[0]}, "
        f"analysis_years={format_year_range_label(analysis_year_range)}"
    )

    if start_phase > 1:
        log_summary(f"Resume mode enabled: starting from phase {start_phase}. Loading cached intermediate tables.")

        annual = _read_cached_csv(tables_dir / "annual_extreme_indices.csv")
        if annual.empty:
            raise FileNotFoundError(f"Required cached file not found or empty: {tables_dir / 'annual_extreme_indices.csv'}")

        qr_summary = _read_cached_csv(tables_dir / "qr_focus_slopes_and_bootstrap_summary.csv")
        feature_table = _read_cached_csv(tables_dir / "clustering_feature_table.csv")
        cluster_df = _read_cached_csv(tables_dir / "cluster_assignments.csv")
        cluster_robustness_df = _read_cached_csv(tables_dir / "cluster_robustness_summary.csv")
        boot_long = _read_cached_csv(tables_dir / "bootstrap_distributions_long.csv")

        if start_phase <= 9 and qr_summary.empty:
            raise FileNotFoundError(f"Required cached file not found or empty: {tables_dir / 'qr_focus_slopes_and_bootstrap_summary.csv'}")
        if start_phase <= 10 and feature_table.empty:
            raise FileNotFoundError(f"Required cached file not found or empty: {tables_dir / 'clustering_feature_table.csv'}")

        artifacts = {}
        if start_phase == 8 and not feature_table.empty:
            _, artifacts = run_clustering(feature_table, cfg)

        advanced_results = {}

        if start_phase == 8:
            log_phase_start(8)
            log_summary("Rendering figures...")
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
            log_summary("Finished rendering figures.")

        if start_phase <= 9:
            log_phase_start(9)
            log_summary("Running advanced publication analyses...")
            advanced_results.update(run_spatial_inference(qr_summary, stations, cfg, outdir, progress_callback=log_detail))
            advanced_results.update(run_method_sensitivity(data, annual, qr_summary, stations, cfg, outdir, progress_callback=log_detail))
            advanced_results.update(run_driver_analysis(feature_table, stations, cfg, outdir, progress_callback=log_detail))
            advanced_results.update(run_regionalization_analysis(feature_table, stations, cfg, outdir, progress_callback=log_detail))
            log_summary("Finished advanced publication analyses.")

            log_summary("Rendering robustness synthesis figure from completed sensitivity outputs...")
            plot_ijoc_robustness_synthesis(figs_dir, cfg)
            log_summary("Saved robustness synthesis figure.")

        log_phase_start(10)
        generate_report(
            annual,
            feature_table,
            cluster_df,
            cfg,
            outdir,
            cluster_robustness_df=cluster_robustness_df,
            advanced_results=advanced_results,
        )
        log_summary("Generated final report.")
        return outdir

    log_phase_start(2)
    log_summary("Running data-quality and homogeneity diagnostics...")
    dq_results = run_data_quality_assessment(data, cfg, outdir)
    log_summary("Saved data-quality and homogeneity diagnostics.")

    log_phase_start(3)
    log_summary("Building annual extreme indices...")
    _, annual = create_extreme_indices(data, cfg, progress_callback=log_detail)
    annual.to_csv(tables_dir / "annual_extreme_indices.csv", index=False)
    log_summary(f"Saved annual indices table with {len(annual)} rows.")

    log_phase_start(4)
    log_summary("Running sensitivity check excluding homogeneity-flagged stations...")
    run_homogeneity_exclusion_sensitivity(annual, dq_results["homogeneity"], cfg, outdir)
    log_summary("Saved homogeneity-exclusion sensitivity diagnostics.")

    metadata = {
        "project_name": cfg["project"]["name"],
        "year_range": [int(annual["year"].min()), int(annual["year"].max())],
        "n_years": int(annual["year"].nunique()),
        "n_stations": int(annual["station_id"].nunique()),
        "bootstrap_reps": int(cfg["bootstrap"]["n_reps"]),
        "focus_quantiles": cfg["quantile_regression"]["focus_quantiles"],
    }
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    log_summary("Wrote run metadata.")

    log_phase_start(5)
    log_summary("Starting quantile regression and bootstrap stage. This can take a while with current settings.")
    qr_all, qr_summary, boot_long = run_station_qr(annual, cfg, progress_callback=log_detail)
    qr_summary = add_sensitivity_check_columns(qr_summary, cfg)
    log_summary("Quantile regression stage completed.")

    log_phase_start(6)
    feature_table = build_feature_table(qr_summary, cfg)
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
    log_summary("Saved quantile, bootstrap, and clustering tables.")

    log_summary("Running alternative-clustering sensitivity analysis...")
    run_alternative_clustering_sensitivity(feature_table, cfg, outdir)
    log_summary("Saved alternative-clustering sensitivity diagnostics.")

    focus_cols = ["index_name", "station_id", "station_name", "n_years", "insufficient_years_for_qr"]
    for tau in get_focus_quantiles(cfg):
        focus_cols.extend([slope_col(tau), ci_low_col(tau), ci_high_col(tau)])
    for delta_name in get_delta_definitions(cfg):
        focus_cols.append(delta_name)
    for tau in get_focus_quantiles(cfg):
        suffix = f"{tau:0.2f}"
        focus_cols.extend([
            boot_mean_col(suffix),
            boot_sd_col(suffix),
            boot_ci_low_col(suffix),
            boot_ci_high_col(suffix),
        ])
    for tau in get_sensitivity_quantiles(cfg):
        suffix = f"{tau:0.2f}"
        focus_cols.extend([
            f"analytic_sig_{suffix}",
            f"bootstrap_sig_{suffix}",
            f"sig_agree_{suffix}",
            f"sensitivity_status_{suffix}",
        ])
    primary_delta = get_primary_delta(cfg)
    focus_cols.extend([
        boot_mean_col(primary_delta),
        boot_sd_col(primary_delta),
        boot_ci_low_col(primary_delta),
        boot_ci_high_col(primary_delta),
        "cluster",
        "cluster_reduced_features",
    ])
    feature_table[[c for c in focus_cols if c in feature_table.columns]].to_csv(tables_dir / "publication_summary_table.csv", index=False)
    log_summary("Saved publication summary table.")

    log_phase_start(7)
    log_summary("Running limited higher-replication bootstrap-depth sensitivity...")
    run_bootstrap_depth_sensitivity(annual, cfg, outdir, progress_callback=log_detail)
    log_summary("Saved bootstrap-depth sensitivity diagnostics.")

    log_phase_start(8)
    log_summary("Rendering figures...")
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
    log_summary("Finished rendering figures.")

    log_phase_start(9)
    log_summary("Running advanced publication analyses...")
    advanced_results = {}
    advanced_results.update(run_spatial_inference(qr_summary, stations, cfg, outdir, progress_callback=log_detail))
    advanced_results.update(run_method_sensitivity(data, annual, qr_summary, stations, cfg, outdir, progress_callback=log_detail))
    advanced_results.update(run_driver_analysis(feature_table, stations, cfg, outdir, progress_callback=log_detail))
    advanced_results.update(run_regionalization_analysis(feature_table, stations, cfg, outdir, progress_callback=log_detail))
    log_summary("Finished advanced publication analyses.")

    log_summary("Rendering robustness synthesis figure from completed sensitivity outputs...")
    plot_ijoc_robustness_synthesis(figs_dir, cfg)
    log_summary("Saved robustness synthesis figure.")

    log_phase_start(10)
    generate_report(
        annual,
        feature_table,
        cluster_df,
        cfg,
        outdir,
        cluster_robustness_df=cluster_robustness_df,
        advanced_results=advanced_results,
    )
    log_summary("Generated final report.")
    return outdir
