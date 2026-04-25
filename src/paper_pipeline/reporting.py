from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config_utils import get_primary_delta, get_sensitivity_quantiles


def generate_report(
    annual: pd.DataFrame,
    features: pd.DataFrame,
    cluster_df: pd.DataFrame,
    cfg: dict,
    outdir: Path,
    cluster_robustness_df: pd.DataFrame | None = None,
    advanced_results: dict | None = None,
):
    n_years_total = annual["year"].nunique()
    station_count = annual["station_id"].nunique()
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures"
    rec = int(cfg["quantile_regression"]["min_years_recommended_for_publication"])
    min_years = int(cfg["quantile_regression"]["min_years_required_to_run"])
    short_run_count = int(features["insufficient_years_for_qr"].sum()) if "insufficient_years_for_qr" in features.columns else 0
    short_pub_count = int(features["publication_warning_short_record"].sum()) if "publication_warning_short_record" in features.columns else 0
    primary_delta = get_primary_delta(cfg)

    lines = [
        "# Quantile Regression + Bootstrap + Clustering Report",
        "",
        "## Data audit",
        f"- Number of stations: **{station_count}**",
        f"- Number of years in uploaded data: **{n_years_total}**",
        f"- Year range: **{annual['year'].min()}-{annual['year'].max()}**",
        f"- Station-index series below QR minimum length ({min_years} years): **{short_run_count}**",
        f"- Station-index series below publication recommendation ({rec}+ years): **{short_pub_count}**",
    ]

    module_checks = [
        ("Data-quality and homogeneity diagnostics", tables_dir / "data_quality_homogeneity_overview.csv"),
        ("Homogeneity-exclusion sensitivity", tables_dir / "homogeneity_flag_exclusion_sensitivity.csv"),
        ("Alternative-clustering sensitivity", tables_dir / "alternative_clustering_sensitivity_summary.csv"),
        ("Bootstrap-depth sensitivity", tables_dir / "bootstrap_depth_sensitivity_summary.csv"),
        ("Advanced publication analyses", tables_dir / "station_significance_fdr.csv"),
        ("Climate-change signal analysis", tables_dir / "climate_fingerprint_component_scores.csv"),
        ("Köppen-Geiger climate-regime analysis", tables_dir / "climate_regime_fingerprint_summary.csv"),
    ]
    lines += ["", "## Pipeline-Integrated Robustness Modules"]
    for label, path in module_checks:
        status = "available" if path.exists() else "not found"
        lines.append(f"- {label}: **{status}**")
    if n_years_total < rec:
        lines += [
            f"- WARNING: record length is below the publication recommendation ({rec}+ years).",
            "- Output is suitable for method validation, not final climatological inference.",
        ]

    lines += ["", f"## Highest {primary_delta} stations by index"]
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        sdf = features.loc[features["index_name"] == idx_name, ["station_name", primary_delta]].dropna().sort_values(primary_delta, ascending=False).head(5)
        lines.append(f"### {idx_cfg['title']}")
        lines += [f"- {r['station_name']}: {primary_delta} = {r[primary_delta]:.3f}" for _, r in sdf.iterrows()] if not sdf.empty else ["- No results"]
        lines.append("")

    sensitivity_taus = get_sensitivity_quantiles(cfg)
    lines.append("## Sensitivity Checks")
    for tau in sensitivity_taus:
        suffix = f"{tau:0.2f}"
        status_col = f"sensitivity_status_{suffix}"
        if status_col not in features.columns:
            continue
        counts = features[status_col].fillna("insufficient").value_counts()
        lines.append(f"### Quantile τ = {tau:.2f}")
        for label in ["agree_significant", "agree_nonsignificant", "analytic_only", "bootstrap_only", "insufficient"]:
            if label in counts.index:
                lines.append(f"- {label}: {int(counts[label])}")
        lines.append("")

    if not cluster_df.empty:
        lines.append("## Cluster sizes")
        for idx_cfg in cfg["indices"]:
            idx_name = idx_cfg["name"]
            counts = cluster_df.loc[cluster_df["index_name"] == idx_name, "cluster"].value_counts().sort_index()
            lines.append(f"### {idx_cfg['title']}")
            lines += [f"- Cluster {int(k)}: {int(v)} stations" for k, v in counts.items()]
            lines.append("")

    if cluster_robustness_df is not None and not cluster_robustness_df.empty:
        lines.append("## Clustering Robustness")
        lines.append("- Comparison between the configured parsimonious clustering baseline and an expanded sensitivity rerun.")
        for idx_cfg in cfg["indices"]:
            idx_name = idx_cfg["name"]
            sdf = cluster_robustness_df.loc[cluster_robustness_df["index_name"] == idx_name]
            if sdf.empty:
                continue
            row = sdf.iloc[0]
            lines.append(f"### {idx_cfg['title']}")
            lines.append(f"- Stations compared: {int(row['n_stations_compared'])}")
            lines.append(f"- Adjusted Rand Index: {float(row['adjusted_rand_index']):.3f}")
            lines.append(f"- Same-label fraction: {float(row['same_label_fraction']):.3f}")
            lines.append("")

    advanced_results = advanced_results or {}

    fdr_df = advanced_results.get("station_significance_fdr")
    moran_df = advanced_results.get("spatial_autocorrelation_moran")
    if isinstance(fdr_df, pd.DataFrame) and not fdr_df.empty:
        lines.append("## Spatial Inference")
        for tau in sorted(pd.to_numeric(fdr_df["tau"], errors="coerce").dropna().unique().tolist()):
            sdf = fdr_df.loc[pd.to_numeric(fdr_df["tau"], errors="coerce") == tau].copy()
            raw_sig = int((pd.to_numeric(sdf["analytic_p"], errors="coerce") < float(cfg["advanced_analyses"]["spatial_inference"]["fdr_alpha"])).sum())
            fdr_sig = int(pd.to_numeric(sdf["fdr_reject"], errors="coerce").fillna(0).sum())
            lines.append(f"- τ = {tau:.2f}: raw significant station-results = {raw_sig}, FDR-retained = {fdr_sig}")
        lines.append("")
    if isinstance(moran_df, pd.DataFrame) and not moran_df.empty:
        alpha = float(cfg["advanced_analyses"]["spatial_inference"]["fdr_alpha"])
        sig_moran = moran_df.loc[pd.to_numeric(moran_df["moran_p_perm"], errors="coerce") < alpha]
        lines.append(f"- Moran's I tests with p < {alpha:.2f}: {len(sig_moran)}/{len(moran_df)} station-slope fields")
        lines.append("")

    ref_sens = advanced_results.get("reference_period_sensitivity_summary")
    boot_sens = advanced_results.get("bootstrap_method_sensitivity_summary")
    interp_sens = advanced_results.get("interpolation_method_sensitivity_summary")
    if isinstance(ref_sens, pd.DataFrame) and not ref_sens.empty:
        lines.append("## Method Sensitivity")
        lines.append("- Reference-period sensitivity:")
        best = ref_sens.sort_values("mean_abs_diff").head(4)
        for _, row in best.iterrows():
            lines.append(f"- {row['index_name']} | {row['metric']}: mean abs diff = {float(row['mean_abs_diff']):.3f}")
        lines.append("")
    if isinstance(boot_sens, pd.DataFrame) and not boot_sens.empty:
        lines.append("- Bootstrap-method sensitivity:")
        best = boot_sens.sort_values("mean_abs_diff").head(4)
        for _, row in best.iterrows():
            lines.append(f"- {row['index_name']} | {row['metric']}: mean abs diff = {float(row['mean_abs_diff']):.3f}")
        lines.append("")
    if isinstance(interp_sens, pd.DataFrame) and not interp_sens.empty:
        top_corr = interp_sens.sort_values("surface_correlation", ascending=False).head(4)
        lines.append("- Interpolation-method agreement (top correlations):")
        for _, row in top_corr.iterrows():
            lines.append(
                f"- {row['index_name']} | τ={float(row['tau']):.2f} | {row['method_left']} vs {row['method_right']}: "
                f"r = {float(row['surface_correlation']):.3f}, RMSE = {float(row['surface_rmse']):.3f}"
            )
        lines.append("")

    driver_df = advanced_results.get("driver_analysis_summary")
    if isinstance(driver_df, pd.DataFrame) and not driver_df.empty:
        lines.append("## Driver Analysis")
        delta_drivers = driver_df.loc[driver_df["metric"] == primary_delta].copy()
        if not delta_drivers.empty:
            top = delta_drivers.reindex(delta_drivers["std_beta"].abs().sort_values(ascending=False).index).head(6)
            for _, row in top.iterrows():
                lines.append(
                    f"- {row['index_name']} | {primary_delta} ~ {row['predictor']}: "
                    f"std beta = {float(row['std_beta']):+.3f}, p = {float(row['std_beta_pvalue']):.3f}, "
                    f"rho = {float(row['spearman_rho']):+.3f}"
                )
            lines.append("")

    regional_df = advanced_results.get("regional_cluster_composites")
    if isinstance(regional_df, pd.DataFrame) and not regional_df.empty:
        lines.append("## Regional Composites")
        delta_region = regional_df.loc[regional_df["metric"] == primary_delta].copy()
        if not delta_region.empty:
            top = delta_region.sort_values("median", ascending=False).head(6)
            for _, row in top.iterrows():
                lines.append(
                    f"- {row['index_name']} | cluster {int(row['cluster'])}: "
                    f"median {primary_delta} = {float(row['median']):+.3f} (n={int(row['n_stations'])})"
                )
            lines.append("")

    fingerprint_df = advanced_results.get("climate_fingerprint_component_scores")
    fixed_df = advanced_results.get("fixed_baseline_period_change_summary")
    warming_df = advanced_results.get("warming_link_network_quantile_response")
    temp_summary = advanced_results.get("regional_temperature_anomaly_summary")
    emergence_df = advanced_results.get("climate_signal_emergence_summary")
    network_fingerprint = advanced_results.get("climate_fingerprint_network_components")
    if isinstance(fingerprint_df, pd.DataFrame) and not fingerprint_df.empty:
        lines.append("## Climate-Change Signal Analysis")
        overall = float(fingerprint_df["overall_fingerprint_score"].iloc[0])
        lines.append(f"- Integrated climate-fingerprint score: **{overall:.3f}**")
        if isinstance(network_fingerprint, pd.DataFrame) and not network_fingerprint.empty:
            obs = float(network_fingerprint["observed_fingerprint_score"].iloc[0])
            p_val = float(network_fingerprint["permutation_p_value"].iloc[0])
            lines.append(f"- Network quantile-fingerprint permutation test: observed score = {obs:.3f}, p = {p_val:.3f}")
        if isinstance(temp_summary, pd.DataFrame) and not temp_summary.empty:
            row = temp_summary.iloc[0]
            lines.append(
                f"- Regional station-temperature anomaly trend: {float(row['temperature_trend_c_per_decade']):+.3f} "
                f"deg C per decade"
            )
        if isinstance(fixed_df, pd.DataFrame) and not fixed_df.empty:
            lines.append("- Fixed-baseline late-minus-early shifts:")
            for _, row in fixed_df.iterrows():
                lines.append(
                    f"- {row['index_name']}: mean shift = {float(row['mean_late_minus_baseline_days']):+.2f} days/year, "
                    f"direction-consistent stations = {int(row['direction_consistent_stations'])}/{int(row['n_stations'])}"
                )
        if isinstance(warming_df, pd.DataFrame) and not warming_df.empty:
            q90 = warming_df.loc[(warming_df["model"] == "QR") & (pd.to_numeric(warming_df["tau"], errors="coerce") == 0.90)]
            if not q90.empty:
                lines.append("- Network q0.90 response to regional warming:")
                for _, row in q90.iterrows():
                    lines.append(f"- {row['index_name']}: {float(row['slope_per_c']):+.2f} days/year per deg C")
        if isinstance(emergence_df, pd.DataFrame) and not emergence_df.empty:
            q50 = emergence_df.loc[emergence_df["metric"] == "0.50"]
            if not q50.empty:
                lines.append("- q0.50 signal emergence counts:")
                for _, row in q50.iterrows():
                    lines.append(f"- {row['index_name']}: {int(row['emerged_stations'])}/{int(row['n_stations'])}")
        lines.append("")

    regime_summary = advanced_results.get("koppen_geiger_regime_summary")
    regime_quantile = advanced_results.get("climate_regime_quantile_summary")
    regime_fingerprint = advanced_results.get("climate_regime_fingerprint_summary")
    regime_tests = advanced_results.get("climate_regime_difference_tests")
    if isinstance(regime_summary, pd.DataFrame) and not regime_summary.empty:
        lines.append("## Köppen-Geiger Climate-Regime Analysis")
        lines.append("- Station counts by present-day Köppen-Geiger analysis regime:")
        for _, row in regime_summary.iterrows():
            lines.append(
                f"- {row['climate_regime_label']}: n = {int(row['n_stations'])}, "
                f"mean elevation = {float(row['mean_elevation_m']):.0f} m, "
                f"mean confidence = {float(row['mean_kg_confidence_pct']):.1f}%"
            )
        if isinstance(regime_fingerprint, pd.DataFrame) and not regime_fingerprint.empty:
            lines.append("- Overall climate-regime fingerprint scores:")
            reg = (
                regime_fingerprint[["climate_regime_label", "overall_regime_fingerprint_score"]]
                .drop_duplicates()
                .sort_values("overall_regime_fingerprint_score", ascending=False)
            )
            for _, row in reg.iterrows():
                lines.append(f"- {row['climate_regime_label']}: {float(row['overall_regime_fingerprint_score']):.3f}")
        if isinstance(regime_quantile, pd.DataFrame) and not regime_quantile.empty:
            warm = regime_quantile.loc[regime_quantile["index_name"] == "warm_days"].copy()
            cool = regime_quantile.loc[regime_quantile["index_name"] == "cool_days"].copy()
            if not warm.empty:
                top = warm.sort_values("mean_slope_0.90", ascending=False).iloc[0]
                lines.append(
                    f"- Strongest warm-day q0.90 mean trend: {top['climate_regime_label']} "
                    f"({float(top['mean_slope_0.90']):+.2f} days/year per decade)"
                )
            if not cool.empty:
                low = cool.sort_values("mean_slope_0.90").iloc[0]
                lines.append(
                    f"- Strongest cool-day q0.90 contraction: {low['climate_regime_label']} "
                    f"({float(low['mean_slope_0.90']):+.2f} days/year per decade)"
                )
        if isinstance(regime_tests, pd.DataFrame) and not regime_tests.empty:
            sig = regime_tests.loc[pd.to_numeric(regime_tests["fdr_q_value"], errors="coerce") < 0.05]
            lines.append(f"- FDR-retained between-regime permutation contrasts: {len(sig)}/{len(regime_tests)}")
        lines.append("")

    (outdir / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
