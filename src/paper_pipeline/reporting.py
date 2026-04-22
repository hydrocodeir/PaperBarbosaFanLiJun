from __future__ import annotations

from pathlib import Path

import pandas as pd


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
    rec = int(cfg["quantile_regression"]["min_years_recommended_for_publication"])
    min_years = int(cfg["quantile_regression"]["min_years_required_to_run"])
    short_run_count = int(features["insufficient_years_for_qr"].sum()) if "insufficient_years_for_qr" in features.columns else 0
    short_pub_count = int(features["publication_warning_short_record"].sum()) if "publication_warning_short_record" in features.columns else 0

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
    if n_years_total < rec:
        lines += [
            f"- WARNING: record length is below the publication recommendation ({rec}+ years).",
            "- Output is suitable for method validation, not final climatological inference.",
        ]

    lines += ["", "## Highest Delta1 stations by index"]
    for idx_cfg in cfg["indices"]:
        idx_name = idx_cfg["name"]
        sdf = features.loc[features["index_name"] == idx_name, ["station_name", "Delta1"]].dropna().sort_values("Delta1", ascending=False).head(5)
        lines.append(f"### {idx_cfg['title']}")
        lines += [f"- {r['station_name']}: Delta1 = {r['Delta1']:.3f}" for _, r in sdf.iterrows()] if not sdf.empty else ["- No results"]
        lines.append("")

    sensitivity_taus = [float(x) for x in cfg["quantile_regression"].get("sensitivity_check_quantiles", [0.05, 0.95])]
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
        lines.append("- Comparison between the configured clustering setup and a reduced-feature rerun.")
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
            raw_sig = int((pd.to_numeric(sdf["analytic_p"], errors="coerce") < 0.05).sum())
            fdr_sig = int(pd.to_numeric(sdf["fdr_reject"], errors="coerce").fillna(0).sum())
            lines.append(f"- τ = {tau:.2f}: raw significant station-results = {raw_sig}, FDR-retained = {fdr_sig}")
        lines.append("")
    if isinstance(moran_df, pd.DataFrame) and not moran_df.empty:
        sig_moran = moran_df.loc[pd.to_numeric(moran_df["moran_p_perm"], errors="coerce") < 0.05]
        lines.append(f"- Moran's I tests with p < 0.05: {len(sig_moran)}/{len(moran_df)} station-slope fields")
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
        delta_drivers = driver_df.loc[driver_df["metric"] == "Delta1"].copy()
        if not delta_drivers.empty:
            top = delta_drivers.reindex(delta_drivers["std_beta"].abs().sort_values(ascending=False).index).head(6)
            for _, row in top.iterrows():
                lines.append(
                    f"- {row['index_name']} | Delta1 ~ {row['predictor']}: "
                    f"std beta = {float(row['std_beta']):+.3f}, p = {float(row['std_beta_pvalue']):.3f}, "
                    f"rho = {float(row['spearman_rho']):+.3f}"
                )
            lines.append("")

    regional_df = advanced_results.get("regional_cluster_composites")
    if isinstance(regional_df, pd.DataFrame) and not regional_df.empty:
        lines.append("## Regional Composites")
        delta_region = regional_df.loc[regional_df["metric"] == "Delta1"].copy()
        if not delta_region.empty:
            top = delta_region.sort_values("median", ascending=False).head(6)
            for _, row in top.iterrows():
                lines.append(
                    f"- {row['index_name']} | cluster {int(row['cluster'])}: "
                    f"median Delta1 = {float(row['median']):+.3f} (n={int(row['n_stations'])})"
                )
            lines.append("")

    (outdir / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
