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

    (outdir / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
