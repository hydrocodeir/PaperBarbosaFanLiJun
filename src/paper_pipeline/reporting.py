from __future__ import annotations

from pathlib import Path

import pandas as pd


def generate_report(annual: pd.DataFrame, features: pd.DataFrame, cluster_df: pd.DataFrame, cfg: dict, outdir: Path):
    n_years_total = annual["year"].nunique()
    station_count = annual["station_id"].nunique()
    rec = int(cfg["quantile_regression"]["min_years_recommended_for_publication"])

    lines = [
        "# Quantile Regression + Bootstrap + Clustering Report",
        "",
        "## Data audit",
        f"- Number of stations: **{station_count}**",
        f"- Number of years in uploaded data: **{n_years_total}**",
        f"- Year range: **{annual['year'].min()}-{annual['year'].max()}**",
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

    if not cluster_df.empty:
        lines.append("## Cluster sizes")
        for idx_cfg in cfg["indices"]:
            idx_name = idx_cfg["name"]
            counts = cluster_df.loc[cluster_df["index_name"] == idx_name, "cluster"].value_counts().sort_index()
            lines.append(f"### {idx_cfg['title']}")
            lines += [f"- Cluster {int(k)}: {int(v)} stations" for k, v in counts.items()]
            lines.append("")

    (outdir / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
