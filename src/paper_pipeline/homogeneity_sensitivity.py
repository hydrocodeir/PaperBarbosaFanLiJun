from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config_utils import (
    compute_defined_deltas,
    get_plot_dpi,
    get_focus_quantiles,
    get_primary_delta,
    get_time_scale_years,
    metric_label,
    slope_col,
)
from .quantile import fit_ols_line, fit_quantile_line


def run_homogeneity_exclusion_sensitivity(
    annual: pd.DataFrame,
    homogeneity: pd.DataFrame,
    cfg: dict,
    outdir: Path,
) -> pd.DataFrame:
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    flagged = set(
        homogeneity.loc[
            homogeneity["any_detrended_homogeneity_flag_p_lt_alpha"].fillna(False),
            "station_name",
        ].astype(str)
    )
    keep_names = sorted(set(annual["station_name"].astype(str)) - flagged)

    focus_quantiles = get_focus_quantiles(cfg)
    max_iter = int(cfg["quantile_regression"]["max_iter"])
    time_scale_years = get_time_scale_years(cfg)
    primary_delta = get_primary_delta(cfg)

    records = []
    for subset_name, frame in {
        "all_stations": annual.copy(),
        "exclude_flagged": annual.loc[annual["station_name"].isin(keep_names)].copy(),
    }.items():
        mean_df = frame.groupby("year", as_index=False)[[idx["name"] for idx in cfg["indices"]]].mean()
        period_years = mean_df["year"].to_numpy(dtype=float)
        for idx_cfg in cfg["indices"]:
            vals = mean_df[idx_cfg["name"]].to_numpy(dtype=float)
            ols = fit_ols_line(period_years, vals, time_scale_years=time_scale_years)
            row = {
                "subset": subset_name,
                "index_name": idx_cfg["name"],
                "n_stations": int(frame["station_name"].nunique()),
                "ols_slope": float(ols["slope"]),
            }
            for tau in focus_quantiles:
                qr = fit_quantile_line(period_years, vals, tau=tau, time_scale_years=time_scale_years, max_iter=max_iter)
                row[slope_col(tau)] = float(qr["slope"])
            records.append(compute_defined_deltas(pd.DataFrame([row]), cfg).iloc[0].to_dict())

    summary = pd.DataFrame(records)
    wide = summary.pivot(index="index_name", columns="subset")
    flat = wide.copy()
    flat.columns = [f"{a}_{b}" for a, b in flat.columns]
    flat = flat.reset_index()
    metrics = ["ols_slope"] + [slope_col(tau) for tau in focus_quantiles] + [primary_delta]
    for metric in metrics:
        flat[f"{metric}_difference_exclude_minus_all"] = (
            flat[f"{metric}_exclude_flagged"] - flat[f"{metric}_all_stations"]
        )

    flat.to_csv(tables_dir / "homogeneity_flag_exclusion_sensitivity.csv", index=False)

    plot_cols = ["index_name"] + [f"{metric}_difference_exclude_minus_all" for metric in metrics]
    plot_df = flat[plot_cols].set_index("index_name")
    plot_df = plot_df.rename(columns={f"{metric}_difference_exclude_minus_all": ("OLS" if metric == "ols_slope" else metric_label(metric)) for metric in metrics})

    fig, ax = plt.subplots(figsize=(8.5, 4.6), dpi=get_plot_dpi(cfg))
    im = ax.imshow(plot_df.to_numpy(), cmap="coolwarm", aspect="auto")
    ax.set_xticks(range(plot_df.shape[1]))
    ax.set_xticklabels(plot_df.columns)
    ax.set_yticks(range(plot_df.shape[0]))
    ax.set_yticklabels([s.replace("_", " ") for s in plot_df.index])
    for i in range(plot_df.shape[0]):
        for j in range(plot_df.shape[1]):
            val = plot_df.iloc[i, j]
            ax.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=8.5, color="#1f2937")
    cbar = fig.colorbar(im, ax=ax, shrink=0.9)
    cbar.set_label("Slope difference after excluding flagged stations")
    ax.set_title("Sensitivity of regional summaries to exclusion of homogeneity-flagged stations")
    fig.tight_layout()
    fig.savefig(figs_dir / "ijoc_homogeneity_sensitivity.png", dpi=get_plot_dpi(cfg), bbox_inches="tight")
    plt.close(fig)

    return flat
