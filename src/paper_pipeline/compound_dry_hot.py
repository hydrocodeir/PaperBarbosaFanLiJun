from __future__ import annotations

from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import anderson_ksamp, cramervonmises_2samp, kendalltau, ks_2samp, linregress

from .config_utils import get_plot_dpi
from .plotting import (
    _draw_boundary,
    _format_geo_axis,
    _get_boundary_path_from_cfg,
    _get_plot_extent,
    _load_boundary_geometry,
    apply_publication_theme,
)


def _cfg(cfg: dict) -> dict:
    return cfg.get("advanced_analyses", {}).get("compound_dry_hot", {})


def _resolve_data_col(cfg: dict, alias: str) -> str:
    data_cfg = cfg.get("data", {})
    key = f"{alias}_col"
    return str(data_cfg.get(key, alias))


def _threshold_col(threshold: int | float) -> str:
    return f"event_rp_ge_{float(threshold):g}"


def _period_label(year: int, baseline_years: tuple[int, int], comparison_years: tuple[int, int]) -> str:
    if baseline_years[0] <= int(year) <= baseline_years[1]:
        return f"{baseline_years[0]}-{baseline_years[1]}"
    if comparison_years[0] <= int(year) <= comparison_years[1]:
        return f"{comparison_years[0]}-{comparison_years[1]}"
    return "outside_comparison"


def _periods_from_cfg(cfg: dict, years: pd.Series) -> tuple[tuple[int, int], tuple[int, int]]:
    module_cfg = _cfg(cfg)
    baseline = module_cfg.get("baseline_years")
    comparison = module_cfg.get("comparison_years")
    if baseline and comparison:
        return (int(baseline[0]), int(baseline[1])), (int(comparison[0]), int(comparison[1]))

    unique_years = sorted(pd.to_numeric(years, errors="coerce").dropna().astype(int).unique().tolist())
    if not unique_years:
        return (np.nan, np.nan), (np.nan, np.nan)
    midpoint = len(unique_years) // 2
    return (unique_years[0], unique_years[midpoint - 1]), (unique_years[midpoint], unique_years[-1])


def _seasonal_definition_label(definition: dict) -> str:
    months = definition.get("months")
    if not months:
        return "annual"
    return "M" + "-".join(f"{int(m):02d}" for m in months)


def _build_definition_aggregates(data: pd.DataFrame, cfg: dict, definition: dict) -> pd.DataFrame:
    station_id_col = _resolve_data_col(cfg, "station_id")
    station_name_col = _resolve_data_col(cfg, "station_name")
    year_col = _resolve_data_col(cfg, "year")
    month_col = _resolve_data_col(cfg, "month")
    precip_col = _resolve_data_col(cfg, "precip")
    temp_col = _resolve_data_col(cfg, str(definition.get("temperature_col", "tmean")))

    local = data.copy()
    months = definition.get("months")
    if months:
        local = local.loc[pd.to_numeric(local[month_col], errors="coerce").isin([int(m) for m in months])].copy()

    coverage_min = float(definition.get("coverage_threshold_pct", _cfg(cfg).get("coverage_threshold_pct", 80.0))) / 100.0
    grouped = (
        local.groupby([station_id_col, station_name_col, year_col], dropna=False)
        .agg(
            n_days=(precip_col, "size"),
            precip_valid_days=(precip_col, "count"),
            temperature_valid_days=(temp_col, "count"),
            precip_total=(precip_col, "sum"),
            temperature_mean=(temp_col, "mean"),
        )
        .reset_index()
        .rename(
            columns={
                station_id_col: "station_id",
                station_name_col: "station_name",
                year_col: "year",
            }
        )
    )
    grouped["definition"] = str(definition["name"])
    grouped["definition_title"] = str(definition.get("title", definition["name"]))
    grouped["season_label"] = _seasonal_definition_label(definition)
    grouped["temperature_variable"] = str(definition.get("temperature_col", "tmean"))
    grouped["precip_coverage_pct"] = 100.0 * grouped["precip_valid_days"] / grouped["n_days"].replace(0, np.nan)
    grouped["temperature_coverage_pct"] = 100.0 * grouped["temperature_valid_days"] / grouped["n_days"].replace(0, np.nan)
    grouped["valid_for_compound"] = (
        (grouped["precip_valid_days"] / grouped["n_days"].replace(0, np.nan) >= coverage_min)
        & (grouped["temperature_valid_days"] / grouped["n_days"].replace(0, np.nan) >= coverage_min)
        & np.isfinite(pd.to_numeric(grouped["precip_total"], errors="coerce"))
        & np.isfinite(pd.to_numeric(grouped["temperature_mean"], errors="coerce"))
    )
    return grouped


def _add_empirical_return_periods(aggregates: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    module_cfg = _cfg(cfg)
    min_years = int(module_cfg.get("min_station_years", 25))
    thresholds = [float(t) for t in module_cfg.get("return_period_thresholds", [5, 10, 20])]
    rows = []

    for (definition, station_id), group in aggregates.loc[aggregates["valid_for_compound"]].groupby(["definition", "station_id"]):
        local = group.sort_values("year").copy()
        if len(local) < min_years:
            continue

        precip = pd.to_numeric(local["precip_total"], errors="coerce").to_numpy(dtype=float)
        temp = pd.to_numeric(local["temperature_mean"], errors="coerce").to_numpy(dtype=float)
        n = len(local)

        joint_counts = np.array([np.sum((precip <= p) & (temp >= t)) for p, t in zip(precip, temp)], dtype=float)
        dry_counts = np.array([np.sum(precip <= p) for p in precip], dtype=float)
        hot_counts = np.array([np.sum(temp >= t) for t in temp], dtype=float)

        local["n_years_for_empirical_copula"] = n
        local["joint_probability_and"] = joint_counts / n
        local["joint_return_period_years"] = 1.0 / local["joint_probability_and"].replace(0, np.nan)
        local["dry_probability"] = dry_counts / n
        local["hot_probability"] = hot_counts / n
        local["dry_return_period_years"] = 1.0 / local["dry_probability"].replace(0, np.nan)
        local["hot_return_period_years"] = 1.0 / local["hot_probability"].replace(0, np.nan)
        local["dominant_driver"] = np.select(
            [
                local["dry_return_period_years"] > local["hot_return_period_years"],
                local["hot_return_period_years"] > local["dry_return_period_years"],
            ],
            ["dry", "hot"],
            default="co-dominant",
        )
        for threshold in thresholds:
            local[_threshold_col(threshold)] = (local["joint_return_period_years"] >= threshold).astype(int)
        rows.append(local)

    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def _make_yearly_extent(station_year: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    thresholds = [float(t) for t in _cfg(cfg).get("return_period_thresholds", [5, 10, 20])]
    rows = []
    for (definition, title, year), group in station_year.groupby(["definition", "definition_title", "year"]):
        n_valid = int(group["station_id"].nunique())
        for threshold in thresholds:
            col = _threshold_col(threshold)
            affected = int(pd.to_numeric(group[col], errors="coerce").fillna(0).sum())
            rows.append(
                {
                    "definition": definition,
                    "definition_title": title,
                    "year": int(year),
                    "return_period_threshold_years": threshold,
                    "n_valid_stations": n_valid,
                    "affected_stations": affected,
                    "affected_fraction_pct": 100.0 * affected / n_valid if n_valid else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _trend_summary(yearly_extent: pd.DataFrame, cfg: dict, baseline_years: tuple[int, int], comparison_years: tuple[int, int]) -> pd.DataFrame:
    rows = []
    for (definition, title, threshold), group in yearly_extent.groupby(["definition", "definition_title", "return_period_threshold_years"]):
        local = group.sort_values("year").copy()
        years = pd.to_numeric(local["year"], errors="coerce").to_numpy(dtype=float)
        counts = pd.to_numeric(local["affected_stations"], errors="coerce").to_numpy(dtype=float)
        fractions = pd.to_numeric(local["affected_fraction_pct"], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(years) & np.isfinite(counts)
        if mask.sum() < 4:
            continue
        ktau, kp = kendalltau(years[mask], counts[mask])
        lr = linregress(years[mask], counts[mask])
        frac_lr = linregress(years[mask], fractions[mask])
        early = local.loc[local["year"].between(*baseline_years)]
        late = local.loc[local["year"].between(*comparison_years)]
        rows.append(
            {
                "definition": definition,
                "definition_title": title,
                "return_period_threshold_years": threshold,
                "n_years": int(mask.sum()),
                "kendall_tau": float(ktau) if np.isfinite(ktau) else np.nan,
                "kendall_p_value": float(kp) if np.isfinite(kp) else np.nan,
                "linear_slope_stations_per_decade": float(lr.slope * 10.0),
                "linear_p_value": float(lr.pvalue),
                "linear_slope_fraction_pct_per_decade": float(frac_lr.slope * 10.0),
                "baseline_mean_affected_stations": float(early["affected_stations"].mean()),
                "comparison_mean_affected_stations": float(late["affected_stations"].mean()),
                "comparison_minus_baseline_stations": float(late["affected_stations"].mean() - early["affected_stations"].mean()),
                "baseline_mean_affected_fraction_pct": float(early["affected_fraction_pct"].mean()),
                "comparison_mean_affected_fraction_pct": float(late["affected_fraction_pct"].mean()),
                "comparison_minus_baseline_fraction_pct": float(late["affected_fraction_pct"].mean() - early["affected_fraction_pct"].mean()),
                "max_affected_year": int(local.loc[local["affected_stations"].idxmax(), "year"]),
                "max_affected_stations": int(local["affected_stations"].max()),
            }
        )
    return pd.DataFrame(rows)


def _distribution_shift_tests(yearly_extent: pd.DataFrame, baseline_years: tuple[int, int], comparison_years: tuple[int, int]) -> pd.DataFrame:
    rows = []
    for (definition, title, threshold), group in yearly_extent.groupby(["definition", "definition_title", "return_period_threshold_years"]):
        early = pd.to_numeric(group.loc[group["year"].between(*baseline_years), "affected_fraction_pct"], errors="coerce").dropna()
        late = pd.to_numeric(group.loc[group["year"].between(*comparison_years), "affected_fraction_pct"], errors="coerce").dropna()
        if len(early) < 4 or len(late) < 4:
            continue
        ks = ks_2samp(early, late, alternative="two-sided", mode="auto")
        cvm = cramervonmises_2samp(early, late)
        try:
            samples = [early.to_numpy(dtype=float), late.to_numpy(dtype=float)]
            try:
                ad = anderson_ksamp(samples, variant="midrank")
            except TypeError:
                ad = anderson_ksamp(samples)
            ad_stat = float(ad.statistic)
            ad_p = float(ad.pvalue)
        except Exception:
            ad_stat = np.nan
            ad_p = np.nan
        rows.append(
            {
                "definition": definition,
                "definition_title": title,
                "return_period_threshold_years": threshold,
                "baseline_period": f"{baseline_years[0]}-{baseline_years[1]}",
                "comparison_period": f"{comparison_years[0]}-{comparison_years[1]}",
                "ks_statistic": float(ks.statistic),
                "ks_p_value": float(ks.pvalue),
                "cramer_von_mises_statistic": float(cvm.statistic),
                "cramer_von_mises_p_value": float(cvm.pvalue),
                "anderson_darling_statistic": ad_stat,
                "anderson_darling_p_value": ad_p,
                "baseline_mean_fraction_pct": float(early.mean()),
                "comparison_mean_fraction_pct": float(late.mean()),
            }
        )
    return pd.DataFrame(rows)


def _station_frequency(station_year: pd.DataFrame, stations: pd.DataFrame, cfg: dict, baseline_years: tuple[int, int], comparison_years: tuple[int, int]) -> pd.DataFrame:
    thresholds = [float(t) for t in _cfg(cfg).get("return_period_thresholds", [5, 10, 20])]
    rows = []
    for (definition, title, station_id, station_name), group in station_year.groupby(["definition", "definition_title", "station_id", "station_name"]):
        for threshold in thresholds:
            col = _threshold_col(threshold)
            events = pd.to_numeric(group[col], errors="coerce").fillna(0)
            event_years = group.loc[events.astype(bool), "year"].astype(int).tolist()
            early = group.loc[group["year"].between(*baseline_years)]
            late = group.loc[group["year"].between(*comparison_years)]
            early_freq = 100.0 * pd.to_numeric(early[col], errors="coerce").fillna(0).mean() if not early.empty else np.nan
            late_freq = 100.0 * pd.to_numeric(late[col], errors="coerce").fillna(0).mean() if not late.empty else np.nan
            driver_counts = group.loc[events.astype(bool), "dominant_driver"].value_counts()
            rows.append(
                {
                    "definition": definition,
                    "definition_title": title,
                    "station_id": station_id,
                    "station_name": station_name,
                    "return_period_threshold_years": threshold,
                    "n_valid_years": int(group["year"].nunique()),
                    "event_years": int(events.sum()),
                    "event_frequency_pct": float(100.0 * events.mean()),
                    "baseline_event_frequency_pct": float(early_freq) if np.isfinite(early_freq) else np.nan,
                    "comparison_event_frequency_pct": float(late_freq) if np.isfinite(late_freq) else np.nan,
                    "comparison_minus_baseline_frequency_pct": float(late_freq - early_freq) if np.isfinite(early_freq) and np.isfinite(late_freq) else np.nan,
                    "first_event_year": min(event_years) if event_years else np.nan,
                    "last_event_year": max(event_years) if event_years else np.nan,
                    "dry_dominant_events": int(driver_counts.get("dry", 0)),
                    "hot_dominant_events": int(driver_counts.get("hot", 0)),
                    "co_dominant_events": int(driver_counts.get("co-dominant", 0)),
                }
            )
    freq = pd.DataFrame(rows)
    if freq.empty:
        return freq
    return freq.merge(stations, on="station_id", how="left", suffixes=("", "_station_meta"))


def _driver_summary(station_year: pd.DataFrame, cfg: dict, baseline_years: tuple[int, int], comparison_years: tuple[int, int]) -> pd.DataFrame:
    thresholds = [float(t) for t in _cfg(cfg).get("return_period_thresholds", [5, 10, 20])]
    rows = []
    local = station_year.copy()
    local["period"] = local["year"].apply(lambda y: _period_label(int(y), baseline_years, comparison_years))
    local = local.loc[local["period"] != "outside_comparison"].copy()
    for (definition, title), group in local.groupby(["definition", "definition_title"]):
        for threshold in thresholds:
            events = group.loc[pd.to_numeric(group[_threshold_col(threshold)], errors="coerce").fillna(0).astype(bool)]
            for period, period_df in events.groupby("period"):
                total = len(period_df)
                counts = period_df["dominant_driver"].value_counts()
                rows.append(
                    {
                        "definition": definition,
                        "definition_title": title,
                        "return_period_threshold_years": threshold,
                        "period": period,
                        "total_events": int(total),
                        "dry_dominant_events": int(counts.get("dry", 0)),
                        "hot_dominant_events": int(counts.get("hot", 0)),
                        "co_dominant_events": int(counts.get("co-dominant", 0)),
                        "hot_dominant_fraction_pct": float(100.0 * counts.get("hot", 0) / total) if total else np.nan,
                        "dry_dominant_fraction_pct": float(100.0 * counts.get("dry", 0) / total) if total else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def _knn_weights(coords: np.ndarray, k_neighbors: int) -> np.ndarray:
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    if n < 2:
        return np.zeros((n, n), dtype=float)
    dmat = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(dmat, np.inf)
    k = max(1, min(int(k_neighbors), n - 1))
    weights = np.zeros((n, n), dtype=float)
    for i in range(n):
        nn_idx = np.argsort(dmat[i])[:k]
        inv_d = 1.0 / np.maximum(dmat[i, nn_idx], 1e-12)
        weights[i, nn_idx] = inv_d / inv_d.sum()
    return weights


def _moran_i(values: np.ndarray, coords: np.ndarray, k_neighbors: int, permutations: int, rng: np.random.Generator) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    coords = np.asarray(coords, dtype=float)
    mask = np.isfinite(values) & np.all(np.isfinite(coords), axis=1)
    values = values[mask]
    coords = coords[mask]
    if len(values) < 4:
        return np.nan, np.nan
    weights = _knn_weights(coords, k_neighbors)
    w_sum = weights.sum()
    centered = values - values.mean()
    denom = np.sum(centered**2)
    if w_sum <= 0 or denom <= 0:
        return np.nan, np.nan
    obs = float((len(values) / w_sum) * ((weights * np.outer(centered, centered)).sum() / denom))
    permuted = []
    for _ in range(int(permutations)):
        shuffled = rng.permutation(centered)
        permuted.append(float((len(values) / w_sum) * ((weights * np.outer(shuffled, shuffled)).sum() / denom)))
    permuted_arr = np.asarray(permuted, dtype=float)
    p_value = float((1 + np.sum(np.abs(permuted_arr) >= abs(obs))) / (len(permuted_arr) + 1))
    return obs, p_value


def _moran_yearly(station_year: pd.DataFrame, stations: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    module_cfg = _cfg(cfg)
    thresholds = [float(t) for t in module_cfg.get("return_period_thresholds", [5, 10, 20])]
    k_neighbors = int(module_cfg.get("moran_k_neighbors", cfg.get("advanced_analyses", {}).get("spatial_inference", {}).get("moran_k_neighbors", 5)))
    permutations = int(module_cfg.get("moran_permutations", cfg.get("advanced_analyses", {}).get("spatial_inference", {}).get("moran_permutations", 499)))
    rng = np.random.default_rng(int(cfg.get("project", {}).get("random_seed", 42)) + 719)
    station_meta = stations[["station_id", "latitude", "longitude"]].copy()
    merged = station_year.merge(station_meta, on="station_id", how="left")
    rows = []
    for (definition, title, year), group in merged.groupby(["definition", "definition_title", "year"]):
        coords = group[["longitude", "latitude"]].to_numpy(dtype=float)
        for threshold in thresholds:
            values = pd.to_numeric(group[_threshold_col(threshold)], errors="coerce").to_numpy(dtype=float)
            moran, p_value = _moran_i(values, coords, k_neighbors, permutations, rng)
            rows.append(
                {
                    "definition": definition,
                    "definition_title": title,
                    "year": int(year),
                    "return_period_threshold_years": threshold,
                    "moran_i": moran,
                    "moran_p_perm": p_value,
                    "n_stations": int(np.isfinite(values).sum()),
                }
            )
    moran_df = pd.DataFrame(rows)
    trend_rows = []
    for (definition, title, threshold), group in moran_df.dropna(subset=["moran_i"]).groupby(["definition", "definition_title", "return_period_threshold_years"]):
        if len(group) < 4:
            continue
        lr = linregress(group["year"], group["moran_i"])
        ktau, kp = kendalltau(group["year"], group["moran_i"])
        trend_rows.append(
            {
                "definition": definition,
                "definition_title": title,
                "return_period_threshold_years": threshold,
                "n_years_with_defined_moran": int(len(group)),
                "mean_moran_i": float(group["moran_i"].mean()),
                "linear_slope_moran_per_decade": float(lr.slope * 10.0),
                "linear_p_value": float(lr.pvalue),
                "kendall_tau": float(ktau) if np.isfinite(ktau) else np.nan,
                "kendall_p_value": float(kp) if np.isfinite(kp) else np.nan,
            }
        )
    return moran_df, pd.DataFrame(trend_rows)


def _plot_extent_timeseries(yearly_extent: pd.DataFrame, outpath: Path, cfg: dict) -> None:
    if yearly_extent.empty:
        return
    apply_publication_theme()
    definitions = yearly_extent[["definition", "definition_title"]].drop_duplicates().to_dict("records")
    thresholds = sorted(yearly_extent["return_period_threshold_years"].dropna().unique().tolist())
    fig, axes = plt.subplots(len(definitions), 1, figsize=(9.6, 3.6 * len(definitions)), sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    colors = {threshold: color for threshold, color in zip(thresholds, ["#35618f", "#b35d2e", "#7b4ea3", "#2f7f5f"])}
    for ax, definition in zip(axes, definitions):
        sdf = yearly_extent.loc[yearly_extent["definition"] == definition["definition"]].sort_values("year")
        for threshold in thresholds:
            line_df = sdf.loc[sdf["return_period_threshold_years"] == threshold]
            ax.plot(
                line_df["year"],
                line_df["affected_fraction_pct"],
                marker="o",
                linewidth=1.6,
                markersize=3.8,
                color=colors.get(threshold),
                label=f"RP >= {threshold:g} yr",
            )
        ax.set_title(definition["definition_title"], loc="left")
        ax.set_ylabel("Affected stations (%)")
        ax.grid(True, axis="y", alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, ncols=len(thresholds), fontsize=8.5)
    axes[-1].set_xlabel("Year")
    fig.savefig(outpath, dpi=get_plot_dpi(cfg))
    plt.close(fig)


def _plot_station_frequency_maps(station_frequency: pd.DataFrame, outpath: Path, cfg: dict, value_col: str, title_suffix: str, cmap: str) -> None:
    if station_frequency.empty:
        return
    module_cfg = _cfg(cfg)
    selected_thresholds = [float(t) for t in module_cfg.get("map_thresholds", [10, 20])]
    definitions = station_frequency[["definition", "definition_title"]].drop_duplicates().to_dict("records")
    boundary_geom = _load_boundary_geometry(_get_boundary_path_from_cfg(cfg))
    apply_publication_theme()
    fig, axes = plt.subplots(len(definitions), len(selected_thresholds), figsize=(6.0 * len(selected_thresholds), 4.8 * len(definitions)), constrained_layout=True)
    axes = np.asarray(axes).reshape(len(definitions), len(selected_thresholds))
    vals = pd.to_numeric(station_frequency[value_col], errors="coerce").replace([np.inf, -np.inf], np.nan)
    if value_col == "comparison_minus_baseline_frequency_pct":
        vmax = max(5.0, float(np.nanmax(np.abs(vals))) if vals.notna().any() else 5.0)
        vmin = -vmax
    else:
        vmin = 0.0
        vmax = max(10.0, float(np.nanmax(vals)) if vals.notna().any() else 10.0)
    last_scatter = None
    for r, definition in enumerate(definitions):
        for c, threshold in enumerate(selected_thresholds):
            ax = axes[r, c]
            sdf = station_frequency.loc[
                (station_frequency["definition"] == definition["definition"])
                & (station_frequency["return_period_threshold_years"] == threshold)
            ].copy()
            _draw_boundary(ax, boundary_geom, linewidth=0.8, zorder=2)
            if not sdf.empty:
                last_scatter = ax.scatter(
                    sdf["longitude"],
                    sdf["latitude"],
                    c=pd.to_numeric(sdf[value_col], errors="coerce"),
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    s=34,
                    edgecolor="black",
                    linewidth=0.35,
                    zorder=4,
                )
            xmin, xmax, ymin, ymax = _get_plot_extent(sdf if not sdf.empty else station_frequency, boundary_geom)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_title(f"{definition['definition_title']} | RP >= {threshold:g} yr", loc="left")
            _format_geo_axis(ax, show_x=(r == len(definitions) - 1), show_y=(c == 0))
            ax.set_aspect("equal", adjustable="box")
    if last_scatter is not None:
        cbar = fig.colorbar(last_scatter, ax=axes.ravel().tolist(), shrink=0.78, pad=0.01)
        cbar.set_label(title_suffix)
    fig.savefig(outpath, dpi=get_plot_dpi(cfg))
    plt.close(fig)


def _plot_driver_shift(driver_summary: pd.DataFrame, outpath: Path, cfg: dict) -> None:
    if driver_summary.empty:
        return
    selected = driver_summary.loc[driver_summary["return_period_threshold_years"].isin([10.0, 20.0])].copy()
    if selected.empty:
        return
    apply_publication_theme()
    selected["label"] = (
        selected["definition_title"]
        + " | RP >= "
        + selected["return_period_threshold_years"].map(lambda v: f"{float(v):g}")
        + " | "
        + selected["period"]
    )
    fig, ax = plt.subplots(figsize=(11.5, 5.8), constrained_layout=True)
    x = np.arange(len(selected))
    dry = selected["dry_dominant_fraction_pct"].to_numpy(dtype=float)
    hot = selected["hot_dominant_fraction_pct"].to_numpy(dtype=float)
    co = 100.0 - dry - hot
    ax.bar(x, dry, label="Dry-dominant", color="#8f6b33")
    ax.bar(x, hot, bottom=dry, label="Hot-dominant", color="#b54432")
    ax.bar(x, co, bottom=dry + hot, label="Co-dominant", color="#6f7378")
    ax.set_ylabel("Compound-event share (%)")
    ax.set_xticks(x)
    ax.set_xticklabels(selected["label"], rotation=35, ha="right")
    ax.set_ylim(0, 100)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False, ncols=3)
    fig.savefig(outpath, dpi=get_plot_dpi(cfg))
    plt.close(fig)


def _plot_moran_timeseries(moran_yearly: pd.DataFrame, outpath: Path, cfg: dict) -> None:
    if moran_yearly.empty:
        return
    selected = moran_yearly.loc[moran_yearly["return_period_threshold_years"].isin([10.0, 20.0])].dropna(subset=["moran_i"]).copy()
    if selected.empty:
        return
    apply_publication_theme()
    definitions = selected[["definition", "definition_title"]].drop_duplicates().to_dict("records")
    fig, axes = plt.subplots(len(definitions), 1, figsize=(9.6, 3.5 * len(definitions)), sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    colors = {10.0: "#35618f", 20.0: "#b35d2e"}
    for ax, definition in zip(axes, definitions):
        sdf = selected.loc[selected["definition"] == definition["definition"]]
        for threshold, line_df in sdf.groupby("return_period_threshold_years"):
            line_df = line_df.sort_values("year")
            ax.plot(line_df["year"], line_df["moran_i"], marker="o", markersize=3.5, linewidth=1.5, color=colors.get(float(threshold)), label=f"RP >= {float(threshold):g} yr")
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_title(definition["definition_title"], loc="left")
        ax.set_ylabel("Moran's I")
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(frameon=False, ncols=2)
    axes[-1].set_xlabel("Year")
    fig.savefig(outpath, dpi=get_plot_dpi(cfg))
    plt.close(fig)


def run_compound_dry_hot_analysis(
    data: pd.DataFrame,
    stations: pd.DataFrame,
    cfg: dict,
    outdir: Path,
    progress_callback: Callable[[str], None] | None = None,
) -> dict[str, pd.DataFrame]:
    module_cfg = _cfg(cfg)
    if not module_cfg.get("enabled", False):
        return {}
    definitions = module_cfg.get("definitions", [])
    if not definitions:
        return {}

    if progress_callback:
        progress_callback("Running compound dry-hot empirical-copula analysis...")

    compound_dir = outdir / str(module_cfg.get("output_subdir", "compound_dry_hot"))
    tables_dir = compound_dir / "tables"
    figs_dir = compound_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    aggregates = pd.concat([_build_definition_aggregates(data, cfg, definition) for definition in definitions], ignore_index=True)
    station_year = _add_empirical_return_periods(aggregates, cfg)
    if station_year.empty:
        return {}

    baseline_years, comparison_years = _periods_from_cfg(cfg, station_year["year"])
    yearly_extent = _make_yearly_extent(station_year, cfg)
    trend = _trend_summary(yearly_extent, cfg, baseline_years, comparison_years)
    distribution_tests = _distribution_shift_tests(yearly_extent, baseline_years, comparison_years)
    frequency = _station_frequency(station_year, stations, cfg, baseline_years, comparison_years)
    driver = _driver_summary(station_year, cfg, baseline_years, comparison_years)
    moran_yearly, moran_trend = _moran_yearly(station_year, stations, cfg)

    station_year.to_csv(tables_dir / "compound_dry_hot_station_year.csv", index=False)
    yearly_extent.to_csv(tables_dir / "compound_dry_hot_yearly_extent.csv", index=False)
    trend.to_csv(tables_dir / "compound_dry_hot_trend_summary.csv", index=False)
    distribution_tests.to_csv(tables_dir / "compound_dry_hot_distribution_shift_tests.csv", index=False)
    frequency.to_csv(tables_dir / "compound_dry_hot_station_frequency.csv", index=False)
    driver.to_csv(tables_dir / "compound_dry_hot_driver_summary.csv", index=False)
    moran_yearly.to_csv(tables_dir / "compound_dry_hot_moran_yearly.csv", index=False)
    moran_trend.to_csv(tables_dir / "compound_dry_hot_moran_trend_summary.csv", index=False)

    _plot_extent_timeseries(yearly_extent, figs_dir / f"compound_dry_hot_extent_timeseries.{cfg['plots']['save_format']}", cfg)
    _plot_station_frequency_maps(
        frequency,
        figs_dir / f"compound_dry_hot_station_frequency_maps.{cfg['plots']['save_format']}",
        cfg,
        "event_frequency_pct",
        "Event frequency over valid years (%)",
        "magma_r",
    )
    _plot_station_frequency_maps(
        frequency,
        figs_dir / f"compound_dry_hot_frequency_change_maps.{cfg['plots']['save_format']}",
        cfg,
        "comparison_minus_baseline_frequency_pct",
        "Late-minus-early event frequency (percentage points)",
        "RdBu_r",
    )
    _plot_driver_shift(driver, figs_dir / f"compound_dry_hot_driver_shift.{cfg['plots']['save_format']}", cfg)
    _plot_moran_timeseries(moran_yearly, figs_dir / f"compound_dry_hot_moran_connectedness.{cfg['plots']['save_format']}", cfg)

    if progress_callback:
        progress_callback(f"Saved compound dry-hot tables and figures under {compound_dir}.")

    return {
        "compound_dry_hot_station_year": station_year,
        "compound_dry_hot_yearly_extent": yearly_extent,
        "compound_dry_hot_trend_summary": trend,
        "compound_dry_hot_distribution_shift_tests": distribution_tests,
        "compound_dry_hot_station_frequency": frequency,
        "compound_dry_hot_driver_summary": driver,
        "compound_dry_hot_moran_yearly": moran_yearly,
        "compound_dry_hot_moran_trend_summary": moran_trend,
    }
