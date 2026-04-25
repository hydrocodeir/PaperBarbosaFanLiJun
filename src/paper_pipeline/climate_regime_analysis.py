from __future__ import annotations

import re
from pathlib import Path
from typing import Callable, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

from .climate_change_signal import _expected_sign, _is_expected_direction
from .config_utils import get_focus_quantiles, get_primary_delta, slope_col
from .plotting import (
    _draw_boundary,
    _format_geo_axis,
    _get_boundary_path_from_cfg,
    _get_plot_extent,
    _load_boundary_geometry,
    _panel_label,
    apply_publication_theme,
)
from .progress_utils import ProgressTracker


REGIME_ORDER = [
    "BWh_hot_desert",
    "BWk_cold_desert",
    "BSh_hot_steppe",
    "BSk_cold_steppe",
    "C_temperate",
    "Dsa_cold_dry_summer",
]

REGIME_LABELS = {
    "BWh_hot_desert": "BWh hot desert",
    "BWk_cold_desert": "BWk cold desert",
    "BSh_hot_steppe": "BSh hot steppe",
    "BSk_cold_steppe": "BSk cold steppe",
    "C_temperate": "C temperate",
    "Dsa_cold_dry_summer": "Dsa cold dry-summer",
}

REGIME_COLORS = {
    "BWh_hot_desert": "#b2182b",
    "BWk_cold_desert": "#ef8a62",
    "BSh_hot_steppe": "#fdb863",
    "BSk_cold_steppe": "#fee08b",
    "C_temperate": "#66bd63",
    "Dsa_cold_dry_summer": "#5e4fa2",
}


def _default_regime_cfg(cfg: dict) -> dict:
    user_cfg = cfg.get("advanced_analyses", {}).get("climate_regime_analysis", {})
    return {
        "enabled": bool(user_cfg.get("enabled", False)),
        "koppen_raster": str(
            user_cfg.get(
                "koppen_raster",
                "data/external/koppen_geiger/Beck_KG_V1_present_0p0083.tif",
            )
        ),
        "koppen_confidence_raster": str(
            user_cfg.get(
                "koppen_confidence_raster",
                "data/external/koppen_geiger/Beck_KG_V1_present_conf_0p0083.tif",
            )
        ),
        "legend_path": str(user_cfg.get("legend_path", "data/external/koppen_geiger/legend.txt")),
        "min_regime_stations": int(user_cfg.get("min_regime_stations", 8)),
        "permutation_tests": int(user_cfg.get("permutation_tests", 999)),
        "nearest_valid_radius_pixels": int(user_cfg.get("nearest_valid_radius_pixels", 10)),
    }


def _read_legend(path: Path) -> dict[int, dict[str, str]]:
    if not path.exists():
        return {}
    pattern = re.compile(r"^\s*(\d+):\s+([A-Za-z]+)\s+(.+?)\s+\[")
    out: dict[int, dict[str, str]] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = pattern.match(line)
        if not m:
            continue
        out[int(m.group(1))] = {"kg_class": m.group(2), "kg_description": m.group(3).strip()}
    return out


def _regime_from_class(kg_class: str) -> tuple[str, str]:
    if kg_class == "BWh":
        return "BWh_hot_desert", REGIME_LABELS["BWh_hot_desert"]
    if kg_class == "BWk":
        return "BWk_cold_desert", REGIME_LABELS["BWk_cold_desert"]
    if kg_class == "BSh":
        return "BSh_hot_steppe", REGIME_LABELS["BSh_hot_steppe"]
    if kg_class == "BSk":
        return "BSk_cold_steppe", REGIME_LABELS["BSk_cold_steppe"]
    if kg_class.startswith("C"):
        return "C_temperate", REGIME_LABELS["C_temperate"]
    if kg_class == "Dsa":
        return "Dsa_cold_dry_summer", REGIME_LABELS["Dsa_cold_dry_summer"]
    return f"{kg_class}_other", kg_class


def _nearest_valid_pixel(arr: np.ndarray, row: int, col: int, max_radius: int) -> tuple[int, int, int] | None:
    nrows, ncols = arr.shape
    if 0 <= row < nrows and 0 <= col < ncols and int(arr[row, col]) > 0:
        return row, col, 0
    for radius in range(1, max_radius + 1):
        r0, r1 = max(0, row - radius), min(nrows, row + radius + 1)
        c0, c1 = max(0, col - radius), min(ncols, col + radius + 1)
        window = arr[r0:r1, c0:c1]
        rr, cc = np.where(window > 0)
        if len(rr) == 0:
            continue
        rr_abs = rr + r0
        cc_abs = cc + c0
        dist2 = (rr_abs - row) ** 2 + (cc_abs - col) ** 2
        best = int(np.argmin(dist2))
        return int(rr_abs[best]), int(cc_abs[best]), radius
    return None


def _assign_koppen_geiger(stations: pd.DataFrame, cfg: dict, regime_cfg: dict, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures" / "advanced_climate_regimes"
    figs_dir.mkdir(parents=True, exist_ok=True)
    legend = _read_legend(Path(regime_cfg["legend_path"]))
    raster_path = Path(regime_cfg["koppen_raster"])
    confidence_path = Path(regime_cfg["koppen_confidence_raster"])
    if not raster_path.exists():
        raise FileNotFoundError(f"Köppen-Geiger raster not found: {raster_path}")

    conf_arr = None
    conf_ds = None
    if confidence_path.exists():
        conf_ds = rasterio.open(confidence_path)
        conf_arr = conf_ds.read(1)

    rows = []
    with rasterio.open(raster_path) as ds:
        arr = ds.read(1)
        for _, station in stations.iterrows():
            lon = float(station["longitude"])
            lat = float(station["latitude"])
            row, col = ds.index(lon, lat)
            nearest = _nearest_valid_pixel(arr, row, col, int(regime_cfg["nearest_valid_radius_pixels"]))
            if nearest is None:
                kg_code, sample_row, sample_col, radius = 0, row, col, np.nan
            else:
                sample_row, sample_col, radius = nearest
                kg_code = int(arr[sample_row, sample_col])
            info = legend.get(kg_code, {"kg_class": "Unclassified", "kg_description": "Unclassified or water"})
            regime_code, regime_label = _regime_from_class(info["kg_class"])
            confidence = np.nan
            if conf_arr is not None and 0 <= sample_row < conf_arr.shape[0] and 0 <= sample_col < conf_arr.shape[1]:
                confidence = float(conf_arr[sample_row, sample_col])
            rows.append(
                {
                    "station_id": station["station_id"],
                    "station_name": station["station_name"],
                    "latitude": lat,
                    "longitude": lon,
                    "elevation": float(station["elevation"]) if "elevation" in station else np.nan,
                    "kg_code": kg_code,
                    "kg_class": info["kg_class"],
                    "kg_description": info["kg_description"],
                    "kg_major_group": str(info["kg_class"])[0] if info["kg_class"] != "Unclassified" else "U",
                    "kg_confidence_pct": confidence,
                    "assigned_by_nearest_valid": bool(radius and radius > 0),
                    "nearest_valid_radius_pixels": radius,
                    "climate_regime": regime_code,
                    "climate_regime_label": regime_label,
                }
            )
    if conf_ds is not None:
        conf_ds.close()

    station_df = pd.DataFrame(rows)
    station_df["climate_regime"] = pd.Categorical(station_df["climate_regime"], categories=REGIME_ORDER, ordered=True)
    station_df = station_df.sort_values(["climate_regime", "station_name"])
    station_df.to_csv(tables_dir / "koppen_geiger_station_assignments.csv", index=False)

    summary_rows = []
    for regime, sdf in station_df.groupby("climate_regime", observed=True):
        classes = ", ".join(f"{k}:{v}" for k, v in sdf["kg_class"].value_counts().sort_index().items())
        summary_rows.append(
            {
                "climate_regime": regime,
                "climate_regime_label": REGIME_LABELS.get(str(regime), str(regime)),
                "n_stations": int(sdf["station_id"].nunique()),
                "kg_classes": classes,
                "mean_latitude": float(sdf["latitude"].mean()),
                "mean_longitude": float(sdf["longitude"].mean()),
                "mean_elevation_m": float(sdf["elevation"].mean()),
                "median_elevation_m": float(sdf["elevation"].median()),
                "mean_kg_confidence_pct": float(pd.to_numeric(sdf["kg_confidence_pct"], errors="coerce").mean()),
                "nearest_valid_assignments": int(sdf["assigned_by_nearest_valid"].sum()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(tables_dir / "koppen_geiger_regime_summary.csv", index=False)
    _plot_regime_station_map(station_df, cfg, figs_dir / f"koppen_geiger_station_regimes.{cfg['plots']['save_format']}")
    return station_df, summary


def _merge_regime(df: pd.DataFrame, station_regimes: pd.DataFrame) -> pd.DataFrame:
    cols = ["station_id", "station_name", "climate_regime", "climate_regime_label"]
    out = df.merge(station_regimes[cols], on=["station_id", "station_name"], how="left")
    out["climate_regime"] = pd.Categorical(out["climate_regime"], categories=REGIME_ORDER, ordered=True)
    return out


def _summarize_quantile_by_regime(
    qr_summary: pd.DataFrame,
    station_regimes: pd.DataFrame,
    cfg: dict,
    outdir: Path,
) -> pd.DataFrame:
    tables_dir = outdir / "tables"
    focus = get_focus_quantiles(cfg)
    primary_delta = get_primary_delta(cfg)
    qdf = _merge_regime(qr_summary, station_regimes)
    fdr_path = tables_dir / "station_significance_fdr.csv"
    fdr = pd.read_csv(fdr_path) if fdr_path.exists() else pd.DataFrame()
    if not fdr.empty:
        fdr = _merge_regime(fdr, station_regimes)
        fdr = fdr.loc[np.isclose(pd.to_numeric(fdr["tau"], errors="coerce"), 0.50)].copy()
    rows = []
    for (regime, idx), sdf in qdf.groupby(["climate_regime", "index_name"], observed=True):
        row: dict[str, object] = {
            "climate_regime": regime,
            "climate_regime_label": REGIME_LABELS.get(str(regime), str(regime)),
            "index_name": idx,
            "n_stations": int(sdf["station_id"].nunique()),
            "mean_ols_slope": float(pd.to_numeric(sdf["ols_slope"], errors="coerce").mean()),
            "median_ols_slope": float(pd.to_numeric(sdf["ols_slope"], errors="coerce").median()),
            f"mean_{primary_delta}": float(pd.to_numeric(sdf[primary_delta], errors="coerce").mean()),
            f"median_{primary_delta}": float(pd.to_numeric(sdf[primary_delta], errors="coerce").median()),
            f"{primary_delta}_expected_direction_fraction": float(
                np.mean([_is_expected_direction(idx, v) for v in pd.to_numeric(sdf[primary_delta], errors="coerce").dropna()])
            ),
        }
        for tau in focus:
            col = slope_col(tau)
            vals = pd.to_numeric(sdf[col], errors="coerce").dropna()
            row[f"mean_{col}"] = float(vals.mean())
            row[f"median_{col}"] = float(vals.median())
            row[f"{col}_expected_direction_fraction"] = float(np.mean([_is_expected_direction(idx, v) for v in vals])) if len(vals) else np.nan
        if not fdr.empty:
            fs = fdr.loc[(fdr["climate_regime"] == regime) & (fdr["index_name"] == idx)]
            retained = fs.loc[pd.to_numeric(fs["fdr_reject"], errors="coerce").fillna(0) > 0]
            row["q50_fdr_retained_stations"] = int(retained["station_id"].nunique())
            row["q50_fdr_retained_fraction"] = float(retained["station_id"].nunique() / fs["station_id"].nunique()) if fs["station_id"].nunique() else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "climate_regime_quantile_summary.csv", index=False)
    _plot_regime_quantile_profiles(out, cfg, outdir / "figures" / "advanced_climate_regimes" / f"climate_regime_quantile_profiles.{cfg['plots']['save_format']}")
    return out


def _summarize_fixed_baseline_by_regime(station_regimes: pd.DataFrame, outdir: Path, cfg: dict) -> pd.DataFrame:
    tables_dir = outdir / "tables"
    path = tables_dir / "fixed_baseline_period_change_station_level.csv"
    if not path.exists():
        return pd.DataFrame()
    df = _merge_regime(pd.read_csv(path), station_regimes)
    rows = []
    for (regime, idx), sdf in df.groupby(["climate_regime", "index_name"], observed=True):
        vals = pd.to_numeric(sdf["late_minus_baseline_days"], errors="coerce").dropna()
        rows.append(
            {
                "climate_regime": regime,
                "climate_regime_label": REGIME_LABELS.get(str(regime), str(regime)),
                "index_name": idx,
                "n_stations": int(len(vals)),
                "mean_late_minus_baseline_days": float(vals.mean()),
                "median_late_minus_baseline_days": float(vals.median()),
                "direction_consistent_stations": int(sdf["direction_consistent"].sum()),
                "direction_consistency_fraction": float(sdf["direction_consistent"].mean()),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "climate_regime_fixed_baseline_summary.csv", index=False)
    _plot_regime_fixed_baseline(out, cfg, outdir / "figures" / "advanced_climate_regimes" / f"climate_regime_fixed_baseline_shifts.{cfg['plots']['save_format']}")
    return out


def _summarize_warming_by_regime(station_regimes: pd.DataFrame, outdir: Path, cfg: dict) -> pd.DataFrame:
    tables_dir = outdir / "tables"
    path = tables_dir / "warming_link_station_quantile_response.csv"
    if not path.exists():
        return pd.DataFrame()
    df = _merge_regime(pd.read_csv(path), station_regimes)
    focus = get_focus_quantiles(cfg)
    rows = []
    for (regime, idx), sdf in df.groupby(["climate_regime", "index_name"], observed=True):
        row: dict[str, object] = {
            "climate_regime": regime,
            "climate_regime_label": REGIME_LABELS.get(str(regime), str(regime)),
            "index_name": idx,
            "n_stations": int(sdf["station_id"].nunique()),
            "mean_ols_slope_per_c": float(pd.to_numeric(sdf["ols_slope_per_c"], errors="coerce").mean()),
            "ols_expected_direction_fraction": float(sdf["ols_expected_direction"].mean()),
        }
        for tau in focus:
            suffix = f"{tau:0.2f}"
            vals = pd.to_numeric(sdf[f"slope_per_c_{suffix}"], errors="coerce").dropna()
            row[f"mean_slope_per_c_{suffix}"] = float(vals.mean())
            row[f"median_slope_per_c_{suffix}"] = float(vals.median())
            row[f"expected_direction_fraction_{suffix}"] = float(sdf[f"expected_direction_{suffix}"].mean())
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "climate_regime_warming_response_summary.csv", index=False)
    return out


def _summarize_emergence_by_regime(station_regimes: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    tables_dir = outdir / "tables"
    path = tables_dir / "climate_signal_emergence_station_level.csv"
    if not path.exists():
        return pd.DataFrame()
    df = _merge_regime(pd.read_csv(path), station_regimes)
    rows = []
    for (regime, idx, metric), sdf in df.groupby(["climate_regime", "index_name", "metric"], observed=True):
        valid = sdf.dropna(subset=["signal_to_noise"])
        rows.append(
            {
                "climate_regime": regime,
                "climate_regime_label": REGIME_LABELS.get(str(regime), str(regime)),
                "index_name": idx,
                "metric": metric,
                "n_stations": int(valid["station_id"].nunique()),
                "expected_direction_stations": int(valid["expected_direction"].sum()),
                "emerged_stations": int(valid["emerged"].sum()),
                "expected_direction_fraction": float(valid["expected_direction"].mean()) if len(valid) else np.nan,
                "emergence_fraction": float(valid["emerged"].mean()) if len(valid) else np.nan,
                "median_signal_to_noise": float(valid["signal_to_noise"].median()) if len(valid) else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "climate_regime_emergence_summary.csv", index=False)
    return out


def _fingerprint_by_regime(
    quantile_summary: pd.DataFrame,
    fixed_summary: pd.DataFrame,
    warming_summary: pd.DataFrame,
    emergence_summary: pd.DataFrame,
    cfg: dict,
    outdir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tables_dir = outdir / "tables"
    primary_delta = get_primary_delta(cfg)
    rows = []
    for _, row in quantile_summary.iterrows():
        base = {
            "climate_regime": row["climate_regime"],
            "climate_regime_label": row["climate_regime_label"],
            "index_name": row["index_name"],
        }
        rows.append({**base, "component": "q50_trend_direction", "component_value": row.get("slope_0.50_expected_direction_fraction", np.nan)})
        rows.append({**base, "component": "tail_asymmetry_direction", "component_value": row.get(f"{primary_delta}_expected_direction_fraction", np.nan)})
        rows.append({**base, "component": "q50_fdr_field_significance", "component_value": row.get("q50_fdr_retained_fraction", np.nan)})
    for _, row in fixed_summary.iterrows():
        rows.append(
            {
                "climate_regime": row["climate_regime"],
                "climate_regime_label": row["climate_regime_label"],
                "index_name": row["index_name"],
                "component": "fixed_baseline_period_shift",
                "component_value": row["direction_consistency_fraction"],
            }
        )
    for _, row in warming_summary.iterrows():
        rows.append(
            {
                "climate_regime": row["climate_regime"],
                "climate_regime_label": row["climate_regime_label"],
                "index_name": row["index_name"],
                "component": "warming_link_q50_direction",
                "component_value": row.get("expected_direction_fraction_0.50", np.nan),
            }
        )
    if not emergence_summary.empty:
        es = emergence_summary.loc[emergence_summary["metric"].astype(str) == "0.50"]
        for _, row in es.iterrows():
            rows.append(
                {
                    "climate_regime": row["climate_regime"],
                    "climate_regime_label": row["climate_regime_label"],
                    "index_name": row["index_name"],
                    "component": "q50_signal_emergence",
                    "component_value": row["emergence_fraction"],
                }
            )
    components = pd.DataFrame(rows)
    components["component_value"] = pd.to_numeric(components["component_value"], errors="coerce")
    components.to_csv(tables_dir / "climate_regime_fingerprint_components.csv", index=False)

    by_index = (
        components.groupby(["climate_regime", "climate_regime_label", "index_name"], observed=True, as_index=False)
        .agg(regime_index_fingerprint_score=("component_value", "mean"))
    )
    by_regime = (
        components.groupby(["climate_regime", "climate_regime_label"], observed=True, as_index=False)
        .agg(overall_regime_fingerprint_score=("component_value", "mean"))
    )
    summary = by_index.merge(by_regime, on=["climate_regime", "climate_regime_label"], how="left")
    summary.to_csv(tables_dir / "climate_regime_fingerprint_summary.csv", index=False)
    _plot_regime_fingerprint_heatmap(summary, cfg, outdir / "figures" / "advanced_climate_regimes" / f"climate_regime_fingerprint_scores.{cfg['plots']['save_format']}")
    return components, summary


def _bh_qvalues(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(q, index=pvalues.index)
    idx = np.where(valid)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.clip(ranked, 0, 1)
    return pd.Series(q, index=pvalues.index)


def _permutation_difference_tests(
    station_regimes: pd.DataFrame,
    qr_summary: pd.DataFrame,
    regime_cfg: dict,
    cfg: dict,
    outdir: Path,
) -> pd.DataFrame:
    rng = np.random.default_rng(int(cfg["project"]["random_seed"]) + 311)
    min_n = int(regime_cfg["min_regime_stations"])
    n_perm = int(regime_cfg["permutation_tests"])
    primary_delta = get_primary_delta(cfg)
    datasets: list[tuple[str, pd.DataFrame, list[str]]] = [
        ("trend", _merge_regime(qr_summary, station_regimes), [slope_col(t) for t in get_focus_quantiles(cfg)] + [primary_delta]),
    ]
    tables_dir = outdir / "tables"
    fixed_path = tables_dir / "fixed_baseline_period_change_station_level.csv"
    if fixed_path.exists():
        datasets.append(("fixed_baseline", _merge_regime(pd.read_csv(fixed_path), station_regimes), ["late_minus_baseline_days"]))
    warming_path = tables_dir / "warming_link_station_quantile_response.csv"
    if warming_path.exists():
        datasets.append(("warming_link", _merge_regime(pd.read_csv(warming_path), station_regimes), [f"slope_per_c_{t:0.2f}" for t in get_focus_quantiles(cfg)]))

    rows = []
    for source, df, metrics in datasets:
        for idx, idf in df.groupby("index_name"):
            for metric in metrics:
                if metric not in idf.columns:
                    continue
                sdf = idf[["station_id", "climate_regime", metric]].dropna().copy()
                counts = sdf.groupby("climate_regime", observed=True)["station_id"].nunique()
                keep = counts[counts >= min_n].index
                sdf = sdf.loc[sdf["climate_regime"].isin(keep)]
                if sdf["climate_regime"].nunique() < 2:
                    continue
                values = pd.to_numeric(sdf[metric], errors="coerce").to_numpy(dtype=float)
                groups = sdf["climate_regime"].astype(str).to_numpy()
                obs_means = pd.Series(values).groupby(groups).mean()
                observed_range = float(obs_means.max() - obs_means.min())
                perm_ranges = []
                for _ in range(n_perm):
                    shuffled = rng.permutation(groups)
                    means = pd.Series(values).groupby(shuffled).mean()
                    perm_ranges.append(float(means.max() - means.min()))
                perm_ranges_arr = np.asarray(perm_ranges, dtype=float)
                p_value = float((1 + np.sum(perm_ranges_arr >= observed_range)) / (len(perm_ranges_arr) + 1))
                rows.append(
                    {
                        "metric_source": source,
                        "index_name": idx,
                        "metric": metric,
                        "n_regimes": int(sdf["climate_regime"].nunique()),
                        "n_stations": int(sdf["station_id"].nunique()),
                        "observed_regime_mean_range": observed_range,
                        "max_regime": str(obs_means.idxmax()),
                        "max_regime_mean": float(obs_means.max()),
                        "min_regime": str(obs_means.idxmin()),
                        "min_regime_mean": float(obs_means.min()),
                        "permutation_p_value": p_value,
                        "n_permutations": n_perm,
                    }
                )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_q_value"] = _bh_qvalues(out["permutation_p_value"])
        out.to_csv(tables_dir / "climate_regime_difference_tests.csv", index=False)
    return out


def _plot_regime_station_map(station_df: pd.DataFrame, cfg: dict, outpath: Path) -> None:
    apply_publication_theme(cfg)
    boundary_geom = _load_boundary_geometry(_get_boundary_path_from_cfg(cfg))
    fig, ax = plt.subplots(figsize=(8.2, 8.4), constrained_layout=True)
    _draw_boundary(ax, boundary_geom, linewidth=0.95, zorder=2)
    for regime in REGIME_ORDER:
        sdf = station_df.loc[station_df["climate_regime"].astype(str) == regime]
        if sdf.empty:
            continue
        ax.scatter(
            sdf["longitude"],
            sdf["latitude"],
            s=42,
            color=REGIME_COLORS.get(regime, "#777777"),
            edgecolor="black",
            linewidth=0.35,
            label=f"{REGIME_LABELS[regime]} (n={len(sdf)})",
            zorder=3,
        )
    xmin, xmax, ymin, ymax = _get_plot_extent(station_df, boundary_geom=boundary_geom, pad_deg=0.25)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    _format_geo_axis(ax, show_x=True, show_y=True)
    ax.set_title("Köppen-Geiger Climate Regimes of Station Network", loc="left")
    ax.legend(frameon=True, loc="lower left", fontsize=8)
    ax.grid(False)
    fig.savefig(outpath, dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def _plot_regime_quantile_profiles(summary: pd.DataFrame, cfg: dict, outpath: Path) -> None:
    apply_publication_theme(cfg)
    focus = get_focus_quantiles(cfg)
    regimes = [r for r in REGIME_ORDER if r in set(summary["climate_regime"].astype(str))]
    x = np.arange(len(regimes))
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.8), constrained_layout=True)
    colors = {0.10: "#2166ac", 0.50: "#202020", 0.90: "#b2182b"}
    for i, idx_cfg in enumerate(cfg["indices"]):
        ax = axes.flat[i]
        sdf = summary.loc[summary["index_name"] == idx_cfg["name"]].set_index("climate_regime")
        for tau in focus:
            col = f"mean_{slope_col(tau)}"
            vals = [float(sdf.loc[r, col]) if r in sdf.index else np.nan for r in regimes]
            ax.plot(x, vals, marker="o", linewidth=1.8, color=colors.get(tau, "#555555"), label=f"q{tau:0.2f}")
        ax.axhline(0, color="#666666", linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([REGIME_LABELS[r].replace(" ", "\n", 1) for r in regimes], rotation=0, fontsize=8)
        ax.set_title(f"{_panel_label(i)} {idx_cfg['title']}", loc="left")
        ax.set_ylabel("Mean trend (days/year per decade)")
        if i == 0:
            ax.legend(frameon=False)
    fig.savefig(outpath, dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def _plot_regime_fixed_baseline(summary: pd.DataFrame, cfg: dict, outpath: Path) -> None:
    if summary.empty:
        return
    apply_publication_theme(cfg)
    regimes = [r for r in REGIME_ORDER if r in set(summary["climate_regime"].astype(str))]
    x = np.arange(len(regimes))
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.8), constrained_layout=True)
    for i, idx_cfg in enumerate(cfg["indices"]):
        ax = axes.flat[i]
        sdf = summary.loc[summary["index_name"] == idx_cfg["name"]].set_index("climate_regime")
        vals = [float(sdf.loc[r, "mean_late_minus_baseline_days"]) if r in sdf.index else np.nan for r in regimes]
        color = "#b2182b" if idx_cfg["name"].startswith("warm") else "#2166ac"
        ax.bar(x, vals, color=color, edgecolor="black", linewidth=0.35)
        ax.axhline(0, color="#666666", linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([REGIME_LABELS[r].replace(" ", "\n", 1) for r in regimes], rotation=0, fontsize=8)
        ax.set_title(f"{_panel_label(i)} {idx_cfg['title']}", loc="left")
        ax.set_ylabel("Late minus baseline (days/year)")
    fig.savefig(outpath, dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def _plot_regime_fingerprint_heatmap(summary: pd.DataFrame, cfg: dict, outpath: Path) -> None:
    if summary.empty:
        return
    apply_publication_theme(cfg)
    regimes = [r for r in REGIME_ORDER if r in set(summary["climate_regime"].astype(str))]
    indices = [idx["name"] for idx in cfg["indices"]]
    mat = np.full((len(regimes), len(indices)), np.nan, dtype=float)
    for i, regime in enumerate(regimes):
        for j, idx in enumerate(indices):
            row = summary.loc[(summary["climate_regime"].astype(str) == regime) & (summary["index_name"] == idx)]
            if not row.empty:
                mat[i, j] = float(row.iloc[0]["regime_index_fingerprint_score"])
    fig, ax = plt.subplots(figsize=(8.8, 5.8), constrained_layout=True)
    im = ax.imshow(mat, vmin=0, vmax=1, cmap="YlGnBu")
    ax.set_xticks(np.arange(len(indices)))
    ax.set_xticklabels([x.replace("_", " ").title() for x in indices], rotation=20, ha="right")
    ax.set_yticks(np.arange(len(regimes)))
    ax.set_yticklabels([REGIME_LABELS[r] for r in regimes])
    for i in range(len(regimes)):
        for j in range(len(indices)):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=8.5)
    ax.set_title("Climate-Regime Fingerprint Scores", loc="left")
    cbar = fig.colorbar(im, ax=ax, shrink=0.82)
    cbar.set_label("Fingerprint score")
    fig.savefig(outpath, dpi=int(cfg["plots"]["dpi"]))
    plt.close(fig)


def run_climate_regime_analysis(
    data: pd.DataFrame,
    annual: pd.DataFrame,
    qr_summary: pd.DataFrame,
    feature_table: pd.DataFrame,
    stations: pd.DataFrame,
    cfg: dict,
    outdir: Path,
    progress_callback: Callable[[str], None] | None = None,
) -> Dict[str, pd.DataFrame]:
    regime_cfg = _default_regime_cfg(cfg)
    if not regime_cfg["enabled"]:
        return {}

    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    (outdir / "figures" / "advanced_climate_regimes").mkdir(parents=True, exist_ok=True)
    tasks = [
        "station Köppen-Geiger assignment",
        "quantile summaries by regime",
        "fixed-baseline summaries by regime",
        "warming and emergence summaries by regime",
        "regime fingerprint and permutation tests",
    ]
    tracker = ProgressTracker("Climate-regime analysis", len(tasks), progress_callback)

    tracker.emit(1, detail="station Köppen-Geiger assignment")
    station_regimes, regime_summary = _assign_koppen_geiger(stations, cfg, regime_cfg, outdir)

    tracker.emit(2, detail="quantile summaries by regime")
    quantile_summary = _summarize_quantile_by_regime(qr_summary, station_regimes, cfg, outdir)

    tracker.emit(3, detail="fixed-baseline summaries by regime")
    fixed_summary = _summarize_fixed_baseline_by_regime(station_regimes, outdir, cfg)

    tracker.emit(4, detail="warming and emergence summaries by regime")
    warming_summary = _summarize_warming_by_regime(station_regimes, outdir, cfg)
    emergence_summary = _summarize_emergence_by_regime(station_regimes, outdir)

    tracker.emit(5, detail="regime fingerprint and permutation tests")
    fingerprint_components, fingerprint_summary = _fingerprint_by_regime(
        quantile_summary,
        fixed_summary,
        warming_summary,
        emergence_summary,
        cfg,
        outdir,
    )
    difference_tests = _permutation_difference_tests(station_regimes, qr_summary, regime_cfg, cfg, outdir)

    return {
        "koppen_geiger_station_assignments": station_regimes,
        "koppen_geiger_regime_summary": regime_summary,
        "climate_regime_quantile_summary": quantile_summary,
        "climate_regime_fixed_baseline_summary": fixed_summary,
        "climate_regime_warming_response_summary": warming_summary,
        "climate_regime_emergence_summary": emergence_summary,
        "climate_regime_fingerprint_components": fingerprint_components,
        "climate_regime_fingerprint_summary": fingerprint_summary,
        "climate_regime_difference_tests": difference_tests,
    }
