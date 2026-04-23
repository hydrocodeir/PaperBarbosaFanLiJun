from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler

from .config_utils import compute_defined_deltas


def build_feature_table(summary: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    return compute_defined_deltas(summary.copy(), cfg)


def _resolve_feature_cols(cfg: dict, feature_cols: list[str] | None = None) -> list[str]:
    if feature_cols is not None:
        return list(feature_cols)
    mode = cfg["clustering"]["feature_mode"]
    if mode == "simple":
        return list(cfg["clustering"]["simple_features"])
    if mode == "uncertainty":
        return list(cfg["clustering"]["uncertainty_features"])
    raise ValueError(f"Unsupported feature mode: {mode}")


def _screen_feature_matrix(
    X: pd.DataFrame,
    feature_cols: list[str],
    cfg: dict,
    *,
    feature_set_label: str,
    index_name: str,
) -> tuple[list[str], pd.DataFrame]:
    screen_cfg = cfg.get("clustering", {}).get("collinearity_screen", {})
    enabled = bool(screen_cfg.get("enabled", False))
    threshold = float(screen_cfg.get("abs_correlation_threshold", 0.95))

    rows: list[dict] = []
    if not enabled:
        for order, feature in enumerate(feature_cols, start=1):
            rows.append(
                {
                    "index_name": index_name,
                    "feature_set": feature_set_label,
                    "feature": feature,
                    "status": "kept_screening_disabled",
                    "reference_feature": "",
                    "abs_correlation_to_reference": np.nan,
                    "screen_order": order,
                }
            )
        return list(feature_cols), pd.DataFrame(rows)

    kept: list[str] = []
    for order, feature in enumerate(feature_cols, start=1):
        series = pd.to_numeric(X[feature], errors="coerce")
        if series.dropna().empty:
            rows.append(
                {
                    "index_name": index_name,
                    "feature_set": feature_set_label,
                    "feature": feature,
                    "status": "dropped_all_missing",
                    "reference_feature": "",
                    "abs_correlation_to_reference": np.nan,
                    "screen_order": order,
                }
            )
            continue
        if float(series.std(ddof=0)) <= 1e-12:
            rows.append(
                {
                    "index_name": index_name,
                    "feature_set": feature_set_label,
                    "feature": feature,
                    "status": "dropped_constant",
                    "reference_feature": "",
                    "abs_correlation_to_reference": np.nan,
                    "screen_order": order,
                }
            )
            continue
        if not kept:
            kept.append(feature)
            rows.append(
                {
                    "index_name": index_name,
                    "feature_set": feature_set_label,
                    "feature": feature,
                    "status": "kept",
                    "reference_feature": "",
                    "abs_correlation_to_reference": np.nan,
                    "screen_order": order,
                }
            )
            continue

        corr_pairs = []
        for ref in kept:
            corr = series.corr(pd.to_numeric(X[ref], errors="coerce"))
            corr_pairs.append((ref, abs(float(corr)) if pd.notna(corr) else 0.0))
        ref_feature, max_abs_corr = max(corr_pairs, key=lambda item: item[1])
        if max_abs_corr >= threshold:
            rows.append(
                {
                    "index_name": index_name,
                    "feature_set": feature_set_label,
                    "feature": feature,
                    "status": "dropped_high_correlation",
                    "reference_feature": ref_feature,
                    "abs_correlation_to_reference": max_abs_corr,
                    "screen_order": order,
                }
            )
            continue

        kept.append(feature)
        rows.append(
            {
                "index_name": index_name,
                "feature_set": feature_set_label,
                "feature": feature,
                "status": "kept",
                "reference_feature": ref_feature,
                "abs_correlation_to_reference": max_abs_corr,
                "screen_order": order,
            }
        )

    return kept, pd.DataFrame(rows)


def screen_clustering_features(
    features: pd.DataFrame,
    cfg: dict,
    feature_cols: list[str] | None = None,
    *,
    feature_set_label: str = "configured",
) -> pd.DataFrame:
    if features.empty or not cfg["clustering"]["enabled"]:
        return pd.DataFrame()

    requested_cols = _resolve_feature_cols(cfg, feature_cols)
    parts: list[pd.DataFrame] = []
    for idx_name, sdf in features.groupby("index_name"):
        X = sdf[requested_cols].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
        X = X.apply(lambda col: col.fillna(col.median()), axis=0)
        _, screen_df = _screen_feature_matrix(
            X,
            requested_cols,
            cfg,
            feature_set_label=feature_set_label,
            index_name=idx_name,
        )
        if not screen_df.empty:
            parts.append(screen_df)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def run_clustering(
    features: pd.DataFrame,
    cfg: dict,
    feature_cols: list[str] | None = None,
    label_col: str = "cluster",
    algorithm: str | None = None,
    linkage_method: str | None = None,
    metric: str | None = None,
) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    if features.empty or not cfg["clustering"]["enabled"]:
        return pd.DataFrame(), {}

    feature_cols = _resolve_feature_cols(cfg, feature_cols)
    feature_set_label = label_col
    algorithm = algorithm or cfg["clustering"]["algorithm"]
    linkage_method = linkage_method or cfg["clustering"]["linkage"]
    metric = metric or cfg["clustering"]["metric"]

    cluster_rows, artifacts = [], {}

    for idx_name, sdf in features.groupby("index_name"):
        if sdf.empty:
            continue
        X = sdf[feature_cols].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
        X = X.apply(lambda col: col.fillna(col.median()), axis=0)
        screened_cols, _ = _screen_feature_matrix(
            X,
            feature_cols,
            cfg,
            feature_set_label=feature_set_label,
            index_name=idx_name,
        )
        use_cols = screened_cols if screened_cols else list(feature_cols)
        X_use = X[use_cols].copy()
        Xv = StandardScaler().fit_transform(X_use) if cfg["clustering"]["standardize"] else X_use.to_numpy()
        n_samples = len(sdf)
        n_clusters = max(1, min(int(cfg["clustering"]["n_clusters"]), n_samples))

        if n_samples == 1:
            labels = np.array([1], dtype=int)
            artifacts[idx_name] = np.empty((0, 4), dtype=float)
        elif algorithm == "hierarchical":
            effective_metric = "euclidean" if linkage_method == "ward" else metric
            Z = linkage(Xv, method=linkage_method, metric=effective_metric)
            labels = fcluster(Z, t=n_clusters, criterion="maxclust")
            artifacts[idx_name] = Z
        elif algorithm == "kmeans":
            km = KMeans(n_clusters=n_clusters, random_state=int(cfg["project"]["random_seed"]), n_init=20)
            labels = km.fit_predict(Xv) + 1
            artifacts[idx_name] = km.cluster_centers_
        else:
            raise ValueError(f"Unsupported clustering algorithm: {algorithm}")

        cdf = sdf[["index_name", "station_id", "station_name"]].copy()
        cdf[label_col] = labels
        cluster_rows.append(cdf)

    if not cluster_rows:
        return pd.DataFrame(), artifacts
    return pd.concat(cluster_rows, ignore_index=True), artifacts


def compare_clusterings(
    primary_df: pd.DataFrame,
    secondary_df: pd.DataFrame,
    primary_label: str = "cluster",
    secondary_label: str = "cluster_reduced_features",
) -> pd.DataFrame:
    if primary_df.empty or secondary_df.empty:
        return pd.DataFrame()

    merged = primary_df.merge(
        secondary_df,
        on=["index_name", "station_id", "station_name"],
        how="inner",
    )
    if merged.empty:
        return pd.DataFrame()

    rows = []
    for idx_name, sdf in merged.groupby("index_name"):
        valid = sdf.dropna(subset=[primary_label, secondary_label]).copy()
        if valid.empty:
            continue
        rows.append(
            {
                "index_name": idx_name,
                "n_stations_compared": int(len(valid)),
                "adjusted_rand_index": float(
                    adjusted_rand_score(
                        valid[primary_label].to_numpy(dtype=int),
                        valid[secondary_label].to_numpy(dtype=int),
                    )
                ),
                "same_label_fraction": float(
                    np.mean(
                        valid[primary_label].to_numpy(dtype=int) == valid[secondary_label].to_numpy(dtype=int)
                    )
                ),
            }
        )
    return pd.DataFrame(rows)
