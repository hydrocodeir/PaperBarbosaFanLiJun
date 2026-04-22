from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler


def build_feature_table(summary: pd.DataFrame) -> pd.DataFrame:
    features = summary.copy()
    features["Delta1"] = features["slope_0.95"] - features["slope_0.05"]
    features["Delta2"] = features["slope_0.95"] - features["slope_0.50"]
    features["Delta3"] = features["slope_0.50"] - features["slope_0.05"]
    return features


def _resolve_feature_cols(cfg: dict, feature_cols: list[str] | None = None) -> list[str]:
    if feature_cols is not None:
        return list(feature_cols)
    mode = cfg["clustering"]["feature_mode"]
    if mode == "simple":
        return list(cfg["clustering"]["simple_features"])
    if mode == "uncertainty":
        return list(cfg["clustering"]["uncertainty_features"])
    raise ValueError(f"Unsupported feature mode: {mode}")


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
    algorithm = algorithm or cfg["clustering"]["algorithm"]
    linkage_method = linkage_method or cfg["clustering"]["linkage"]
    metric = metric or cfg["clustering"]["metric"]

    cluster_rows, artifacts = [], {}

    for idx_name, sdf in features.groupby("index_name"):
        if sdf.empty:
            continue
        X = sdf[feature_cols].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
        X = X.apply(lambda col: col.fillna(col.median()), axis=0)
        Xv = StandardScaler().fit_transform(X) if cfg["clustering"]["standardize"] else X.to_numpy()
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
