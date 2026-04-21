from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


def build_feature_table(summary: pd.DataFrame) -> pd.DataFrame:
    features = summary.copy()
    features["Delta1"] = features["slope_0.95"] - features["slope_0.05"]
    features["Delta2"] = features["slope_0.95"] - features["slope_0.50"]
    features["Delta3"] = features["slope_0.50"] - features["slope_0.05"]
    return features


def run_clustering(features: pd.DataFrame, cfg: dict) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    if features.empty or not cfg["clustering"]["enabled"]:
        return pd.DataFrame(), {}

    mode = cfg["clustering"]["feature_mode"]
    if mode == "simple":
        feature_cols = cfg["clustering"]["simple_features"]
    elif mode == "uncertainty":
        feature_cols = cfg["clustering"]["uncertainty_features"]
    else:
        raise ValueError(f"Unsupported feature mode: {mode}")

    cluster_rows, artifacts = [], {}

    for idx_name, sdf in features.groupby("index_name"):
        X = sdf[feature_cols].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
        X = X.apply(lambda col: col.fillna(col.median()), axis=0)
        Xv = StandardScaler().fit_transform(X) if cfg["clustering"]["standardize"] else X.to_numpy()

        if cfg["clustering"]["algorithm"] == "hierarchical":
            Z = linkage(Xv, method=cfg["clustering"]["linkage"], metric=cfg["clustering"]["metric"])
            labels = fcluster(Z, t=int(cfg["clustering"]["n_clusters"]), criterion="maxclust")
            artifacts[idx_name] = Z
        elif cfg["clustering"]["algorithm"] == "kmeans":
            km = KMeans(n_clusters=int(cfg["clustering"]["n_clusters"]), random_state=int(cfg["project"]["random_seed"]), n_init=20)
            labels = km.fit_predict(Xv) + 1
            artifacts[idx_name] = km.cluster_centers_
        else:
            raise ValueError(f"Unsupported clustering algorithm: {cfg['clustering']['algorithm']}")

        cdf = sdf[["index_name", "station_id", "station_name"]].copy()
        cdf["cluster"] = labels
        cluster_rows.append(cdf)

    return pd.concat(cluster_rows, ignore_index=True), artifacts
