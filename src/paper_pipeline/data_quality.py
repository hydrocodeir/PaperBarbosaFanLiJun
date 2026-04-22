from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _pettitt_test(x: np.ndarray) -> tuple[float, float, int | None]:
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 8:
        return np.nan, np.nan, None
    U = np.zeros(n)
    for k in range(n):
        left = x[: k + 1]
        right = x[k + 1 :]
        if right.size == 0:
            U[k] = 0
        else:
            diffs = left[:, None] - right[None, :]
            U[k] = np.sign(diffs).sum()
    K = float(np.max(np.abs(U)))
    k_idx = int(np.argmax(np.abs(U)))
    p = float(2 * np.exp((-6 * K**2) / (n**3 + n**2)))
    p = min(max(p, 0.0), 1.0)
    return K, p, k_idx


def _snht_test(x: np.ndarray) -> tuple[float, int | None]:
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 8 or np.std(x, ddof=1) == 0:
        return np.nan, None
    z = (x - x.mean()) / np.std(x, ddof=1)
    scores = []
    for k in range(1, n):
        z1 = z[:k].mean()
        z2 = z[k:].mean()
        t = k * z1**2 + (n - k) * z2**2
        scores.append(t)
    arr = np.asarray(scores, dtype=float)
    idx = int(np.argmax(arr))
    return float(arr[idx]), idx + 1


def _buishand_range_test(x: np.ndarray) -> tuple[float, int | None]:
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 8 or np.std(x, ddof=1) == 0:
        return np.nan, None
    centered = x - x.mean()
    s = np.cumsum(centered)
    rng = (s.max() - s.min()) / np.std(x, ddof=1)
    idx = int(np.argmax(np.abs(s)))
    return float(rng), idx


def _permutation_pvalue(
    x: np.ndarray,
    stat_fn,
    n_perm: int = 499,
    seed: int = 42,
) -> float:
    x = np.asarray(x, dtype=float)
    if x.size < 8:
        return np.nan
    rng = np.random.default_rng(seed)
    observed = stat_fn(x)[0]
    if not np.isfinite(observed):
        return np.nan
    exceed = 0
    for _ in range(n_perm):
        shuffled = rng.permutation(x)
        stat = stat_fn(shuffled)[0]
        if np.isfinite(stat) and stat >= observed:
            exceed += 1
    return float((exceed + 1) / (n_perm + 1))


def run_data_quality_assessment(
    data: pd.DataFrame,
    cfg: dict,
    outdir: Path,
) -> dict[str, pd.DataFrame]:
    tables_dir = outdir / "tables"
    figs_dir = outdir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    station_id_col = cfg["data"]["station_id_col"]
    station_name_col = cfg["data"]["station_name_col"]
    year_col = cfg["data"]["year_col"]
    month_col = cfg["data"]["month_col"]
    day_col = cfg["data"]["day_col"]
    tmin_col = cfg["data"]["tmin_col"]
    tmax_col = cfg["data"]["tmax_col"]
    tmean_col = cfg["data"]["tmean_col"]

    df = data.copy()
    df["date"] = pd.to_datetime(
        dict(year=df[year_col], month=df[month_col], day=df[day_col]),
        errors="coerce",
    )
    df["tmean_candidate"] = df[tmean_col]
    fill_mask = df["tmean_candidate"].isna() & df[tmin_col].notna() & df[tmax_col].notna()
    df.loc[fill_mask, "tmean_candidate"] = (df.loc[fill_mask, tmin_col] + df.loc[fill_mask, tmax_col]) / 2.0

    def summarize_station(g: pd.DataFrame) -> pd.Series:
        n_rows = len(g)
        dup_dates = int(g["date"].duplicated().sum())
        valid_tmin = g[tmin_col].notna().sum()
        valid_tmax = g[tmax_col].notna().sum()
        valid_tmean = g["tmean_candidate"].notna().sum()
        inconsistent_tmin_tmax = int(((g[tmin_col] > g[tmax_col]) & g[tmin_col].notna() & g[tmax_col].notna()).sum())
        inconsistent_tmean = int(
            (
                (
                    (g["tmean_candidate"] < g[tmin_col]) | (g["tmean_candidate"] > g[tmax_col])
                )
                & g["tmean_candidate"].notna()
                & g[tmin_col].notna()
                & g[tmax_col].notna()
            ).sum()
        )
        return pd.Series(
            {
                "n_rows": n_rows,
                "duplicate_dates": dup_dates,
                "tmin_completeness_pct": 100.0 * valid_tmin / n_rows if n_rows else np.nan,
                "tmax_completeness_pct": 100.0 * valid_tmax / n_rows if n_rows else np.nan,
                "tmean_completeness_pct": 100.0 * valid_tmean / n_rows if n_rows else np.nan,
                "inconsistent_tmin_gt_tmax": inconsistent_tmin_tmax,
                "inconsistent_tmean_outside_range": inconsistent_tmean,
            }
        )

    station_qc = (
        df.groupby([station_id_col, station_name_col], as_index=False)
        .apply(summarize_station, include_groups=False)
        .reset_index()
    )
    if "level_2" in station_qc.columns:
        station_qc = station_qc.drop(columns=["level_2"])

    annual = (
        df.groupby([station_id_col, station_name_col, year_col], as_index=False)
        .agg(
            annual_mean_tmean=("tmean_candidate", "mean"),
            valid_days=("tmean_candidate", lambda s: int(s.notna().sum())),
        )
    )
    annual["coverage_pct"] = annual["valid_days"] / 365.0 * 100.0
    annual = annual.loc[annual["coverage_pct"] >= 80.0].copy()

    homogeneity_rows = []
    for (sid, sname), g in annual.groupby([station_id_col, station_name_col]):
        series = g.sort_values(year_col)["annual_mean_tmean"].to_numpy(dtype=float)
        years = g.sort_values(year_col)[year_col].to_numpy(dtype=int)
        t = np.arange(series.size, dtype=float)
        if series.size >= 3:
            slope, intercept = np.polyfit(t, series, 1)
            detrended = series - (slope * t + intercept)
        else:
            slope = np.nan
            detrended = series.copy()
        pettitt_stat, pettitt_p, pettitt_idx = _pettitt_test(series)
        snht_stat, snht_idx = _snht_test(series)
        buishand_stat, buishand_idx = _buishand_range_test(series)
        snht_p = _permutation_pvalue(series, _snht_test)
        buishand_p = _permutation_pvalue(series, _buishand_range_test)
        pettitt_dt_stat, pettitt_dt_p, pettitt_dt_idx = _pettitt_test(detrended)
        snht_dt_stat, snht_dt_idx = _snht_test(detrended)
        buishand_dt_stat, buishand_dt_idx = _buishand_range_test(detrended)
        snht_dt_p = _permutation_pvalue(detrended, _snht_test)
        buishand_dt_p = _permutation_pvalue(detrended, _buishand_range_test)
        homogeneity_rows.append(
            {
                station_id_col: sid,
                station_name_col: sname,
                "n_years_homogeneity_test": int(series.size),
                "mean_annual_coverage_pct": float(g["coverage_pct"].mean()),
                "linear_trend_degC_per_year": float(slope),
                "pettitt_stat": pettitt_stat,
                "pettitt_pvalue": pettitt_p,
                "pettitt_change_year": int(years[pettitt_idx]) if pettitt_idx is not None else np.nan,
                "snht_stat": snht_stat,
                "snht_pvalue": snht_p,
                "snht_change_year": int(years[snht_idx]) if snht_idx is not None and snht_idx < years.size else np.nan,
                "buishand_range_stat": buishand_stat,
                "buishand_pvalue": buishand_p,
                "buishand_change_year": int(years[buishand_idx]) if buishand_idx is not None else np.nan,
                "pettitt_detrended_stat": pettitt_dt_stat,
                "pettitt_detrended_pvalue": pettitt_dt_p,
                "pettitt_detrended_change_year": int(years[pettitt_dt_idx]) if pettitt_dt_idx is not None else np.nan,
                "snht_detrended_stat": snht_dt_stat,
                "snht_detrended_pvalue": snht_dt_p,
                "snht_detrended_change_year": int(years[snht_dt_idx]) if snht_dt_idx is not None and snht_dt_idx < years.size else np.nan,
                "buishand_detrended_stat": buishand_dt_stat,
                "buishand_detrended_pvalue": buishand_dt_p,
                "buishand_detrended_change_year": int(years[buishand_dt_idx]) if buishand_dt_idx is not None else np.nan,
                "any_homogeneity_flag_p_lt_0_05": bool(
                    (np.isfinite(pettitt_p) and pettitt_p < 0.05)
                    or (np.isfinite(snht_p) and snht_p < 0.05)
                    or (np.isfinite(buishand_p) and buishand_p < 0.05)
                ),
                "any_detrended_homogeneity_flag_p_lt_0_05": bool(
                    (np.isfinite(pettitt_dt_p) and pettitt_dt_p < 0.05)
                    or (np.isfinite(snht_dt_p) and snht_dt_p < 0.05)
                    or (np.isfinite(buishand_dt_p) and buishand_dt_p < 0.05)
                ),
            }
        )

    homogeneity = pd.DataFrame(homogeneity_rows)
    summary = pd.DataFrame(
        {
            "metric": [
                "n_stations",
                "median_tmin_completeness_pct",
                "median_tmax_completeness_pct",
                "median_tmean_completeness_pct",
                "stations_with_duplicate_dates",
                "stations_with_tmin_gt_tmax",
                "stations_with_tmean_outside_range",
                "stations_flagged_by_any_raw_homogeneity_test_p_lt_0_05",
                "stations_flagged_by_any_detrended_homogeneity_test_p_lt_0_05",
            ],
            "value": [
                int(station_qc.shape[0]),
                float(station_qc["tmin_completeness_pct"].median()),
                float(station_qc["tmax_completeness_pct"].median()),
                float(station_qc["tmean_completeness_pct"].median()),
                int((station_qc["duplicate_dates"] > 0).sum()),
                int((station_qc["inconsistent_tmin_gt_tmax"] > 0).sum()),
                int((station_qc["inconsistent_tmean_outside_range"] > 0).sum()),
                int(homogeneity["any_homogeneity_flag_p_lt_0_05"].sum()),
                int(homogeneity["any_detrended_homogeneity_flag_p_lt_0_05"].sum()),
            ],
        }
    )

    station_qc.to_csv(tables_dir / "data_quality_station_summary.csv", index=False)
    homogeneity.to_csv(tables_dir / "data_homogeneity_tests_station_summary.csv", index=False)
    summary.to_csv(tables_dir / "data_quality_homogeneity_overview.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), dpi=300)
    fig.patch.set_facecolor("white")

    completeness = station_qc[["tmin_completeness_pct", "tmax_completeness_pct", "tmean_completeness_pct"]]
    axes[0].boxplot(
        [completeness[c].dropna().to_numpy() for c in completeness.columns],
        labels=["Tmin", "Tmax", "Tmean*"],
        patch_artist=True,
        boxprops=dict(facecolor="#dceaf7", color="#476582"),
        medianprops=dict(color="#b54708", linewidth=1.5),
        whiskerprops=dict(color="#476582"),
        capprops=dict(color="#476582"),
    )
    axes[0].set_ylim(0, 100)
    axes[0].set_ylabel("Completeness (%)")
    axes[0].set_title("(a) Station completeness")
    axes[0].grid(axis="y", alpha=0.25)

    test_counts = [
        int((homogeneity["pettitt_detrended_pvalue"] < 0.05).sum()),
        int((homogeneity["snht_detrended_pvalue"] < 0.05).sum()),
        int((homogeneity["buishand_detrended_pvalue"] < 0.05).sum()),
        int(homogeneity["any_detrended_homogeneity_flag_p_lt_0_05"].sum()),
    ]
    axes[1].bar(
        ["Pettitt", "SNHT", "Buishand", "Any"],
        test_counts,
        color=["#7aa6d1", "#97c79c", "#f0b775", "#d98f8f"],
        edgecolor="#44546a",
    )
    axes[1].set_ylabel("Stations flagged (p < 0.05)")
    axes[1].set_title("(b) Detrended homogeneity-test summary")
    axes[1].grid(axis="y", alpha=0.25)

    fig.suptitle("Data completeness and homogeneity diagnostics", fontsize=12.5, fontweight="bold")
    fig.text(
        0.5,
        0.01,
        "* Tmean completeness uses observed Tmean where available and otherwise the average of Tmin and Tmax when both exist. Homogeneity bars summarize detrended annual mean temperature tests.",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#4a5568",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(figs_dir / "ijoc_data_quality_homogeneity.png", bbox_inches="tight")
    plt.close(fig)

    return {
        "station_qc": station_qc,
        "homogeneity": homogeneity,
        "summary": summary,
    }
