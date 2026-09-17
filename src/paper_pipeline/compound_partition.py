"""Fixed-threshold compound frequencies and an exact descriptive partition.

The excess-joint term is indicator covariance, NOT a copula-only dependence
effect or physical attribution. All probabilities are computed station first.
Year blocks are synchronized across stations and both variables, and sampled
separately within the two periods. Baseline thresholds are refitted by default.
"""
from __future__ import annotations

import calendar
import hashlib
import importlib.metadata
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

COMPONENTS = ("joint_change", "dry_marginal", "hot_marginal", "excess_joint")


def partition(d0, h0, j0, d1, h1, j1):
    """Exact symmetric product partition, in probability units."""
    dry = (d1 - d0) * (h0 + h1) / 2
    hot = (h1 - h0) * (d0 + d1) / 2
    excess = (j1 - d1 * h1) - (j0 - d0 * h0)
    return np.stack((j1 - j0, dry, hot, excess), axis=-1)


def circular_blocks(n, length, reps, rng):
    starts = rng.integers(0, n, size=(reps, int(np.ceil(n / length))))
    return ((starts[..., None] + np.arange(length)) % n).reshape(reps, -1)[:, :n]


def aggregate_daily(raw, years):
    raw = raw.loc[raw.year.between(*years)].copy()
    date = pd.to_datetime(raw[["year", "month", "day"]], errors="raise")
    if raw.assign(date=date).duplicated(["station_id", "date"]).any():
        raise ValueError("Duplicate station-days must be resolved before analysis")
    conflict = raw.tmin.notna() & raw.tmax.notna() & (raw.tmin > raw.tmax)
    outside = (~conflict & raw.tmean.notna() & raw.tmin.notna() & raw.tmax.notna()
               & ((raw.tmean < raw.tmin) | (raw.tmean > raw.tmax)))
    audit = {"rows": len(raw), "stations": int(raw.station_id.nunique()),
             "tmin_gt_tmax": int(conflict.sum()), "tmean_outside_nonoverlap": int(outside.sum()),
             "negative_precipitation": int((raw.precip < 0).sum())}
    raw.loc[conflict, ["tmin", "tmax", "tmean"]] = np.nan
    raw.loc[outside, "tmean"] = np.nan
    raw.loc[raw.precip < 0, "precip"] = np.nan
    output = []
    for definition, months, temperature in [
        ("annual", list(range(1, 13)), "tmean"),
        ("warm_season", [6, 7, 8, 9], "tmax"),
    ]:
        local = raw.loc[raw.month.isin(months)]
        agg = local.groupby(["station_id", "year"]).agg(
            precip=("precip", lambda s: s.sum(min_count=1)),
            temperature=(temperature, "mean"), precip_days=("precip", "count"),
            temperature_days=(temperature, "count"), observed_rows=("day", "size")).reset_index()
        # Calendar denominators also detect entirely absent daily rows.
        agg["expected_days"] = agg.year.map(lambda y: sum(calendar.monthrange(int(y), m)[1] for m in months))
        agg["coverage"] = agg[["precip_days", "temperature_days"]].min(axis=1) / agg.expected_days
        agg["definition"] = definition
        output.append(agg)
    return pd.concat(output, ignore_index=True), audit


def event_rates(p, t, pcut, tcut, inclusive=False):
    dry = p <= pcut if inclusive else p < pcut
    hot = t >= tcut if inclusive else t > tcut
    return dry.mean(axis=-2), hot.mean(axis=-2), (dry & hot).mean(axis=-2)


def analyze_scenario(aggregates, cfg, scenario, definition, return_station_draws=False):
    years = np.arange(cfg["analysis_years"][0], cfg["analysis_years"][1] + 1)
    eligible = aggregates.loc[(aggregates.definition == definition) & (aggregates.coverage >= scenario["coverage"])]
    p = eligible.pivot(index="year", columns="station_id", values="precip").reindex(years)
    t = eligible.pivot(index="year", columns="station_id", values="temperature").reindex(years)
    keep = p.notna().all() & t.notna().all()
    if scenario.get("common_network", False):
        other = aggregates.loc[aggregates.coverage >= scenario["coverage"]]
        counts = other.groupby(["definition", "station_id"]).year.nunique().unstack(0)
        common_ids = counts.index[(counts == len(years)).all(axis=1)]
        keep &= p.columns.isin(common_ids)
    ids = p.columns[keep].to_numpy()
    p, t = p.loc[:, keep].to_numpy(), t.loc[:, keep].to_numpy()
    if len(ids) < 2:
        raise ValueError(f"Insufficient balanced stations for {definition}/{scenario['name']}")
    early = (years >= cfg["baseline_years"][0]) & (years <= cfg["baseline_years"][1])
    late = (years >= cfg["comparison_years"][0]) & (years <= cfg["comparison_years"][1])
    p0, p1, t0, t1 = p[early], p[late], t[early], t[late]
    qd, qh = scenario["dry_quantile"], scenario["hot_quantile"]
    pc, tc = np.quantile(p0, qd, axis=0), np.quantile(t0, qh, axis=0)
    if scenario.get("exclude_zero_threshold", False):
        selected = pc > 0
        ids, pc, tc = ids[selected], pc[selected], tc[selected]
        p, t = p[:, selected], t[:, selected]
        p0, p1, t0, t1 = p0[:, selected], p1[:, selected], t0[:, selected], t1[:, selected]
    a = event_rates(p0, t0, pc, tc, scenario["inclusive"])
    b = event_rates(p1, t1, pc, tc, scenario["inclusive"])
    point = partition(*a, *b) * 100
    assert np.allclose(point[:, 0], point[:, 1:].sum(axis=1), atol=1e-12)
    station = pd.DataFrame(point, columns=COMPONENTS)
    station["station_id"] = ids
    station["precip_threshold"] = pc
    station["temperature_threshold"] = tc
    station["zero_precip_threshold"] = pc == 0
    for label, rates in [("early", a), ("late", b)]:
        for metric, values in zip(["dry", "hot", "joint"], rates):
            station[f"{metric}_{label}_pct"] = values * 100
    # Seeds shared among sensitivity choices, reducing simulation noise in comparisons.
    rng = np.random.default_rng(cfg["seed"] + (0 if definition == "annual" else 10000))
    reps = int(cfg["bootstrap_reps"])
    i0 = circular_blocks(len(p0), scenario["block_length"], reps, rng)
    i1 = circular_blocks(len(p1), scenario["block_length"], reps, rng)
    draws = np.empty((reps, len(ids), 4))
    for start in range(0, reps, 100):
        end = min(start + 100, reps)
        bp0, bt0, bp1, bt1 = p0[i0[start:end]], t0[i0[start:end]], p1[i1[start:end]], t1[i1[start:end]]
        if scenario["refit_thresholds"]:
            bpc = np.quantile(bp0, qd, axis=1)[:, None, :]
            btc = np.quantile(bt0, qh, axis=1)[:, None, :]
        else:
            bpc, btc = pc, tc
        aa = event_rates(bp0, bt0, bpc, btc, scenario["inclusive"])
        bb = event_rates(bp1, bt1, bpc, btc, scenario["inclusive"])
        draws[start:end] = 100 * partition(*aa, *bb)
    assert np.allclose(draws[..., 0], draws[..., 1:].sum(axis=-1), atol=1e-12)
    network = draws.mean(axis=1)
    summary = []
    for k, component in enumerate(COMPONENTS):
        lo, hi = np.quantile(network[:, k], [0.025, 0.975])
        # Eight primary contrasts: four components x two definitions.
        flo, fhi = np.quantile(network[:, k], [0.05 / 16, 1 - 0.05 / 16])
        summary.append(dict(definition=definition, scenario=scenario["name"], component=component,
                            estimate_pp=point[:, k].mean(), ci_low_pp=lo, ci_high_pp=hi,
                            family_ci_low_pp=flo, family_ci_high_pp=fhi, n_stations=len(ids),
                            joint_early_pct=100*a[2].mean(), joint_late_pct=100*b[2].mean(),
                            zero_precip_threshold_stations=int((pc == 0).sum()), **{k: v for k, v in scenario.items() if k != "name"}))
        station[f"{component}_ci_low"] = np.quantile(draws[:, :, k], .025, axis=0)
        station[f"{component}_ci_high"] = np.quantile(draws[:, :, k], .975, axis=0)
    dry = p <= pc if scenario["inclusive"] else p < pc
    hot = t >= tc if scenario["inclusive"] else t > tc
    series = pd.DataFrame({"year": years, "dry_pct": 100*dry.mean(axis=1),
                           "hot_pct": 100*hot.mean(axis=1), "joint_pct": 100*(dry & hot).mean(axis=1),
                           "n_stations": len(ids), "definition": definition})
    station["definition"] = definition
    station["scenario"] = scenario["name"]
    result = (pd.DataFrame(summary), station, series, network)
    return result + (draws,) if return_station_draws else result


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_partition(cfg, root, config_path=None):
    years = list(range(cfg["analysis_years"][0], cfg["analysis_years"][1] + 1))
    periods = list(range(cfg["baseline_years"][0], cfg["baseline_years"][1] + 1)) + list(range(cfg["comparison_years"][0], cfg["comparison_years"][1] + 1))
    if years != periods or cfg["bootstrap_reps"] < 999:
        raise ValueError("Periods must partition analysis years; use at least 999 replicates")
    names = [s["name"] for s in cfg["scenarios"]]
    if len(set(names)) != len(names) or names.count("primary") != 1:
        raise ValueError("Scenario names must be unique and include primary")
    for scenario in cfg["scenarios"]:
        if not (0 < scenario["dry_quantile"] < .5 < scenario["hot_quantile"] < 1):
            raise ValueError("Dry/hot quantiles must straddle the median")
        if not (0 < scenario["coverage"] <= 1 and 1 <= scenario["block_length"] <= min(len(years)//2, len(years))):
            raise ValueError("Invalid coverage or block length")
    out = root / cfg["output_dir"]
    tables = out / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(root / cfg["data_csv"])
    aggregates, audit = aggregate_daily(raw, cfg["analysis_years"])
    aggregates.to_csv(tables / "screened_annual_seasonal_aggregates.csv", index=False)
    summaries, stations, series, all_stations = [], [], [], []
    for scenario in cfg["scenarios"]:
        for definition in ["annual", "warm_season"]:
            print(f"Partition: {scenario['name']} / {definition}", flush=True)
            summary, station, ts, draws = analyze_scenario(aggregates, cfg, scenario, definition)
            summaries.append(summary)
            all_stations.append(station)
            if scenario["name"] == "primary":
                stations.append(station)
                series.append(ts)
                pd.DataFrame(draws, columns=COMPONENTS).to_csv(tables / f"bootstrap_network_{definition}.csv", index=False)
    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(tables / "compound_partition_sensitivity.csv", index=False)
    summary.loc[summary.scenario == "primary"].to_csv(tables / "compound_partition_primary.csv", index=False)
    pd.concat(stations, ignore_index=True).to_csv(tables / "compound_partition_stations.csv", index=False)
    pd.concat(all_stations, ignore_index=True).to_csv(tables / "compound_partition_all_scenario_stations.csv", index=False)
    pd.concat(series, ignore_index=True).to_csv(tables / "compound_fixed_threshold_extent.csv", index=False)
    inputs = [cfg["data_csv"], cfg["station_csv"], str(config_path or "publication_config.yaml"),
              "src/paper_pipeline/compound_partition.py", "run_publication.py"]
    metadata = {"created_utc": datetime.now(timezone.utc).isoformat(), "config": cfg, "raw_data_audit": audit, "input_sha256": {p: file_hash(root/p) for p in inputs},
                "python": platform.python_version(), "platform": platform.platform(),
                "packages": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scipy", "matplotlib", "statsmodels", "geopandas", "PyYAML"]},
                "scope": "New compound partition rebuilt from raw data; thermal figures use separately verified historical derived tables.",
                "uncertainty": "Circular year blocks synchronized across stations/variables, separately within periods; percentile intervals; thresholds refitted except named sensitivity.",
                "interpretation": "Excess-joint contribution is changing indicator covariance, not isolated copula change or causal attribution."}
    (out / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (out / "environment_versions.txt").write_text("\n".join(f"{d.metadata['Name']}=={d.version}" for d in sorted(importlib.metadata.distributions(), key=lambda d: d.metadata['Name'].lower())), encoding="utf-8")
    return out
