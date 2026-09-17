"""Scientific regression checks for the new partition and historical inputs."""
from pathlib import Path
import json
import sys
import copy
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from paper_pipeline.compound_partition import aggregate_daily, partition, circular_blocks, event_rates
from paper_pipeline.indices import create_extreme_indices


def main():
    checks = {}
    # Known independent marginals: exact product change .06 -> .20.
    independent = partition(.2, .3, .06, .4, .5, .20)
    assert np.allclose(independent, [.14, .08, .06, 0])
    checks["known_independent_probability_partition"] = "PASS"
    rng = np.random.default_rng(42)
    for _ in range(100):
        p0, p1 = rng.dirichlet(np.ones(4), size=2)
        values = partition(p0[1]+p0[3], p0[2]+p0[3], p0[3], p1[1]+p1[3], p1[2]+p1[3], p1[3])
        assert abs(values[0] - values[1:].sum()) < 1e-12
    checks["random_valid_joint_distributions_close_exactly"] = "PASS"
    idx = circular_blocks(17, 4, 100, rng)
    assert idx.shape == (100, 17) and idx.min() >= 0 and idx.max() < 17
    assert np.all((np.diff(idx[:, :4], axis=1) % 17) == 1)
    checks["circular_block_lengths_and_wrap"] = "PASS"
    dry, hot, joint = event_rates(np.zeros((17, 2)), np.ones((17, 2)), np.zeros(2), np.zeros(2))
    assert np.all(dry == 0) and np.all(hot == 1) and np.all(joint == 0)
    dry, _, joint = event_rates(np.zeros((17, 2)), np.ones((17, 2)), np.zeros(2), np.zeros(2), inclusive=True)
    assert np.all(dry == 1) and np.all(joint == 1)
    checks["zero_rainfall_strict_and_inclusive_thresholds"] = "PASS"
    # Entirely absent rows must still count against calendar coverage.
    tiny = pd.DataFrame(dict(station_id=[1], year=[2000], month=[6], day=[1],
                             tmin=[10.], tmax=[30.], tmean=[20.], precip=[1.]))
    agg, _ = aggregate_daily(tiny, [2000, 2000])
    assert np.isclose(agg.loc[agg.definition == "annual", "coverage"].iloc[0], 1/366)
    assert np.isclose(agg.loc[agg.definition == "warm_season", "coverage"].iloc[0], 1/122)
    checks["calendar_denominator_for_absent_days"] = "PASS"
    table = ROOT / "outputs/publication_v2/tables"
    stations = pd.read_csv(table / "compound_partition_stations.csv")
    assert np.allclose(stations.joint_change, stations[["dry_marginal", "hot_marginal", "excess_joint"]].sum(axis=1))
    assert ((stations.joint_late_pct <= stations.dry_late_pct+1e-10) & (stations.joint_late_pct <= stations.hot_late_pct+1e-10)).all()
    checks["real_data_joint_event_subset_and_partition"] = "PASS"
    # Independently reconstruct all historical temperature indices from raw data.
    cfg = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    raw = pd.read_csv(ROOT / "data/data.csv")
    raw = raw.loc[raw.year.between(*cfg["data"]["analysis_years"])]
    for tag, reference, filename in [("full_record", None, "annual_extreme_indices.csv"),
                                     ("fixed_baseline", [1991, 2007], "fixed_baseline_annual_extreme_indices.csv")]:
        local_cfg = copy.deepcopy(cfg)
        local_cfg["index_construction"]["reference_years"] = reference
        print(f"Rebuild temperature indices: {tag}", flush=True)
        _, rebuilt = create_extreme_indices(raw, local_cfg)
        historical = pd.read_csv(ROOT / "outputs/tables" / filename)
        cols = ["warm_days", "warm_nights", "cool_days", "cool_nights"]
        left = rebuilt.set_index(["station_id", "year"])[cols].sort_index()
        right = historical.set_index(["station_id", "year"])[cols].sort_index()
        pd.testing.assert_index_equal(left.index, right.index)
        assert np.allclose(left.to_numpy(), right.to_numpy(), equal_nan=True), tag
        checks[f"raw_rebuild_{tag}_all_station_year_indices"] = "PASS"
    # Full-population and station-level estimates are not interchangeable.
    profiles = pd.read_csv(table / "network_quantile_profiles_recomputed.csv")
    baseline = pd.read_csv(ROOT / "outputs/tables/homogeneity_flag_exclusion_sensitivity.csv").set_index("index_name")
    for _, row in profiles.loc[profiles.tau.round(2).isin([.1, .5, .9])].iterrows():
        assert abs(row.network_slope - baseline.loc[row.index_name, f"slope_{row.tau:.2f}_all_stations"]) < .01
    checks["recomputed_network_focal_quantiles_match_manuscript"] = "PASS"
    result = {"checks": checks, "scope": "Complete raw reconstruction of both thermal index sets; new compound extension and figures. Historical station bootstrap/clustering full pipeline NOT rerun."}
    (ROOT / "outputs/publication_v2/validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)

if __name__ == "__main__":
    main()
