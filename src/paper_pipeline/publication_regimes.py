"""Climate-regime summaries with synchronized, station-first uncertainty."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from .compound_partition import analyze_scenario, COMPONENTS, file_hash

REGIMES = ["BWh_hot_desert", "BWk_cold_desert", "BSh_hot_steppe", "BSk_cold_steppe", "C_temperate", "Dsa_cold_dry_summer"]
LABELS = ["BWh hot desert", "BWk cold desert", "BSh hot steppe", "BSk cold steppe", "Csa/Cfa temperate", "Dsa cold dry-summer"]
COLORS = ["#C96239", "#CBAD65", "#CC7BA6", "#7666AA", "#258D76", "#327CAC"]


def build_regime_tables(root, out, cfg):
    assignments = pd.read_csv(root / "outputs/tables/koppen_geiger_station_assignments.csv")
    assert len(assignments) == 124 and assignments.station_id.is_unique
    # Exact CSV float round trips preserve strict comparisons at tied thresholds.
    aggregates = pd.read_csv(out / "tables/screened_annual_seasonal_aggregates.csv", float_precision="round_trip")
    primary = next(s for s in cfg["scenarios"] if s["name"] == "primary")
    summaries, draw_rows = [], []
    validation = {}
    for definition in ["annual", "warm_season"]:
        _, stations, _, network, draws = analyze_scenario(aggregates, cfg, primary, definition, return_station_draws=True)
        existing = pd.read_csv(out / f"tables/bootstrap_network_{definition}.csv", float_precision="round_trip").to_numpy()
        assert np.allclose(network, existing, atol=1e-10)
        validation[f"{definition}_bootstrap_matches_prior_network"] = True
        regimes = assignments.set_index("station_id").loc[stations.station_id, "climate_regime"].to_numpy()
        for regime in REGIMES:
            mask = regimes == regime
            local = stations.loc[mask]
            ensemble = draws[:, mask, :].mean(axis=1)
            for k, component in enumerate(COMPONENTS):
                lo, hi = np.quantile(ensemble[:, k], [.025, .975])
                summaries.append(dict(definition=definition, climate_regime=regime, component=component,
                    n_stations=int(mask.sum()), n_zero_dry_threshold=int(local.zero_precip_threshold.sum()),
                    joint_early_pct=local.joint_early_pct.mean(), joint_late_pct=local.joint_late_pct.mean(),
                    estimate_pp=local[component].mean(), ci_low_pp=lo, ci_high_pp=hi))
            draw_rows.append(pd.DataFrame(ensemble, columns=COMPONENTS).assign(definition=definition, climate_regime=regime, replicate=np.arange(len(ensemble))))
    summary = pd.DataFrame(summaries)
    summary.to_csv(out / "tables/climate_regime_compound_partition.csv", index=False)
    pd.concat(draw_rows, ignore_index=True).to_csv(out / "tables/climate_regime_compound_bootstrap.csv", index=False)
    climate = pd.read_csv(root / "outputs/tables/climate_regime_quantile_summary.csv")
    meta = pd.read_csv(root / "outputs/tables/koppen_geiger_regime_summary.csv").set_index("climate_regime")
    fixed = pd.read_csv(root / "outputs/tables/climate_regime_fixed_baseline_summary.csv")
    q90 = climate.pivot(index="climate_regime", columns="index_name", values="mean_slope_0.90")
    thermal = meta[["n_stations", "mean_elevation_m", "kg_classes"]].join(q90).loc[REGIMES]
    thermal.to_csv(out / "tables/table02_climate_regime_thermal.csv")
    summer = summary.loc[summary.definition == "warm_season"]
    # Independently reproduce group means from station slopes.
    station_qr = pd.read_csv(root / "outputs/tables/qr_focus_slopes_and_bootstrap_summary.csv").merge(assignments[["station_id", "climate_regime"]], on="station_id", validate="many_to_one")
    check = station_qr.groupby(["climate_regime", "index_name"])["slope_0.90"].mean().unstack()
    assert np.allclose(check.loc[REGIMES, q90.columns], q90.loc[REGIMES])
    validation["thermal_regime_means_reproduced_from_station_slopes"] = True
    validation["scope"] = "Within-regime intervals are exploratory; no between-regime significance claim or causal attribution. Synchronized primary draws, with threshold re-estimation."
    validation["source_hashes"] = {str(p.relative_to(root)):file_hash(p) for p in [root / "src/paper_pipeline/publication_regimes.py",root / "outputs/tables/koppen_geiger_station_assignments.csv",root / "outputs/tables/qr_focus_slopes_and_bootstrap_summary.csv",out / "tables/screened_annual_seasonal_aggregates.csv"]}
    (out / "regime_validation.json").write_text(json.dumps(validation,indent=2),encoding="utf-8")
    print(summary.loc[summary.component == "joint_change"].round(3).to_string(index=False),flush=True)


def plot_regimes(root, out, countries, records, save, map_base):
    stations = pd.read_csv(root / "outputs/tables/koppen_geiger_station_assignments.csv")
    thermal = pd.read_csv(out / "tables/table02_climate_regime_thermal.csv").set_index("climate_regime").loc[REGIMES]
    fig = plt.figure(figsize=(7.1, 5.7), layout="constrained")
    gs = fig.add_gridspec(2,2,height_ratios=[1.8,1])
    ax = fig.add_subplot(gs[0,:])
    map_base(ax,countries,"(a) Station climate classification · Beck et al. (2018)")
    for regime,label,color in zip(REGIMES,LABELS,COLORS):
        local=stations.loc[stations.climate_regime == regime]
        ax.scatter(local.longitude,local.latitude,c=color,s=16,edgecolor="white",linewidth=.25,zorder=4,label=f"{label} (n={len(local)})")
    ax.legend(loc="center left",bbox_to_anchor=(1.02,.5),frameon=False,fontsize=6.5)
    heat=fig.add_subplot(gs[1,:])
    values=thermal[["warm_days","warm_nights","cool_days","cool_nights"]].to_numpy()
    im=heat.imshow(values,aspect="auto",cmap="RdBu_r",norm=TwoSlopeNorm(vmin=-25,vcenter=0,vmax=25))
    heat.set(xticks=range(4),xticklabels=["Warm days","Warm nights","Cool days","Cool nights"],yticks=range(6),yticklabels=LABELS,title="(b) Mean station upper-quantile slopes")
    for y in range(6):
        for x in range(4):
            heat.text(x,y,f"{values[y,x]:.1f}",ha="center",va="center",fontsize=7,color="white" if abs(values[y,x])>14 else "#202528")
    fig.colorbar(im,ax=heat,fraction=.026,pad=.02,label="Days per decade")
    save(fig,out,"fig05_climate_regimes",records)
    summary=pd.read_csv(out / "tables/climate_regime_compound_partition.csv")
    fig,axes=plt.subplots(1,2,figsize=(7.1,3.9),sharex=True,sharey=True,layout="constrained")
    for k,(definition,ax) in enumerate(zip(["annual","warm_season"],axes)):
        local=summary.loc[(summary.definition==definition)&(summary.component=="joint_change")].set_index("climate_regime").loc[REGIMES]
        for y,((_,row),color) in enumerate(zip(local.iterrows(),COLORS)):
            ax.plot([row.ci_low_pp,row.ci_high_pp],[y,y],color=color,lw=1.6)
            ax.scatter(row.estimate_pp,y,color=color,s=23)
            ax.text(.98,y,f"n={row.n_stations:d}" if isinstance(row.n_stations,int) else f"n={int(row.n_stations)}",transform=ax.get_yaxis_transform(),ha="right",va="center",fontsize=6,color="#50575B")
        ax.axvline(0,color="gray",lw=.6)
        ax.set_xlim(-3,70)
        ax.set_xticks([0,20,40,60])
        ax.set(yticks=range(6),yticklabels=LABELS,title=f"({'ab'[k]}) {'Annual' if k==0 else 'June–September'}")
    axes[0].invert_yaxis()
    fig.supxlabel("Joint-frequency change (percentage points); 95% within-regime intervals",fontsize=8)
    save(fig,out,"fig09_compound_climate_regimes",records)
