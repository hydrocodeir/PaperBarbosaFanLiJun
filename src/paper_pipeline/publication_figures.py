"""Journal-sized, vector-first figures; station maps never interpolate data."""
from pathlib import Path
import json

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import statsmodels.api as sm

from .compound_partition import COMPONENTS, file_hash

INDEXES = ["warm_days", "warm_nights", "cool_days", "cool_nights"]
NAMES = ["Warm days", "Warm nights", "Cool days", "Cool nights"]
COLORS = ["#B74335", "#D18B26", "#277CA6", "#465AA3"]
COMP_NAMES = ["Joint-frequency change", "Dry-frequency term", "Hot-frequency term", "Excess-joint term"]
COMP_COLORS = ["#242A30", "#AD7B27", "#C4483B", "#277C87"]


def theme():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
        "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 7, "axes.linewidth": .6,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "axes.unicode_minus": True})


def save(fig, out, name, records):
    folder = out / "figures"
    folder.mkdir(exist_ok=True)
    for fmt in ["png", "pdf", "svg", "tiff"]:
        kwargs = {"pil_kwargs": {"compression": "tiff_lzw"}} if fmt == "tiff" else {}
        fig.savefig(folder / f"{name}.{fmt}", dpi=600 if fmt == "tiff" else 350,
                    bbox_inches="tight", pad_inches=.06, **kwargs)
    records.append(dict(figure=name, width_inches=fig.get_figwidth(), height_inches=fig.get_figheight(),
                        formats="PDF; SVG; PNG 350 dpi; TIFF 600 dpi LZW"))
    plt.close(fig)


def map_base(ax, countries, title):
    ax.set_facecolor("#EEF5F7")
    countries.plot(ax=ax, color="#F1EFE9", edgecolor="#9CA0A0", linewidth=.35, zorder=1)
    countries.loc[countries.cca2.str.lower() == "ir"].plot(ax=ax, color="#FFFFFF", edgecolor="#454B4F", linewidth=.7, zorder=2)
    ax.set(xlim=(43.2, 64), ylim=(24, 40.5), title=title)
    ax.set_xticks([45, 50, 55, 60])
    ax.set_yticks([25, 30, 35, 40])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.0f}°E"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.0f}°N"))
    ax.set_aspect(1 / np.cos(np.deg2rad(32)))
    ax.grid(color="#D6DEDF", linewidth=.3, linestyle=":", zorder=0)


def create_figures(root: Path, out: Path):
    theme()
    records = []
    table = out / "tables"
    historical = root / "outputs/tables"
    stations = pd.read_csv(root / "data/stationsInfo.csv")
    countries = gpd.read_file(root / "data/Iran_Sea_Ne.geojson").to_crs(4326)
    annual = pd.read_csv(historical / "annual_extreme_indices.csv")
    qr = pd.read_csv(historical / "qr_focus_slopes_and_bootstrap_summary.csv")
    station_partition = pd.read_csv(table / "compound_partition_stations.csv").merge(stations, on="station_id", validate="many_to_one")
    summary = pd.read_csv(table / "compound_partition_primary.csv")

    fig, ax = plt.subplots(figsize=(7.1, 5.4), layout="constrained")
    map_base(ax, countries, "124-station observing network · 1991–2024")
    dots = ax.scatter(stations.longitude, stations.latitude, c=stations.elevation, cmap="cividis",
                      s=23, edgecolor="white", linewidth=.3, zorder=4)
    for name, xy in [("Caspian Sea", (51, 38.5)), ("Persian Gulf", (51.2, 26.3)),
                     ("Gulf of Oman", (59, 24.5)), ("Iraq", (44.5, 33)),
                     ("Afghanistan", (62, 34)), ("Turkmenistan", (58, 39))]:
        ax.text(*xy, name, fontsize=7, color="#4D5C61", ha="center", zorder=5)
    fig.colorbar(dots, ax=ax, fraction=.034, pad=.025, label="Elevation (m)")
    save(fig, out, "fig01_station_network", records)

    # Explicitly distinguish quantiles of the network mean from mean station slopes.
    network = annual.groupby("year")[INDEXES].mean()
    x = (network.index.to_numpy() - network.index.min()) / 10
    design = sm.add_constant(x)
    all_qr = pd.read_csv(historical / "qr_all_quantiles_long.csv")
    taus = sorted(all_qr.tau.unique())
    numerical = []
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.4), sharex=True, sharey=True, layout="constrained")
    for k, (idx, name, color, ax) in enumerate(zip(INDEXES, NAMES, COLORS, axes.flat)):
        slopes = [float(sm.QuantReg(network[idx].to_numpy(), design).fit(q=tau, max_iter=50000).params[1]) for tau in taus]
        mean_slope = float(sm.OLS(network[idx].to_numpy(), design).fit().params[1])
        station_profile = all_qr.loc[all_qr.index_name == idx].groupby("tau").slope
        lo, hi = station_profile.quantile(.25).reindex(taus), station_profile.quantile(.75).reindex(taus)
        ax.fill_between(taus, lo, hi, color=color, alpha=.13, linewidth=0, label="Station-slope IQR")
        ax.plot(taus, slopes, color=color, lw=1.8, label="Network-mean QR")
        ax.axhline(mean_slope, color="#4D5559", lw=1, ls="--", label="Network-mean OLS")
        ax.axhline(0, color="#8B9498", lw=.6)
        ax.set(title=f"({chr(97+k)}) {name}", xticks=[.1, .3, .5, .7, .9])
        for tau, slope in zip(taus, slopes):
            numerical.append(dict(index_name=idx, tau=tau, network_slope=slope, ols_slope=mean_slope))
    axes[0, 0].legend(loc="upper left", frameon=False)
    fig.supxlabel("Conditional quantile of annual event counts", fontsize=8)
    fig.supylabel("Trend (days per decade)", fontsize=8)
    pd.DataFrame(numerical).to_csv(table / "network_quantile_profiles_recomputed.csv", index=False)
    save(fig, out, "fig02_quantile_profiles", records)

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 6.4), layout="constrained")
    norm = TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60)
    for k, (idx, name, ax) in enumerate(zip(INDEXES, NAMES, axes.flat)):
        local = qr.loc[qr.index_name == idx].merge(stations, on="station_id")
        map_base(ax, countries, f"({chr(97+k)}) {name}")
        dots = ax.scatter(local.longitude, local.latitude, c=local.Delta1, s=17, cmap="RdBu_r", norm=norm,
                          edgecolor="#4C5054", linewidth=.2, zorder=4)
    fig.colorbar(dots, ax=list(axes.flat), orientation="horizontal", fraction=.045, pad=.02,
                 extend="both", label="Upper minus lower quantile slope, Δ₁ (days per decade)")
    save(fig, out, "fig03_thermal_asymmetry_maps", records)

    fixed = pd.read_csv(historical / "fixed_baseline_period_change_station_level.csv")
    fig, ax = plt.subplots(figsize=(7.1, 3.1), layout="constrained")
    rng = np.random.default_rng(42)
    for k, (idx, color) in enumerate(zip(INDEXES, COLORS)):
        values = fixed.loc[fixed.index_name == idx, "late_minus_baseline_days"].dropna().to_numpy()
        ax.scatter(values, k + rng.uniform(-.18, .18, len(values)), s=8, color=color, alpha=.4, linewidths=0)
        q1, median, q3 = np.quantile(values, [.25, .5, .75])
        ax.plot([q1, q3], [k, k], color="#242A30", lw=2)
        ax.scatter([median], [k], marker="|", s=90, color="#242A30", zorder=4)
        ax.scatter([values.mean()], [k], marker="D", s=25, facecolor="white", edgecolor="black", zorder=5)
    ax.axvline(0, lw=.7, color="#7A8285")
    ax.set(yticks=range(4), yticklabels=NAMES, xlabel="Late minus early annual count (days per year)")
    ax.invert_yaxis()
    ax.legend(handles=[Line2D([], [], marker="D", markerfacecolor="white", color="black", linestyle="", label="Network mean"),
                       Line2D([], [], color="black", lw=2, marker="|", label="Station median and IQR")], frameon=False, loc="lower right")
    save(fig, out, "fig04_fixed_baseline_changes", records)

    ts = pd.read_csv(table / "compound_fixed_threshold_extent.csv")
    fig, axes = plt.subplots(2, 1, figsize=(7.1, 4.8), sharex=True, sharey=True, layout="constrained")
    for k, (definition, ax) in enumerate(zip(["annual", "warm_season"], axes)):
        local = ts.loc[ts.definition == definition]
        ax.axvspan(1991, 2007.5, color="#ECEDE9", zorder=0)
        for col, name, color, ls in [("dry_pct", "Dry", COMP_COLORS[1], ":"),
                                     ("hot_pct", "Hot", COMP_COLORS[2], "--"),
                                     ("joint_pct", "Dry AND hot", COMP_COLORS[0], "-")]:
            ax.plot(local.year, local[col], lw=1.3, label=name, color=color, ls=ls)
        ax.set(title=f"({chr(97+k)}) {'Annual' if definition == 'annual' else 'June–September'} · n = {local.n_stations.iloc[0]}",
               ylim=(0, 102), ylabel="Stations (%)", xlim=(1991, 2024))
    axes[0].legend(frameon=False, ncol=3, loc="upper left")
    axes[-1].set_xlabel("Year")
    save(fig, out, "fig06_fixed_compound_extent", records)

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.5), sharex=True, sharey=True, layout="constrained")
    for k, (definition, ax) in enumerate(zip(["annual", "warm_season"], axes)):
        local = summary.loc[summary.definition == definition].set_index("component").loc[list(COMPONENTS)]
        for y, (_, row) in enumerate(local.iterrows()):
            ax.plot([row.family_ci_low_pp, row.family_ci_high_pp], [y, y], color=COMP_COLORS[y], lw=.6)
            ax.plot([row.ci_low_pp, row.ci_high_pp], [y, y], color=COMP_COLORS[y], lw=2.6)
            ax.scatter(row.estimate_pp, y, color=COMP_COLORS[y], s=30, zorder=3)
        ax.axvline(0, color="#8F9698", lw=.7)
        ax.set(title=f"({chr(97+k)}) {'Annual' if k == 0 else 'June–September'}", yticks=range(4),
               yticklabels=COMP_NAMES, xlabel="Frequency change (percentage points)")
    axes[0].invert_yaxis()
    save(fig, out, "fig07_compound_partition", records)

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 6.4), layout="constrained")
    local = station_partition.loc[station_partition.definition == "warm_season"]
    norm = TwoSlopeNorm(vmin=-40, vcenter=0, vmax=40)
    for k, (component, name, ax) in enumerate(zip(COMPONENTS, COMP_NAMES, axes.flat)):
        map_base(ax, countries, f"({chr(97+k)}) {name}")
        eligible = local.loc[~local.zero_precip_threshold]
        zero = local.loc[local.zero_precip_threshold]
        dots = ax.scatter(eligible.longitude, eligible.latitude, c=eligible[component], s=20, norm=norm, cmap="RdBu_r",
                          edgecolor="#454B50", linewidth=.25, zorder=4)
        ax.scatter(zero.longitude, zero.latitude, marker="x", s=13, color="#787F82", linewidth=.6, zorder=4,
                   label=f"Zero dry threshold (n = {len(zero)})")
        if k == 0:
            ax.legend(loc="lower left", frameon=True, facecolor="white", edgecolor="none", fontsize=6)
    fig.colorbar(dots, ax=list(axes.flat), orientation="horizontal", fraction=.045, pad=.02,
                 extend="both", label="June–September: late minus early frequency (percentage points)")
    save(fig, out, "fig08_compound_partition_maps", records)

    sensitivity = pd.read_csv(table / "compound_partition_sensitivity.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 4.1), sharex=True, sharey=True, layout="constrained")
    for k, (definition, ax) in enumerate(zip(["annual", "warm_season"], axes)):
        local = sensitivity.loc[(sensitivity.definition == definition) & (sensitivity.component == "joint_change")]
        for y, (_, row) in enumerate(local.iterrows()):
            ax.plot([row.ci_low_pp, row.ci_high_pp], [y, y], color="#277C87", lw=1.3)
            ax.scatter(row.estimate_pp, y, color="#277C87", s=18)
        ax.axvline(0, color="gray", lw=.6)
        ax.set(yticks=range(len(local)), yticklabels=[x.replace("_", " ") for x in local.scenario],
               title=f"({'ab'[k]}) {'Annual' if k == 0 else 'June–September'}", xlabel="Joint-frequency change (pp)")
    axes[0].invert_yaxis()
    save(fig, out, "fig10_partition_sensitivity", records)

    from .publication_regimes import plot_regimes
    plot_regimes(root, out, countries, records, save, map_base)
    records.sort(key=lambda item: item["figure"])
    pd.DataFrame(records).to_csv(out / "figure_manifest.csv", index=False)
    sources = ["outputs/tables/annual_extreme_indices.csv", "outputs/tables/qr_all_quantiles_long.csv",
               "outputs/tables/qr_focus_slopes_and_bootstrap_summary.csv", "outputs/tables/fixed_baseline_period_change_station_level.csv",
               "outputs/tables/climate_regime_quantile_summary.csv", "data/Iran_Sea_Ne.geojson",
               "src/paper_pipeline/publication_figures.py", "src/paper_pipeline/publication_regimes.py",
               "outputs/tables/koppen_geiger_station_assignments.csv"]
    sources += [str(p.relative_to(root)) for p in sorted(table.glob("*.csv"))]
    (out / "figure_source_hashes.json").write_text(json.dumps({p: file_hash(root/p) for p in sources}, indent=2), encoding="utf-8")
    print(f"Created {len(records)} figures in four formats", flush=True)
