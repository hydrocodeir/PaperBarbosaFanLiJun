"""Build two uncluttered 60 x 50 mm graphical-abstract concepts."""

from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "submission_IJC_2026"
TABLES = ROOT / "outputs" / "publication_v2" / "tables"

NAVY = "#15364D"
DEEP = "#102B3D"
INK = "#182A35"
MUTED = "#63747E"
IVORY = "#F7F5F0"
WHITE = "#FFFFFF"
WARM = "#D84E3F"
WARM2 = "#EC9A32"
COOL = "#2F7FA3"
COOL2 = "#5368AD"
DRY = "#C88A24"
EXCESS = "#20808A"
PALE = "#DFE7EA"


def panel(fig, xywh, fc, ec="none", radius=0.018, lw=0.5, z=-3):
    x, y, w, h = xywh
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.005,rounding_size={radius}",
        transform=fig.transFigure, facecolor=fc, edgecolor=ec,
        linewidth=lw, zorder=z,
    )
    fig.add_artist(patch)
    return patch


def read_data():
    stations = pd.read_csv(ROOT / "data" / "stationsInfo.csv")
    countries = gpd.read_file(ROOT / "data" / "Iran_Sea_Ne.geojson").to_crs(4326)
    profiles = pd.read_csv(TABLES / "thermal_network_profiles.csv")
    primary = pd.read_csv(TABLES / "compound_partition_primary.csv")
    station_parts = pd.read_csv(TABLES / "compound_partition_stations.csv")
    summer_map = station_parts.loc[
        station_parts.definition.eq("warm_season") & station_parts.scenario.eq("primary")
    ].merge(stations, on="station_id", validate="one_to_one")
    return stations, countries, profiles, primary, summer_map


def thermal_plot(fig, profiles, dark=False):
    fg = WHITE if dark else INK
    grid = "#446172" if dark else "#D9E2E5"
    muted = "#B9CAD3" if dark else MUTED
    ax = fig.add_axes([0.645, 0.595, 0.315, 0.125], facecolor="none")
    order = ["warm_days", "warm_nights", "cool_days", "cool_nights"]
    names = ["Warm days", "Warm nights", "Cool days", "Cool nights"]
    colors = [WARM, WARM2, COOL, COOL2]
    ypos = np.arange(4)[::-1]
    for idx, color, y in zip(order, colors, ypos):
        loc = profiles.loc[(profiles.index_name == idx) & profiles.network.eq("common_fixed")]
        vals = {round(float(r.tau), 2): float(r.slope) for r in loc.itertuples()}
        lo, mid, hi = vals[0.1], vals[0.5], vals[0.9]
        ax.plot([min(lo, hi), max(lo, hi)], [y, y], color=color, lw=2.4,
                solid_capstyle="round", zorder=2)
        ax.scatter([lo, hi], [y, y], s=5, color=color, zorder=3)
        ax.scatter([mid], [y], s=16, color=color, edgecolor=WHITE,
                   linewidth=0.35, zorder=4)
    ax.axvline(0, color=muted, lw=0.6)
    ax.set_xlim(-20, 25)
    ax.set_ylim(-0.5, 3.5)
    ax.set_yticks(ypos, names)
    ax.set_xticks([-20, 0, 20])
    ax.tick_params(axis="both", labelsize=3.75, length=1.6, pad=1.1, colors=muted)
    ax.grid(axis="x", color=grid, lw=0.35, linestyle=":")
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(grid)
    ax.set_xlabel("q10–q90 slopes  (days decade⁻¹)", fontsize=3.7, color=muted, labelpad=1)
    for tick, color in zip(ax.get_yticklabels(), colors):
        tick.set_color(color)
        tick.set_fontweight("bold")


def map_plot(fig, countries, summer_map, dark=True):
    ax = fig.add_axes([0.055, 0.235, 0.405, 0.425], facecolor="none")
    context = "#284455" if dark else "#E4E7E6"
    iran_fill = "#F5F3EC" if dark else WHITE
    edge = "#AFC0C9" if dark else "#68777F"
    countries.plot(ax=ax, color=context, edgecolor=context, linewidth=0.22, zorder=1)
    countries.loc[countries.cca2.str.lower().eq("ir")].plot(
        ax=ax, color=iran_fill, edgecolor=edge, linewidth=0.55, zorder=2
    )
    cmap = LinearSegmentedColormap.from_list("joint", ["#277DA1", "#F3EFE3", "#D74D3F"])
    norm = TwoSlopeNorm(vmin=-40, vcenter=0, vmax=40)
    ok = summer_map.loc[~summer_map.zero_precip_threshold]
    zero = summer_map.loc[summer_map.zero_precip_threshold]
    dots = ax.scatter(ok.longitude, ok.latitude, c=ok.joint_change,
                      cmap=cmap, norm=norm, s=8.5, edgecolor=DEEP,
                      linewidth=0.18, zorder=4)
    ax.scatter(zero.longitude, zero.latitude, marker="x", s=8,
               color="#70818A", linewidth=0.55, zorder=5)
    ax.set_xlim(43.2, 64)
    ax.set_ylim(24, 40.5)
    ax.set_aspect(1 / np.cos(np.deg2rad(32)))
    ax.axis("off")

    cax = fig.add_axes([0.095, 0.210, 0.325, 0.012])
    cb = fig.colorbar(dots, cax=cax, orientation="horizontal", ticks=[-40, 0, 40])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=3.2, length=1.5, pad=0.5,
                      colors="#C6D4DA" if dark else MUTED)
    cb.ax.set_xticklabels(["−40", "0", "+40"])


def save_all(fig, stem):
    svg = OUT / f"{stem}.svg"
    png = OUT / f"{stem}_60x50mm_300dpi.png"
    tif300 = OUT / f"{stem}_60x50mm_300dpi.tiff"
    tif600 = OUT / f"{stem}_60x50mm_600dpi.tiff"
    fig.savefig(svg, format="svg", dpi=300)
    fig.savefig(png, format="png", dpi=300)
    fig.savefig(tif300, format="tiff", dpi=300,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(tif600, format="tiff", dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return svg, png, tif300, tif600


def build_variant(stem, dark_right=False):
    _, countries, profiles, primary, summer_map = read_data()

    width_in, height_in = 60 / 25.4, 50 / 25.4
    bg = DEEP if dark_right else IVORY
    fig = plt.figure(figsize=(width_in, height_in), dpi=300, facecolor=bg)

    if dark_right:
        title_color, body, muted, rule = WHITE, WHITE, "#BDD0D9", "#3B5666"
        map_fc = "#17384C"
        accent_fc = "#F2C46D"
        accent_text = DEEP
    else:
        title_color, body, muted, rule = NAVY, INK, MUTED, "#D8E0E3"
        map_fc = DEEP
        accent_fc = "#F1C56E"
        accent_text = DEEP

    # Header: title and authors are deliberately isolated from the figure.
    fig.text(0.04, 0.970,
             "Distributional thermal change and the\n"
             "components of increasing dry–hot\n"
             "concurrence across Iran, 1991–2024",
             ha="left", va="top", color=title_color, fontsize=5.05,
             fontweight="bold", linespacing=1.03)
    fig.text(0.04, 0.823, "[AUTHOR NAMES; *corresponding author]",
             ha="left", va="top", color=muted, fontsize=3.7)
    fig.add_artist(Rectangle((0.04, 0.795), 0.92, 0.0025,
                             transform=fig.transFigure, facecolor=rule, edgecolor="none"))

    # Dominant scientific visual: station-level summer change across Iran.
    panel(fig, (0.03, 0.17, 0.445, 0.595), map_fc, radius=0.022)
    fig.text(0.055, 0.738, "SUMMER DRY–HOT CHANGE",
             color=WHITE, fontsize=4.65, fontweight="bold", va="top")
    fig.text(0.055, 0.707, "103 stations  •  change (pp)  •  × zero-rain cutoff",
             color="#C4D3DA", fontsize=3.35, va="top")
    map_plot(fig, countries, summer_map, dark=True)

    # Right-hand narrative: three analytical steps, no nested cards.
    x0 = 0.515
    fig.text(x0, 0.758, "1", color=WARM, fontsize=7.2, fontweight="bold", va="top")
    fig.text(x0 + 0.035, 0.756, "THERMAL COUNTS", color=body,
             fontsize=4.8, fontweight="bold", va="top")
    thermal_plot(fig, profiles, dark=dark_right)

    fig.add_artist(Rectangle((x0, 0.555), 0.445, 0.0018,
                             transform=fig.transFigure, facecolor=rule, edgecolor="none"))
    fig.text(x0, 0.530, "2", color=WARM, fontsize=7.2, fontweight="bold", va="top")
    fig.text(x0 + 0.035, 0.528, "DRY–HOT YEARS", color=body,
             fontsize=4.8, fontweight="bold", va="top")

    annual = primary.loc[primary.definition.eq("annual") & primary.component.eq("joint_change")].iloc[0]
    summer = primary.loc[primary.definition.eq("warm_season") & primary.component.eq("joint_change")].iloc[0]
    fig.text(x0, 0.465, "JUN–SEP", color=muted, fontsize=3.8, fontweight="bold", va="center")
    fig.text(x0 + 0.105, 0.465,
             f"{summer.joint_early_pct:.1f}%  →  {summer.joint_late_pct:.1f}%",
             color=WARM, fontsize=6.35, fontweight="bold", va="center")
    fig.text(x0, 0.415, "ANNUAL", color=muted, fontsize=3.8, fontweight="bold", va="center")
    fig.text(x0 + 0.105, 0.415,
             f"{annual.joint_early_pct:.1f}%  →  {annual.joint_late_pct:.1f}%",
             color=body, fontsize=5.5, fontweight="bold", va="center")

    fig.add_artist(Rectangle((x0, 0.372), 0.445, 0.0018,
                             transform=fig.transFigure, facecolor=rule, edgecolor="none"))
    fig.text(x0, 0.347, "3", color=WARM, fontsize=7.2, fontweight="bold", va="top")
    fig.text(x0 + 0.035, 0.345, "SUMMER PARTITION", color=body,
             fontsize=4.8, fontweight="bold", va="top")

    parts = primary.loc[primary.definition.eq("warm_season")].set_index("component")
    hot = float(parts.loc["hot_marginal", "estimate_pp"])
    dry = float(parts.loc["dry_marginal", "estimate_pp"])
    excess = float(parts.loc["excess_joint", "estimate_pp"])
    ax_bar = fig.add_axes([x0, 0.277, 0.445, 0.035], facecolor="none")
    ax_bar.barh([0], [hot], color=WARM, height=0.58)
    ax_bar.barh([0], [dry], left=[hot], color=DRY, height=0.58)
    ax_bar.barh([0], [excess], left=[hot + dry], color=EXCESS, height=0.58)
    ax_bar.set_xlim(0, hot + dry + excess)
    ax_bar.axis("off")
    fig.text(x0, 0.247, "HOT  +9.91 pp", color=WARM, fontsize=3.9,
             fontweight="bold", va="center")
    fig.text(x0 + 0.20, 0.247, "DRY  +2.89 pp", color=DRY, fontsize=3.9,
             fontweight="bold", va="center")
    fig.text(x0, 0.215, "EXCESS-JOINT  +0.06 pp*", color=EXCESS, fontsize=3.45,
             fontweight="bold", va="center")

    # Single, reviewer-safe conclusion.
    panel(fig, (0.03, 0.035, 0.93, 0.10), accent_fc, radius=0.017)
    fig.text(0.055, 0.108,
             "MORE DRY–HOT YEARS—MAINLY FROM MORE HEAT.",
             color=accent_text, fontsize=4.8, fontweight="bold", va="top")
    fig.text(0.055, 0.068,
             "*Thermal asymmetry and excess-joint change were not resolved.",
             color=accent_text, fontsize=3.75, va="top")

    return save_all(fig, stem)


if __name__ == "__main__":
    for args in [("Graphical_Abstract", False),
                 ("Graphical_Abstract_Alternative_Dark", True)]:
        paths = build_variant(*args)
        for path in paths:
            print(f"Wrote {path}")
