"""Build the data-faithful 50 x 60 mm graphical abstract for IJC submission."""

from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "submission_IJC_2026"

NAVY = "#15324B"
INK = "#1D2A35"
MUTED = "#5F6E78"
PAPER = "#F3F6F8"
WHITE = "#FFFFFF"
GRID = "#DCE4E8"
WARM = "#D05243"
WARM2 = "#E6962F"
COOL = "#2883A8"
COOL2 = "#5167B1"
DRY = "#C18A2A"
EXCESS = "#21818A"


def rounded_panel(fig, xywh, radius=0.018, fc=WHITE, ec="#D9E1E5", lw=0.55):
    x, y, w, h = xywh
    panel = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.005,rounding_size={radius}",
        transform=fig.transFigure,
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        zorder=-2,
    )
    fig.add_artist(panel)
    return panel


def build():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 5,
            "axes.linewidth": 0.45,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": PAPER,
        }
    )

    stations = pd.read_csv(ROOT / "data" / "stationsInfo.csv")
    countries = gpd.read_file(ROOT / "data" / "Iran_Sea_Ne.geojson").to_crs(4326)
    profiles = pd.read_csv(OUT.parent / "outputs" / "publication_v2" / "tables" / "thermal_network_profiles.csv")
    partition = pd.read_csv(OUT.parent / "outputs" / "publication_v2" / "tables" / "compound_partition_primary.csv")

    # Exact journal display footprint: 50 mm x 60 mm.
    width_in, height_in = 50 / 25.4, 60 / 25.4
    fig = plt.figure(figsize=(width_in, height_in), dpi=300, facecolor=PAPER)

    # Header
    fig.add_artist(Rectangle((0, 0.835), 1, 0.165, transform=fig.transFigure,
                             facecolor=NAVY, edgecolor="none", zorder=-3))
    fig.text(
        0.5, 0.970,
        "Distributional thermal change and the\n"
        "components of increasing dry–hot\n"
        "concurrence across Iran, 1991–2024",
        ha="center", va="top", color=WHITE, fontsize=4.9, fontweight="bold",
        linespacing=1.08,
    )
    fig.text(
        0.5, 0.846, "[AUTHOR NAMES; *corresponding author]",
        ha="center", va="center", color="#D9E7EF", fontsize=3.8,
    )

    # Panel 1: observation network map
    rounded_panel(fig, (0.045, 0.565, 0.43, 0.245))
    fig.text(0.065, 0.787, "OBSERVATIONS", color=NAVY, fontsize=5.0,
             fontweight="bold", va="top")
    fig.text(0.065, 0.762, "124 stations  •  1991–2024", color=MUTED,
             fontsize=4.35, va="top")
    ax_map = fig.add_axes([0.064, 0.582, 0.392, 0.164])
    ax_map.set_facecolor("#EAF2F5")
    countries.plot(ax=ax_map, color="#ECE9E1", edgecolor="#A6AFB2", linewidth=0.24, zorder=1)
    iran = countries.loc[countries.cca2.str.lower() == "ir"]
    iran.plot(ax=ax_map, color=WHITE, edgecolor="#47555D", linewidth=0.52, zorder=2)
    ax_map.scatter(
        stations.longitude, stations.latitude,
        c=stations.elevation, cmap="cividis", s=4.2,
        edgecolor=WHITE, linewidth=0.13, zorder=3,
    )
    ax_map.set_xlim(43.2, 64)
    ax_map.set_ylim(24.0, 40.5)
    ax_map.set_aspect(1 / np.cos(np.deg2rad(32)))
    ax_map.axis("off")
    ax_map.text(53.6, 31.8, "IRAN", ha="center", va="center", color="#87939A",
                fontsize=4.0, fontweight="bold", zorder=2)

    # Panel 2: distributional thermal signals
    rounded_panel(fig, (0.50, 0.565, 0.455, 0.245))
    fig.text(0.52, 0.787, "THERMAL DISTRIBUTION", color=NAVY, fontsize=4.45,
             fontweight="bold", va="top")
    fig.text(0.52, 0.762, "q10–q90 slopes  •  fixed n = 108", color=MUTED,
             fontsize=3.8, va="top")
    ax_q = fig.add_axes([0.64, 0.592, 0.288, 0.142])
    order = ["warm_days", "warm_nights", "cool_days", "cool_nights"]
    names = ["Warm days", "Warm nights", "Cool days", "Cool nights"]
    colors = [WARM, WARM2, COOL, COOL2]
    ypos = np.arange(4)[::-1]
    for idx, name, color, y in zip(order, names, colors, ypos):
        loc = profiles.loc[(profiles.index_name == idx) & (profiles.network == "common_fixed")]
        vals = {round(float(row.tau), 2): float(row.slope) for row in loc.itertuples()}
        lo, mid, hi = vals[0.1], vals[0.5], vals[0.9]
        ax_q.plot([min(lo, hi), max(lo, hi)], [y, y], color=color, lw=2.0,
                  solid_capstyle="round", zorder=2)
        ax_q.scatter([lo, hi], [y, y], s=3.5, color=color, zorder=3)
        ax_q.scatter([mid], [y], s=13, color=color, edgecolor=WHITE, linewidth=0.35, zorder=4)
    ax_q.axvline(0, color="#95A2A9", lw=0.5, zorder=0)
    ax_q.set_xlim(-20, 25)
    ax_q.set_ylim(-0.55, 3.55)
    ax_q.set_yticks(ypos, names)
    ax_q.set_xticks([-20, 0, 20])
    ax_q.tick_params(axis="both", labelsize=3.8, length=1.8, pad=1.2, colors=MUTED)
    ax_q.grid(axis="x", color=GRID, lw=0.35, linestyle=":")
    ax_q.spines[["top", "right", "left"]].set_visible(False)
    ax_q.spines["bottom"].set_color("#9AA6AC")
    ax_q.set_xlabel("trend (days decade⁻¹)", fontsize=3.75, color=MUTED, labelpad=1.2)
    for label, color in zip(ax_q.get_yticklabels(), colors):
        label.set_color(color)
        label.set_fontweight("bold")

    # Panel 3: observed dry-hot expansion
    fig.text(0.05, 0.535, "DRY–HOT FREQUENCY", color=NAVY, fontsize=5.0,
             fontweight="bold", va="top")
    annual = partition.loc[(partition.definition == "annual") & (partition.component == "joint_change")].iloc[0]
    summer = partition.loc[(partition.definition == "warm_season") & (partition.component == "joint_change")].iloc[0]
    for x, label, row, n in [
        (0.045, "ANNUAL", annual, 104),
        (0.505, "JUNE–SEPTEMBER", summer, 103),
    ]:
        rounded_panel(fig, (x, 0.405, 0.45, 0.107), radius=0.015)
        fig.text(x + 0.02, 0.487, label, color=MUTED, fontsize=4.2,
                 fontweight="bold", va="top")
        fig.text(
            x + 0.225, 0.453,
            f"{row.joint_early_pct:.2f}%  →  {row.joint_late_pct:.2f}%",
            ha="center", va="center", color=INK, fontsize=5.75, fontweight="bold",
        )
        fig.text(x + 0.225, 0.419, f"+{row.estimate_pp:.2f} pp  •  n = {n}",
                 ha="center", va="center", color=WARM, fontsize=4.35, fontweight="bold")

    # Panel 4: summer partition
    rounded_panel(fig, (0.045, 0.215, 0.91, 0.158))
    fig.text(0.065, 0.350, "WHAT DROVE SUMMER CHANGE?", color=NAVY,
             fontsize=4.45, fontweight="bold", va="top")
    fig.text(0.935, 0.350, "+12.85 pp", color=WARM, fontsize=4.45,
             fontweight="bold", va="top", ha="right")
    parts = partition.loc[partition.definition.eq("warm_season")].set_index("component")
    hot = float(parts.loc["hot_marginal", "estimate_pp"])
    dry = float(parts.loc["dry_marginal", "estimate_pp"])
    excess = float(parts.loc["excess_joint", "estimate_pp"])
    ax_bar = fig.add_axes([0.07, 0.277, 0.86, 0.042])
    ax_bar.barh([0], [hot], left=[0], color=WARM, height=0.48)
    ax_bar.barh([0], [dry], left=[hot], color=DRY, height=0.48)
    ax_bar.barh([0], [excess], left=[hot + dry], color=EXCESS, height=0.48)
    ax_bar.set_xlim(0, hot + dry + excess)
    ax_bar.set_ylim(-0.6, 0.6)
    ax_bar.axis("off")
    ax_bar.text(hot / 2, 0, "+9.91", color=WHITE, fontsize=4.6,
                fontweight="bold", ha="center", va="center")
    ax_bar.text(hot + dry / 2, 0, "+2.89", color=WHITE, fontsize=4.2,
                fontweight="bold", ha="center", va="center")
    legend_y = 0.244
    fig.text(0.07, legend_y, "●", color=WARM, fontsize=5.0, va="center")
    fig.text(0.092, legend_y, "Hot frequency  +9.91 pp", color=INK, fontsize=3.75, va="center")
    fig.text(0.51, legend_y, "●", color=DRY, fontsize=5.0, va="center")
    fig.text(0.532, legend_y, "Dry frequency  +2.89 pp", color=INK, fontsize=3.75, va="center")
    fig.text(0.07, 0.222, "●", color=EXCESS, fontsize=5.0, va="center")
    fig.text(0.092, 0.222, "Excess joint  +0.06 pp  (95% interval crosses zero)*",
             color=INK, fontsize=3.55, va="center")

    # Take-home message: two sentences, calibrated to the evidence.
    rounded_panel(fig, (0.045, 0.045, 0.91, 0.135), radius=0.018, fc=NAVY, ec=NAVY)
    fig.text(0.07, 0.145, "TAKE-HOME", color="#8FD0D2", fontsize=4.55,
             fontweight="bold", va="top")
    fig.text(0.07, 0.116,
             "Warming-consistent thermal counts and\n"
             "dry–hot concurrence increased.",
             color=WHITE, fontsize=4.55, fontweight="bold", va="top", linespacing=1.05)
    fig.text(0.07, 0.067,
             "Asymmetry and excess-joint change remain unresolved.*",
             color="#D9E7EF", fontsize=3.95, va="top")

    svg_path = OUT / "Graphical_Abstract.svg"
    png_path = OUT / "Graphical_Abstract_50x60mm.png"
    png_300_path = OUT / "Graphical_Abstract_50x60mm_300dpi.png"
    tiff_300_path = OUT / "Graphical_Abstract_50x60mm_300dpi.tiff"
    tiff_600_path = OUT / "Graphical_Abstract_50x60mm_600dpi.tiff"

    fig.savefig(svg_path, format="svg", dpi=300, facecolor=PAPER)
    fig.savefig(png_path, format="png", dpi=300, facecolor=PAPER)
    fig.savefig(png_300_path, format="png", dpi=300, facecolor=PAPER)
    fig.savefig(tiff_300_path, format="tiff", dpi=300, facecolor=PAPER,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(tiff_600_path, format="tiff", dpi=600, facecolor=PAPER,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    print(f"Wrote {svg_path}")
    print(f"Wrote {png_300_path}")
    print(f"Wrote {tiff_300_path}")
    print(f"Wrote {tiff_600_path}")


if __name__ == "__main__":
    build()
