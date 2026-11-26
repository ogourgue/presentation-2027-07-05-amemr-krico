"""
FA plot: where the recruits are when the sea ice returns.

Reads data/aggregated.nc and data/bathymetry.nc and produces
recruit_map.png: one polar-stereographic map of the recruits of both
classes together, by their position on the day sea ice advances over them,
with the outline of the slope where the eggs are spawned and the two
regions where almost no egg gives a recruit.

The recruits are called "larvae ready for winter" on the figure, as in the
abstract.
"""

from __future__ import annotations

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.lines import Line2D


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

mpl.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9,
    "figure.figsize": (7.0, 4.6),
})

INK = "#222222"       # titles and the two regions the abstract explains
INK_2 = "#555555"     # other labels
FRAME = "#999999"
LAND = "#d9d9d9"
COAST = "#8c8c8c"
SLOPE = "#333333"     # outline of the spawning slope


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
AGGREGATED = HERE / "data" / "aggregated.nc"
BATHYMETRY = HERE / "data" / "bathymetry.nc"
OUTPUT = HERE / "recruit_map.png"

# Map extent (lon_min, lon_max, lat_min, lat_max). Narrower than the model
# domain: the western and eastern ends hold almost no recruit.
EXTENT = (-100.0, 28.0, -77.5, -46.5)

# Cartopy projection, as in the figures of krico-paper1.
PROJECTION = ccrs.SouthPolarStereo(central_longitude=-37.5)
DATA_TRANSFORM = ccrs.PlateCarree()

# Spawning slope: sea floor between these depths (m), inside the CCAMLR
# subareas of the project (data/bathymetry.nc is masked to them), where
# particles are released (76°S to 41°S, 114°W to 39°E).
SLOPE_DEPTHS = (1000.0, 2000.0)
RELEASE_LAT = (-76.0, -41.0)
RELEASE_LON = (-114.0, 39.0)

# Colour scale: clipped at this percentile of the non-empty cells, which
# brightens the bulk of the signal and saturates a few hotspot cells.
COLOR_PERCENTILE = 99.0

# Labels: text, lon, lat.
PLACE_LABELS = (
    ("South Orkney\nIslands", -50.5, -56.3),
    ("South Sandwich\nIslands", -14.5, -55.6),
    ("Bellingshausen\nSea", -98.0, -64.5),
)
# Labels with a leader line: text, (lon, lat) of the text, (lon, lat) of
# the point, emphasised or not.
CALLOUTS = (
    ("Antarctic\nPeninsula", (-80.0, -60.5), (-62.5, -65.5), False),
    ("Weddell Sea\ntoo much ice", (-38.0, -68.5), (-52.0, -72.6), True),
    ("South Georgia\ntoo little ice", (-57.0, -50.6), (-38.6, -54.0), True),
)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load_recruits():
    """
    Share of the recruits of both classes in each cell (%), with the cell
    edges, the off-shelf/slope share of the recruits and the share of
    recruits among all particles.
    """
    with xr.open_dataset(AGGREGATED) as ds:
        counts = ds["counts"].values.astype(float)
        classes = [str(c) for c in ds["recruit_class"].values]
        total = ds["total"].values.astype(float)
        lon = ds["lon"].values.astype(float)
        lat = ds["lat"].values.astype(float)
        dataset = ds.attrs.get("dataset", "")
        step = float(ds.attrs["grid_resolution_deg"])

    both = counts.sum(axis=0)
    off = counts[classes.index("recruit_off_shelf_slope")].sum()
    share = np.ma.masked_where(both == 0, 100.0 * both / both.sum())
    lon_edges = np.append(lon - step / 2, lon[-1] + step / 2)
    lat_edges = np.append(lat - step / 2, lat[-1] + step / 2)

    print(f"Input: {dataset}")
    print(f"Recruits, both classes: {both.sum() / total.sum():.2%} of the "
          f"particles")
    print(f"Off-shelf/slope share of the recruits: {off / both.sum():.2%}")
    return share, lon_edges, lat_edges, off / both.sum()


def load_slope():
    """Mask of the spawning slope on the model grid, in map coordinates."""
    with xr.open_dataset(BATHYMETRY) as ds:
        bathymetry = ds["bathymetry_masked"].values.astype(float)
        lon = ds["nav_lon"].values.astype(float)
        lat = ds["nav_lat"].values.astype(float)
    with np.errstate(invalid="ignore"):
        slope = ((bathymetry >= SLOPE_DEPTHS[0]) & (bathymetry <= SLOPE_DEPTHS[1])
                 & (lat >= RELEASE_LAT[0]) & (lat <= RELEASE_LAT[1])
                 & (lon >= RELEASE_LON[0]) & (lon <= RELEASE_LON[1]))
    xy = PROJECTION.transform_points(DATA_TRANSFORM, lon, lat)
    return xy[..., 0], xy[..., 1], slope.astype(float)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def main():
    share, lon_edges, lat_edges, off_share = load_recruits()
    slope_x, slope_y, slope = load_slope()
    vmax = float(np.percentile(share.compressed(), COLOR_PERCENTILE))
    print(f"Colour scale: 0 to {vmax:.3f} % (P{COLOR_PERCENTILE:.0f} of the "
          f"non-empty cells)")

    fig = plt.figure(dpi=300)
    ax = fig.add_axes([0.02, 0.145, 0.96, 0.715], projection=PROJECTION)
    ax.set_extent(EXTENT, crs=DATA_TRANSFORM)
    ax.spines["geo"].set_edgecolor(FRAME)
    ax.spines["geo"].set_linewidth(0.6)

    ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "50m"),
                   facecolor=LAND, edgecolor=COAST, linewidth=0.3, zorder=1)

    mesh = ax.pcolormesh(lon_edges, lat_edges, share, cmap="Greens",
                         vmin=0.0, vmax=vmax, transform=DATA_TRANSFORM,
                         zorder=2, rasterized=True)

    ax.contour(slope_x, slope_y, slope, levels=[0.5], colors=[SLOPE],
               linewidths=0.35, zorder=3)

    for text, lon0, lat0 in PLACE_LABELS:
        ax.text(lon0, lat0, text, transform=DATA_TRANSFORM, ha="center",
                va="center", fontsize=7.5, color=INK_2, linespacing=1.15,
                zorder=6)
    geodetic = DATA_TRANSFORM._as_mpl_transform(ax)
    for text, xy_text, xy_point, emphasised in CALLOUTS:
        color = INK if emphasised else INK_2
        ax.annotate(
            text, xy=xy_point, xytext=xy_text, xycoords=geodetic,
            textcoords=geodetic, ha="center", va="center",
            fontsize=8.5 if emphasised else 7.5, color=color,
            fontweight="bold" if emphasised else "normal", linespacing=1.15,
            zorder=7,
            arrowprops=dict(arrowstyle="-", color=color, linewidth=0.6,
                            shrinkA=2, shrinkB=1),
        )

    # Title and subtitle.
    fig.text(0.02, 0.955,
             "Where krill larvae ready for winter are when the sea ice returns",
             fontsize=11.5, fontweight="bold", color=INK, ha="left",
             va="center")
    fig.text(0.02, 0.895,
             f"Eggs are spawned on the slope. Of the larvae ready for winter, "
             f"{100 * off_share:.0f}% are then over deep water, off the slope.",
             fontsize=8.5, color=INK_2, ha="left", va="center")

    # Colour bar, key of the slope outline and period.
    fig.text(0.02, 0.118, "Larvae ready for winter, % of their total per 1° cell",
             fontsize=7.5, color=INK_2, ha="left", va="center")
    cax = fig.add_axes([0.02, 0.075, 0.30, 0.022])
    cbar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="max")
    cbar.outline.set_linewidth(0.4)
    cbar.outline.set_edgecolor(FRAME)
    cbar.ax.tick_params(labelsize=7.5, length=2, width=0.4, colors=INK_2)
    cbar.ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=3))
    fig.legend(
        handles=[Line2D([0], [0], color=SLOPE, linewidth=0.9)],
        labels=["Outline of the slope where eggs are spawned (1000 to 2000 m)"],
        loc="lower left", bbox_to_anchor=(0.36, 0.052), frameon=False,
        fontsize=7.5, labelcolor=INK_2, handlelength=1.6,
    )
    fig.text(0.98, 0.118, "32 seasons, 1993/94 to 2024/25", fontsize=7.5,
             color=INK_2, ha="right", va="center")

    fig.savefig(OUTPUT, dpi=300, facecolor="white")
    print(f"Written {OUTPUT}")


if __name__ == "__main__":
    main()
