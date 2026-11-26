"""
FA plot: where the recruits come from, and where they meet the ice.

Reads data/aggregated.nc, data/bathymetry.nc and the CCAMLR shapefile, and
produces recruit_map.png: one polar-stereographic map with two readings.

  - Numbers: where the recruits come from. For each region, the share of
    the eggs spawned there that give a recruit, as quoted in the abstract.
  - Colour: where the recruits meet the ice. Recruits of both classes
    together, by their position on the day sea ice advances over them. The
    colour scale carries no numbers, so that every number on the figure is
    a share of eggs.

The recruits are called "larvae ready for winter" on the figure, as in the
abstract. The regions are the CCAMLR subareas of the project, with 48.6
split at 60°S.
"""

from __future__ import annotations

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, TextArea, VPacker
from matplotlib.patheffects import withStroke
from shapely.geometry import LineString
from shapely.ops import unary_union


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "figure.figsize": (7.0, 4.6),
})

INK = "#222222"       # titles and numbers
INK_2 = "#555555"     # other text
FRAME = "#999999"
LAND = "#d9d9d9"
COAST = "#8c8c8c"
SLOPE = "#333333"     # outline of the spawning slope
SUBAREA = "#b3b3b3"   # outlines of the regions


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
AGGREGATED = HERE / "data" / "aggregated.nc"
BATHYMETRY = HERE / "data" / "bathymetry.nc"
SHAPEFILE = HERE.parent / "shared" / "ccamlr-data" / "CCAMLR_ASD_EPSG4326.shp"
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

# Regions outlined on the map: the CCAMLR subareas of the project, with
# 48.6 split at 60°S into Bouvet (north) and Dronning Maud Land (south).
SHAPE_FIELD = "GAR_Long_L"
SUBAREA_CODES = ("48.1", "48.2", "48.3", "48.4", "48.5", "48.6", "88.3")
SPLIT_LAT_486 = -60.0
N_DIVIDER_POINTS = 200

# Colour scale: clipped at this percentile of the non-empty cells, which
# brightens the bulk of the signal and saturates a few hotspot cells.
COLOR_PERCENTILE = 99.0

# Region labels: subarea, name, reason given in the abstract (or None),
# (lon, lat) of the label, (lon, lat) of the point its line leads to (or
# None for no line). The number is computed from the data.
REGION_LABELS = (
    ("88.3", "Bellingshausen Sea", None, (-99.5, -63.2), None),
    ("48.1", "Antarctic Peninsula", None, (-80.0, -60.0), (-62.5, -65.5)),
    ("48.2", "South Orkney Islands", None, (-53.5, -56.4), None),
    ("48.4", "South Sandwich Islands", None, (-15.5, -54.6), None),
    ("48.5", "Weddell Sea", "too much ice", (-38.0, -68.3), (-52.0, -72.6)),
    ("48.3", "South Georgia", "too little ice", (-58.0, -50.4), (-38.6, -54.0)),
    ("48.6N", "Bouvet", None, (11.5, -56.3), None),
    ("48.6S", "Dronning Maud Land", None, (14.0, -64.5), None),
)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load_recruits():
    """
    Share of the recruits of both classes in each cell, with the cell
    edges; the off-shelf/slope share of the recruits; and, per subarea of
    release, the share of the particles that are recruits.
    """
    with xr.open_dataset(AGGREGATED) as ds:
        counts = ds["counts"].values.astype(float)
        classes = [str(c) for c in ds["recruit_class"].values]
        total = ds["total"].values.astype(float)
        by_subarea = ds["counts_release"].values.astype(float)
        subareas = [str(s) for s in ds["subarea"].values]
        outcomes = [str(o) for o in ds["outcome"].values]
        lon = ds["lon"].values.astype(float)
        lat = ds["lat"].values.astype(float)
        dataset = ds.attrs.get("dataset", "")
        step = float(ds.attrs["grid_resolution_deg"])

    both = counts.sum(axis=0)
    off = counts[classes.index("recruit_off_shelf_slope")].sum()
    share = np.ma.masked_where(both == 0, both / both.sum())
    lon_edges = np.append(lon - step / 2, lon[-1] + step / 2)
    lat_edges = np.append(lat - step / 2, lat[-1] + step / 2)

    recruit = [outcomes.index(c) for c in classes]
    potential = {}
    print(f"Input: {dataset}")
    print(f"Recruits, both classes: {both.sum() / total.sum():.2%} of the "
          f"particles")
    print(f"Off-shelf/slope share of the recruits: {off / both.sum():.2%}")
    print("Share of the particles released in each subarea that are recruits:")
    for k, name in enumerate(subareas):
        released = by_subarea[k].sum()
        if released:
            potential[name] = by_subarea[k, recruit].sum() / released
            print(f"  {name:8s} {potential[name]:7.2%}")
    return share, lon_edges, lat_edges, off / both.sum(), potential


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


def split_486(ccamlr):
    """The 60°S parallel inside subarea 48.6, densified to follow the map."""
    area = unary_union(ccamlr[ccamlr[SHAPE_FIELD] == "48.6"].geometry.values)
    lon_min, _, lon_max, _ = area.bounds
    lons = np.linspace(lon_min, lon_max, N_DIVIDER_POINTS)
    return LineString([(lon, SPLIT_LAT_486) for lon in lons]).intersection(area)


def format_share(value: float) -> str:
    """A share of eggs as in the abstract: whole percent, or 'below 1%'."""
    return "below 1%" if value < 0.01 else f"{100 * value:.0f}%"


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def lead_and_text(fig, x, y, lead, text, size):
    """One line of text whose first words are bold."""
    first = fig.text(x, y, lead, fontsize=size, fontweight="bold", color=INK,
                     ha="left", va="center")
    fig.canvas.draw()
    box = first.get_window_extent().transformed(fig.transFigure.inverted())
    fig.text(box.x1 + 0.006, y, text, fontsize=size, color=INK_2, ha="left",
             va="center")


def region_label(ax, name, number, reason, xy_label, xy_point):
    """Name, number and reason of a region, with an optional leader line."""
    halo = [withStroke(linewidth=2.0, foreground="white")]
    lines = [
        TextArea(name, textprops=dict(fontsize=7.5, color=INK_2,
                                      path_effects=halo)),
        TextArea(number, textprops=dict(fontsize=10.5, color=INK,
                                        fontweight="bold", path_effects=halo)),
    ]
    if reason:
        lines.append(TextArea(reason, textprops=dict(
            fontsize=8.5, color=INK, fontweight="bold", path_effects=halo)))
    box = VPacker(children=lines, align="center", pad=0, sep=1.5)
    geodetic = DATA_TRANSFORM._as_mpl_transform(ax)
    arrow = None
    if xy_point is not None:
        arrow = dict(arrowstyle="-", color=INK, linewidth=0.6, shrinkA=2,
                     shrinkB=1)
    ax.add_artist(AnnotationBbox(
        box, xy=xy_point if xy_point is not None else xy_label,
        xybox=xy_label, xycoords=geodetic, boxcoords=geodetic,
        frameon=False, pad=0.0, arrowprops=arrow, zorder=7))


def main():
    share, lon_edges, lat_edges, off_share, potential = load_recruits()
    slope_x, slope_y, slope = load_slope()
    vmax = float(np.percentile(share.compressed(), COLOR_PERCENTILE))

    fig = plt.figure(dpi=300)
    ax = fig.add_axes([0.02, 0.105, 0.96, 0.715], projection=PROJECTION)
    ax.set_extent(EXTENT, crs=DATA_TRANSFORM)
    ax.spines["geo"].set_edgecolor(FRAME)
    ax.spines["geo"].set_linewidth(0.6)

    ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "50m"),
                   facecolor=LAND, edgecolor=COAST, linewidth=0.3, zorder=1)

    mesh = ax.pcolormesh(lon_edges, lat_edges, share, cmap="Greens",
                         vmin=0.0, vmax=vmax, transform=DATA_TRANSFORM,
                         zorder=2, rasterized=True)

    ccamlr = gpd.read_file(SHAPEFILE)
    ccamlr = ccamlr[ccamlr[SHAPE_FIELD].isin(SUBAREA_CODES)]
    ax.add_geometries(ccamlr.geometry, crs=DATA_TRANSFORM, facecolor="none",
                      edgecolor=SUBAREA, linewidth=0.5, zorder=3)
    ax.add_geometries([split_486(ccamlr)], crs=DATA_TRANSFORM,
                      facecolor="none", edgecolor=SUBAREA, linewidth=0.5,
                      zorder=3)

    ax.contour(slope_x, slope_y, slope, levels=[0.5], colors=[SLOPE],
               linewidths=0.35, zorder=4)

    for subarea, name, reason, xy_label, xy_point in REGION_LABELS:
        region_label(ax, name, format_share(potential[subarea]), reason,
                     xy_label, xy_point)

    # Title and the two readings of the map.
    fig.text(0.02, 0.955,
             "Krill larvae ready for winter: where they come from, "
             "where they meet the ice",
             fontsize=11.5, fontweight="bold", color=INK, ha="left",
             va="center")
    lead_and_text(fig, 0.02, 0.895, "Numbers: where they come from.",
                  "Share of the eggs spawned in each region that give a larva "
                  "ready for winter.", 8.5)
    lead_and_text(fig, 0.02, 0.853, "Colour: where they meet the ice.",
                  f"{100 * off_share:.0f}% of them are then over deep water, "
                  f"off the slope where eggs are spawned.", 8.5)

    # Colour bar without numbers, keys of the outlines and period.
    cax = fig.add_axes([0.068, 0.048, 0.15, 0.022])
    cbar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="max")
    cbar.outline.set_linewidth(0.4)
    cbar.outline.set_edgecolor(FRAME)
    cbar.set_ticks([])
    fig.text(0.061, 0.059, "fewer", fontsize=7.5, color=INK_2, ha="right",
             va="center")
    fig.text(0.232, 0.059, "more larvae", fontsize=7.5, color=INK_2,
             ha="left", va="center")
    fig.legend(
        handles=[Line2D([0], [0], color=SLOPE, linewidth=0.9),
                 Line2D([0], [0], color=SUBAREA, linewidth=1.2)],
        labels=["Slope where eggs are spawned (1000–2000 m)",
                "Regions (CCAMLR subareas)"],
        loc="center left", bbox_to_anchor=(0.345, 0.059), frameon=False,
        fontsize=7.5, labelcolor=INK_2, handlelength=1.6, ncols=2,
        columnspacing=1.6, borderpad=0, borderaxespad=0,
    )
    ax.text(0.012, 0.03, "32 seasons, 1993/94 to 2024/25", fontsize=7.5,
            color=INK_2, ha="left", va="center", transform=ax.transAxes,
            zorder=6)

    fig.savefig(OUTPUT, dpi=300, facecolor="white")
    print(f"Written {OUTPUT}")


if __name__ == "__main__":
    main()
