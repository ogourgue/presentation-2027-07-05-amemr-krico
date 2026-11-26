"""
FA plot: where the recruits come from, and where they meet the ice.

Reads data/aggregated.nc, data/bathymetry.nc and the CCAMLR shapefile, and
produces recruit_map.png: one polar-stereographic map with two readings.

  - Numbers: where the recruits come from. For each region, the share of
    the eggs spawned there that give a recruit, as quoted in the abstract.
  - Colour: where the recruits meet the ice. Recruits of both classes
    together, by their position on the day sea ice advances over them. The
    colour scale carries no numbers, so that every number on the figure is
    a share of eggs. Its low end is transparent, so that the shelf and
    the slope, in two shades of grey, show through where the recruits are
    few. The slope is outlined again above the colours, and the land, in a
    third grey, is drawn above everything.

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
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, TextArea, VPacker
from matplotlib.patches import Circle, Patch
from matplotlib.patheffects import withStroke
from shapely.geometry import LineString
from shapely.ops import linemerge, unary_union


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
LAND = "#b8b8b8"
SHELF = "#e4e4e4"     # sea floor shallower than the slope
SLOPE = "#7d7d7d"     # slope where the eggs are spawned
SUBAREA = "#9a9a9a"   # outlines of the regions


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
# subareas of the project (data/bathymetry.nc is masked to them). Particles
# are released over all of it except 96 model cells south of 76°S (0.3% of
# its cells), which are shaded like the rest. The shelf is the sea floor
# shallower than the slope, inside the same subareas.
SLOPE_DEPTHS = (1000.0, 2000.0)

# Regions outlined on the map: the CCAMLR subareas of the project, with
# 48.6 split at 60°S into Bouvet (north) and Dronning Maud Land (south).
SHAPE_FIELD = "GAR_Long_L"
SUBAREA_CODES = ("48.1", "48.2", "48.3", "48.4", "48.5", "48.6", "88.3")
SPLIT_LAT_486 = -60.0
N_DIVIDER_POINTS = 200

# A segment of a region outline is a border, and not coast, when it runs
# along a meridian or a parallel: one coordinate constant within
# BORDER_TOLERANCE (degrees) over at least BORDER_MIN_LENGTH (degrees).
BORDER_TOLERANCE = 1e-7
BORDER_MIN_LENGTH = 0.01

# Colour scale: clipped at this percentile of the non-empty cells, which
# brightens the bulk of the signal and saturates a few hotspot cells.
COLOR_PERCENTILE = 99.0

# Colour map: the greens of matplotlib, transparent at the low end and
# opaque from OPAQUE_FROM (fraction of the scale) upward.
OPAQUE_FROM = 0.40

# Width in pixels of the image of the recruits (see draw_recruits); above
# the width of the map in the figure, so that the cell edges stay sharp.
IMAGE_WIDTH = 4000

# Region labels: subarea, name, reason given in the abstract (or None),
# (lon, lat) of the label, and (lon, lat) of the point its line leads to.
# A label sits inside its region and has no line (None). Where the region
# is too crowded the label sits outside, and its line ends on a border
# that the region shares with no other. The number is computed from the
# data.
REGION_LABELS = (
    ("88.3", "Bellingshausen Sea", None, (-99.5, -63.2), None),
    ("48.1", "Antarctic Peninsula", None, (-73.5, -57.3), (-66.0, -60.0)),
    ("48.2", "South Orkney Islands", None, (-57.5, -56.3), (-50.0, -58.4)),
    ("48.4", "South\nSandwich\nIslands", None, (-25.0, -52.9), None),
    ("48.5", "Weddell Sea", "too much ice", (-38.0, -68.3), None),
    ("48.3", "South Georgia", "too little ice", (-59.0, -50.2), (-50.0, -51.6)),
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


def load_bathymetry():
    """
    Masks of the shelf and of the spawning slope on the model grid, in map
    coordinates, and the radius of the southern edge of the grid. Both
    masks exist inside the regions only.
    """
    with xr.open_dataset(BATHYMETRY) as ds:
        bathymetry = ds["bathymetry_masked"].values.astype(float)
        lon = ds["nav_lon"].values.astype(float)
        lat = ds["nav_lat"].values.astype(float)
    with np.errstate(invalid="ignore"):
        shelf = bathymetry < SLOPE_DEPTHS[0]
        slope = (bathymetry >= SLOPE_DEPTHS[0]) & (bathymetry <= SLOPE_DEPTHS[1])
    xy = PROJECTION.transform_points(DATA_TRANSFORM, lon, lat)
    # Distance from the pole, in map coordinates, of the second southernmost
    # row of the grid: the shelf shade south of the grid overlaps it by one
    # row, so that no gap shows between the two.
    south_radius = float(np.hypot(*PROJECTION.transform_point(
        0.0, float(lat[1].max()), DATA_TRANSFORM)))
    return (xy[..., 0], xy[..., 1], shelf.astype(float), slope.astype(float),
            south_radius)


def green_colormap():
    """
    Greens with a transparent low end.

    Over white the colours are exactly those of the matplotlib map: each
    one is replaced by a more saturated green with an opacity below 1 that
    gives the same colour once blended with white. Over the greys the low
    end therefore tints what lies beneath in place of hiding it. The
    opacity rises linearly to 1 at OPAQUE_FROM, and is never lower than
    what the colour needs.
    """
    x = np.linspace(0.0, 1.0, 256)
    target = mpl.colormaps["Greens"](x)[:, :3]
    alpha = np.maximum(np.clip(x / OPAQUE_FROM, 0.0, 1.0),
                       1.0 - target.min(axis=1))
    colors = np.ones((x.size, 4))
    colors[:, :3] = 1.0 - (1.0 - target) / alpha[:, None]
    colors[:, 3] = alpha
    return ListedColormap(np.clip(colors, 0.0, 1.0), name="greens_transparent")


def region_borders(ccamlr):
    """
    The borders of the regions that run along a meridian or a parallel.

    The CCAMLR polygons end at the coast of Antarctica, so their full
    outline would also trace that coast, and no other. Only the straight
    borders between and around the regions are kept.
    """
    lines = []
    for geometry in ccamlr.geometry:
        polygons = getattr(geometry, "geoms", [geometry])
        for polygon in polygons:
            xy = np.asarray(polygon.exterior.coords)
            step = np.abs(np.diff(xy, axis=0))
            straight = (step.min(axis=1) < BORDER_TOLERANCE) \
                & (step.max(axis=1) > BORDER_MIN_LENGTH)
            lines += [LineString([xy[k], xy[k + 1]])
                      for k in np.flatnonzero(straight)]
    return linemerge(lines)


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

def draw_recruits(ax, share, lon_edges, lat_edges, vmax):
    """
    Draw the share of recruits per cell as one image in map coordinates,
    and return the mappable for the colour bar.

    One image and not one patch per cell: with transparent colours, patches
    leave their edges visible as a grid. Each pixel takes the colour of the
    lon/lat cell that holds its centre.
    """
    cmap = green_colormap()
    norm = mpl.colors.Normalize(vmin=0.0, vmax=vmax)
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    nx = IMAGE_WIDTH
    ny = int(round(nx * (y1 - y0) / (x1 - x0)))
    x = x0 + (np.arange(nx) + 0.5) * (x1 - x0) / nx
    y = y0 + (np.arange(ny) + 0.5) * (y1 - y0) / ny
    lonlat = DATA_TRANSFORM.transform_points(PROJECTION, *np.meshgrid(x, y))
    j = np.searchsorted(lon_edges, lonlat[..., 0], side="right") - 1
    i = np.searchsorted(lat_edges, lonlat[..., 1], side="right") - 1
    inside = ((i >= 0) & (i < share.shape[0])
              & (j >= 0) & (j < share.shape[1]))
    values = np.ma.masked_all(inside.shape)
    values[inside] = share[i[inside], j[inside]]
    rgba = cmap(norm(values))
    rgba[np.ma.getmaskarray(values)] = 0.0
    ax.imshow(rgba, extent=(x0, x1, y0, y1), origin="lower",
              interpolation="nearest", zorder=2)
    return mpl.cm.ScalarMappable(norm=norm, cmap=cmap)


def lead_and_text(fig, x, y, lead, text, size):
    """One line of text whose first words are bold."""
    first = fig.text(x, y, lead, fontsize=size, fontweight="bold", color=INK,
                     ha="left", va="center")
    fig.canvas.draw()
    box = first.get_window_extent().transformed(fig.transFigure.inverted())
    fig.text(box.x1 + 0.006, y, text, fontsize=size, color=INK_2, ha="left",
             va="center")


def region_label(ax, name, number, reason, xy_label, xy_point):
    """
    Name, number and reason of a region; with a line to xy_point, a point
    on the border of the region, when the label sits outside it.
    """
    halo = [withStroke(linewidth=2.0, foreground="white")]
    lines = [
        TextArea(name, textprops=dict(fontsize=7.5, color=INK_2,
                                      path_effects=halo, linespacing=1.1,
                                      multialignment="center")),
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
                     shrinkB=0)
    ax.add_artist(AnnotationBbox(
        box, xy=xy_point if xy_point is not None else xy_label,
        xybox=xy_label, xycoords=geodetic, boxcoords=geodetic,
        frameon=False, pad=0.0, arrowprops=arrow, zorder=7))


def main():
    share, lon_edges, lat_edges, off_share, potential = load_recruits()
    bathy_x, bathy_y, shelf, slope, south_radius = load_bathymetry()
    vmax = float(np.percentile(share.compressed(), COLOR_PERCENTILE))

    fig = plt.figure(dpi=300)
    ax = fig.add_axes([0.02, 0.105, 0.96, 0.715], projection=PROJECTION)
    ax.set_extent(EXTENT, crs=DATA_TRANSFORM)
    ax.spines["geo"].set_edgecolor(FRAME)
    ax.spines["geo"].set_linewidth(0.6)

    # Three shades of grey: the shelf and the slope under the recruits, the
    # land above them, so that the islands stay visible where the recruits
    # are many. South of the bathymetry grid everything is drawn as shelf.
    ax.add_patch(Circle((0.0, 0.0), south_radius, facecolor=SHELF,
                        edgecolor="none", transform=ax.transData, zorder=0.9))
    ax.contourf(bathy_x, bathy_y, shelf, levels=[0.5, 1.5], colors=[SHELF],
                zorder=1.0)
    ax.contourf(bathy_x, bathy_y, slope, levels=[0.5, 1.5], colors=[SLOPE],
                zorder=1.2)
    mesh = draw_recruits(ax, share, lon_edges, lat_edges, vmax)

    # The slope again, as a thin outline above the recruits: it keeps the
    # relief readable where the colours are opaque.
    ax.contour(bathy_x, bathy_y, slope, levels=[0.5], colors=[SLOPE],
               linewidths=0.3, zorder=2.2)

    ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "50m"),
                   facecolor=LAND, edgecolor=LAND, linewidth=0.2, zorder=2.5)

    ccamlr = gpd.read_file(SHAPEFILE)
    ccamlr = ccamlr[ccamlr[SHAPE_FIELD].isin(SUBAREA_CODES)]
    ax.add_geometries([region_borders(ccamlr)], crs=DATA_TRANSFORM,
                      facecolor="none", edgecolor=SUBAREA, linewidth=0.5,
                      zorder=3)
    ax.add_geometries([split_486(ccamlr)], crs=DATA_TRANSFORM,
                      facecolor="none", edgecolor=SUBAREA, linewidth=0.5,
                      zorder=3)

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
        handles=[Patch(facecolor=SLOPE, edgecolor="none"),
                 Patch(facecolor=SHELF, edgecolor="none"),
                 Patch(facecolor=LAND, edgecolor="none"),
                 Line2D([0], [0], color=SUBAREA, linewidth=1.0)],
        labels=["Slope, the assumed spawning zone (1000-2000 m)", "Shelf",
                "Land", "Regions"],
        loc="center left", bbox_to_anchor=(0.325, 0.059), frameon=False,
        fontsize=7.5, labelcolor=INK_2, handlelength=1.3, handleheight=0.8,
        ncols=4, columnspacing=1.3, handletextpad=0.5, borderpad=0,
        borderaxespad=0,
    )
    ax.text(0.012, 0.03, "32 seasons, 1993/94 to 2024/25", fontsize=7.5,
            color=INK_2, ha="left", va="center", transform=ax.transAxes,
            zorder=6)

    fig.savefig(OUTPUT, dpi=300, facecolor="white")
    print(f"Written {OUTPUT}")


if __name__ == "__main__":
    main()
