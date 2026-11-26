"""
FA aggregation: where the recruits are when the sea ice returns.

Reads the per-particle recruitment outcomes of every cohort file of a
folder. Counts the particles of each outcome, and bins the recruits of the
two classes (on-shelf/slope and off-shelf/slope) on a regular 1° × 1°
lon/lat grid by their position on the day sea ice advances over them
(final_lon / final_lat, the position at the moment fate is determined).
Also counts the particles of each outcome per CCAMLR subarea of release
(48.6 split at 60°S), which gives the share of the eggs spawned in a
subarea that give a recruit.

Subarea assignment uses a 0.05° lookup raster of the subarea polygons, as
the per-subarea aggregates of the paper do (FS1_ice_product/aggregate.py
in krico-paper2 @ 3637d65).

Input:   DATA_DIR/YYYY_MM_DD.nc, the cohort files of version 3.0.0 of the
         KRICO Recruitment Outcomes Dataset (1994-2025), classified with
         GLORYS12 sea ice (data_v3.0_glorys12/ of krico-post-production)
         ../shared/ccamlr-data/CCAMLR_ASD_EPSG4326.shp
Output:  data/aggregated.nc

Usage:   python aggregate.py DATA_DIR

Every cohort file is checked against that version before it is counted; the
script stops at the first file that does not match.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import geopandas as gpd
import netCDF4
import numpy as np
import shapely
import xarray as xr
from shapely.ops import unary_union


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Input dataset. Version 3.0.0 is produced by krico-post-production at tag
# v3.0.0-rc1 (https://github.com/ogourgue/krico-post-production).
DATASET_NAME = "KRICO Recruitment Outcomes Dataset (1994-2025)"
DATASET_VERSION = "3.0.0"
SIC_SOURCE = "GLORYS12"

# Outcome integer flags of version 3.0.0 (CF flag values). Earlier versions
# number the outcomes differently, so the labels of every file are compared
# with these (see check_dataset_version).
OUTCOME_NAMES = (
    "recruit_on_shelf_slope",    # 0
    "recruit_off_shelf_slope",   # 1
    "M1_spawning_under_ice",     # 2
    "M3_cold_sensitivity",       # 3
    "M4_calyptopis_starvation",  # 4
    "M5_overwintering_stage",    # 5
    "M6_overwintering_habitat",  # 6
    "censored",                  # 7
    "exited_domain",             # 8
)
N_OUTCOMES = len(OUTCOME_NAMES)

# The two classes of recruits, by outcome code.
RECRUIT_CLASSES = ("recruit_on_shelf_slope", "recruit_off_shelf_slope")
RECRUIT_CODES = tuple(OUTCOME_NAMES.index(name) for name in RECRUIT_CLASSES)

# 1° × 1° lon/lat grid covering the KRICO domain (115°W to 40°E, 78°S to
# 40°S), as in F5_destination_maps of krico-paper1.
LON_MIN, LON_MAX = -115.0, 40.0
LAT_MIN, LAT_MAX = -78.0, -40.0
DLON = DLAT = 1.0

# CCAMLR subareas of the project, as named in the shapefile, and the order
# of the output. Subarea 48.6 is reported as two halves split at 60°S;
# "outside" collects the positions in none of them.
SHAPEFILE = (Path(__file__).resolve().parent.parent / "shared" / "ccamlr-data"
             / "CCAMLR_ASD_EPSG4326.shp")
SHAPE_FIELD = "GAR_Long_L"
SUBAREA_CODES = ("48.1", "48.2", "48.3", "48.4", "48.5", "48.6", "88.3")
SPLIT_LAT_486 = -60.0
SUBAREAS = ("88.3", "48.1", "48.2", "48.3", "48.4", "48.5", "48.6N", "48.6S",
            "outside")
N_SUBAREAS = len(SUBAREAS)
OUTSIDE = SUBAREAS.index("outside")

# Lookup raster of the subareas: the model domain with a small margin.
LOOKUP_RES = 0.05
LOOKUP_LON = (-116.0, 41.0)
LOOKUP_LAT = (-79.0, -39.0)

OUTPUT = Path(__file__).resolve().parent / "data" / "aggregated.nc"


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

COHORT_FILENAME_RE = re.compile(r"^(\d{4})_(\d{2})_(\d{2})\.nc$")


def find_cohort_files(data_dir: Path) -> list[Path]:
    return [path for path in sorted(data_dir.iterdir())
            if COHORT_FILENAME_RE.match(path.name)]


# ---------------------------------------------------------------------------
# Subarea lookup raster
# ---------------------------------------------------------------------------

def build_lookup(shp_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Rasterise the subarea polygons. Returns (labels[n_lat, n_lon] int8,
    lon_edges, lat_edges). Cell value = index into SUBAREAS; OUTSIDE where
    no polygon contains the cell centre.
    """
    ccamlr = gpd.read_file(shp_path)
    ccamlr = ccamlr[ccamlr[SHAPE_FIELD].isin(SUBAREA_CODES)]

    lon_edges = np.arange(LOOKUP_LON[0], LOOKUP_LON[1] + LOOKUP_RES / 2, LOOKUP_RES)
    lat_edges = np.arange(LOOKUP_LAT[0], LOOKUP_LAT[1] + LOOKUP_RES / 2, LOOKUP_RES)
    lon_c = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    lat_c = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    LON, LAT = np.meshgrid(lon_c, lat_c)

    labels = np.full(LON.shape, OUTSIDE, dtype=np.int8)
    for code in SUBAREA_CODES:
        rows = ccamlr[ccamlr[SHAPE_FIELD] == code]
        if rows.empty:
            sys.exit(f"ERROR: subarea {code} is not in {shp_path}")
        geom = unary_union(rows.geometry.values)
        inside = shapely.contains_xy(geom, LON, LAT) & (labels == OUTSIDE)
        if code == "48.6":
            labels[inside & (LAT >= SPLIT_LAT_486)] = SUBAREAS.index("48.6N")
            labels[inside & (LAT < SPLIT_LAT_486)] = SUBAREAS.index("48.6S")
        else:
            labels[inside] = SUBAREAS.index(code)
    return labels, lon_edges, lat_edges


def assign(labels: np.ndarray, lon_edges: np.ndarray, lat_edges: np.ndarray,
           lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Subarea index of each position; NaN or off the raster gives OUTSIDE."""
    out = np.full(lon.shape, OUTSIDE, dtype=np.int8)
    ok = np.isfinite(lon) & np.isfinite(lat)
    i = np.searchsorted(lat_edges, lat[ok], side="right") - 1
    j = np.searchsorted(lon_edges, lon[ok], side="right") - 1
    inb = (i >= 0) & (i < labels.shape[0]) & (j >= 0) & (j < labels.shape[1])
    vals = np.full(int(ok.sum()), OUTSIDE, dtype=np.int8)
    vals[inb] = labels[i[inb], j[inb]]
    out[ok] = vals
    return out


# ---------------------------------------------------------------------------
# Dataset version check
# ---------------------------------------------------------------------------

def check_dataset_version(nc: netCDF4.Dataset, path: Path) -> None:
    """
    Stop unless a cohort file comes from version 3.0.0 of the dataset,
    classified with GLORYS12 sea ice.

    Version 3.0.0 files carry the global attributes dataset_version and
    sic_source; earlier versions carry no dataset_version. The outcome
    labels are compared as well, because the codes of the recruits (0 and 1)
    mean something else in earlier versions.
    """
    attributes = set(nc.ncattrs())
    meanings = None
    if "outcome" in nc.variables:
        meanings = getattr(nc["outcome"], "flag_meanings", None)

    if "dataset_version" not in attributes:
        reason = ("it has no global attribute dataset_version, as files of "
                  "versions 1.0.0 and 2.0.0 do")
    elif str(nc.dataset_version) != DATASET_VERSION:
        reason = f"dataset_version is {nc.dataset_version}"
    elif "sic_source" not in attributes:
        reason = "it has no global attribute sic_source"
    elif str(nc.sic_source) != SIC_SOURCE:
        reason = (f"sic_source is {nc.sic_source}: the file was classified "
                  f"with another sea-ice product")
    elif meanings != " ".join(OUTCOME_NAMES):
        reason = f"the outcome labels are not those expected: {meanings}"
    else:
        return

    sys.exit(
        f"ERROR: {path} is not from version {DATASET_VERSION} of the "
        f"{DATASET_NAME}, classified with {SIC_SOURCE} sea ice:\n"
        f"  {reason}.\n"
        f"This script builds data/aggregated.nc from that version only."
    )


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate(data_dir: Path) -> xr.Dataset:
    """
    Walk all cohort files; count the particles of each outcome, over the
    domain and per subarea of release, and the recruits of each class on
    the lon/lat grid by final_lon / final_lat.
    """
    lon_edges = np.arange(LON_MIN, LON_MAX + DLON, DLON)
    lat_edges = np.arange(LAT_MIN, LAT_MAX + DLAT, DLAT)
    lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2.0
    lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2.0

    counts = np.zeros((len(RECRUIT_CODES), len(lat_centers), len(lon_centers)),
                      dtype=np.int64)
    total = np.zeros(N_OUTCOMES, dtype=np.int64)
    by_subarea = np.zeros((N_SUBAREAS, N_OUTCOMES), dtype=np.int64)
    labels, raster_lon, raster_lat = build_lookup(SHAPEFILE)
    n_particles = 0
    n_not_binned = 0
    code_versions = set()
    code_commits = set()

    files = find_cohort_files(data_dir)
    if not files:
        sys.exit(f"ERROR: no cohort file (YYYY_MM_DD.nc) in {data_dir}")

    print(f"Found {len(files)} cohort files in {data_dir}")
    print(f"Grid: lon {LON_MIN}..{LON_MAX} ({len(lon_centers)} cells), "
          f"lat {LAT_MIN}..{LAT_MAX} ({len(lat_centers)} cells)")

    for n, path in enumerate(files, start=1):
        with netCDF4.Dataset(path) as nc:
            check_dataset_version(nc, path)
            nc.set_auto_maskandscale(False)
            outcome = np.asarray(nc["outcome"][:])
            release_lon = np.asarray(nc["release_lon"][:], dtype=np.float64)
            release_lat = np.asarray(nc["release_lat"][:], dtype=np.float64)
            final_lon = np.asarray(nc["final_lon"][:], dtype=np.float64)
            final_lat = np.asarray(nc["final_lat"][:], dtype=np.float64)
            code_versions.add(str(getattr(nc, "code_version", "unknown")))
            code_commits.add(str(getattr(nc, "code_commit", "unknown")))

        if outcome.min() < 0 or outcome.max() >= N_OUTCOMES:
            sys.exit(f"ERROR: {path} holds outcome codes outside 0.."
                     f"{N_OUTCOMES - 1}")
        total += np.bincount(outcome, minlength=N_OUTCOMES)
        n_particles += outcome.size

        subarea = assign(labels, raster_lon, raster_lat, release_lon, release_lat)
        by_subarea += np.bincount(
            subarea.astype(np.int64) * N_OUTCOMES + outcome,
            minlength=N_SUBAREAS * N_OUTCOMES).reshape(N_SUBAREAS, N_OUTCOMES)

        # Recruits with no position, or with one outside the grid, are
        # counted in total but cannot be placed on the map.
        for k, code in enumerate(RECRUIT_CODES):
            mask = outcome == code
            h, _, _ = np.histogram2d(final_lat[mask], final_lon[mask],
                                     bins=[lat_edges, lon_edges])
            counts[k] += h.astype(np.int64)
            n_not_binned += int(mask.sum()) - int(h.sum())

        if n % 200 == 0:
            print(f"  processed {n} / {len(files)}")

    print("Done.")
    n_recruits = int(total[list(RECRUIT_CODES)].sum())
    print(f"  particles: {n_particles:,}")
    for name, value in zip(OUTCOME_NAMES, total):
        print(f"  {name:26s} {int(value):>14,}  {value / n_particles:7.2%}")
    print(f"  recruits, both classes: {n_recruits / n_particles:.2%} of the "
          f"particles")
    if n_recruits:
        print(f"  off-shelf/slope share of the recruits: "
              f"{total[RECRUIT_CODES[1]] / n_recruits:.2%}")
    print(f"  recruits not placed on the grid: {n_not_binned:,}")
    print("  per subarea of release: particles, share of the particles, "
          "recruits of both classes")
    released = by_subarea.sum(axis=1)
    recruits = by_subarea[:, list(RECRUIT_CODES)].sum(axis=1)
    for name, n_rel, n_rec in zip(SUBAREAS, released, recruits):
        rate = f"{n_rec / n_rel:7.2%}" if n_rel else "      -"
        print(f"    {name:8s} {int(n_rel):>14,}  {n_rel / n_particles:7.2%}  {rate}")
    print(f"  code version(s): {', '.join(sorted(code_versions))}")

    return xr.Dataset(
        data_vars={
            "counts": (
                ("recruit_class", "lat", "lon"),
                counts,
                {"long_name": "recruits per class per 1° × 1° cell, by "
                              "position on the day sea ice advances over "
                              "them (final_lon / final_lat)"},
            ),
            "total": (
                ("outcome",),
                total,
                {"long_name": "particles of each outcome, over all cohorts"},
            ),
            "counts_release": (
                ("subarea", "outcome"),
                by_subarea,
                {"long_name": "particles of each outcome per CCAMLR subarea "
                              "of release (release_lon / release_lat)",
                 "comment": f"subarea from a {LOOKUP_RES}° lookup raster of "
                            f"the CCAMLR shapefile; 48.6 split at 60°S"},
            ),
        },
        coords={
            "recruit_class": ("recruit_class", np.array(RECRUIT_CLASSES)),
            "outcome": ("outcome", np.array(OUTCOME_NAMES)),
            "subarea": ("subarea", np.array(SUBAREAS)),
            "lat": ("lat", lat_centers.astype(np.float32),
                    {"units": "degrees_north",
                     "long_name": "latitude (cell center)"}),
            "lon": ("lon", lon_centers.astype(np.float32),
                    {"units": "degrees_east",
                     "long_name": "longitude (cell center)"}),
        },
        attrs={
            "title": "FA recruit map: climatological aggregation",
            "description": (
                "Recruits of the two classes on a 1° × 1° lon/lat grid, by "
                "their position on the day sea ice advances over them, and "
                "particle counts per outcome, over the domain and per CCAMLR "
                "subarea of release, over all cohort files."
            ),
            "dataset": f"{DATASET_NAME}, version {DATASET_VERSION}",
            "sic_source": SIC_SOURCE,
            "shapefile": SHAPEFILE.name,
            "code_version": ", ".join(sorted(code_versions)),
            "code_commit": ", ".join(sorted(code_commits)),
            "n_files": np.int64(len(files)),
            "n_particles": np.int64(n_particles),
            "n_recruits_not_binned": np.int64(n_not_binned),
            "grid_lon_min": LON_MIN,
            "grid_lon_max": LON_MAX,
            "grid_lat_min": LAT_MIN,
            "grid_lat_max": LAT_MAX,
            "grid_resolution_deg": DLON,
        },
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python aggregate.py DATA_DIR\n"
                 "DATA_DIR is the folder of the version 3.0.0 cohort files, "
                 "e.g. /path/to/krico-post-production/recruitment/"
                 "data_v3.0_glorys12")
    data_dir = Path(sys.argv[1])
    if not data_dir.is_dir():
        sys.exit(f"ERROR: {data_dir} is not a folder")

    ds = aggregate(data_dir)
    OUTPUT.parent.mkdir(exist_ok=True)
    ds.to_netcdf(OUTPUT)
    print(f"Written {OUTPUT}")


if __name__ == "__main__":
    main()
