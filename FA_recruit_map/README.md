# FA — Map of the recruits when the sea ice returns

Image submitted with the abstract. One map: where the recruits are on the day sea ice advances over them, with the outline of the slope where the eggs are spawned.

A recruit is a larva that passes the rules acting before winter, has reached furcilia IV and is under sea ice when the ice advances. The figure calls the recruits "larvae ready for winter", as the abstract does. They come in two classes, by where they are at that moment: over the shelf or slope (sea floor shallower than 2000 m) or beyond it. The map shows both classes together.

Adapted from `F5_destination_maps` in [krico-paper1](https://github.com/ogourgue/krico-paper1) @ `d883793`, whose panels "Recruitment success" and "Off-shelf (M5b)" are these two classes under their earlier names. `aggregate.py` and `plot.py` are rewritten for one map and for version 3.0.0 of the dataset.

## What it shows

- All 32 spawning years (1994–2025) and all release dates together.
- Per-cell quantity: the share of all recruits whose position, on the day sea ice advances over them, fell in the cell (%). It sums to 100 over the map.
- Grid: 1° × 1° lon/lat. Colours clipped at the 99th percentile of the non-empty cells.
- Outline: the slope where particles are released, sea floor between 1000 and 2000 m inside the CCAMLR subareas of the project, from 76°S to 41°S.
- Projection: South Polar Stereographic, central longitude 37.5°W. The map is cropped to 100°W–28°E and 77.5°S–46.5°S; the model domain is wider.
- The subtitle gives the off-shelf/slope share of the recruits, computed from the data.
- The two emphasised labels mark the regions where almost no egg gives a recruit, with the reason given in the abstract.

## Input data

`aggregate.py` reads the per-particle recruitment outcomes of version 3.0.0 of the KRICO Recruitment Outcomes Dataset (1994–2025), the set classified with GLORYS12 sea ice (`data_v3.0_glorys12/`), from the folder given on the command line. It checks every file it reads and stops if one comes from another version or another sea-ice product. The main README says how the version is recognised.

`data/bathymetry.nc` is an unchanged copy of `F1_domain_map/data/aggregated.nc` in [krico-paper1](https://github.com/ogourgue/krico-paper1) @ `d883793`: the GLORYS12v1 bathymetry on the model grid, masked to the CCAMLR subareas of the project. The file also holds the sea-ice climatology of that figure, which is not used here.

## Files

- `aggregate.py` — counts the particles of each outcome and bins the recruits of each class on the 1° grid, writes `data/aggregated.nc`.
- `plot.py` — reads `data/aggregated.nc` and `data/bathymetry.nc`, writes `recruit_map.png`.
- `data/bathymetry.nc` — bathymetry on the model grid, for the slope outline.
- `data/aggregated.nc` — recruits per class on the 1° grid, and particle counts per outcome. Not committed yet.
- `recruit_map.png` — the rendered figure. Not committed yet.

## Run

```bash
python aggregate.py /path/to/krico-post-production/recruitment/data_v3.0_glorys12
python plot.py
```
