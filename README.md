# Presentation: KRICO at AMEMR 2027, 5-8 July 2027

Figures and scripts for a contribution on the first paper of the KRICO project (Antarctic Krill Connectivity) to AMEMR 2027 (Advances in Marine Ecosystem Modelling Research), University of Plymouth, UK, 5-8 July 2027.

Author: Olivier Gourgue (RBINS)

Abstract: [AMEMR 2027 abstract KRICO](https://docs.google.com/document/d/17gDisDJ9yEJo6pr7_cLbUosD3myWUPuhntwLWtXJMV0/edit?usp=sharing) (Google Docs)

For now the repository holds one figure, the image submitted with the abstract. It is adapted from [krico-paper1](https://github.com/ogourgue/krico-paper1) and shows the work as submitted; later versions of the paper may differ.

Related repositories:

* __[krico-paper1](https://github.com/ogourgue/krico-paper1)__ — Figures and analysis for the paper
* __[krico-post-production](https://github.com/ogourgue/krico-post-production)__ — Post-processing and recruitment classification pipeline
* __[krico-templates](https://github.com/ogourgue/krico-templates)__ — Simulation templates (Parcels + GLORYS12v1)

## Reproducing the figure

Each figure has its own folder, with `aggregate.py`, `plot.py`, the aggregated data in `data/` and the rendered figure. Once the aggregated data are committed, redrawing a figure needs only:

```bash
cd <folder>
python plot.py
```

`FA_recruit_map/data/aggregated.nc` and the rendered figure are not committed yet: they are built from version 3.0.0 of the dataset (see below).

## Re-running the aggregation

`aggregate.py` rebuilds `data/aggregated.nc` from the data upstream. Each folder README lists what its script reads.

FA reads the per-particle recruitment outcomes of **version 3.0.0** of the KRICO Recruitment Outcomes Dataset (1994–2025), the set classified with GLORYS12 sea ice. That version is produced by [krico-post-production](https://github.com/ogourgue/krico-post-production) at tag `v3.0.0-rc1` and is not on Zenodo yet. Earlier versions number the outcomes differently, so the script checks every file: a version 3.0.0 file has the global attributes `dataset_version` (`3.0.0`) and `sic_source` (`GLORYS12`). The script stops at the first file that does not match.

The folder that holds the cohort files is given on the command line:

```bash
cd FA_recruit_map
python aggregate.py /path/to/krico-post-production/recruitment/data_v3.0_glorys12
python plot.py
```

## Folders

| Folder | Use | Source |
|---|---|---|
| `FA_recruit_map` | Image submitted with the abstract | Recruitment and off-shelf panels of `F5_destination_maps` in [krico-paper1](https://github.com/ogourgue/krico-paper1) @ `d883793`, merged into one map |
