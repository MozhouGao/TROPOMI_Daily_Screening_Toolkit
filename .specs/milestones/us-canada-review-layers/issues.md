# Issues

| # | Slice | Status | Notes |
|---|-------|--------|-------|
| 1 | Shared file cache and unit conversions | Done | `backend/sources/file_cache.py` |
| 2 | GHGSat landfill catalog adapter | Done | Zenodo CSV, skip non-detects |
| 3 | Cusworth 2024 landfill campaign adapter | Done | GHGRP ID + parsed time |
| 4 | US and Canada GHGRP landfill / facility points | Done | Latest year, bbox, kg/hr; landfills included via subpart / NAICS |
| 5 | Nesser 2024 landfill posterior adapter | Done | Gg/yr → kg/hr |
| 6 | Canada NIR inventory overlay | Done | Four sectors, ≥ 1 kg/hr, 0.1° cells |
| 7 | Register sources in FastAPI | Done | `PLUME_SOURCES` |
| 8 | Frontend registry, i18n, swatches, copy | Done | kinds + count labels |
| 9 | README and cache gitignore | Done | licenses and archival date note |
| 10 | EPA gridded GHGI inventory overlay | Planned | Maasakkers v3 NetCDF, 0.1° CONUS; same cache + bbox + cap as Canada NIR. Do not implement until this slice is picked up. |
| 11 | MethaneSAT (EDF / GEE) layer listing | Done | Dashboard layer only. Empty catalog + attribution; no raw L3/L4 from GEE. |
