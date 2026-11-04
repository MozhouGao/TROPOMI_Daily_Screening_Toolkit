# US and Canada review layers

## Goal

Add the public, map-ready US and Canada datasets from Vollrath, Xing,
Hugenholtz et al. 2025 (ERL 20 083003) as toggleable layers on the existing
dashboard, so a drawn region can show satellite/aircraft plumes next to
self-reported facilities and US/Canada inventory grids.

## Planned layers

| Layer | Source | Provides | Status |
|-------|--------|----------|--------|
| GHGSat landfill catalog | Dogniaux et al.; Zenodo [16641834](https://doi.org/10.5281/zenodo.16641834) | Point detections at landfills (`Q_t_per_hr` → kg/hr). Catalog layer. | Done |
| US GHGRP landfill / facility points | EPA Envirofacts CH4 join, latest reporting year | Self-reported facility CH4 (includes landfills via subparts). tCO2e / GWP 25 → kg/hr. | Done |
| Canada GHGRP landfill / facility points | ECCC Data Mart GHGRP CSV, latest reference year | Self-reported facility CH4 (includes waste NAICS). Reported CH4 tonnes → kg/hr. | Done |
| Canada NIR overlay | Scarpelli et al. 2018 gridded NIR; [10.7910/DVN/CC3KLO](https://doi.org/10.7910/DVN/CC3KLO) | 0.1° cells (solid waste, wastewater, residential combustion, oil/gas leakage) ≥ 1 kg/hr. Inventory overlay. | Done |
| EPA GHGI overlay | Maasakkers et al. gridded US GHGI v3; Zenodo [16782734](https://doi.org/10.5281/zenodo.16782734) | 0.1° CONUS annual CH4 (2012–2022, 34 source categories) as inventory cells, same pattern as Canada NIR. Complements GHGRP facility points with a spatial bottom-up grid. | Planned (not implemented) |
| MethaneSAT (EDF, Google Earth Engine) | Publisher catalog `edf-methanesat-ee` (L4 points / L4 area / L3 XCH4; request-gated) | Candidate plume/overlay source alongside Carbon Mapper, MARS, TROPOMI, GHGSat, Cusworth, GHGRP, Nesser, Canada NIR. Feasibility lookup already in progress (GEE Python API or GCP download; license forbids redistributing raw L3/L4). | Evaluating — not confirmed; do not implement until access + license are decided |

The four review-paper products above that were previously called out as deferred (GHGRP landfill points, GHGSat landfill catalog, EPA GHGI, Canada NIR) are now in this milestone: the first three rows plus Canada NIR are already shipped in this branch; EPA GHGI remains the open implementation slice.

## Scope

In:

- GHGSat landfill plumes (Dogniaux et al.; Zenodo 16641834) as a catalog layer.
- Cusworth et al. 2024 US landfill campaign detections (Zenodo 10642570) as a
  catalog layer with GHGRP IDs.
- US EPA GHGRP facility CH4 (Envirofacts, latest reporting year) as a facility
  layer. Convert reported tCO2e with AR4 GWP 25 to kg/hr. Landfills are included
  via reported subparts, not a separate source id.
- Canada GHGRP facility CH4 (ECCC Data Mart CSV, latest reference year) as a
  facility layer. Use reported CH4 tonnes. Waste/landfill facilities included
  via NAICS / sector.
- Nesser et al. 2024 TROPOMI inversion landfill posteriors (GitHub
  `ghgrp_processed.csv`) as a facility layer. Convert Gg/yr to kg/hr.
- Canada NIR 2018 gridded inventory (Scarpelli et al.; solid waste, wastewater,
  residential combustion, oil/gas leakage) as an inventory overlay of 0.1°
  cells above 1 kg/hr.
- EPA gridded GHGI NetCDF inventory overlay (Maasakkers et al. v3): 0.1° CONUS
  cells, archival, bbox-only, cap 500, default off. Plan only until this slice
  is picked up — do not treat as shipped.
- Layer registry, i18n, attribution, local file cache, backend tests and
  frontend registry tests.

Out:

- MethaneSAT / MethaneAIR on Earth Engine until evaluation is confirmed
  (request form, EE or GCP auth, no raw L3/L4 redistribution). Tracked as
  Evaluating, not an implementation slice.
- Nesser city polygons (`summ_cities.csv` has no coordinates).
- Street-level leak campaign maps, tower networks, EDGAR.
- New API shapes; every source still returns `PlumeCollection` via
  `/api/plumes/{source}`.

## Approach

- Reuse the Carbon Mapper / MARS adapter pattern: cache a public file, parse
  to `Plume`, filter by bbox, cap at 500 rows with `truncated`.
- Archival layers ignore the request date range. The default UI window is
  “last N days”; these products are 2016–2024 vintages. Freshness copy states
  the year.
- `kind` on the registry: `catalog` | `facility` | `inventory` | `screening`.
  All non-screening kinds use `/api/plumes/{source}`.
- Default on: existing three layers plus `ghgsat`, `ghgrp_us`, `ghgrp_ca`.
  Campaign, inversion, and inventory layers default off.
- Decisions: see `../../decisions/0001-typescript-frontend-python-backend.md`.

## Acceptance criteria

1. `GET /api/plumes/ghgsat` with a bbox returns detections with `lat`/`lon`
   not both zero; `Q_t_per_hr` is converted to kg/hr; unknown source is still
   404.
2. `GET /api/plumes/cusworth` returns landfill campaign points with GHGRP ID
   in `details` and a parsed observation time from the plume candidate ID.
3. `GET /api/plumes/ghgrp_us` and `ghgrp_ca` return one point per facility for
   the latest year in the cached file, CH4 > 0, filtered by bbox.
4. `GET /api/plumes/nesser` returns landfill posteriors with `post_mean` as
   kg/hr and GHGRP / site name in details.
5. `GET /api/plumes/canada_nir` returns 0.1° cell centroids for the four
   sectors inside the bbox with emission ≥ 1 kg/hr, newest-sector-agnostic
   sort by rate, truncated at 500.
6. Frontend registry lists all new layers; catalog/facility/inventory layers
   fetch when enabled; screening stays TROPOMI-only.
7. README names each new source, cache directory, license, and that archival
   layers ignore the time window.
8. `GET /api/plumes/epa_ghgi` (planned): 0.1° CONUS GHGI cells inside the bbox
   with emission ≥ 1 kg/hr, newest-year vintage in the cached NetCDF, truncated
   at 500. Not an acceptance gate until the Planned slice is implemented.

## Risks and open questions

| Item | Status |
|------|--------|
| Default last-N-days window would hide 2021–2022 catalogs | Accepted: archival sources ignore dates |
| EPA GHGRP reports CH4 as tCO2e; GWP 25 (AR4) | Accepted, shown in details |
| Canada GHGRP CSV is ~10 MB; first fetch is slow | Accepted, same cache pattern as MARS |
| Canada NIR arrays are 0.1° global; must subset by bbox | Accepted |
| Cusworth overlaps live Carbon Mapper | Default off |
| EPA gridded GHGI NetCDF is large; same bbox + cap as Canada NIR | Planned; not implemented |
| MethaneSAT GEE is request-gated; license forbids serving raw L3/L4 | Evaluating; not in build scope |

## References

- Vollrath et al. 2025: https://doi.org/10.1088/1748-9326/ade9e1
- GHGSat landfills: https://doi.org/10.5281/zenodo.16641834
- Cusworth landfills: https://doi.org/10.5281/zenodo.10642570
- Nesser inversion: https://doi.org/10.5281/zenodo.10946769
- Canada NIR: https://doi.org/10.7910/DVN/CC3KLO
- US GHGRP Envirofacts: https://www.epa.gov/enviro/envirofacts-data-service-api
- Canada GHGRP: https://open.canada.ca/data/en/dataset/a8ba14b7-7f23-462a-bdbb-83b0ef629823
- EPA gridded GHGI v3: https://doi.org/10.5281/zenodo.16782734
- MethaneSAT EE catalog: https://developers.google.com/earth-engine/datasets/publisher/edf-methanesat-ee
- MethaneSAT access form: https://forms.gle/jqw4Mvr63dsV1fUF8
