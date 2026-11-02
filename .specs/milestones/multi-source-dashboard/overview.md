# Multi-source methane dashboard

## Goal

Turn the TROPOMI-only Dash app into a dashboard where one map shows several
methane plume sources as toggleable layers over a shared region and time
window: TROPOMI patch screening (existing algorithm), Carbon Mapper, and UNEP
IMEO MARS. A TypeScript frontend talks to a FastAPI backend that owns provider
access and normalizes every source into one plume shape.

## Scope

In:

- `backend/` FastAPI service: `/api/plumes/{source}` for `carbonmapper` and
  `mars`, `/api/tropomi/*` wrapping `TROPOMI_toolkit.py`.
- `frontend/` Vite + React app: MapLibre map, rectangle draw, layer registry,
  time window, plume list, detail panel, credential manager.
- Attribution and license text for every external source.
- Removal of the Dash app once the new UI reaches parity.

Out:

- CO2 plumes, Carbon Mapper scene rasters or vector tiles, MARS source
  clustering, user accounts, hosted deployment.

## Approach

- One shared plume shape crosses the API boundary (`backend/schemas.py`):
  `id, source, observedAt, lon, lat, emissionKgHr, uncertaintyKgHr, instrument,
  sector, imageUrl, attribution, details`. Adapters normalize; the frontend
  never sees provider field names.
- Carbon Mapper: public catalog, no key. Repeated `bbox` keys in W, S, E, N
  order (comma form returns 422). 200-row cap with `truncated` flag.
- MARS: no API. Public zip export cached under `MARS_data/`, refreshed weekly,
  filtered in memory per request. Stale cache is kept when the download fails.
- Decisions: see `../../decisions/0001-typescript-frontend-python-backend.md`.

## Acceptance criteria

1. `GET /api/plumes/carbonmapper?west&south&east&north&start&end` returns the
   shared shape; a bbox count above the cap sets `truncated: true`.
2. `GET /api/plumes/mars` with the same query answers from the local cache
   after the first call and never redownloads within seven days.
3. Unknown source is 404, inverted bbox or dates are 422, provider failure is
   502 with the provider message.
4. Frontend shows each enabled layer's plumes on the map and in the list for
   the drawn rectangle and chosen time window; selecting a plume shows its
   detail with the layer's attribution line.
5. TROPOMI screening runs from the new UI with the same inputs as the Dash app
   and draws the XCH4 grid and hotspots on the map.
6. README documents install, run, and the non-commercial terms of both
   external sources.

## Risks and open questions

| Item | Status |
|------|--------|
| Carbon Mapper image URLs are signed and expire; display immediately, never cache | Accepted |
| MARS `tile_background` is a product ID, not an image; no MARS imagery in v1 | Accepted |
| MARS export is 11 MB; first request per week is slow (about 4 s) | Accepted, loading state in UI |
| Both external sources are non-commercial licenses | README states it |

## References

- Carbon Mapper product guide: https://carbonmapper.org/articles/product-guide
- Carbon Mapper terms: https://carbonmapper.org/terms
- MARS data dictionary: https://methanedata.unep.org/dict-mars-plumes
- UI pattern reference (ideas only, AGPL code not reused):
  https://github.com/koala73/worldmonitor
