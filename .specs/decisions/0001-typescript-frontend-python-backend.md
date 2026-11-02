# ADR 0001: TypeScript frontend over a Python backend

Date: 2026-11-01. Status: accepted.

## Context

The toolkit was a single Dash app. The goal is a multi-source dashboard in the
style of a layer registry (one map, toggleable sources, shared time window).
Dash can draw markers but makes a registry-driven UI, per-layer panels, and
30k-point layers awkward. The TROPOMI algorithm (numpy, netCDF4, cartopy) must
stay in Python.

## Decision

- Frontend: Vite, React, TypeScript strict, antd, `@tanstack/react-query`,
  jotai, luxon, MapLibre GL JS via `react-map-gl`, Vitest. OSM raster tiles.
- Backend: FastAPI + uvicorn, Python 3.9 compatible, `urllib` for provider
  HTTP. Pydantic models with camelCase aliases define the API contract.
- The backend owns all provider access; the frontend only calls `/api/*`.
- Design ideas come from worldmonitor; no code is copied (AGPL-3.0 vs MIT).

## Consequences

- Two processes to run locally (uvicorn and vite dev server); README covers it.
- The Dash app is removed in the last slice of the milestone.
- Existing `tests/` keep running with `unittest`; backend tests live in
  `backend/tests/` and run with `python -m unittest discover -s backend/tests -t .`.
