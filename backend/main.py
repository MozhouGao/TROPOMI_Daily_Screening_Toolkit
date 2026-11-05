"""FastAPI backend for the methane data screening toolkit."""
from datetime import date
from typing import Callable, Dict

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware

from backend.schemas import (
    BBox,
    CredentialStatus,
    CredentialWrite,
    DateRange,
    DateWindow,
    NewsCollection,
    PlumeCollection,
    ScreenRequest,
    SourceError,
    TropomiScreenResult,
)
from backend.sources import (
    canada_nir,
    carbonmapper,
    cusworth,
    ghgrp_ca,
    ghgrp_us,
    ghgsat,
    mars,
    methanesat,
    nesser,
    news,
    tropomi,
)

PlumeFetcher = Callable[[BBox, DateRange], PlumeCollection]

PLUME_SOURCES: Dict[str, PlumeFetcher] = {
    carbonmapper.SOURCE_ID: carbonmapper.fetch_plumes,
    mars.SOURCE_ID: mars.fetch_plumes,
    ghgsat.SOURCE_ID: ghgsat.fetch_plumes,
    cusworth.SOURCE_ID: cusworth.fetch_plumes,
    ghgrp_us.SOURCE_ID: ghgrp_us.fetch_plumes,
    ghgrp_ca.SOURCE_ID: ghgrp_ca.fetch_plumes,
    nesser.SOURCE_ID: nesser.fetch_plumes,
    canada_nir.SOURCE_ID: canada_nir.fetch_plumes,
    methanesat.SOURCE_ID: methanesat.fetch_plumes,
}

app = FastAPI(title="Methane Data Screening Toolkit API", version="2.0.0")
app.add_middleware(
    CORSMiddleware,
    allow_origin_regex=r"https?://(localhost|127\.0\.0\.1)(:\d+)?",
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/api/health")
def health() -> Dict[str, str]:
    return {"status": "ok"}


@app.get("/api/plumes/{source}", response_model=PlumeCollection)
def list_plumes(
    source: str,
    west: float = Query(..., ge=-180, le=180),
    south: float = Query(..., ge=-90, le=90),
    east: float = Query(..., ge=-180, le=180),
    north: float = Query(..., ge=-90, le=90),
    start: date = Query(...),
    end: date = Query(...),
) -> PlumeCollection:
    fetcher = PLUME_SOURCES.get(source)
    if fetcher is None:
        raise HTTPException(status_code=404, detail=f"Unknown plume source '{source}'")
    bbox = _bbox_or_422(west, south, east, north)
    dates = _dates_or_422(start, end)
    try:
        return fetcher(bbox, dates)
    except SourceError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc


@app.get("/api/news", response_model=NewsCollection)
def list_news(
    west: float = Query(..., ge=-180, le=180),
    south: float = Query(..., ge=-90, le=90),
    east: float = Query(..., ge=-180, le=180),
    north: float = Query(..., ge=-90, le=90),
) -> NewsCollection:
    bbox = _bbox_or_422(west, south, east, north)
    try:
        return news.fetch_news(bbox)
    except SourceError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc


def _bbox_or_422(west: float, south: float, east: float, north: float) -> BBox:
    if west > east or south > north:
        raise HTTPException(status_code=422, detail="Bounding box edges are inverted")
    return BBox(west, south, east, north)


def _dates_or_422(start: date, end: date) -> DateRange:
    if start > end:
        raise HTTPException(status_code=422, detail="start must be on or before end")
    return DateRange(start, end)


@app.get("/api/tropomi/credentials", response_model=CredentialStatus)
def get_credentials() -> CredentialStatus:
    return tropomi.credential_status()


@app.put("/api/tropomi/credentials", response_model=CredentialStatus)
def put_credentials(payload: CredentialWrite) -> CredentialStatus:
    try:
        return tropomi.save_credentials(payload)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.post("/api/tropomi/download")
def tropomi_download(window: DateWindow) -> Dict[str, str]:
    dates = _dates_or_422(window.start, window.end)
    try:
        return {"message": tropomi.download_observations(dates)}
    except SourceError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc


@app.post("/api/tropomi/screen", response_model=TropomiScreenResult)
def tropomi_screen(payload: ScreenRequest) -> TropomiScreenResult:
    bbox = _bbox_or_422(payload.west, payload.south, payload.east, payload.north)
    dates = _dates_or_422(payload.start, payload.end)
    try:
        return tropomi.screen_region(
            bbox, dates, payload.threshold_delta, payload.min_pixels
        )
    except SourceError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc
