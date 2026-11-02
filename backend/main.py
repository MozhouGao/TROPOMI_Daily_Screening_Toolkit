"""FastAPI backend for the methane data screening toolkit."""
from datetime import date
from typing import Callable, Dict

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware

from backend.schemas import BBox, DateRange, PlumeCollection, SourceError
from backend.sources import carbonmapper, mars

PlumeFetcher = Callable[[BBox, DateRange], PlumeCollection]

PLUME_SOURCES: Dict[str, PlumeFetcher] = {
    carbonmapper.SOURCE_ID: carbonmapper.fetch_plumes,
    mars.SOURCE_ID: mars.fetch_plumes,
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
    if west > east or south > north:
        raise HTTPException(status_code=422, detail="Bounding box edges are inverted")
    if start > end:
        raise HTTPException(status_code=422, detail="start must be on or before end")
    try:
        return fetcher(BBox(west, south, east, north), DateRange(start, end))
    except SourceError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc
