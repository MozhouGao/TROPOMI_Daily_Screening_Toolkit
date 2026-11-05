"""Cusworth et al. 2024 US landfill campaign plumes (Zenodo 10642570)."""
import csv
import io
import os
import re
from datetime import datetime, timedelta, timezone
from typing import Callable, List, Optional, Tuple

from backend.schemas import BBox, DateRange, Plume, PlumeCollection
from backend.sources.file_cache import (
    capped_collection,
    ensure_file,
    in_bbox,
    optional_float,
    row_value,
    year_start,
)
from backend.sources.http import fetch_bytes

SOURCE_ID = "cusworth"
EXPORT_URL = "https://zenodo.org/records/10642570/files/Landfill_plume_detections.csv?download=1"
ATTRIBUTION = "Source: Cusworth et al. 2024 US landfill campaign"
LICENSE_URL = "https://doi.org/10.5281/zenodo.10642570"
CACHE_FILENAME = "Landfill_plume_detections.csv"
CACHE_MAX_AGE = timedelta(days=30)
STAMP = re.compile(r"(\d{8})t(\d{6})", re.IGNORECASE)


def default_cache_dir() -> str:
    return os.path.join(os.getcwd(), "Cusworth_data")


def parse_observed_at(candidate: str) -> Tuple[datetime, str]:
    match = STAMP.search(candidate)
    instrument = "Carbon Mapper"
    lowered = candidate.lower()
    if lowered.startswith("gao"):
        instrument = "GAO"
    elif lowered.startswith("ang"):
        instrument = "AVIRIS-NG"
    if match is None:
        return year_start(2018), instrument
    day, clock = match.group(1), match.group(2)
    observed = datetime.strptime(day + clock, "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc)
    return observed, instrument


def parse_row(row: dict) -> Optional[Plume]:
    candidate = row_value(row, "plume candidate")
    lat = optional_float(row_value(row, "latitude"))
    lon = optional_float(row_value(row, "longitude"))
    rate = optional_float(row_value(row, "instantaneous emission rate", exclude=("uncertainty",)))
    if not candidate or lat is None or lon is None:
        return None
    observed, instrument = parse_observed_at(candidate)
    ghgrp_id = row_value(row, "ghgrp")
    return Plume(
        id=candidate,
        source=SOURCE_ID,
        observed_at=observed,
        lon=lon,
        lat=lat,
        emission_kg_hr=rate,
        uncertainty_kg_hr=optional_float(row_value(row, "uncertainty")),
        instrument=instrument,
        sector="solid waste",
        image_url=None,
        attribution=ATTRIBUTION,
        details={"GHGRP ID": ghgrp_id} if ghgrp_id else {},
    )


def parse_csv(text: str) -> List[Plume]:
    return [plume for plume in (parse_row(row) for row in csv.DictReader(io.StringIO(text))) if plume]


def fetch_plumes(
    bbox: BBox,
    dates: DateRange,
    cache_dir: Optional[str] = None,
    fetch: Callable[[str], bytes] = fetch_bytes,
) -> PlumeCollection:
    del dates
    path = ensure_file(
        EXPORT_URL,
        os.path.join(cache_dir or default_cache_dir(), CACHE_FILENAME),
        CACHE_MAX_AGE,
        fetch=fetch,
    )
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        items = [plume for plume in parse_csv(handle.read()) if in_bbox(plume, bbox)]
    return capped_collection(SOURCE_ID, items, ATTRIBUTION, LICENSE_URL)
