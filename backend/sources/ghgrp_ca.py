"""Canada GHGRP facility methane from the ECCC Data Mart CSV."""
import csv
import io
import os
from datetime import timedelta
from typing import Callable, List, Optional
from urllib.parse import quote

from backend.schemas import BBox, DateRange, Plume, PlumeCollection
from backend.sources.file_cache import (
    capped_collection,
    ensure_file,
    in_bbox,
    optional_float,
    row_value,
    tonnes_per_year_to_kg_hr,
    year_start,
)
from backend.sources.http import fetch_bytes

SOURCE_ID = "ghgrp_ca"
FILE_PATH = (
    "/substances/monitor/greenhouse-gas-reporting-program-ghgrp-facility-greenhouse-gas-ghg-data/"
    "PDGES-GHGRP-GHGEmissionsGES-2004-Present.csv"
)
EXPORT_URL = f"https://data-donnees.az.ec.gc.ca/api/file?path={quote(FILE_PATH)}"
ATTRIBUTION = "Source: Environment and Climate Change Canada GHGRP"
LICENSE_URL = "https://open.canada.ca/data/en/dataset/a8ba14b7-7f23-462a-bdbb-83b0ef629823"
CACHE_FILENAME = "PDGES-GHGRP-GHGEmissionsGES-2004-Present.csv"
CACHE_MAX_AGE = timedelta(days=30)


def default_cache_dir() -> str:
    return os.path.join(os.getcwd(), "GHGRP_data")


def parse_row(row: dict, latest_year: int) -> Optional[Plume]:
    year = int(optional_float(row_value(row, "reference year")) or 0)
    if year != latest_year:
        return None
    ch4_tonnes = optional_float(row_value(row, "ch4 (tonnes)", exclude=("co2e", "éq", "eq"))) or 0.0
    if ch4_tonnes <= 0:
        return None
    lat = optional_float(row_value(row, "latitude"))
    lon = optional_float(row_value(row, "longitude"))
    facility_id = row_value(row, "ghgrp id")
    if lat is None or lon is None or not facility_id:
        return None
    sector = row_value(row, "english facility naics")
    name = row_value(row, "facility name")
    return Plume(
        id=facility_id,
        source=SOURCE_ID,
        observed_at=year_start(year),
        lon=lon,
        lat=lat,
        emission_kg_hr=tonnes_per_year_to_kg_hr(ch4_tonnes),
        instrument="Canada GHGRP",
        sector=sector or None,
        image_url=None,
        attribution=ATTRIBUTION,
        details={key: value for key, value in (("Name", name), ("Year", str(year))) if value},
    )


def parse_csv(text: str) -> List[Plume]:
    rows = list(csv.DictReader(io.StringIO(text)))
    years = [int(optional_float(row_value(row, "reference year")) or 0) for row in rows]
    latest_year = max(years) if years else 0
    return [plume for plume in (parse_row(row, latest_year) for row in rows) if plume]


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
