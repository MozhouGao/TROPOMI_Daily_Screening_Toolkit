"""US EPA GHGRP facility methane (Envirofacts, latest reporting year)."""
import csv
import io
import os
from collections import defaultdict
from datetime import timedelta
from typing import Callable, Dict, List, Optional, Tuple

from backend.schemas import BBox, DateRange, Plume, PlumeCollection, SourceError
from backend.sources.file_cache import (
    GHGRP_CH4_GWP,
    capped_collection,
    co2e_tonnes_to_ch4_kg_hr,
    ensure_file,
    in_bbox,
    optional_float,
    year_start,
)
from backend.sources.http import fetch_bytes

SOURCE_ID = "ghgrp_us"
ATTRIBUTION = "Source: US EPA Greenhouse Gas Reporting Program"
LICENSE_URL = "https://www.epa.gov/ghgreporting/data-sets"
CACHE_FILENAME = "ghgrp_us_ch4.csv"
CACHE_MAX_AGE = timedelta(days=30)
PAGE_SIZE = 10000
REPORT_YEARS = tuple(range(2025, 2017, -1))


def default_cache_dir() -> str:
    return os.path.join(os.getcwd(), "GHGRP_data")


def page_url(year: int, start: int, end: int) -> str:
    return (
        "https://data.epa.gov/efservice/PUB_DIM_FACILITY/YEAR/=/"
        f"{year}/PUB_FACTS_SUBP_GHG_EMISSION/PUB_DIM_GHG/GAS_CODE/=/CH4/"
        f"ROWS/{start}:{end}/CSV"
    )


def _facility(row: dict) -> Optional[Tuple[str, dict]]:
    facility_id = (row.get("facility_id") or "").strip().split(".")[0]
    lat = optional_float(row.get("latitude"))
    lon = optional_float(row.get("longitude"))
    co2e = optional_float(row.get("co2e_emission")) or 0.0
    year = int(optional_float(row.get("year")) or 0)
    if not facility_id or lat is None or lon is None or year == 0:
        return None
    return facility_id, {
        "lat": lat,
        "lon": lon,
        "year": year,
        "name": (row.get("facility_name") or "").strip(),
        "subparts": (row.get("reported_subparts") or "").strip(),
        "parent": (row.get("parent_company") or "").strip(),
        "state": (row.get("state") or "").strip(),
        "co2e": co2e,
    }


def parse_csv(text: str) -> List[Plume]:
    grouped: Dict[str, dict] = {}
    totals: Dict[str, float] = defaultdict(float)
    for row in csv.DictReader(io.StringIO(text)):
        parsed = _facility(row)
        if parsed is None:
            continue
        facility_id, record = parsed
        totals[facility_id] += record["co2e"]
        grouped.setdefault(facility_id, record)
    plumes = []
    for facility_id, record in grouped.items():
        co2e = totals[facility_id]
        if co2e <= 0:
            continue
        plumes.append(
            Plume(
                id=facility_id,
                source=SOURCE_ID,
                observed_at=year_start(record["year"]),
                lon=record["lon"],
                lat=record["lat"],
                emission_kg_hr=co2e_tonnes_to_ch4_kg_hr(co2e),
                instrument="US GHGRP",
                sector=record["subparts"] or None,
                image_url=None,
                attribution=ATTRIBUTION,
                details={
                    key: value
                    for key, value in (
                        ("Name", record["name"]),
                        ("State", record["state"]),
                        ("Parent", record["parent"]),
                        ("Subparts", record["subparts"]),
                        ("GWP", str(GHGRP_CH4_GWP)),
                    )
                    if value
                },
            )
        )
    return plumes


def _download_year(year: int, fetch: Callable[[str], bytes]) -> bytes:
    chunks: List[bytes] = []
    header = b""
    start = 0
    while True:
        body = fetch(page_url(year, start, start + PAGE_SIZE - 1))
        if not body.strip() or b"facility_id" not in body[:400].lower():
            break
        text = body.decode("utf-8-sig")
        lines = text.splitlines()
        if len(lines) <= 1:
            break
        if not header:
            header = (lines[0] + "\n").encode("utf-8")
            chunks.append(body if start == 0 else ("\n".join(lines[1:]) + "\n").encode("utf-8"))
        else:
            chunks.append(("\n".join(lines[1:]) + "\n").encode("utf-8"))
        if len(lines) - 1 < PAGE_SIZE:
            break
        start += PAGE_SIZE
    return b"".join(chunks)


def fetch_export(fetch: Callable[[str], bytes] = fetch_bytes) -> bytes:
    last_error = None
    for year in REPORT_YEARS:
        try:
            body = _download_year(year, fetch)
        except SourceError as exc:
            last_error = exc
            continue
        if body.strip():
            return body
    if last_error:
        raise last_error
    return b""


def fetch_plumes(
    bbox: BBox,
    dates: DateRange,
    cache_dir: Optional[str] = None,
    fetch: Callable[[str], bytes] = fetch_bytes,
) -> PlumeCollection:
    del dates
    path = os.path.join(cache_dir or default_cache_dir(), CACHE_FILENAME)

    def load(_url: str) -> bytes:
        return fetch_export(fetch)

    ensure_file("https://data.epa.gov/efservice/ghgrp_us", path, CACHE_MAX_AGE, fetch=load)
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        items = [plume for plume in parse_csv(handle.read()) if in_bbox(plume, bbox)]
    return capped_collection(SOURCE_ID, items, ATTRIBUTION, LICENSE_URL)
