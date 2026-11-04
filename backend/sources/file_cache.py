"""Local cache for public archive files and CH4 rate unit conversions."""
import os
from datetime import datetime, timedelta, timezone
from typing import Callable, List, Optional

from backend.schemas import BBox, Plume, PlumeCollection, SourceError
from backend.sources.http import fetch_bytes

HOURS_PER_YEAR = 365.25 * 24
SECONDS_PER_HOUR = 3600
GHGRP_CH4_GWP = 25
MAX_ITEMS = 500


def cache_is_fresh(path: str, now: datetime, max_age: timedelta) -> bool:
    if os.path.isfile(path):
        modified = datetime.fromtimestamp(os.path.getmtime(path), tz=timezone.utc)
        return now - modified < max_age
    return False


def ensure_file(
    url: str,
    path: str,
    max_age: timedelta,
    now: Optional[datetime] = None,
    fetch: Callable[[str], bytes] = fetch_bytes,
) -> str:
    now = now or datetime.now(timezone.utc)
    if cache_is_fresh(path, now, max_age):
        return path
    try:
        body = fetch(url)
    except SourceError:
        if os.path.isfile(path):
            return path
        raise
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "wb") as handle:
        handle.write(body)
    return path


def tonnes_per_year_to_kg_hr(tonnes: float) -> float:
    return tonnes * 1000.0 / HOURS_PER_YEAR


def gg_per_year_to_kg_hr(gigagrams: float) -> float:
    return tonnes_per_year_to_kg_hr(gigagrams * 1000.0)


def co2e_tonnes_to_ch4_kg_hr(co2e_tonnes: float, gwp: float = GHGRP_CH4_GWP) -> float:
    return tonnes_per_year_to_kg_hr(co2e_tonnes / gwp)


def flux_kg_m2_s_to_kg_hr(flux: float, area_m2: float) -> float:
    return flux * area_m2 * SECONDS_PER_HOUR


def year_start(year: int) -> datetime:
    return datetime(year, 1, 1, tzinfo=timezone.utc)


def in_bbox(plume: Plume, bbox: BBox) -> bool:
    return bbox.contains(plume.lon, plume.lat)


def capped_collection(
    source: str,
    items: List[Plume],
    attribution: str,
    license_url: str,
    limit: int = MAX_ITEMS,
) -> PlumeCollection:
    items = sorted(items, key=lambda plume: (plume.emission_kg_hr or 0, plume.observed_at), reverse=True)
    return PlumeCollection(
        source=source,
        items=items[:limit],
        total_count=len(items),
        truncated=len(items) > limit,
        attribution=attribution,
        license_url=license_url,
    )


def row_value(row: dict, *needles: str, exclude: tuple = ()) -> str:
    for key, value in row.items():
        lowered = (key or "").lower()
        if any(token.lower() in lowered for token in exclude):
            continue
        if all(needle.lower() in lowered for needle in needles):
            return (value or "").strip()
    return ""


def optional_float(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if text == "":
        return None
    return float(text)
