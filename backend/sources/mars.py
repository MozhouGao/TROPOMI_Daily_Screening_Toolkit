"""UNEP IMEO Methane Alert and Response System (MARS) public plume export.

There is no query API; UNEP publishes one zip with every plume. The zip is cached
locally and refreshed weekly, then filtered in memory per request.
"""
import csv
import io
import os
import zipfile
from datetime import datetime, timedelta, timezone
from functools import lru_cache
from typing import Callable, List, Optional

from backend.schemas import BBox, DateRange, Plume, PlumeCollection, SourceError
from backend.sources.http import fetch_bytes

SOURCE_ID = "mars"
EXPORT_URL = (
    "https://unepazeconomyadlsstorage.blob.core.windows.net/public/"
    "unep_methanedata_detected_plumes_csv.zip"
)
ATTRIBUTION = "Source: UNEP IMEO Eye on Methane (CC BY-NC-SA 4.0)"
LICENSE_URL = "https://methanedata.unep.org/dict-mars-plumes"
PLUMES_CSV_NAME = "unep_methanedata_detected_plumes.csv"
CACHE_DIRNAME = "MARS_data"
CACHE_FILENAME = "unep_methanedata_detected_plumes_csv.zip"
CACHE_MAX_AGE = timedelta(days=7)


def default_cache_dir() -> str:
    return os.path.join(os.getcwd(), CACHE_DIRNAME)


def cache_is_fresh(path: str, now: datetime) -> bool:
    if os.path.isfile(path):
        modified = datetime.fromtimestamp(os.path.getmtime(path), tz=timezone.utc)
        return now - modified < CACHE_MAX_AGE
    return False


def ensure_export(
    cache_dir: str,
    now: Optional[datetime] = None,
    fetch: Callable[[str], bytes] = fetch_bytes,
) -> str:
    """Return the path of a MARS export zip, downloading when missing or older than a week."""
    now = now or datetime.now(timezone.utc)
    path = os.path.join(cache_dir, CACHE_FILENAME)
    if cache_is_fresh(path, now):
        return path
    try:
        body = fetch(EXPORT_URL)
    except SourceError:
        if os.path.isfile(path):
            return path
        raise
    os.makedirs(cache_dir, exist_ok=True)
    with open(path, "wb") as handle:
        handle.write(body)
    return path


def _optional_float(value: str) -> Optional[float]:
    stripped = (value or "").strip()
    return float(stripped) if stripped else None


def parse_row(row: dict) -> Plume:
    details = {
        label: row[column].strip()
        for label, column in (
            ("Country", "country"),
            ("Source", "source_name"),
            ("Detected by", "detection_institution"),
            ("Actionable", "actionable"),
            ("Notified", "notified"),
            ("Satellite product", "tile"),
        )
        if row.get(column, "").strip()
    }
    return Plume(
        id=row["id_plume"],
        source=SOURCE_ID,
        observed_at=datetime.fromisoformat(row["tile_date"]).replace(tzinfo=timezone.utc),
        lon=float(row["lon"]),
        lat=float(row["lat"]),
        emission_kg_hr=_optional_float(row.get("ch4_fluxrate")),
        uncertainty_kg_hr=_optional_float(row.get("ch4_fluxrate_std")),
        instrument=row.get("satellite") or None,
        sector=row.get("sector") or None,
        image_url=None,
        attribution=ATTRIBUTION,
        details=details,
    )


def read_plumes(zip_path: str) -> List[Plume]:
    try:
        with zipfile.ZipFile(zip_path) as archive:
            text = archive.read(PLUMES_CSV_NAME).decode("utf-8-sig")
    except (zipfile.BadZipFile, KeyError, OSError) as exc:
        raise SourceError(f"MARS export unreadable: {exc}") from exc
    return [parse_row(row) for row in csv.DictReader(io.StringIO(text))]


@lru_cache(maxsize=2)
def _cached_plumes(zip_path: str, modified_ns: int) -> List[Plume]:
    return read_plumes(zip_path)


def load_plumes(zip_path: str) -> List[Plume]:
    return _cached_plumes(zip_path, os.stat(zip_path).st_mtime_ns)


def filter_plumes(plumes: List[Plume], bbox: BBox, dates: DateRange) -> List[Plume]:
    return [
        plume
        for plume in plumes
        if bbox.contains(plume.lon, plume.lat) and dates.contains(plume.observed_at)
    ]


def fetch_plumes(
    bbox: BBox,
    dates: DateRange,
    cache_dir: Optional[str] = None,
    fetch: Callable[[str], bytes] = fetch_bytes,
) -> PlumeCollection:
    zip_path = ensure_export(cache_dir or default_cache_dir(), fetch=fetch)
    items = filter_plumes(load_plumes(zip_path), bbox, dates)
    items.sort(key=lambda plume: plume.observed_at, reverse=True)
    return PlumeCollection(
        source=SOURCE_ID,
        items=items,
        total_count=len(items),
        truncated=False,
        attribution=ATTRIBUTION,
        license_url=LICENSE_URL,
    )
