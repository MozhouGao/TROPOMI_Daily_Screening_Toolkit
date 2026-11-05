"""GHGSat landfill plumes from Dogniaux et al. (Zenodo 16641834)."""
import csv
import io
import os
from datetime import datetime, timedelta, timezone
from typing import Callable, List, Optional

from backend.schemas import BBox, DateRange, Plume, PlumeCollection
from backend.sources.file_cache import (
    capped_collection,
    ensure_file,
    in_bbox,
    optional_float,
)
from backend.sources.http import fetch_bytes

SOURCE_ID = "ghgsat"
EXPORT_URL = "https://zenodo.org/records/16641834/files/GHGSat_detected_plumes.csv?download=1"
ATTRIBUTION = "Source: GHGSat landfill survey (Dogniaux et al. 2025)"
LICENSE_URL = "https://doi.org/10.5281/zenodo.16641834"
CACHE_FILENAME = "GHGSat_detected_plumes.csv"
CACHE_MAX_AGE = timedelta(days=30)


def default_cache_dir() -> str:
    return os.path.join(os.getcwd(), "GHGSat_data")


def parse_row(row: dict) -> Optional[Plume]:
    lat = optional_float(row.get("lat"))
    lon = optional_float(row.get("lon"))
    if lat is None or lon is None or (lat == 0.0 and lon == 0.0):
        return None
    rate_t_hr = optional_float(row.get("Q_t_per_hr")) or 0.0
    if rate_t_hr <= 0:
        return None
    observed = datetime.fromisoformat(row["date"].strip().replace(" ", "T")).replace(
        tzinfo=timezone.utc
    )
    sat_id = (row.get("sat_ID") or "").strip()
    return Plume(
        id=(row.get("obs_ID") or "").strip(),
        source=SOURCE_ID,
        observed_at=observed,
        lon=lon,
        lat=lat,
        emission_kg_hr=rate_t_hr * 1000.0,
        uncertainty_kg_hr=(
            None
            if optional_float(row.get("Q_error_t_per_hr")) is None
            else optional_float(row.get("Q_error_t_per_hr")) * 1000.0
        ),
        instrument=f"GHGSat-C{sat_id}" if sat_id else "GHGSat",
        sector="solid waste",
        image_url=None,
        attribution=ATTRIBUTION,
        details={
            key: value
            for key, value in (
                ("Site", (row.get("site_ID") or "").strip()),
                ("Plume file", (row.get("plume_tif_file_name") or "").strip()),
            )
            if value
        },
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
