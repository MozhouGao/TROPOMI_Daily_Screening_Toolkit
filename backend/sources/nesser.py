"""Nesser et al. 2024 TROPOMI inversion landfill posteriors."""
import csv
import io
import os
from datetime import timedelta
from typing import Callable, List, Optional

from backend.schemas import BBox, DateRange, Plume, PlumeCollection
from backend.sources.file_cache import (
    capped_collection,
    ensure_file,
    gg_per_year_to_kg_hr,
    in_bbox,
    optional_float,
    year_start,
)
from backend.sources.http import fetch_bytes

SOURCE_ID = "nesser"
EXPORT_URL = (
    "https://raw.githubusercontent.com/hannahnesser/TROPOMI_inversion/v1.0/"
    "inversion_data/landfills/ghgrp_processed.csv"
)
ATTRIBUTION = "Source: Nesser et al. 2024 TROPOMI inversion (CC BY 4.0)"
LICENSE_URL = "https://doi.org/10.5281/zenodo.10946769"
CACHE_FILENAME = "ghgrp_processed.csv"
CACHE_MAX_AGE = timedelta(days=30)


def default_cache_dir() -> str:
    return os.path.join(os.getcwd(), "Nesser_data")


def parse_row(row: dict) -> Optional[Plume]:
    lat = optional_float(row.get("lat"))
    lon = optional_float(row.get("lon"))
    posterior = optional_float(row.get("post_mean"))
    facility_id = (row.get("id") or "").strip()
    if lat is None or lon is None or posterior is None or not facility_id:
        return None
    name = (row.get("name") or "").strip()
    year = int(optional_float(row.get("year")) or 2019)
    return Plume(
        id=facility_id,
        source=SOURCE_ID,
        observed_at=year_start(year),
        lon=lon,
        lat=lat,
        emission_kg_hr=gg_per_year_to_kg_hr(posterior),
        uncertainty_kg_hr=None,
        instrument="TROPOMI inversion",
        sector="solid waste",
        image_url=None,
        attribution=ATTRIBUTION,
        details={
            key: value
            for key, value in (
                ("Name", name),
                ("GHGRP ID", facility_id),
                ("City", (row.get("city") or "").strip()),
                ("State", (row.get("state") or "").strip()),
                ("GHGRP Gg/yr", (row.get("ghgrp") or "").strip()),
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
