"""Carbon Mapper public plume catalog (no API key required)."""
import json
from datetime import datetime, timezone
from typing import Callable, List, Optional, Tuple
from urllib.parse import urlencode

from backend.schemas import BBox, DateRange, Plume, PlumeCollection, SourceError
from backend.sources.http import fetch_bytes

SOURCE_ID = "carbonmapper"
API_URL = "https://api.carbonmapper.org/api/v1/catalog/plumes/annotated"
ATTRIBUTION = "Source: Carbon Mapper"
LICENSE_URL = "https://carbonmapper.org/terms"
MAX_PLUMES = 200


def build_query(bbox: BBox, dates: DateRange, limit: int = MAX_PLUMES) -> List[Tuple[str, str]]:
    """Repeated bbox keys in west, south, east, north order; the comma form is rejected upstream."""
    return [
        ("plume_gas", "CH4"),
        ("status", "published"),
        ("sort", "desc"),
        ("limit", str(limit)),
        ("bbox", str(bbox.west)),
        ("bbox", str(bbox.south)),
        ("bbox", str(bbox.east)),
        ("bbox", str(bbox.north)),
        ("datetime", dates.as_rfc3339_interval()),
    ]


def build_url(bbox: BBox, dates: DateRange, limit: int = MAX_PLUMES) -> str:
    return f"{API_URL}?{urlencode(build_query(bbox, dates, limit))}"


def _optional_float(value) -> Optional[float]:
    return None if value is None else float(value)


def _point(item: dict) -> Tuple[float, float]:
    geometry = item.get("geometry_json") or {}
    coordinates = geometry.get("coordinates")
    if coordinates and len(coordinates) >= 2:
        return float(coordinates[0]), float(coordinates[1])
    bounds = item.get("plume_bounds")
    if bounds and len(bounds) == 4:
        return (bounds[0] + bounds[2]) / 2, (bounds[1] + bounds[3]) / 2
    raise SourceError(f"Carbon Mapper plume {item.get('plume_id')} has no location")


def parse_plume(item: dict) -> Plume:
    lon, lat = _point(item)
    details = {
        key: str(item[field])
        for key, field in (("Platform", "platform"), ("Scene", "scene_id"), ("Published", "published_at"))
        if item.get(field)
    }
    return Plume(
        id=str(item["plume_id"]),
        source=SOURCE_ID,
        observed_at=datetime.fromisoformat(item["scene_timestamp"].replace("Z", "+00:00")).astimezone(timezone.utc),
        lon=lon,
        lat=lat,
        emission_kg_hr=_optional_float(item.get("emission_auto")),
        uncertainty_kg_hr=_optional_float(item.get("emission_uncertainty_auto")),
        instrument=item.get("instrument"),
        sector=item.get("sector"),
        image_url=item.get("plume_png"),
        attribution=ATTRIBUTION,
        details=details,
    )


def parse_response(body: bytes) -> PlumeCollection:
    try:
        payload = json.loads(body)
    except ValueError as exc:
        raise SourceError("Carbon Mapper returned malformed JSON") from exc
    items = [parse_plume(item) for item in payload.get("items") or []]
    total_count = int(payload.get("bbox_count") or len(items))
    return PlumeCollection(
        source=SOURCE_ID,
        items=items,
        total_count=total_count,
        truncated=total_count > len(items),
        attribution=ATTRIBUTION,
        license_url=LICENSE_URL,
    )


def fetch_plumes(
    bbox: BBox,
    dates: DateRange,
    fetch: Callable[[str], bytes] = fetch_bytes,
) -> PlumeCollection:
    return parse_response(fetch(build_url(bbox, dates)))
