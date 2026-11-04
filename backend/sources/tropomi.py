"""TROPOMI download, credentials, and patch screening."""
from datetime import datetime, timezone
from typing import Callable, Tuple

from TROPOMI_toolkit import (
    AWS_KEYS_EXAMPLE,
    Load_CH4,
    aws_keys_are_configured,
    download_TROPOMI_CH4_L2_data,
    generate_results,
    iter_screening_dates,
    load_aws_keys,
    screening_plumes,
    write_aws_keys,
    write_aws_keys_text,
)
from backend.schemas import (
    BBox,
    CredentialStatus,
    CredentialWrite,
    DateRange,
    GridOverlay,
    Plume,
    SourceError,
    TropomiScreenResult,
)

SOURCE_ID = "tropomi"
ATTRIBUTION = "Source: Copernicus Sentinel-5P TROPOMI"
LICENSE_URL = "https://documentation.dataspace.copernicus.eu/"
DEFAULT_THRESHOLD_DELTA = 15
DEFAULT_MIN_PIXELS = 1

LoadCh4 = Callable[..., Tuple]


def credential_status(keys_path=None) -> CredentialStatus:
    access_key, _, raw_text = load_aws_keys(keys_path)
    return CredentialStatus(
        configured=aws_keys_are_configured(keys_path),
        access_key_id=access_key,
        file_text=raw_text or AWS_KEYS_EXAMPLE,
        example=AWS_KEYS_EXAMPLE,
    )


def save_credentials(payload: CredentialWrite, keys_path=None) -> CredentialStatus:
    if payload.file_text:
        write_aws_keys_text(payload.file_text, keys_path)
    elif payload.access_key_id and payload.secret_access_key:
        write_aws_keys(payload.access_key_id, payload.secret_access_key, keys_path)
    else:
        raise ValueError("Provide access_key_id and secret_access_key, or file_text.")
    return credential_status(keys_path)


def download_observations(dates: DateRange) -> str:
    try:
        return download_TROPOMI_CH4_L2_data(dates.start, dates.end)
    except (FileNotFoundError, ValueError) as exc:
        raise SourceError(str(exc)) from exc


def _as_plume(record: dict, observed_at: datetime) -> Plume:
    return Plume(
        id=record["id"],
        source=SOURCE_ID,
        observed_at=observed_at,
        lon=record["lon"],
        lat=record["lat"],
        emission_kg_hr=record.get("emission_kg_hr"),
        instrument="TROPOMI",
        sector=None,
        image_url=None,
        attribution=ATTRIBUTION,
        details={"Enhancement (ppb)": f"{record['enhancement_ppb']:.1f}"},
    )


def screen_region(
    bbox: BBox,
    dates: DateRange,
    threshold_delta: int = DEFAULT_THRESHOLD_DELTA,
    min_pixels: int = DEFAULT_MIN_PIXELS,
    load_ch4: LoadCh4 = Load_CH4,
    screen: Callable = screening_plumes,
    write_results: Callable = generate_results,
) -> TropomiScreenResult:
    items = []
    overlay = None
    figure_path = None
    screened_days = 0
    plume_days = 0
    for screening_date in iter_screening_dates(dates.start, dates.end):
        try:
            grid_lon, grid_lat, fch4, fwind, fpressure = load_ch4(
                bbox.south, bbox.north, bbox.west, bbox.east, screening_date, qa_pass=0.5
            )
        except FileNotFoundError as exc:
            raise SourceError(str(exc)) from exc
        screened_days += 1
        detected = screen(
            fch4, fwind, fpressure, grid_lon, grid_lat, threshold_delta, min_pixels
        )
        result = write_results(
            grid_lon,
            grid_lat,
            fch4,
            *detected,
            screening_date.strftime("%Y%m%d"),
        )
        figure_path = result["figure_path"]
        observed_at = datetime(screening_date.year, screening_date.month, screening_date.day, tzinfo=timezone.utc)
        if result.get("overlay"):
            overlay = GridOverlay.model_validate(result["overlay"])
        if detected[0]:
            plume_days += 1
            items.extend(_as_plume(record, observed_at) for record in result["hotspots"])

    if screened_days == 0:
        message = "No dates to screen."
    elif plume_days == 0:
        message = "No ultra-emitter event detected in your region."
    else:
        message = (
            f"Suspect plumes found on {plume_days} of {screened_days} day(s). "
            "Maps and CSVs are in ./assets."
        )
    return TropomiScreenResult(
        source=SOURCE_ID,
        items=items,
        total_count=len(items),
        truncated=False,
        attribution=ATTRIBUTION,
        license_url=LICENSE_URL,
        overlay=overlay,
        message=message,
        figure_path=figure_path,
        screened_days=screened_days,
        plume_days=plume_days,
    )
