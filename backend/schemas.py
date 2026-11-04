"""Shared shapes that cross the API boundary for every plume source."""
from datetime import date, datetime
from typing import Dict, List, NamedTuple, Optional

from pydantic import BaseModel, ConfigDict
from pydantic.alias_generators import to_camel


class BBox(NamedTuple):
    west: float
    south: float
    east: float
    north: float

    def contains(self, lon: float, lat: float) -> bool:
        return self.west <= lon <= self.east and self.south <= lat <= self.north


class DateRange(NamedTuple):
    start: date
    end: date

    def contains(self, moment: datetime) -> bool:
        return self.start <= moment.date() <= self.end

    def as_rfc3339_interval(self) -> str:
        return f"{self.start.isoformat()}T00:00:00Z/{self.end.isoformat()}T23:59:59Z"


class CamelModel(BaseModel):
    model_config = ConfigDict(alias_generator=to_camel, populate_by_name=True)


class Plume(CamelModel):
    id: str
    source: str
    observed_at: datetime
    lon: float
    lat: float
    emission_kg_hr: Optional[float] = None
    uncertainty_kg_hr: Optional[float] = None
    instrument: Optional[str] = None
    sector: Optional[str] = None
    image_url: Optional[str] = None
    attribution: str
    details: Dict[str, str] = {}


class PlumeCollection(CamelModel):
    source: str
    items: List[Plume]
    total_count: int
    truncated: bool
    attribution: str
    license_url: str


class GridOverlay(CamelModel):
    image_url: str
    west: float
    south: float
    east: float
    north: float


class TropomiScreenResult(PlumeCollection):
    overlay: Optional[GridOverlay] = None
    message: str
    figure_path: Optional[str] = None
    screened_days: int
    plume_days: int


class CredentialStatus(CamelModel):
    configured: bool
    access_key_id: str = ""
    file_text: str = ""
    example: str = ""


class CredentialWrite(CamelModel):
    access_key_id: Optional[str] = None
    secret_access_key: Optional[str] = None
    file_text: Optional[str] = None


class DateWindow(CamelModel):
    start: date
    end: date


class ScreenRequest(CamelModel):
    west: float
    south: float
    east: float
    north: float
    start: date
    end: date
    threshold_delta: int = 15
    min_pixels: int = 1


class SourceError(Exception):
    """Raised when an upstream provider cannot be read."""
