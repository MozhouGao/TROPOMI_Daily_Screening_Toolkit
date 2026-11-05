"""MethaneSAT layer stub: attribution only, no raw L3/L4 redistribution."""
from typing import List

from backend.schemas import BBox, DateRange, Plume, PlumeCollection
from backend.sources.file_cache import capped_collection

SOURCE_ID = "methanesat"
ATTRIBUTION = (
    "Data from MethaneSAT. Download the most current dataset at Google Earth Engine "
    "and/or Google Cloud. This toolkit does not redistribute raw L3 or L4."
)
LICENSE_URL = (
    "https://www.methanesat.org/sites/default/files/2025-02/"
    "MethaneSAT%20-%20Content%20License%20Terms%20of%20Use%20%28Revised%202-12-2025%29%5B25%5D.pdf"
)


def fetch_plumes(bbox: BBox, dates: DateRange) -> PlumeCollection:
    del bbox, dates
    items: List[Plume] = []
    return capped_collection(SOURCE_ID, items, ATTRIBUTION, LICENSE_URL)
