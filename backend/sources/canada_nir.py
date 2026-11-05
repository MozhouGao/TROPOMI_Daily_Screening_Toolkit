"""Canada 2018 NIR gridded methane (Scarpelli et al.) as 0.1° inventory cells."""
import os
from datetime import timedelta
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from netCDF4 import Dataset

from backend.schemas import BBox, DateRange, Plume, PlumeCollection
from backend.sources.file_cache import (
    capped_collection,
    ensure_file,
    flux_kg_m2_s_to_kg_hr,
    year_start,
)
from backend.sources.http import fetch_bytes

SOURCE_ID = "canada_nir"
ATTRIBUTION = "Source: Scarpelli et al. 2022 gridded Canada NIR (2018)"
LICENSE_URL = "https://doi.org/10.7910/DVN/CC3KLO"
CACHE_MAX_AGE = timedelta(days=365)
MIN_KG_HR = 1.0
MIN_PLAUSIBLE_AREA_M2 = 1.0e6
METERS_PER_DEG_LAT = 110574.0
METERS_PER_DEG_LON = 111320.0
BASE_URL = "http://geoschemdata.wustl.edu/ExtData/HEMCO/CH4/v2022-01/Scarpelli_Canada"
SECTORS: Dict[str, str] = {
    "solid waste": "can_emis_solid_waste_2018.nc",
    "wastewater": "can_emis_wastewater_2018.nc",
    "residential combustion": "can_emis_residential_combustion_2018.nc",
    "oil and gas leakage": "can_emis_oil_gas_leakage_2018.nc",
}
Grid = Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
LoadGrid = Callable[[str], Grid]


def default_cache_dir() -> str:
    return os.path.join(os.getcwd(), "CanadaNIR_data")


def cell_area_m2(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    dlat = abs(float(np.diff(lats).mean())) if lats.size > 1 else 0.1
    dlon = abs(float(np.diff(lons).mean())) if lons.size > 1 else 0.1
    lon_m = METERS_PER_DEG_LON * dlon * np.cos(np.radians(lats))
    return np.outer(lon_m * METERS_PER_DEG_LAT * dlat, np.ones(lons.size))


def resolve_area(lats: np.ndarray, lons: np.ndarray, area: np.ndarray) -> np.ndarray:
    values = np.asarray(area, dtype=float)
    if values.size == 0 or float(np.nanmedian(values)) < MIN_PLAUSIBLE_AREA_M2:
        return cell_area_m2(lats, lons)
    return values


def cells_from_grid(
    lons: np.ndarray,
    lats: np.ndarray,
    flux: np.ndarray,
    area: np.ndarray,
    sector: str,
    bbox: BBox,
    min_kg_hr: float = MIN_KG_HR,
) -> List[Plume]:
    lon_idx = np.where((lons >= bbox.west) & (lons <= bbox.east))[0]
    lat_idx = np.where((lats >= bbox.south) & (lats <= bbox.north))[0]
    if lon_idx.size == 0 or lat_idx.size == 0:
        return []
    flux_slice = np.atleast_2d(flux)[np.ix_(lat_idx, lon_idx)]
    area_slice = np.atleast_2d(area)[np.ix_(lat_idx, lon_idx)]
    items: List[Plume] = []
    for lat_i, lat in enumerate(lats[lat_idx]):
        for lon_i, lon in enumerate(lons[lon_idx]):
            rate = float(flux_kg_m2_s_to_kg_hr(float(flux_slice[lat_i, lon_i]), float(area_slice[lat_i, lon_i])))
            if rate < min_kg_hr:
                continue
            items.append(
                Plume(
                    id=f"{sector}:{lat:.2f}:{lon:.2f}",
                    source=SOURCE_ID,
                    observed_at=year_start(2018),
                    lon=float(lon),
                    lat=float(lat),
                    emission_kg_hr=rate,
                    instrument="Canada NIR 2018",
                    sector=sector,
                    image_url=None,
                    attribution=ATTRIBUTION,
                    details={"Grid": "0.1°"},
                )
            )
    return items


def read_grid(path: str) -> Grid:
    with Dataset(path) as dataset:
        lons = np.asarray(dataset.variables["lon"][:], dtype=float)
        lats = np.asarray(dataset.variables["lat"][:], dtype=float)
        total_name = next(name for name in dataset.variables if name.endswith("_total"))
        flux = np.asarray(dataset.variables[total_name][:], dtype=float)
        area = np.asarray(dataset.variables["area"][:], dtype=float)
    return lons, lats, np.squeeze(flux), resolve_area(lats, lons, np.squeeze(area))


def fetch_plumes(
    bbox: BBox,
    dates: DateRange,
    cache_dir: Optional[str] = None,
    fetch: Callable[[str], bytes] = fetch_bytes,
    load_grid: LoadGrid = read_grid,
) -> PlumeCollection:
    del dates
    folder = cache_dir or default_cache_dir()
    items: List[Plume] = []
    for sector, filename in SECTORS.items():
        path = ensure_file(
            f"{BASE_URL}/{filename}",
            os.path.join(folder, filename),
            CACHE_MAX_AGE,
            fetch=fetch,
        )
        lons, lats, flux, area = load_grid(path)
        items.extend(cells_from_grid(lons, lats, flux, area, sector, bbox))
    return capped_collection(SOURCE_ID, items, ATTRIBUTION, LICENSE_URL)
