import os
import tempfile
import unittest
from datetime import date, datetime, timezone

from backend.schemas import BBox, DateRange
from backend.sources import ghgsat

HEADER = (
    "site_ID,lat,lon,year,month,day,hour,minute,second,Q_t_per_hr,Q_error_t_per_hr,"
    "wind_speed_m_per_s,sat_ID,obs_ID,date,IME_kg,intermediate_results_L_m,"
    "intermediate_results_effective_wind_speed_m_per_s,conversion_ch4_ppb_to_molm2,"
    "manually_pinned_sources,plume_tif_file_name"
)
DETECTION = (
    "1,20.54346163,-103.17285782,2022,10,31,17,13,28,2.0683827242227,1.055682832583701,"
    "1.98,1,6axJ3gW,2022-10-31 17:13:28,545.37,1265.49,1.33,0.000298,"
    "(20.542700;-103.173100),C1_20221031_20221103_6axJ3gW_7367_CH4PL.tif"
)
NON_DETECT = (
    "2,0.0,0.0,2022,11,1,0,0,0,0.0,0.0,0,1,none,2022-11-01 00:00:00,0,0,0,0,,"
)
OUTSIDE = (
    "3,51.04,-114.07,2022,10,31,12,0,0,0.5,0.1,2,2,outside,2022-10-31 12:00:00,"
    "10,100,1,0.0003,(51.04;-114.07),C2_outside.tif"
)
CSV = "\n".join([HEADER, DETECTION, NON_DETECT, OUTSIDE]) + "\n"
BBOX = BBox(west=-104.0, south=20.0, east=-102.0, north=21.0)
DATES = DateRange(date(2026, 10, 1), date(2026, 10, 31))


class GhgsatTests(unittest.TestCase):
    def test_parse_converts_tonnes_per_hour_and_skips_nondetects(self):
        plumes = ghgsat.parse_csv(CSV)
        self.assertEqual([plume.id for plume in plumes], ["6axJ3gW", "outside"])
        plume = plumes[0]
        self.assertEqual(plume.source, "ghgsat")
        self.assertEqual(plume.observed_at, datetime(2022, 10, 31, 17, 13, 28, tzinfo=timezone.utc))
        self.assertAlmostEqual(plume.lon, -103.17285782)
        self.assertAlmostEqual(plume.lat, 20.54346163)
        self.assertAlmostEqual(plume.emission_kg_hr, 2068.3827242227)
        self.assertAlmostEqual(plume.uncertainty_kg_hr, 1055.682832583701)
        self.assertEqual(plume.instrument, "GHGSat-C1")
        self.assertEqual(plume.sector, "solid waste")
        self.assertEqual(plume.details["Site"], "1")

    def test_fetch_filters_bbox_and_ignores_dates(self):
        with tempfile.TemporaryDirectory() as tmp:
            collection = ghgsat.fetch_plumes(
                BBOX,
                DATES,
                cache_dir=tmp,
                fetch=lambda url: CSV.encode("utf-8"),
            )
        self.assertEqual([plume.id for plume in collection.items], ["6axJ3gW"])
        self.assertEqual(collection.source, "ghgsat")
        self.assertEqual(collection.license_url, ghgsat.LICENSE_URL)
        self.assertFalse(collection.truncated)


if __name__ == "__main__":
    unittest.main()
