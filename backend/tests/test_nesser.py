import tempfile
import unittest
from datetime import date

from backend.schemas import BBox, DateRange
from backend.sources import nesser
from backend.sources.file_cache import gg_per_year_to_kg_hr

CSV = (
    ",lat_center,lon_center,ghgrp,id,xa,post_mean,post_min,post_max,dofs,frac,"
    "year,name,lat,lon,city,county,state,company,subparts,diff_abs,diff_rel\n"
    "0,41.0,-85.3125,3.41316,1007520,9.64,44.313277,10.18,14.28,0.28,0.88,2019,"
    "NATIONAL SERV-ALL LANDFILL,41.031121,-85.216174,FORT WAYNE,ALLEN COUNTY,IN,"
    "REPUBLIC SERVICES INC (100%),\"C,HH\",40.9,11.98\n"
    "1,51.04,-114.07,1.0,999,1.0,2.0,1.0,3.0,0.2,0.5,2019,OUTSIDE,51.04,-114.07,"
    "CALGARY,DIVISION NO. 6,AB,TEST,HH,1.0,1.0\n"
)
BBOX = BBox(west=-86.0, south=40.5, east=-85.0, north=41.5)


class NesserTests(unittest.TestCase):
    def test_parse_converts_gg_year_to_kg_hr(self):
        plumes = nesser.parse_csv(CSV)
        self.assertEqual(len(plumes), 2)
        plume = plumes[0]
        self.assertEqual(plume.id, "1007520")
        self.assertEqual(plume.source, "nesser")
        self.assertEqual(plume.observed_at.year, 2019)
        self.assertAlmostEqual(plume.lat, 41.031121)
        self.assertAlmostEqual(plume.lon, -85.216174)
        self.assertAlmostEqual(plume.emission_kg_hr, gg_per_year_to_kg_hr(44.313277))
        self.assertEqual(plume.instrument, "TROPOMI inversion")
        self.assertEqual(plume.sector, "solid waste")
        self.assertEqual(plume.details["Name"], "NATIONAL SERV-ALL LANDFILL")
        self.assertEqual(plume.details["GHGRP ID"], "1007520")

    def test_fetch_filters_bbox(self):
        with tempfile.TemporaryDirectory() as tmp:
            collection = nesser.fetch_plumes(
                BBOX,
                DateRange(date(2026, 1, 1), date(2026, 1, 31)),
                cache_dir=tmp,
                fetch=lambda url: CSV.encode("utf-8"),
            )
        self.assertEqual([plume.id for plume in collection.items], ["1007520"])

    def test_midwest_site_is_outside_western_canada_box(self):
        western = BBox(west=-127.96, south=42.31, east=-100.18, north=58.39)
        midwest_only = "\n".join(CSV.strip().splitlines()[:2]) + "\n"
        with tempfile.TemporaryDirectory() as tmp:
            collection = nesser.fetch_plumes(
                western,
                DateRange(date(2026, 1, 1), date(2026, 1, 31)),
                cache_dir=tmp,
                fetch=lambda url: midwest_only.encode("utf-8"),
            )
        self.assertEqual(collection.items, [])


if __name__ == "__main__":
    unittest.main()
