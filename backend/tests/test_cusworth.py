import tempfile
import unittest
from datetime import date, datetime, timezone

from backend.schemas import BBox, DateRange
from backend.sources import cusworth

CSV = (
    "GHGRP ID,Latitude of Plume Origin,Longitude of Plume Origin,"
    "Plume Candidate ID (UTC time of detection),"
    "Instantaneous Emission Rate (kg/h),Instantaneous Emission Rate Uncertainty (kg/h)\n"
    "1005717,29.935453,-90.262612,GAO20211026t174729p0000-B,10151,5395\n"
    "1007607,34.5948,-117.26764,ang20160911t174834-C,61,16\n"
)
BBOX = BBox(west=-91.0, south=29.0, east=-90.0, north=31.0)


class CusworthTests(unittest.TestCase):
    def test_parse_reads_ghgrp_id_and_timestamp(self):
        plumes = cusworth.parse_csv(CSV)
        self.assertEqual(len(plumes), 2)
        plume = plumes[0]
        self.assertEqual(plume.id, "GAO20211026t174729p0000-B")
        self.assertEqual(plume.source, "cusworth")
        self.assertEqual(plume.observed_at, datetime(2021, 10, 26, 17, 47, 29, tzinfo=timezone.utc))
        self.assertEqual(plume.emission_kg_hr, 10151.0)
        self.assertEqual(plume.uncertainty_kg_hr, 5395.0)
        self.assertEqual(plume.instrument, "GAO")
        self.assertEqual(plume.sector, "solid waste")
        self.assertEqual(plume.details["GHGRP ID"], "1005717")
        self.assertEqual(plumes[1].observed_at, datetime(2016, 9, 11, 17, 48, 34, tzinfo=timezone.utc))
        self.assertEqual(plumes[1].instrument, "AVIRIS-NG")

    def test_fetch_filters_bbox(self):
        with tempfile.TemporaryDirectory() as tmp:
            collection = cusworth.fetch_plumes(
                BBOX,
                DateRange(date(2026, 1, 1), date(2026, 1, 31)),
                cache_dir=tmp,
                fetch=lambda url: CSV.encode("utf-8"),
            )
        self.assertEqual([plume.id for plume in collection.items], ["GAO20211026t174729p0000-B"])

    def test_western_canada_box_is_outside_campaign_extent(self):
        western = BBox(west=-127.96, south=42.31, east=-100.18, north=58.39)
        with tempfile.TemporaryDirectory() as tmp:
            collection = cusworth.fetch_plumes(
                western,
                DateRange(date(2026, 1, 1), date(2026, 1, 31)),
                cache_dir=tmp,
                fetch=lambda url: CSV.encode("utf-8"),
            )
        self.assertEqual(collection.items, [])


if __name__ == "__main__":
    unittest.main()
