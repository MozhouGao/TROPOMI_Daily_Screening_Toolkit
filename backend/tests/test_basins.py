import unittest

from backend.schemas import BBox
from backend.sources import basins

PERMIAN_BOX = BBox(west=-102.8, south=31.6, east=-101.6, north=32.4)
WESTERN_CANADA = BBox(west=-127.96, south=42.31, east=-100.18, north=58.39)
PACIFIC = BBox(west=-160.0, south=0.0, east=-150.0, north=10.0)


class BasinOverlapTests(unittest.TestCase):
    def test_permian_box_selects_permian(self):
        names = {basin.name for basin in basins.overlapping(PERMIAN_BOX)}
        self.assertIn("Permian", names)

    def test_western_canada_box_selects_montney_and_wcsb(self):
        names = {basin.name for basin in basins.overlapping(WESTERN_CANADA)}
        self.assertIn("Montney", names)
        self.assertIn("WCSB", names)

    def test_ocean_box_selects_none(self):
        self.assertEqual(basins.overlapping(PACIFIC), [])
