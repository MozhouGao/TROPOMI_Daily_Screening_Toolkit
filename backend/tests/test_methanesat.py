import unittest
from datetime import date

from backend.schemas import BBox, DateRange
from backend.sources import methanesat

BBOX = BBox(west=-114.2, south=51.0, east=-113.8, north=51.2)


class MethaneSatTests(unittest.TestCase):
    def test_returns_empty_collection_without_raw_points(self):
        collection = methanesat.fetch_plumes(
            BBOX, DateRange(date(2024, 6, 1), date(2025, 6, 1))
        )
        self.assertEqual(collection.source, "methanesat")
        self.assertEqual(collection.items, [])
        self.assertEqual(collection.total_count, 0)
        self.assertIn("does not redistribute raw L3 or L4", collection.attribution)
        self.assertIn("methanesat.org", collection.license_url)
