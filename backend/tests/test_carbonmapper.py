import json
import unittest
from datetime import date, timezone
from unittest import mock
from urllib.error import HTTPError
from urllib.parse import parse_qsl, urlsplit

from backend.schemas import BBox, DateRange, SourceError
from backend.sources import carbonmapper
from backend.sources.http import fetch_bytes

BBOX = BBox(west=-102.8, south=31.6, east=-101.6, north=32.4)
DATES = DateRange(date(2025, 6, 1), date(2025, 6, 30))

SYNTHETIC_ITEM = {
    "plume_id": "tan20250625t181304c33s4001-A",
    "scene_id": "9805d37c",
    "gas": "CH4",
    "scene_timestamp": "2025-06-25T18:13:04.330Z",
    "instrument": "tan",
    "platform": "Planet Tanager 1",
    "sector": "1B2",
    "emission_auto": 554.4959968936353,
    "emission_uncertainty_auto": 105.06958340294109,
    "geometry_json": {"type": "Point", "coordinates": [-101.828, 32.264]},
    "plume_bounds": [-101.849, 32.252, -101.812, 32.283],
    "plume_png": "https://catalog.carbonmapper.org/x/plume.png",
    "published_at": "2025-07-25T18:21:53.401Z",
}


def _payload(items, bbox_count=None):
    return json.dumps({
        "bbox_count": len(items) if bbox_count is None else bbox_count,
        "total_count": 907,
        "limit": 200,
        "offset": 0,
        "items": items,
        "nearby_items": [],
    }).encode()


class QueryTests(unittest.TestCase):
    def test_bbox_is_repeated_in_west_south_east_north_order(self):
        query = carbonmapper.build_query(BBOX, DATES)
        self.assertEqual(
            [value for key, value in query if key == "bbox"],
            ["-102.8", "31.6", "-101.6", "32.4"],
        )

    def test_datetime_is_a_full_day_interval(self):
        query = dict(carbonmapper.build_query(BBOX, DATES))
        self.assertEqual(query["datetime"], "2025-06-01T00:00:00Z/2025-06-30T23:59:59Z")
        self.assertEqual(query["plume_gas"], "CH4")
        self.assertEqual(query["limit"], str(carbonmapper.MAX_PLUMES))

    def test_url_keeps_repeated_keys(self):
        url = carbonmapper.build_url(BBOX, DATES)
        self.assertTrue(url.startswith(carbonmapper.API_URL + "?"))
        pairs = parse_qsl(urlsplit(url).query)
        self.assertEqual(sum(1 for key, _ in pairs if key == "bbox"), 4)
        self.assertNotIn("bbox=-102.8%2C", url)


class ParseTests(unittest.TestCase):
    def test_synthetic_item_maps_to_shared_shape(self):
        plume = carbonmapper.parse_plume(SYNTHETIC_ITEM)
        self.assertEqual(plume.id, "tan20250625t181304c33s4001-A")
        self.assertEqual(plume.source, "carbonmapper")
        self.assertEqual(plume.observed_at.tzinfo, timezone.utc)
        self.assertEqual(plume.observed_at.isoformat(), "2025-06-25T18:13:04.330000+00:00")
        self.assertAlmostEqual(plume.lon, -101.828)
        self.assertAlmostEqual(plume.lat, 32.264)
        self.assertAlmostEqual(plume.emission_kg_hr, 554.496, places=3)
        self.assertAlmostEqual(plume.uncertainty_kg_hr, 105.070, places=3)
        self.assertEqual(plume.instrument, "tan")
        self.assertEqual(plume.sector, "1B2")
        self.assertEqual(plume.image_url, "https://catalog.carbonmapper.org/x/plume.png")
        self.assertEqual(plume.attribution, carbonmapper.ATTRIBUTION)
        self.assertEqual(plume.details["Platform"], "Planet Tanager 1")

    def test_falls_back_to_bounds_centre_without_point(self):
        item = dict(SYNTHETIC_ITEM, geometry_json=None)
        plume = carbonmapper.parse_plume(item)
        self.assertAlmostEqual(plume.lon, (-101.849 + -101.812) / 2)
        self.assertAlmostEqual(plume.lat, (32.252 + 32.283) / 2)

    def test_missing_emission_is_none(self):
        item = dict(SYNTHETIC_ITEM, emission_auto=None, emission_uncertainty_auto=None)
        plume = carbonmapper.parse_plume(item)
        self.assertIsNone(plume.emission_kg_hr)
        self.assertIsNone(plume.uncertainty_kg_hr)

    def test_json_output_uses_camel_case(self):
        data = carbonmapper.parse_plume(SYNTHETIC_ITEM).model_dump(by_alias=True)
        self.assertIn("observedAt", data)
        self.assertIn("emissionKgHr", data)
        self.assertNotIn("emission_kg_hr", data)


class FetchTests(unittest.TestCase):
    def test_collection_marks_truncation_from_bbox_count(self):
        collection = carbonmapper.fetch_plumes(
            BBOX, DATES, fetch=lambda url: _payload([SYNTHETIC_ITEM], bbox_count=350)
        )
        self.assertEqual(len(collection.items), 1)
        self.assertEqual(collection.total_count, 350)
        self.assertTrue(collection.truncated)
        self.assertEqual(collection.license_url, carbonmapper.LICENSE_URL)

    def test_empty_result(self):
        collection = carbonmapper.fetch_plumes(BBOX, DATES, fetch=lambda url: _payload([]))
        self.assertEqual(collection.items, [])
        self.assertFalse(collection.truncated)

    def test_malformed_json_is_a_source_error(self):
        with self.assertRaises(SourceError):
            carbonmapper.fetch_plumes(BBOX, DATES, fetch=lambda url: b"<html>oops</html>")

    def test_http_error_becomes_source_error_string(self):
        error = HTTPError("https://api.carbonmapper.org/x", 422, "Unprocessable", {}, None)
        with mock.patch("backend.sources.http.urlopen", side_effect=error):
            with self.assertRaises(SourceError) as raised:
                fetch_bytes("https://api.carbonmapper.org/x")
        self.assertIn("HTTP 422", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
