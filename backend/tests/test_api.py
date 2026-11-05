import unittest
from datetime import date, datetime, timezone
from unittest import mock

from fastapi.testclient import TestClient

from backend import main
from backend.schemas import BBox, DateRange, Plume, PlumeCollection, SourceError

QUERY = {
    "west": -102.8, "south": 31.6, "east": -101.6, "north": 32.4,
    "start": "2025-06-01", "end": "2025-06-30",
}


def _collection(source):
    plume = Plume(
        id="p1",
        source=source,
        observed_at=datetime(2025, 6, 25, 18, 13, 4, tzinfo=timezone.utc),
        lon=-101.8,
        lat=32.2,
        emission_kg_hr=554.5,
        attribution="Source: test",
    )
    return PlumeCollection(
        source=source, items=[plume], total_count=1, truncated=False,
        attribution="Source: test", license_url="https://example.org/terms",
    )


class PlumeEndpointTests(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(main.app)

    def test_health(self):
        self.assertEqual(self.client.get("/api/health").json(), {"status": "ok"})

    def test_dispatches_to_source_with_bbox_and_dates(self):
        fake = mock.Mock(return_value=_collection("carbonmapper"))
        with mock.patch.dict(main.PLUME_SOURCES, {"carbonmapper": fake}):
            response = self.client.get("/api/plumes/carbonmapper", params=QUERY)
        self.assertEqual(response.status_code, 200)
        fake.assert_called_once_with(
            BBox(-102.8, 31.6, -101.6, 32.4), DateRange(date(2025, 6, 1), date(2025, 6, 30))
        )
        body = response.json()
        self.assertEqual(body["source"], "carbonmapper")
        self.assertEqual(body["items"][0]["emissionKgHr"], 554.5)
        self.assertEqual(body["items"][0]["observedAt"], "2025-06-25T18:13:04Z")
        self.assertEqual(body["licenseUrl"], "https://example.org/terms")

    def test_unknown_source_is_404(self):
        self.assertEqual(self.client.get("/api/plumes/nope", params=QUERY).status_code, 404)

    def test_missing_bbox_is_422(self):
        response = self.client.get("/api/plumes/mars", params={"start": "2025-06-01", "end": "2025-06-30"})
        self.assertEqual(response.status_code, 422)

    def test_inverted_bbox_is_422(self):
        response = self.client.get("/api/plumes/mars", params=dict(QUERY, west=0, east=-10))
        self.assertEqual(response.status_code, 422)

    def test_inverted_dates_are_422(self):
        response = self.client.get(
            "/api/plumes/mars", params=dict(QUERY, start="2025-07-01", end="2025-06-01")
        )
        self.assertEqual(response.status_code, 422)

    def test_source_error_is_502_with_message(self):
        fake = mock.Mock(side_effect=SourceError("upstream returned HTTP 503"))
        with mock.patch.dict(main.PLUME_SOURCES, {"mars": fake}):
            response = self.client.get("/api/plumes/mars", params=QUERY)
        self.assertEqual(response.status_code, 502)
        self.assertIn("HTTP 503", response.json()["detail"])

    def test_sources_are_registered(self):
        self.assertEqual(
            set(main.PLUME_SOURCES),
            {
                "carbonmapper",
                "mars",
                "ghgsat",
                "cusworth",
                "ghgrp_us",
                "ghgrp_ca",
                "nesser",
                "canada_nir",
                "methanesat",
            },
        )


if __name__ == "__main__":
    unittest.main()
