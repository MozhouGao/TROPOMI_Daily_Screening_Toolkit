import os
import tempfile
import unittest
from datetime import date
from unittest import mock

import numpy as np
from fastapi.testclient import TestClient

from TROPOMI_toolkit import grid_overlay_png, hotspot_records
from backend import main
from backend.schemas import BBox, CredentialStatus, CredentialWrite, DateRange, SourceError, TropomiScreenResult
from backend.sources import tropomi


def _grid():
    lons = np.array([[-114.0, -113.95], [-114.0, -113.95]])
    lats = np.array([[51.10, 51.10], [51.05, 51.05]])
    field = np.full((2, 2), 1800.0)
    field[0, 0] = 1900.0
    return lons, lats, field, np.ones((2, 2)), np.full((2, 2), 1000.0)


class HotspotRecordTests(unittest.TestCase):
    def test_records_include_id_location_and_rate(self):
        plume = np.zeros((3, 3))
        plume[1, 1] = 40
        lons = np.array([[10.0, 10.05, 10.1], [10.0, 10.05, 10.1], [10.0, 10.05, 10.1]])
        lats = np.array([[20.1, 20.1, 20.1], [20.05, 20.05, 20.05], [20.0, 20.0, 20.0]])
        records, polygons, enhance, *_ = hotspot_records(
            [plume], [np.ones((3, 3))], [np.full((3, 3), 1000.0)], [lons], [lats], "20260830"
        )
        self.assertEqual(len(records), 1)
        self.assertEqual(len(polygons), 1)
        self.assertEqual(records[0]["id"], "20260830-1")
        self.assertAlmostEqual(records[0]["lon"], 10.05)
        self.assertAlmostEqual(records[0]["lat"], 20.05)
        self.assertGreater(records[0]["emission_kg_hr"], 0)
        self.assertEqual(enhance[0], 40)

    def test_overlay_png_is_a_data_url_with_bounds(self):
        lons, lats, field, *_ = _grid()
        overlay = grid_overlay_png(lons, lats, field)
        self.assertTrue(overlay["image_url"].startswith("data:image/png;base64,"))
        self.assertEqual(overlay["west"], -114.0)
        self.assertEqual(overlay["north"], 51.10)

    def test_empty_grid_has_no_overlay(self):
        empty = np.array([]).reshape(0, 0)
        self.assertIsNone(grid_overlay_png(empty, empty, empty))


class CredentialFileTests(unittest.TestCase):
    def test_save_fields_and_file_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "AWS_Keys.txt")
            status = tropomi.save_credentials(
                CredentialWrite(access_key_id="AKIAUI", secret_access_key="secret-ui"),
                path,
            )
            self.assertTrue(status.configured)
            self.assertEqual(status.access_key_id, "AKIAUI")
            tropomi.save_credentials(
                CredentialWrite(file_text="access_key_id: edited\nsecret_access_key: from-file\n"),
                path,
            )
            reloaded = tropomi.credential_status(path)
            self.assertEqual(reloaded.access_key_id, "edited")
            self.assertIn("from-file", reloaded.file_text)

    def test_missing_file_is_not_configured(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "AWS_Keys.txt")
            status = tropomi.credential_status(path)
            self.assertFalse(status.configured)
            self.assertIn("access_key_id", status.file_text)


class ScreenTests(unittest.TestCase):
    def test_collects_hotspots_and_last_overlay(self):
        lons, lats, field, wind, pressure = _grid()

        def write_results(*_args, **_kwargs):
            return {
                "figure_path": "assets/day.jpg",
                "hotspots": [{
                    "id": "20260601-1",
                    "lon": -114.0,
                    "lat": 51.1,
                    "enhancement_ppb": 100.0,
                    "emission_kg_hr": 2500.0,
                }],
                "overlay": {
                    "image_url": "data:image/png;base64,aaa",
                    "west": -114.0,
                    "south": 51.05,
                    "east": -113.95,
                    "north": 51.1,
                },
            }

        result = tropomi.screen_region(
            BBox(-114.2, 51.0, -113.8, 51.2),
            DateRange(date(2026, 6, 1), date(2026, 6, 1)),
            load_ch4=lambda *_args, **_kwargs: (lons, lats, field, wind, pressure),
            screen=lambda *_args, **_kwargs: ([field], [wind], [pressure], [lons], [lats]),
            write_results=write_results,
        )
        self.assertEqual(result.plume_days, 1)
        self.assertEqual(result.items[0].source, "tropomi")
        self.assertEqual(result.items[0].instrument, "TROPOMI")
        self.assertEqual(result.overlay.image_url, "data:image/png;base64,aaa")
        self.assertIn("Suspect plumes", result.message)

    def test_missing_files_are_source_errors(self):
        def load_ch4(*_args, **_kwargs):
            raise FileNotFoundError("No TROPOMI CH4 files")

        with self.assertRaises(SourceError):
            tropomi.screen_region(
                BBox(-114.2, 51.0, -113.8, 51.2),
                DateRange(date(2026, 6, 1), date(2026, 6, 1)),
                load_ch4=load_ch4,
            )


def _empty_screen():
    return TropomiScreenResult(
        source="tropomi",
        items=[],
        total_count=0,
        truncated=False,
        attribution=tropomi.ATTRIBUTION,
        license_url=tropomi.LICENSE_URL,
        message="No ultra-emitter event detected in your region.",
        screened_days=1,
        plume_days=0,
    )


class TropomiEndpointTests(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(main.app)

    def test_screen_requires_region_fields(self):
        response = self.client.post(
            "/api/tropomi/screen",
            json={"start": "2026-06-01", "end": "2026-06-01"},
        )
        self.assertEqual(response.status_code, 422)

    def test_screen_dispatches(self):
        with mock.patch.object(tropomi, "screen_region", return_value=_empty_screen()) as fake:
            response = self.client.post(
                "/api/tropomi/screen",
                json={
                    "west": -114.2, "south": 51.0, "east": -113.8, "north": 51.2,
                    "start": "2026-06-01", "end": "2026-06-01",
                },
            )
        self.assertEqual(response.status_code, 200)
        fake.assert_called_once()
        self.assertIn("No ultra-emitter", response.json()["message"])

    def test_put_credentials(self):
        with mock.patch.object(
            tropomi,
            "save_credentials",
            return_value=CredentialStatus(
                configured=True,
                access_key_id="AKIAUI",
                file_text="access_key_id: AKIAUI\n",
                example="",
            ),
        ) as save:
            response = self.client.put(
                "/api/tropomi/credentials",
                json={"accessKeyId": "AKIAUI", "secretAccessKey": "secret-ui"},
            )
        self.assertEqual(response.status_code, 200)
        save.assert_called_once()
        self.assertEqual(response.json()["accessKeyId"], "AKIAUI")
        self.assertTrue(response.json()["configured"])
