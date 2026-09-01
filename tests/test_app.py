import os
import tempfile
import unittest

from dash.exceptions import PreventUpdate

import app as tropomi_app


class AppCallbackTests(unittest.TestCase):
    def test_index_renders(self):
        client = tropomi_app.app.server.test_client()
        response = client.get("/")
        self.assertEqual(response.status_code, 200)
        layout = client.get("/_dash-layout")
        self.assertEqual(layout.status_code, 200)
        self.assertIn(b"TROPOMI Methane Plume Daily Screening Toolkit - V1.2", layout.data)

    def test_date_label_single_day(self):
        label = tropomi_app.update_date_range("2026-08-30", None)
        self.assertIn("August 30, 2026", label)
        self.assertIn("single day", label)

    def test_bbox_from_geojson(self):
        geojson = {
            "features": [{
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [10.0, 20.0],
                        [12.0, 20.0],
                        [12.0, 23.0],
                        [10.0, 23.0],
                        [10.0, 20.0],
                    ]],
                }
            }]
        }
        label = tropomi_app.update_polygon_data(geojson)
        self.assertIn("lon [10.00, 12.00]", label)
        self.assertIn("lat [20.00, 23.00]", label)

    def test_download_requires_date(self):
        self.assertEqual(tropomi_app.download_data(1, None, None), "Please select a start date.")

    def test_screening_requires_region(self):
        self.assertEqual(
            tropomi_app.run_analysis(1, "2026-08-30", "2026-08-30", None, None, None),
            "Please define the region.",
        )

    def test_prevent_update_on_empty_click(self):
        with self.assertRaises(PreventUpdate):
            tropomi_app.download_data(0, "2026-08-30", "2026-08-30")

    def test_layout_includes_credential_manager(self):
        client = tropomi_app.app.server.test_client()
        layout = client.get("/_dash-layout")
        self.assertIn(b"Copernicus credential manager", layout.data)
        self.assertIn(b"cred-open", layout.data)

    def test_open_manager_without_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "AWS_Keys.txt")
            style, access, secret, raw, status, banner = tropomi_app.open_credential_manager(path)
            self.assertEqual(style["display"], "block")
            self.assertEqual(access, "")
            self.assertEqual(secret, "")
            self.assertIn("access_key_id", raw)
            self.assertIn("No AWS_Keys.txt yet", status)
            self.assertEqual(banner, "Credentials not set")

    def test_save_and_reload_from_manager(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "AWS_Keys.txt")
            tropomi_app.save_credential_fields("AKIAUI", "secret-ui", path)
            style, access, secret, raw, status, banner = tropomi_app.open_credential_manager(path)
            self.assertEqual(access, "AKIAUI")
            self.assertEqual(secret, "secret-ui")
            self.assertIn("Existing AWS_Keys.txt loaded", status)
            self.assertEqual(banner, "Credentials saved locally")
            tropomi_app.save_credential_file_text(
                "access_key_id: edited\nsecret_access_key: from-file\n",
                path,
            )
            _, access, secret, raw, _, _ = tropomi_app.open_credential_manager(path)
            self.assertEqual(access, "edited")
            self.assertEqual(secret, "from-file")
            self.assertTrue(os.path.isfile(path))


if __name__ == "__main__":
    unittest.main()

