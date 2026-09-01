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


if __name__ == "__main__":
    unittest.main()
