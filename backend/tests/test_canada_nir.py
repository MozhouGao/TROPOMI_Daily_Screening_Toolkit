import tempfile
import unittest
from datetime import date

import numpy as np

from backend.schemas import BBox, DateRange
from backend.sources import canada_nir
from backend.sources.file_cache import flux_kg_m2_s_to_kg_hr

BBOX = BBox(west=-80.2, south=43.6, east=-79.8, north=44.0)


class CanadaNirTests(unittest.TestCase):
    def test_cells_from_grid_converts_flux_and_filters(self):
        lons = np.array([-80.15, -80.05, -79.0])
        lats = np.array([43.65, 43.85])
        flux = np.array([[1e-8, 0.0, 1.0], [0.0, 2e-8, 0.0]])
        area = np.ones_like(flux) * 1.0e8
        plumes = canada_nir.cells_from_grid(
            lons, lats, flux, area, sector="solid waste", bbox=BBOX, min_kg_hr=1.0
        )
        self.assertEqual(len(plumes), 2)
        rates = sorted(plume.emission_kg_hr for plume in plumes)
        self.assertAlmostEqual(rates[0], flux_kg_m2_s_to_kg_hr(1e-8, 1.0e8))
        self.assertAlmostEqual(rates[1], flux_kg_m2_s_to_kg_hr(2e-8, 1.0e8))
        self.assertTrue(all(plume.source == "canada_nir" for plume in plumes))
        self.assertTrue(all(plume.observed_at.year == 2018 for plume in plumes))
        self.assertNotIn(-79.0, [plume.lon for plume in plumes])

    def test_fetch_reads_injected_grids(self):
        lons = np.array([-80.05])
        lats = np.array([43.85])
        flux = np.array([[5e-8]])
        area = np.array([[1.0e8]])

        def loader(path):
            return lons, lats, flux, area

        with tempfile.TemporaryDirectory() as tmp:
            collection = canada_nir.fetch_plumes(
                BBOX,
                DateRange(date(2026, 1, 1), date(2026, 1, 31)),
                cache_dir=tmp,
                fetch=lambda url: b"nc",
                load_grid=loader,
            )
        self.assertEqual(len(collection.items), 4)
        self.assertEqual(collection.source, "canada_nir")
        self.assertEqual({plume.sector for plume in collection.items}, set(canada_nir.SECTORS))


if __name__ == "__main__":
    unittest.main()
