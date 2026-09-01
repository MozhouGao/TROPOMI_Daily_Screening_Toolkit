import os
import unittest
from datetime import date

import numpy as np

from TROPOMI_toolkit import (
    bbox_from_coordinates,
    is_ch4_l2_nc,
    iter_screening_dates,
    local_download_path,
    nanargmax,
    screening_plumes,
    tropomi_sensing_start_date,
    unique_hotspots,
)


class DateWalkTests(unittest.TestCase):
    def test_same_day_is_included(self):
        days = list(iter_screening_dates("2026-08-30", "2026-08-30"))
        self.assertEqual(days, [date(2026, 8, 30)])

    def test_inclusive_end_date(self):
        days = list(iter_screening_dates("2026-08-01", "2026-08-05"))
        self.assertEqual(days[0], date(2026, 8, 1))
        self.assertEqual(days[-1], date(2026, 8, 5))
        self.assertEqual(len(days), 5)

    def test_cross_month(self):
        days = list(iter_screening_dates(date(2026, 8, 28), date(2026, 9, 2)))
        self.assertEqual(len(days), 6)
        self.assertEqual(days[-1], date(2026, 9, 2))

    def test_swaps_inverted_range(self):
        days = list(iter_screening_dates("2026-08-10", "2026-08-08"))
        self.assertEqual(days, [date(2026, 8, 8), date(2026, 8, 9), date(2026, 8, 10)])

    def test_missing_end_is_single_day(self):
        days = list(iter_screening_dates("2026-08-30"))
        self.assertEqual(days, [date(2026, 8, 30)])


class BBoxTests(unittest.TestCase):
    def test_min_max_independent_of_ring_order(self):
        ring = [
            [10.0, 20.0],
            [12.0, 20.0],
            [12.0, 23.0],
            [10.0, 23.0],
            [10.0, 20.0],
        ]
        reversed_ring = list(reversed(ring))
        self.assertEqual(bbox_from_coordinates(ring), (20.0, 23.0, 10.0, 12.0))
        self.assertEqual(bbox_from_coordinates(reversed_ring), (20.0, 23.0, 10.0, 12.0))


class FilenameTests(unittest.TestCase):
    def test_sensing_start_not_end_time(self):
        name = "S5P_OFFL_L2__CH4____20260828T235000_20260829T013000_00000_03_020901_20260830T000000.nc"
        self.assertEqual(tropomi_sensing_start_date(name), "20260828")

    def test_product_type_filter(self):
        offl = "S5P_OFFL_L2__CH4____20260829T011418_20260829T025548_45991_03_020901_20260830T171040.nc"
        nrti = "S5P_NRTI_L2__CH4____20260830T184616_20260830T185116_46015_03_020901_20260830T191951.nc"
        other = "S5P_OFFL_L2__NO2____20260829T011418_20260829T025548_45991_03_020901_20260830T171040.nc"
        self.assertTrue(is_ch4_l2_nc(offl, ("OFFL",)))
        self.assertFalse(is_ch4_l2_nc(nrti, ("OFFL",)))
        self.assertTrue(is_ch4_l2_nc(nrti, ("NRTI",)))
        self.assertFalse(is_ch4_l2_nc(other))


class PathTests(unittest.TestCase):
    def test_download_path_uses_os_join(self):
        target = os.path.join(os.sep, "tmp", "TROPOMI_data")
        path = local_download_path(target, "S5P_OFFL_L2__CH4____x.nc")
        self.assertEqual(path, os.path.join(target, "S5P_OFFL_L2__CH4____x.nc"))
        self.assertFalse("\\" in path.replace(os.sep, "/").replace(target.replace(os.sep, "/"), ""))


class ScreeningTests(unittest.TestCase):
    def test_flat_patch_does_not_divide_by_zero(self):
        field = np.full((20, 20), 1800.0)
        wind = np.ones_like(field)
        pressure = np.full_like(field, 1000.0)
        lons = np.linspace(-114, -113, 20)[None, :].repeat(20, axis=0)
        lats = np.linspace(51, 52, 20)[:, None].repeat(20, axis=1)
        result = screening_plumes(field, wind, pressure, lons, lats, 15, 1)
        self.assertEqual(result[0], [])

    def test_empty_field(self):
        empty = np.array([]).reshape(0, 0)
        result = screening_plumes(empty, empty, empty, empty, empty, 15, 1)
        self.assertEqual(result[0], [])

    def test_detects_synthetic_plume(self):
        field = np.full((20, 20), 1800.0)
        field[8:13, 8:13] = 1900.0
        wind = np.ones_like(field)
        pressure = np.full_like(field, 1000.0)
        lons = np.linspace(-114, -113, 20)[None, :].repeat(20, axis=0)
        lats = np.linspace(51, 52, 20)[:, None].repeat(20, axis=1)
        detected, *_ = screening_plumes(field, wind, pressure, lons, lats, 15, 1)
        self.assertGreater(len(detected), 0)
        self.assertGreater(np.nanmax(detected[0]), 15)

    def test_hotspot_dedupe_uses_coordinate_pairs(self):
        plume = np.zeros((3, 3))
        plume[1, 1] = 40
        lons_a = np.array([[10.0, 10.05, 10.1], [10.0, 10.05, 10.1], [10.0, 10.05, 10.1]])
        lats_a = np.array([[20.1, 20.1, 20.1], [20.05, 20.05, 20.05], [20.0, 20.0, 20.0]])
        lons_b = np.array([[10.0, 10.05, 10.1], [10.0, 10.05, 10.1], [10.0, 10.05, 10.1]])
        lats_b = np.array([[21.1, 21.1, 21.1], [21.05, 21.05, 21.05], [21.0, 21.0, 21.0]])
        wind = np.ones((3, 3))
        pressure = np.full((3, 3), 1000.0)
        polygons, enhance, max_lons, max_lats, *_ = unique_hotspots(
            [plume, plume],
            [wind, wind],
            [pressure, pressure],
            [lons_a, lons_b],
            [lats_a, lats_b],
        )
        self.assertEqual(len(polygons), 2)
        self.assertEqual(len(enhance), 2)
        self.assertNotEqual((max_lons[0], max_lats[0]), (max_lons[1], max_lats[1]))

    def test_nanargmax_ignores_nan(self):
        arr = np.array([[np.nan, 1.0], [3.0, np.nan]])
        self.assertEqual(nanargmax(arr), (1, 0))


if __name__ == "__main__":
    unittest.main()
