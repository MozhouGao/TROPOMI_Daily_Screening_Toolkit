import csv
import io
import os
import tempfile
import unittest
import zipfile
from datetime import date, datetime, timedelta, timezone

from backend.schemas import BBox, DateRange, SourceError
from backend.sources import mars

COLUMNS = [
    "id_plume", "source_name", "satellite", "tile_date", "lat", "lon", "actionable",
    "notified", "country", "sector", "detection_institution", "quantification_institution",
    "tile", "ch4_fluxrate", "ch4_fluxrate_std", "wind_u", "wind_v", "total_emission",
    "total_emission_std", "wind_speed", "last_update", "insert_date", "tile_background",
]


def _row(**overrides):
    row = {column: "" for column in COLUMNS}
    row.update({
        "id_plume": "bd5f3817-ded2-4cc1-9a5d-a0d02968865e",
        "source_name": "ARG_S_002",
        "satellite": "PRISMA - ASI",
        "tile_date": "2020-10-29T14:34:04",
        "lat": "-36.36529",
        "lon": "-67.07484",
        "actionable": "Yes",
        "notified": "True",
        "country": "Argentina",
        "sector": "Oil and Gas",
        "detection_institution": "UNEP IMEO MARS",
        "tile": "S2B_MSIL1C_20221102T102059",
        "ch4_fluxrate": "3148.0",
        "ch4_fluxrate_std": "1244.0",
    })
    row.update(overrides)
    return row


def _zip_bytes(rows):
    text = io.StringIO()
    writer = csv.DictWriter(text, fieldnames=COLUMNS)
    writer.writeheader()
    writer.writerows(rows)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr(mars.PLUMES_CSV_NAME, "\ufeff" + text.getvalue())
        archive.writestr("LICENSE-CC-BY-NC-SA.txt", "CC BY-NC-SA 4.0")
    return buffer.getvalue()


ROWS = [
    _row(),
    _row(id_plume="inside-2024", tile_date="2024-03-05T10:00:00", lat="-36.40", lon="-67.10",
         satellite="EMIT - NASA", ch4_fluxrate="900", ch4_fluxrate_std=""),
    _row(id_plume="outside-bbox", tile_date="2024-03-05T10:00:00", lat="51.0", lon="-114.0"),
]
BBOX = BBox(west=-68.0, south=-37.0, east=-66.0, north=-36.0)


class CacheTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cache_dir = os.path.join(self.tmp.name, "MARS_data")
        self.calls = []

    def tearDown(self):
        self.tmp.cleanup()

    def _fetch(self, url):
        self.calls.append(url)
        return _zip_bytes(ROWS)

    def test_downloads_when_missing(self):
        path = mars.ensure_export(self.cache_dir, fetch=self._fetch)
        self.assertTrue(os.path.isfile(path))
        self.assertEqual(self.calls, [mars.EXPORT_URL])

    def test_reuses_fresh_cache(self):
        mars.ensure_export(self.cache_dir, fetch=self._fetch)
        mars.ensure_export(self.cache_dir, fetch=self._fetch)
        self.assertEqual(len(self.calls), 1)

    def test_refreshes_stale_cache(self):
        path = mars.ensure_export(self.cache_dir, fetch=self._fetch)
        stale = datetime.now(timezone.utc) + mars.CACHE_MAX_AGE + timedelta(hours=1)
        mars.ensure_export(self.cache_dir, now=stale, fetch=self._fetch)
        self.assertEqual(len(self.calls), 2)
        self.assertTrue(os.path.isfile(path))

    def test_keeps_stale_cache_when_download_fails(self):
        path = mars.ensure_export(self.cache_dir, fetch=self._fetch)

        def failing_fetch(url):
            raise SourceError("offline")

        stale = datetime.now(timezone.utc) + mars.CACHE_MAX_AGE + timedelta(hours=1)
        self.assertEqual(mars.ensure_export(self.cache_dir, now=stale, fetch=failing_fetch), path)

    def test_download_failure_without_cache_raises(self):
        def failing_fetch(url):
            raise SourceError("offline")

        with self.assertRaises(SourceError):
            mars.ensure_export(self.cache_dir, fetch=failing_fetch)


class ParseTests(unittest.TestCase):
    def test_row_maps_to_shared_shape(self):
        plume = mars.parse_row(_row())
        self.assertEqual(plume.source, "mars")
        self.assertEqual(plume.observed_at, datetime(2020, 10, 29, 14, 34, 4, tzinfo=timezone.utc))
        self.assertAlmostEqual(plume.lon, -67.07484)
        self.assertAlmostEqual(plume.lat, -36.36529)
        self.assertEqual(plume.emission_kg_hr, 3148.0)
        self.assertEqual(plume.uncertainty_kg_hr, 1244.0)
        self.assertEqual(plume.instrument, "PRISMA - ASI")
        self.assertEqual(plume.sector, "Oil and Gas")
        self.assertIsNone(plume.image_url)
        self.assertEqual(plume.attribution, mars.ATTRIBUTION)
        self.assertEqual(plume.details["Satellite product"], "S2B_MSIL1C_20221102T102059")
        self.assertEqual(plume.details["Country"], "Argentina")
        self.assertEqual(plume.details["Notified"], "True")

    def test_blank_numbers_are_none(self):
        plume = mars.parse_row(_row(ch4_fluxrate="", ch4_fluxrate_std=""))
        self.assertIsNone(plume.emission_kg_hr)
        self.assertIsNone(plume.uncertainty_kg_hr)

    def test_read_plumes_from_zip(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "export.zip")
            with open(path, "wb") as handle:
                handle.write(_zip_bytes(ROWS))
            self.assertEqual(len(mars.read_plumes(path)), 3)

    def test_bad_zip_is_a_source_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "export.zip")
            with open(path, "wb") as handle:
                handle.write(b"not a zip")
            with self.assertRaises(SourceError):
                mars.read_plumes(path)


class FilterTests(unittest.TestCase):
    def test_filters_by_bbox_and_time(self):
        plumes = [mars.parse_row(row) for row in ROWS]
        filtered = mars.filter_plumes(plumes, BBOX, DateRange(date(2024, 1, 1), date(2024, 12, 31)))
        self.assertEqual([plume.id for plume in filtered], ["inside-2024"])

    def test_fetch_sorts_newest_first(self):
        with tempfile.TemporaryDirectory() as tmp:
            collection = mars.fetch_plumes(
                BBOX,
                DateRange(date(2020, 1, 1), date(2024, 12, 31)),
                cache_dir=os.path.join(tmp, "MARS_data"),
                fetch=lambda url: _zip_bytes(ROWS),
            )
        self.assertEqual([plume.id for plume in collection.items], ["inside-2024", ROWS[0]["id_plume"]])
        self.assertEqual(collection.total_count, 2)
        self.assertFalse(collection.truncated)
        self.assertEqual(collection.license_url, mars.LICENSE_URL)


if __name__ == "__main__":
    unittest.main()
