import os
import tempfile
import unittest
from datetime import datetime, timedelta, timezone

from backend.schemas import SourceError
from backend.sources import file_cache


class FileCacheTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "export.csv")
        self.calls = []

    def tearDown(self):
        self.tmp.cleanup()

    def _fetch(self, url):
        self.calls.append(url)
        return b"id,lat,lon\n1,0,0\n"

    def test_downloads_when_missing(self):
        path = file_cache.ensure_file(
            "https://example.org/export.csv",
            self.path,
            timedelta(days=7),
            fetch=self._fetch,
        )
        self.assertEqual(path, self.path)
        self.assertTrue(os.path.isfile(path))
        self.assertEqual(self.calls, ["https://example.org/export.csv"])

    def test_reuses_fresh_cache(self):
        file_cache.ensure_file(
            "https://example.org/export.csv",
            self.path,
            timedelta(days=7),
            fetch=self._fetch,
        )
        file_cache.ensure_file(
            "https://example.org/export.csv",
            self.path,
            timedelta(days=7),
            fetch=self._fetch,
        )
        self.assertEqual(len(self.calls), 1)

    def test_keeps_stale_cache_when_download_fails(self):
        path = file_cache.ensure_file(
            "https://example.org/export.csv",
            self.path,
            timedelta(days=7),
            fetch=self._fetch,
        )

        def failing_fetch(url):
            raise SourceError("offline")

        stale = datetime.now(timezone.utc) + timedelta(days=8)
        self.assertEqual(
            file_cache.ensure_file(
                "https://example.org/export.csv",
                self.path,
                timedelta(days=7),
                now=stale,
                fetch=failing_fetch,
            ),
            path,
        )

    def test_download_failure_without_cache_raises(self):
        def failing_fetch(url):
            raise SourceError("offline")

        with self.assertRaises(SourceError):
            file_cache.ensure_file(
                "https://example.org/export.csv",
                self.path,
                timedelta(days=7),
                fetch=failing_fetch,
            )


class RateConversionTests(unittest.TestCase):
    def test_tonnes_per_year_to_kg_hr(self):
        self.assertAlmostEqual(file_cache.tonnes_per_year_to_kg_hr(8.766), 1.0, places=5)

    def test_gg_per_year_to_kg_hr(self):
        self.assertAlmostEqual(file_cache.gg_per_year_to_kg_hr(0.008766), 1.0, places=5)

    def test_co2e_tonnes_to_ch4_kg_hr(self):
        self.assertAlmostEqual(
            file_cache.co2e_tonnes_to_ch4_kg_hr(25 * 8.766, gwp=25), 1.0, places=5
        )


if __name__ == "__main__":
    unittest.main()
