import tempfile
import unittest

from backend.schemas import BBox, SourceError
from backend.sources import news

PERMIAN_BOX = BBox(west=-102.8, south=31.6, east=-101.6, north=32.4)
WESTERN_CANADA = BBox(west=-127.96, south=42.31, east=-100.18, north=58.39)

PERMIAN_RSS = """<?xml version="1.0" encoding="UTF-8"?>
<rss version="2.0"><channel>
  <item>
    <title>Permian methane leak found - Reuters</title>
    <link>https://example.org/permian-leak</link>
    <pubDate>Wed, 04 Nov 2026 12:00:00 GMT</pubDate>
    <description>Operators in the Permian Basin reported a methane leak.</description>
    <source url="https://reuters.com">Reuters</source>
    <guid>permian-1</guid>
  </item>
  <item>
    <title>Global methane pledge update - AP</title>
    <link>https://example.org/global</link>
    <pubDate>Tue, 03 Nov 2026 08:00:00 GMT</pubDate>
    <description>Countries discuss methane emissions.</description>
    <source url="https://apnews.com">AP</source>
    <guid>global-1</guid>
  </item>
</channel></rss>
"""

CANADA_RSS = """<?xml version="1.0" encoding="UTF-8"?>
<rss version="2.0"><channel>
  <item>
    <title>Montney producers cut methane - CBC</title>
    <link>https://example.org/montney</link>
    <pubDate>Mon, 02 Nov 2026 15:00:00 GMT</pubDate>
    <description>WCSB operators in the Montney play reported lower leak rates.</description>
    <source url="https://cbc.ca">CBC</source>
    <guid>montney-1</guid>
  </item>
</channel></rss>
"""


def _item(index: int) -> str:
    return (
        "<item>"
        f"<title>Methane story {index} - Wire</title>"
        f"<link>https://example.org/story-{index}</link>"
        "<pubDate>Wed, 04 Nov 2026 12:00:00 GMT</pubDate>"
        "<description>methane emissions</description>"
        "<source>Wire</source>"
        f"<guid>story-{index}</guid>"
        "</item>"
    )


MANY_RSS = (
    '<?xml version="1.0" encoding="UTF-8"?><rss version="2.0"><channel>'
    + "".join(_item(index) for index in range(31))
    + "</channel></rss>"
)


class NewsAdapterTests(unittest.TestCase):
    def test_permian_fixture_tags_permian_and_sorts_tagged_first(self):
        with tempfile.TemporaryDirectory() as tmp:
            collection = news.fetch_news(
                PERMIAN_BOX, cache_dir=tmp, fetch=lambda url: PERMIAN_RSS.encode()
            )
        self.assertEqual([item.url for item in collection.items], [
            "https://example.org/permian-leak",
            "https://example.org/global",
        ])
        self.assertEqual(collection.items[0].basins, ["Permian"])
        self.assertEqual(collection.items[1].basins, [])
        self.assertEqual(collection.items[0].source, "Reuters")
        self.assertIn("Google News RSS", collection.attribution)

    def test_western_canada_fixture_tags_montney_and_wcsb(self):
        with tempfile.TemporaryDirectory() as tmp:
            collection = news.fetch_news(
                WESTERN_CANADA, cache_dir=tmp, fetch=lambda url: CANADA_RSS.encode()
            )
        self.assertEqual(set(collection.items[0].basins), {"Montney", "WCSB"})

    def test_cap_30_sets_truncated(self):
        with tempfile.TemporaryDirectory() as tmp:
            collection = news.fetch_news(
                PERMIAN_BOX, cache_dir=tmp, fetch=lambda url: MANY_RSS.encode()
            )
        self.assertEqual(len(collection.items), 30)
        self.assertEqual(collection.total_count, 31)
        self.assertTrue(collection.truncated)

    def test_dedupes_by_url_across_feeds(self):
        with tempfile.TemporaryDirectory() as tmp:
            collection = news.fetch_news(
                PERMIAN_BOX, cache_dir=tmp, fetch=lambda url: PERMIAN_RSS.encode()
            )
        urls = [item.url for item in collection.items]
        self.assertEqual(len(urls), len(set(urls)))

    def test_stale_cache_is_kept_when_fetch_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            news.fetch_news(PERMIAN_BOX, cache_dir=tmp, fetch=lambda url: PERMIAN_RSS.encode())

            def fail(_url: str) -> bytes:
                raise SourceError("Google News returned HTTP 503")

            collection = news.fetch_news(PERMIAN_BOX, cache_dir=tmp, fetch=fail)
        self.assertEqual(collection.items[0].url, "https://example.org/permian-leak")

    def test_cold_fetch_failure_raises_source_error(self):
        def fail(_url: str) -> bytes:
            raise SourceError("Google News returned HTTP 503")

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SourceError):
                news.fetch_news(PERMIAN_BOX, cache_dir=tmp, fetch=fail)

    def test_queries_global_and_overlapping_basins(self):
        urls = []

        def capture(url: str) -> bytes:
            urls.append(url)
            return PERMIAN_RSS.encode()

        with tempfile.TemporaryDirectory() as tmp:
            news.fetch_news(PERMIAN_BOX, cache_dir=tmp, fetch=capture)
        joined = " ".join(urls)
        self.assertIn("news.google.com/rss/search", joined)
        self.assertTrue(any("Permian" in url for url in urls))
        self.assertTrue(any("methane" in url.lower() and "emissions" in url.lower() for url in urls))
