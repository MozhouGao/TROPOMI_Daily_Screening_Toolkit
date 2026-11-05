"""Google News RSS headlines tagged with overlapping oil-and-gas basins."""
import hashlib
import os
import re
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from email.utils import parsedate_to_datetime
from typing import Callable, Dict, List, Optional
from urllib.parse import urlencode

from backend.schemas import BBox, NewsArticle, NewsCollection
from backend.sources.basins import BASINS, overlapping
from backend.sources.file_cache import ensure_file
from backend.sources.http import fetch_bytes

ATTRIBUTION = "Headlines via Google News RSS. Link through to the publisher."
LICENSE_URL = "https://news.google.com"
CACHE_DIRNAME = "NEWS_data"
CACHE_MAX_AGE = timedelta(minutes=30)
MAX_ITEMS = 30
GLOBAL_QUERY = "methane emissions"
BASIN_QUERY = "methane (emissions OR leak)"


def default_cache_dir() -> str:
    return os.path.join(os.getcwd(), CACHE_DIRNAME)


def rss_url(query: str) -> str:
    return "https://news.google.com/rss/search?" + urlencode(
        {"q": query, "hl": "en-US", "gl": "US", "ceid": "US:en"}
    )


def cache_path(cache_dir: str, query: str) -> str:
    digest = hashlib.sha256(query.encode()).hexdigest()[:16]
    return os.path.join(cache_dir, f"{digest}.xml")


def strip_html(value: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", value)).strip()


def parse_published(value: Optional[str]) -> Optional[datetime]:
    if not value:
        return None
    try:
        moment = parsedate_to_datetime(value)
    except (TypeError, ValueError, IndexError):
        return None
    return moment


def tag_basins(text: str) -> List[str]:
    names = []
    for basin in BASINS:
        if any(re.search(r"\b" + re.escape(alias) + r"\b", text, re.I) for alias in basin.aliases):
            names.append(basin.name)
    return names


def parse_items(body: bytes) -> List[NewsArticle]:
    root = ET.fromstring(body)
    articles = []
    for item in root.findall(".//item"):
        title = (item.findtext("title") or "").strip()
        url = (item.findtext("link") or "").strip()
        if title == "" or url == "":
            continue
        source = (item.findtext("source") or "").strip()
        snippet = strip_html(item.findtext("description") or "")
        guid = (item.findtext("guid") or url).strip()
        articles.append(
            NewsArticle(
                id=guid,
                title=title,
                url=url,
                source=source or "Google News",
                published_at=parse_published(item.findtext("pubDate")),
                snippet=snippet,
                basins=tag_basins(f"{title} {snippet}"),
                attribution=ATTRIBUTION,
            )
        )
    return articles


def queries_for(bbox: BBox) -> List[str]:
    return [GLOBAL_QUERY] + [f"{BASIN_QUERY} {basin.name}" for basin in overlapping(bbox)]


def fetch_news(
    bbox: BBox,
    cache_dir: Optional[str] = None,
    now: Optional[datetime] = None,
    fetch: Callable[[str], bytes] = fetch_bytes,
) -> NewsCollection:
    cache_dir = cache_dir or default_cache_dir()
    seen: Dict[str, NewsArticle] = {}
    for query in queries_for(bbox):
        path = ensure_file(rss_url(query), cache_path(cache_dir, query), CACHE_MAX_AGE, now, fetch)
        with open(path, "rb") as handle:
            for article in parse_items(handle.read()):
                existing = seen.get(article.url)
                if existing is None or len(article.basins) > len(existing.basins):
                    seen[article.url] = article
    items = sorted(
        seen.values(),
        key=lambda article: (
            len(article.basins) == 0,
            -(article.published_at.timestamp() if article.published_at else 0),
        ),
    )
    return NewsCollection(
        items=items[:MAX_ITEMS],
        total_count=len(items),
        truncated=len(items) > MAX_ITEMS,
        attribution=ATTRIBUTION,
        license_url=LICENSE_URL,
    )
