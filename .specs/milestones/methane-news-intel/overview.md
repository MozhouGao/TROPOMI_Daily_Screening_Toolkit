# Methane news intel v1

## Goal

Show recent methane-emissions news for the drawn map region, tagged with
intersecting basin names, so screening sits next to public reporting without
a paid news API.

## Scope

In:

- `GET /api/news?west&south&east&north` returning a `NewsCollection` (not
  `PlumeCollection`).
- Google News RSS proxy in the FastAPI backend (no API key). Query
  `methane (emissions OR leak)` plus each overlapping basin name, and one
  global `methane emissions` feed.
- Static basin list with axis-aligned bboxes (Permian, Eagle Ford,
  Haynesville, Marcellus, Utica, Bakken, DJ, Anadarko, Uinta, Permian Delaware
  folded into Permian, Montney, Duvernay, WCSB, Appalachian). Overlap with
  the drawn rectangle selects which basin queries to run.
- Tag an article if title or snippet contains a basin alias (case-insensitive).
  Tagged items sort first. Cap 30, `truncated` if more.
- Cache RSS XML under `NEWS_data/` for ~30 minutes via `file_cache` (same
  stale-on-failure pattern as MARS). Gitignore the cache dir.
- `NewsPanel` in the right sidebar under the plume list. Headline, source,
  date, basin tags, open-in-new-tab. Empty state until a region is drawn.
  i18n + tests.

Out:

- GDELT, NewsAPI.org, GNews, Guardian keys.
- Map markers for articles, user accounts, full EIA shapefiles.
- Filtering by the plume time window (RSS is “recent”; pubDate is display-only).
- MethaneSAT L3/L4 fetch (stays listing-only).

## Approach

- New schemas: `NewsArticle` (`id, title, url, source, publishedAt, snippet,
  basins, attribution`) and `NewsCollection` (`items, totalCount, truncated,
  attribution, licenseUrl`). CamelCase on the wire like plumes.
- Adapter `backend/sources/news.py`: pick basins via `bbox` overlap, fetch RSS
  with injectable `fetch`, parse XML, dedupe by URL, tag, cap.
- Basins live in `backend/sources/basins.py` as data, not a GIS stack. A
  rectangle that overlaps several basins runs one RSS query per basin plus
  the global query; results are merged.
- Frontend: React Query `useNews(bbox)`, no layer registry entry (not a plume
  source).
- Attribution: “Headlines via Google News RSS. Link through to the publisher.”
- Decisions: `../../decisions/0001-typescript-frontend-python-backend.md`.

## Acceptance criteria

1. `GET /api/news` without a bbox is 422. Inverted bbox is 422. Unknown path
   still 404.
2. With a Permian-ish bbox and a fixture RSS, items include a tagged
   “Permian” article; a western-Canada bbox tags Montney/WCSB when those
   strings appear; cap 30 sets `truncated`.
3. RSS fetch failure with a warm cache still 200; with no cache is 502.
4. Drawing a region loads the panel; clearing the region shows the empty
   state. Links use `target="_blank"` and `rel="noreferrer"`.
5. README names the endpoint, cache dir, and that this is headlines only.

## Risks and open questions

| Item | Status |
|------|--------|
| Google News RSS has no SLA and may throttle | Accepted: cache 30 min, 502 if cold miss |
| RSS items have no reliable lat/lon | Accepted: basin name overlap, not true geocode |
| Headlines are not detections | Copy states they are news, not inventory |
| Google ToS / robots | Link-out only; do not scrape article bodies |

## References

- Google News RSS search: `https://news.google.com/rss/search?q=...&hl=en-US&gl=US&ceid=US:en`
- Earlier probe: RSS search for `methane+emissions` returned a feed (~237 KB)
