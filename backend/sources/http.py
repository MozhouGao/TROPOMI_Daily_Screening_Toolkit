from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from backend.schemas import SourceError

USER_AGENT = "TROPOMI-Daily-Screening-Toolkit"
DEFAULT_TIMEOUT_SECONDS = 120


def fetch_bytes(url: str, timeout: float = DEFAULT_TIMEOUT_SECONDS) -> bytes:
    request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
    try:
        with urlopen(request, timeout=timeout) as response:
            return response.read()
    except HTTPError as exc:
        raise SourceError(f"{url} returned HTTP {exc.code}") from exc
    except URLError as exc:
        raise SourceError(f"{url} unreachable: {exc.reason}") from exc
