"""Local file cache for HTTP downloads. Re-runs are offline once the cache is warm."""
from __future__ import annotations
import hashlib, json, re
from pathlib import Path
from typing import Any, Optional
import time

import requests



RETRIES = 3
RETRY_WAIT_S = 2.0

class Cache:
    def __init__(self, root: Path, timeout: int = 300) -> None:
        self.root = Path(root)
        self.timeout = timeout
        self.root.mkdir(parents=True, exist_ok=True)

    def _path(self, url: str, suffix: str, method: str, data: Optional[str]) -> Path:
        key = hashlib.sha1(f"{method} {url} {data or ''}".encode()).hexdigest()[:16]
        host = re.sub(r"^https?://", "", url).split("/")[0].replace(".", "_")
        slug = re.sub(r"[^A-Za-z0-9._-]+", "_", re.sub(r"^https?://[^/]+/", "", url))[:80].strip("_")
        d = self.root / host
        d.mkdir(parents=True, exist_ok=True)
        return d / f"{slug}-{key}{suffix}"

    def fetch_bytes(self, url: str, suffix: str = ".bin", method: str = "GET", data: Optional[str] = None, force: bool = False) -> Path:
        p = self._path(url, suffix, method, data)
        if p.exists() and not force:
            return p
        headers = {"Content-Type": "application/x-www-form-urlencoded"} if data is not None else None
        # transient failures are retried: transport errors (connection reset, response cut short) and server errors (HTTP 5xx, e.g. the PDBe
        # API's intermittent 500s); client errors (4xx) are not. The last attempt's error is raised.
        for attempt in range(RETRIES):
            last = attempt == RETRIES - 1
            try:
                r = requests.request(method, url, data=data, headers=headers, timeout=self.timeout)
            except (requests.ConnectionError, requests.exceptions.ChunkedEncodingError, requests.Timeout):
                if last:
                    raise
                time.sleep(RETRY_WAIT_S * (attempt + 1))
                continue
            if getattr(r, "status_code", 200) >= 500 and not last:
                time.sleep(RETRY_WAIT_S * (attempt + 1))
                continue
            break
        r.raise_for_status()
        tmp = p.with_suffix(p.suffix + ".part")
        tmp.write_bytes(r.content)
        tmp.replace(p)
        return p

    def fetch_text(self, url: str, suffix: str = ".txt", method: str = "GET", data: Optional[str] = None, force: bool = False) -> str:
        return self.fetch_bytes(url, suffix, method, data, force).read_text()

    def fetch_json(self, url: str, method: str = "GET", data: Optional[str] = None, force: bool = False) -> Any:
        p = self.fetch_bytes(url, ".json", method, data, force)
        try:
            return json.loads(p.read_text())
        except json.JSONDecodeError:
            p.unlink(missing_ok=True)
            raise
