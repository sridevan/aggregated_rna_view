"""NCBI Taxonomy lookups (E-utilities) to resolve strain/subspecies tax ids to species-level names. Cached; <= 3 requests/s without an API key."""
from __future__ import annotations
import re, time
from typing import Any, Optional
from .cache import Cache

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_last = [0.0]


def _throttle() -> None:
    now = time.time()
    wait = 0.35 - (now - _last[0])
    if wait > 0:
        time.sleep(wait)
    _last[0] = time.time()


def taxon(cache: Cache, tax_id: int | str) -> Optional[dict[str, Any]]:
    """Return {'tax_id', 'scientific_name', 'rank', 'species_name', 'species_tax_id'} or None if the lookup fails.
    species_name is the taxon itself when its rank is 'species', otherwise the species-rank ancestor from LineageEx (None if there is none)."""
    try:
        _throttle()
        xml = cache.fetch_text(f"{EUTILS}/efetch.fcgi?db=taxonomy&id={tax_id}&retmode=xml", suffix=".taxonomy.xml")
    except Exception:   # network or HTTP error: caller falls back to the source name
        return None
    m_name = re.search(r"<ScientificName>(.*?)</ScientificName>", xml)
    m_rank = re.search(r"<Rank>(.*?)</Rank>", xml)
    if not m_name:
        return None
    name, rank = m_name.group(1), (m_rank.group(1) if m_rank else "")
    species_name, species_tax = (name, str(tax_id)) if rank == "species" else (None, None)
    if species_name is None:
        lin = re.search(r"<LineageEx>(.*?)</LineageEx>", xml, re.S)
        if lin:
            for m in re.finditer(r"<Taxon>\s*<TaxId>(\d+)</TaxId>\s*<ScientificName>(.*?)</ScientificName>\s*<Rank>(.*?)</Rank>", lin.group(1)):
                if m.group(3) == "species":
                    species_tax, species_name = m.group(1), m.group(2)
    return {"tax_id": str(tax_id), "scientific_name": name, "rank": rank, "species_name": species_name, "species_tax_id": species_tax}
