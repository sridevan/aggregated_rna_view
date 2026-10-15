"""Rfam retrieval: family metadata, covariance model, release, PDB structure hits (with bit scores), clan membership."""
from __future__ import annotations
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Optional
from .cache import Cache
from .config import RFAM_BASE, RFAM_FTP


@dataclass
class CMInfo:
    acc: str
    name: str
    description: str
    clen: int
    checksum: str
    build_date: str
    infernal_format: str
    path: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def get_family(cache: Cache, acc: str) -> dict[str, Any]:
    d = cache.fetch_json(f"{RFAM_BASE}/family/{acc}?content-type=application/json")
    return d.get("rfam", d)


def get_cm(cache: Cache, acc: str) -> CMInfo:
    p = cache.fetch_bytes(f"{RFAM_BASE}/family/{acc}/cm", suffix=".cm")
    head: dict[str, str] = {}
    fmt = ""
    with open(p) as f:
        for line in f:
            if line.startswith("INFERNAL"):
                fmt = line.strip()
            m = re.match(r"^(NAME|ACC|DESC|CLEN|CKSUM|DATE)\s+(.+?)\s*$", line)
            if m and m.group(1) not in head:
                head[m.group(1)] = m.group(2)
            if line.startswith("CM"):
                break
    if head.get("ACC") != acc:
        raise ValueError(f"CM accession mismatch: {head.get('ACC')} != {acc}")
    return CMInfo(acc, head["NAME"], head.get("DESC", ""), int(head["CLEN"]), head["CKSUM"], head.get("DATE", ""), fmt, str(p))


def get_release(cache: Cache) -> str:
    txt = cache.fetch_text(f"{RFAM_FTP}/README", suffix=".README.txt")
    m = re.search(r"Release\s+([\d.]+)", txt)
    return m.group(1) if m else "unknown"


def get_structures(cache: Cache, acc: str) -> list[dict[str, Any]]:
    """Rfam 'structures' endpoint: one record per PDB chain hit with bit_score, evalue_score, cm_start/cm_end, pdb_start/pdb_end."""
    d = cache.fetch_json(f"{RFAM_BASE}/family/{acc}/structures?content-type=application/json")
    return list(d.get("mapping", []))


def get_clan_members(cache: Cache, acc: str) -> list[str]:
    """Accessions of the families in the same clan as `acc` (empty if the family is not in a clan)."""
    fam = get_family(cache, acc)
    clan: Optional[str] = None
    for key in ("clan", "clan_acc"):
        v = fam.get(key)
        if isinstance(v, dict):
            v = v.get("clan_acc") or v.get("acc")
        if isinstance(v, str) and v.startswith("CL"):
            clan = v
    if not clan:
        return []
    d = cache.fetch_json(f"{RFAM_BASE}/clan/{clan}?content-type=application/json")
    c = d.get("clan", d)
    members = [m.get("acc") if isinstance(m, dict) else m for m in c.get("members", [])]
    return [m for m in members if isinstance(m, str)]


def structures_by_chain(records: list[dict[str, Any]]) -> dict[tuple[str, str], dict[str, Any]]:
    return {(r["pdb_id"].lower(), r["chain"]): r for r in records}
