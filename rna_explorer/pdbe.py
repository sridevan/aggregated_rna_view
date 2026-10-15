"""PDBe retrieval: Rfam nucleic mappings, molecules (batch), residue listings, search API metadata, updated mmCIF files."""
from __future__ import annotations
from pathlib import Path
from typing import Any, Iterable, Optional
from urllib.parse import urlencode
import requests
from .cache import Cache
from .config import PDBE_API, PDBE_FILES, PDBE_SEARCH

SEARCH_FIELDS = ["pdb_id", "entity_id", "chain_id", "struct_asym_id", "organism_scientific_name", "tax_id", "molecule_name",
                 "experimental_method", "resolution", "release_year", "title", "rfam_accession"]


def get_nucleic_mappings(cache: Cache, acc: str) -> dict[str, list[dict[str, Any]]]:
    d = cache.fetch_json(f"{PDBE_API}/v2/nucleic_mappings/{acc}")
    return d.get(acc, {}).get("PDB", {})


def get_entry_nucleic_mappings(cache: Cache, pdb_id: str) -> dict[str, Any]:
    return cache.fetch_json(f"{PDBE_API}/v2/nucleic_mappings/{pdb_id}").get(pdb_id, {})


def get_molecules(cache: Cache, pdb_ids: Iterable[str], chunk: int = 40) -> dict[str, list[dict[str, Any]]]:
    ids = sorted({p.lower() for p in pdb_ids})
    out: dict[str, list[dict[str, Any]]] = {}
    for i in range(0, len(ids), chunk):
        part = ids[i:i + chunk]
        try:
            out.update(cache.fetch_json(f"{PDBE_API}/pdb/entry/molecules/", method="POST", data=",".join(part)))
        except requests.HTTPError:   # batch endpoint refused; fall back to per-entry GET below
            pass
        for p in part:
            if p not in out:
                out.update(cache.fetch_json(f"{PDBE_API}/pdb/entry/molecules/{p}"))
    return out


def get_residue_listing(cache: Cache, pdb_id: str) -> dict[str, Any]:
    return cache.fetch_json(f"{PDBE_API}/pdb/entry/residue_listing/{pdb_id}").get(pdb_id, {})


def search_rfam_entities(cache: Cache, acc: str) -> dict[tuple[str, int], dict[str, Any]]:
    """Entity-level metadata (organism, method, resolution, ...) for every entity annotated with the Rfam accession."""
    q = urlencode({"q": f"rfam_accession:{acc}", "fl": ",".join(SEARCH_FIELDS), "rows": 10000, "wt": "json"})
    docs = cache.fetch_json(f"{PDBE_SEARCH}?{q}")["response"]["docs"]
    return {(d["pdb_id"].lower(), int(d["entity_id"])): d for d in docs}


def get_updated_cif(cache: Cache, pdb_id: str) -> Path:
    return cache.fetch_bytes(f"{PDBE_FILES}/{pdb_id.lower()}_updated.cif", suffix=".cif")


PDBE_STATIC = "https://www.ebi.ac.uk/pdbe/static/entry"


def pdbe_2d_available(cache: Cache, pdb_id: str, entity_id: int | str, auth_asym_id: str) -> Optional[bool]:
    """True when PDBe serves both the RNA 2D topology and the FR3D base-pair file for the chain (inputs of the PDBe RNA viewer); None if unknown."""
    urls = [f"{PDBE_STATIC}/{pdb_id.lower()}_{entity_id}_{auth_asym_id.upper()}.json", f"{PDBE_STATIC}/{pdb_id.lower()}_{auth_asym_id.upper()}_basepair.json"]
    try:
        for u in urls:
            p = cache._path(u, ".head", "HEAD", None)
            if p.exists():
                status = p.read_text().strip()
            else:
                status = str(requests.head(u, timeout=cache.timeout, allow_redirects=True).status_code)
                p.write_text(status)
            if status != "200":
                return False
        return True
    except requests.RequestException:
        return None
