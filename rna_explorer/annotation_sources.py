"""Retrieval and parsing of external annotation sources (Milestone 3B).

Sources
  * BGSU nucleotide annotations: family reference file ``data/nucleotide_annotations/<ACC>_annotations.txt`` and the per-chain
    transfer endpoint ``correspondence/nucleotide_annotation?chain=PDB|model|chain``.
  * BGSU RNA 3D Motif Atlas: loop membership ``rna3dhub/loops/download/<PDB>`` (CSV) and the per-entry motifs page
    ``rna3dhub/pdb/<PDB>/motifs`` (HTML; the only place exposing loop -> motif-group assignments and the 'problematic' status).
    The current Motif Atlas release id is read from ``rna3dhub/motifs/release/<type>/current``.
  * PDBe bound molecules: ``api/v2/pdb/entry/binding_sites/<pdb>/<entity>``.

Every response is cached through :class:`rna_explorer.cache.Cache`. A ``raw_dir`` (e.g. ``annotation_assessment/cache``) can seed
the cache with files downloaded during Milestone 3A so that builds are reproducible offline. Each parsed result carries a
``SourceMeta`` (url, retrieval date = cache file mtime, status) so that provenance can be written to the outputs.
"""
from __future__ import annotations
import csv
import datetime as _dt
import html as _html
import io
import json
import re
import shutil
from dataclasses import replace, dataclass, asdict
from pathlib import Path
from typing import Any, Optional
import requests
from .cache import Cache

BGSU = "https://rna.bgsu.edu"
PDBE_API_V2 = "https://www.ebi.ac.uk/pdbe/api/v2"
LOOP_TYPES = ("HL", "IL", "J3", "J4")


@dataclass
class SourceMeta:
    url: str
    status: str            # ok | empty | http_error | parse_error | not_available (the source states it has no data)
    retrieved: Optional[str]
    detail: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class AnnotationFetcher:
    """Cache-backed fetcher with optional seeding from a directory of previously downloaded raw files."""

    def __init__(self, cache: Cache, raw_dir: Optional[Path] = None, offline: bool = False) -> None:
        self.cache = cache
        self.raw_dir = Path(raw_dir) if raw_dir else None
        self.offline = offline

    def _seed(self, path: Path, seed_names: list[str]) -> None:
        if path.exists() or not self.raw_dir:
            return
        for name in seed_names:
            for sub in ("", "bgsu", "pdbe"):
                src = self.raw_dir / sub / name
                if src.is_file():
                    shutil.copy2(src, path)
                    return

    def text(self, url: str, suffix: str, seed_name: str | list[str] | None = None) -> tuple[Optional[str], SourceMeta]:
        path = self.cache._path(url, suffix, "GET", None)
        self._seed(path, [seed_name] if isinstance(seed_name, str) else list(seed_name or []))
        if not path.exists() and self.offline:
            return None, SourceMeta(url, "http_error", None, "offline and not cached")
        try:
            txt = self.cache.fetch_text(url, suffix)
        except requests.HTTPError as e:
            return None, SourceMeta(url, "http_error", None, str(e))
        except requests.RequestException as e:
            return None, SourceMeta(url, "http_error", None, f"request failed: {e}")
        retrieved = _dt.datetime.fromtimestamp(path.stat().st_mtime, _dt.timezone.utc).date().isoformat()
        return txt, SourceMeta(url, "ok" if txt.strip() else "empty", retrieved)


# ---------------------------------------------------------------- BGSU nucleotide annotations
@dataclass
class ReferenceAnnotation:
    unit_id: str
    labels: list[str]          # all labels listed for the nucleotide, in source order (first = preferred)


def _split_labels(label: str) -> list[str]:
    return [x.strip() for x in re.split(r";|\|\|", label) if x.strip()]


def parse_reference_annotations(text: str) -> tuple[list[ReferenceAnnotation], list[str]]:
    """Family file: ``# comment`` lines, then ``unit_id, label`` (comma or tab separated)."""
    rows: list[ReferenceAnnotation] = []
    comments: list[str] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            comments.append(line.lstrip("# ").strip())
            continue
        if "\t" in line:
            unit, label = line.split("\t", 1)
        elif "," in line:
            unit, label = line.split(",", 1)
        else:
            # not a ``unit_id, label`` row (e.g. BGSU's plain-text "Could not find annotation file ..." reply): never a label
            continue
        rows.append(ReferenceAnnotation(unit.strip(), _split_labels(label)))
    return rows, comments


def fetch_family_reference(fetcher: AnnotationFetcher, acc: str) -> tuple[list[ReferenceAnnotation], list[str], SourceMeta]:
    url = f"{BGSU}/data/nucleotide_annotations/{acc}_annotations.txt"
    txt, meta = fetcher.text(url, ".txt", f"{acc}_annotations.txt")
    if txt is None:
        return [], [], meta
    rows, comments = parse_reference_annotations(txt)
    return rows, comments, meta


def fetch_chain_annotations(fetcher: AnnotationFetcher, pdb_id: str, auth_asym_id: str, model: int = 1) -> tuple[list[ReferenceAnnotation], SourceMeta]:
    pdb = pdb_id.upper()
    url = f"{BGSU}/correspondence/nucleotide_annotation?chain={pdb}%7C{model}%7C{auth_asym_id}"
    txt, meta = fetcher.text(url, ".txt", f"nucleotide_annotation_{pdb}_{model}_{auth_asym_id}.txt")
    if txt is None or not txt.strip():
        return [], meta
    rows, _ = parse_reference_annotations(txt)
    if not rows:
        # BGSU answers HTTP 200 with a plain-text message when it has no annotation file for the family (e.g. RF00162): the source
        # is unavailable for this chain, which is recorded as such (never as "no structural elements")
        first = txt.strip().splitlines()[0][:200]
        status = "not_available" if first.lower().startswith("could not find annotation file") else "parse_error"
        return [], replace(meta, status=status, detail=first)
    return rows, meta


# ---------------------------------------------------------------- BGSU Motif Atlas
@dataclass
class LoopMembership:
    loop_id: str
    units: list[str]

    @property
    def loop_type(self) -> str:
        return self.loop_id.split("_")[0]


def parse_loops_csv(text: str) -> list[LoopMembership]:
    out: list[LoopMembership] = []
    for row in csv.reader(io.StringIO(text)):
        if not row or not row[0].strip():
            continue
        units = [u.strip() for u in row[1].split(",") if u.strip()] if len(row) > 1 else []
        out.append(LoopMembership(row[0].strip(), units))
    return out


def fetch_entry_loops(fetcher: AnnotationFetcher, pdb_id: str) -> tuple[list[LoopMembership], SourceMeta]:
    pdb = pdb_id.upper()
    txt, meta = fetcher.text(f"{BGSU}/rna3dhub/loops/download/{pdb}", ".csv", f"loops_{pdb}.csv")
    return (parse_loops_csv(txt) if txt else []), meta


@dataclass
class MotifAssignment:
    loop_id: str
    atlas_status: str                  # valid | problematic
    location: Optional[str]
    annotation_text: Optional[str]     # None when the page says 'No text annotation'
    motif_group: Optional[str]
    problem: Optional[str]
    assignment_type: Optional[str] = None     # member (loop belongs to the release) | homologous | geometric (mapped through another loop)
    mapped_via_loop_id: Optional[str] = None  # the release loop through which a homologous/geometric assignment was made


_TAG = re.compile(r"<[^>]+>")
_ROW = re.compile(r"<tr>(.*?)</tr>", re.S)
_CELL = re.compile(r"<td[^>]*>(.*?)</td>", re.S)
_LOOP_LINK = re.compile(r"loops/view/((?:HL|IL|J\d)_[0-9A-Z]{4}_\d{3})")
_GROUP_LINK = re.compile(r"motif/view/((?:HL|IL|J\d)_\d+\.\d+)")


def _clean(s: str) -> str:
    return _html.unescape(_TAG.sub(" ", s.replace("<br>", "\n").replace("<br/>", "\n"))).strip()


def parse_motifs_html(text: str) -> list[MotifAssignment]:
    """Parse the per-entry motifs page: 'Valid loops' tables (loop id, location, annotation + motif group) and
    'Problematic loops' tables (loop id, problem, annotation)."""
    out: list[MotifAssignment] = []
    # split the page at the h3 headers so that each table knows its status
    parts = re.split(r"<h3>\s*(Valid loops|Problematic loops)\s*</h3>", text)
    for i in range(1, len(parts) - 1, 2):
        status = "valid" if parts[i].startswith("Valid") else "problematic"
        for row in _ROW.findall(parts[i + 1]):
            m = _LOOP_LINK.search(row)
            if not m:
                continue
            cells = [c for c in _CELL.findall(row)]
            if status == "valid":
                location = _clean(cells[2]) if len(cells) > 2 else None
                ann_cell = cells[3] if len(cells) > 3 else ""
                map_cell = cells[4] if len(cells) > 4 else ""
                g = _GROUP_LINK.search(ann_cell)
                ann_text = _clean(_TAG.sub(" ", ann_cell.split("<br>")[0]))
                if g and ann_text == g.group(1):
                    ann_text = ""
                assignment, via, group = None, None, None
                if g:
                    assignment, group = "member", g.group(1)
                else:
                    gm = _GROUP_LINK.search(map_cell)
                    if gm:
                        group = gm.group(1)
                        words = _clean(_TAG.sub(" ", map_cell)).replace(group, "").split()
                        via = next((w for w in words if _LOOP_LINK.pattern and re.fullmatch(r"(?:HL|IL|J\d)_[0-9A-Z]{4}_\d{3}", w)), None)
                        assignment = next((w.lower() for w in words if w != via), "mapped")
                out.append(MotifAssignment(m.group(1), status, location or None, None if ann_text in ("", "No text annotation") else ann_text,
                                           group, None, assignment, via))
            else:
                problem = _clean(cells[2]) if len(cells) > 2 else None
                ann_text = _clean(cells[3]) if len(cells) > 3 else ""
                out.append(MotifAssignment(m.group(1), status, None, None if ann_text in ("", "No text annotation") else ann_text, None, problem or None))
    return out


def fetch_entry_motifs(fetcher: AnnotationFetcher, pdb_id: str) -> tuple[list[MotifAssignment], SourceMeta]:
    pdb = pdb_id.upper()
    txt, meta = fetcher.text(f"{BGSU}/rna3dhub/pdb/{pdb}/motifs", ".html", f"motifs_{pdb}.html")
    if txt is None:
        return [], meta
    try:
        return parse_motifs_html(txt), meta
    except Exception as e:   # pragma: no cover - defensive
        meta.status, meta.detail = "parse_error", str(e)
        return [], meta


def fetch_motif_atlas_release(fetcher: AnnotationFetcher, loop_type: str) -> tuple[Optional[str], SourceMeta]:
    txt, meta = fetcher.text(f"{BGSU}/rna3dhub/motifs/release/{loop_type.lower()}/current", ".html", f"release_{loop_type.lower()}_current.html")
    if txt is None:
        return None, meta
    m = re.search(r"Release\s+(\d+\.\d+)", txt)
    return (m.group(1) if m else None), meta


# ---------------------------------------------------------------- PDBe binding sites
def fetch_binding_sites(fetcher: AnnotationFetcher, pdb_id: str, entity_id: int) -> tuple[Optional[dict[str, Any]], SourceMeta]:
    pdb = pdb_id.lower()
    seeds = [f"binding_sites_{pdb}_{entity_id}.json"] + ([f"binding_sites_{pdb}.json"] if int(entity_id) == 1 else [])
    txt, meta = fetcher.text(f"{PDBE_API_V2}/pdb/entry/binding_sites/{pdb}/{entity_id}", ".json", seeds)
    if txt is None:
        return None, meta
    try:
        d = json.loads(txt)
    except json.JSONDecodeError as e:
        meta.status, meta.detail = "parse_error", str(e)
        return None, meta
    ent = d.get(pdb)
    if not ent:
        meta.status = "empty"
        return None, meta
    return ent, meta
