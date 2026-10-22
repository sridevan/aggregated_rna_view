"""Nonstandard nucleotide components of the deposited RNA polymers of a family, mapped to Rfam consensus positions.

Method (validated in experiments/modifications_RF00167/REPORT.md):
  * Primary source: the deposited polymer sequence of each RNA chain instance in the cached PDBe updated mmCIF. The component at every
    label_seq_id of the instance's label_asym_id is read from _pdbx_poly_seq_scheme and checked against _entity_poly_seq, _atom_site
    (model 1, coordinates are not required) and the residue_name of the validated family mapping. Polymer membership therefore never
    depends on the CCD component type, and non-polymer ligands with the same CCD id are never read.
  * Standard = A, C, G, U. Any other reliably identified component of the RNA polymer is nonstandard. No biological classification.
  * Identity is "missing" (not counted) when it cannot be established: microheterogeneity, unknown components (N, UNK, X),
    disagreement between sequence records, coordinates or the mapping, or a polymer position absent from the mapping.
  * Consensus positions come only from the validated mapping (match state). Insert, local-end and flank residues keep a null position.
  * Occurrence per consensus position = instances with a nonstandard component / instances with a determined component (positive +
    negative); one count per instance; instances without a residue there, with missing identity, or whose mmCIF is unavailable are not counted.
  * Cross-validation (optional): PDBe modified_AA_or_NA records. Absence of a PDBe record is never used as a negative.
  * Enrichment (optional): CCD parent (wwPDB CCD), NAKB/BGSU standard_base and change locations. Missing enrichment never blocks the core output.
"""
from __future__ import annotations
import datetime as _dt
import hashlib
import json
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional
from . import __version__
from .annotation_mapping import InstanceIndex
from .cache import Cache
from .pdbe import family_uses_chain_files
from .export import validate_json

log = logging.getLogger("rna_explorer.modifications")
SCHEMA_VERSION = "1.0"
STANDARD = frozenset({"A", "C", "G", "U"})
UNKNOWN = frozenset({"N", "UNK", "X", "?", "."})
PDBE_FILES = "https://www.ebi.ac.uk/pdbe/entry-files/download"
PDBE_MODS_URL = "https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/modified_AA_or_NA"
NAKB_URL = "https://rna.bgsu.edu/modified/modified_to_change_data.json"
NAKB_PAGE = "https://www.nakb.org/modifiednt.html"
#: wwPDB Chemical Component Dictionary definition: PDBe first, RCSB only when PDBe cannot provide it (same atoms, leaving flags and parents for
#: all 41 components used so far, checked 2026-10-22)
CCD_URLS = ("https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{ccd}.cif", "https://files.rcsb.org/ligands/view/{ccd}.cif")
MISSING_REASONS = ("microheterogeneity", "unknown_component", "sequence_record_conflict", "coordinate_identity_conflict", "mapping_identity_conflict", "not_in_mapping")
METHOD = {
    "primary_source": "deposited RNA polymer sequence in the PDBe updated mmCIF (_pdbx_poly_seq_scheme.mon_id at label_asym_id + label_seq_id)",
    "consistency_checks": ["_entity_poly_seq", "_atom_site (model 1, where coordinates exist)", "residue_name of the validated family mapping"],
    "standard_components": sorted(STANDARD),
    "nonstandard_definition": "any other reliably identified component of the RNA polymer; components are not classified (natural, synthetic, engineered, post-transcriptional)",
    "missing_identity_reasons": list(MISSING_REASONS),
    "consensus_positions": "validated family mapping only (match state); insert, local-end and outside-mapped-region residues have no consensus position",
    "cross_validation": "PDBe modified_AA_or_NA (an absent PDBe record is not evidence of a standard nucleotide)",
    "enrichment": "CCD parent (mon_nstd_parent_comp_id, wwPDB CCD); NAKB/BGSU standard_base and change locations (chemistry reference only)",
}


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def _mtime(p: Path) -> str:
    return _dt.datetime.fromtimestamp(p.stat().st_mtime, _dt.timezone.utc).isoformat(timespec="seconds")


# ----------------------------------------------------------------------------- mmCIF extraction
def _rows(block, cat: str, cols: list[str]) -> list[list[str]]:
    import gemmi
    if not block.find_mmcif_category(cat):
        return []
    t = block.find(cat, cols)
    return [[gemmi.cif.as_string(r[i]) for i in range(len(cols))] for r in t]


def read_polymer(block, label_asym_id: str, entity_id: str) -> dict[str, Any]:
    """Deposited component records of one polymer chain: {label_seq_id: {...}} plus entity/chem_comp information."""
    scheme: dict[int, list[dict[str, str]]] = defaultdict(list)
    for a, e, sid, mon, pnum, ins, het in _rows(block, "_pdbx_poly_seq_scheme.", ["asym_id", "entity_id", "seq_id", "mon_id", "pdb_seq_num", "pdb_ins_code", "hetero"]):
        if a == label_asym_id:
            if e != entity_id:
                raise ValueError(f"asym {label_asym_id} belongs to entity {e}, expected {entity_id}")
            scheme[int(sid)].append({"mon": mon, "auth_seq_id": pnum, "ins": "" if ins in (".", "?") else ins, "hetero": het})
    if not scheme:
        raise ValueError(f"no polymer sequence scheme for asym {label_asym_id}")
    eps: dict[int, list[tuple[str, str]]] = defaultdict(list)
    for e, num, mon, het in _rows(block, "_entity_poly_seq.", ["entity_id", "num", "mon_id", "hetero"]):
        if e == entity_id:
            eps[int(num)].append((mon, het))
    coords: dict[int, set[str]] = defaultdict(set)
    first_model = None
    for a, sid, comp, model in _rows(block, "_atom_site.", ["label_asym_id", "label_seq_id", "label_comp_id", "pdbx_PDB_model_num"]):
        first_model = first_model or model
        if a == label_asym_id and model == first_model and sid not in (".", "?"):
            coords[int(sid)].add(comp)
    poly_type = {e: t for e, t in _rows(block, "_entity_poly.", ["entity_id", "type"])}.get(entity_id)
    chem = {r[0]: {"name": r[2], "type": r[1]} for r in _rows(block, "_chem_comp.", ["id", "type", "name"])}
    return {"scheme": dict(scheme), "entity_poly_seq": dict(eps), "coords": dict(coords), "entity_poly_type": poly_type, "chem_comp": chem}


def classify_position(seq_id: int, poly: dict[str, Any], mapping_residue: Optional[dict[str, Any]]) -> dict[str, Any]:
    """Status of one deposited polymer position: positive (nonstandard), negative (A/C/G/U) or missing (with reason)."""
    s = poly["scheme"][seq_id]
    mons = sorted({x["mon"] for x in s}); eps = sorted({m for m, _ in poly["entity_poly_seq"].get(seq_id, [])})
    coord = sorted(poly["coords"].get(seq_id, set()))
    comp = mons[0] if len(mons) == 1 else None
    rec = {"label_seq_id": seq_id, "component": comp, "components": mons, "has_coordinates": bool(coord)}
    hetero = len(mons) > 1 or len(eps) > 1 or any(x["hetero"].lower() in ("y", "yes") for x in s) or any(h.lower() in ("y", "yes") for _, h in poly["entity_poly_seq"].get(seq_id, []))
    reason = None
    if hetero:
        reason = "microheterogeneity"
    elif comp in UNKNOWN:
        reason = "unknown_component"
    elif eps and eps != mons:
        reason = "sequence_record_conflict"
    elif coord and coord != mons:
        reason = "coordinate_identity_conflict"
    elif mapping_residue is None:
        reason = "not_in_mapping"
    elif mapping_residue["residue_name"] != comp:
        reason = "mapping_identity_conflict"
    if reason:
        return rec | {"status": "missing", "missing_reason": reason}
    return rec | {"status": "negative" if comp in STANDARD else "positive", "missing_reason": None}


def extract_instance(block, inst: dict[str, Any], mapping: dict[str, Any]) -> dict[str, Any]:
    """Per-instance record: deposited components (index = label_seq_id - 1), nonstandard observations and positions with missing identity."""
    poly = read_polymer(block, inst["label_asym_id"], str(inst["entity_id"]))
    mres = {int(r["polymer_seq_position"]): r for r in mapping["residues"]}
    ids = sorted(poly["scheme"])
    if ids != list(range(1, len(ids) + 1)):
        raise ValueError("label_seq_id values are not contiguous from 1")
    first, last = ids[0], ids[-1]
    classes = {sid: classify_position(sid, poly, mres.get(sid)) for sid in ids}
    observations, missing = [], []
    for sid, c in classes.items():
        s0 = poly["scheme"][sid][0]; m = mres.get(sid) or {}
        base = {"label_seq_id": sid, "auth_seq_id": s0["auth_seq_id"], "auth_ins_code": s0["ins"], "has_coordinates": c["has_coordinates"],
                "alignment_state": m.get("alignment_state"), "rfam_consensus_pos": m.get("rfam_consensus_pos") if m.get("alignment_state") == "match" else None,
                "insertion_anchor": m.get("insertion_anchor")}
        if c["status"] == "positive":
            ccd = c["component"]; cc = poly["chem_comp"].get(ccd, {})
            observations.append({"observation_id": f"{inst['seq_name']}:{sid}", "seq_name": inst["seq_name"], "pdb_id": inst["pdb_id"], "entity_id": int(inst["entity_id"]),
                                 "label_asym_id": inst["label_asym_id"], "auth_asym_id": inst["auth_asym_id"], **base,
                                 "ccd_id": ccd, "chem_comp_name": cc.get("name"), "chem_comp_type": cc.get("type"),
                                 "on_consensus_axis": base["alignment_state"] == "match", "outside_mapped_region": base["alignment_state"] == "not_submitted_flank",
                                 "is_5prime_terminal": sid == first, "is_3prime_terminal": sid == last,
                                 "mapping_nucleotide_parent": m.get("nucleotide_parent")})
        elif c["status"] == "missing":
            missing.append({**base, "components": c["components"], "reason": c["missing_reason"]})
    counts = Counter(c["status"] for c in classes.values())
    return {"seq_name": inst["seq_name"], "pdb_id": inst["pdb_id"], "entity_id": int(inst["entity_id"]), "label_asym_id": inst["label_asym_id"], "auth_asym_id": inst["auth_asym_id"],
            "status": "ok", "entity_poly_type": poly["entity_poly_type"], "n_polymer_residues": len(ids),
            "n_standard": counts.get("negative", 0), "n_nonstandard": counts.get("positive", 0), "n_missing_identity": counts.get("missing", 0),
            "n_without_coordinates": sum(not c["has_coordinates"] for c in classes.values()),
            "deposited_components": [classes[s]["component"] if classes[s]["component"] is not None else "|".join(classes[s]["components"]) for s in ids],
            "label_seq_ids_without_coordinates": [s for s in ids if not classes[s]["has_coordinates"]],
            "observations": observations, "missing_identity": missing, "_status_by_seq": {s: c["status"] for s, c in classes.items()},
            "_missing_reason_by_seq": {s: c["missing_reason"] for s, c in classes.items() if c["missing_reason"]}}


# ----------------------------------------------------------------------------- aggregation
def aggregate(records: list[dict[str, Any]], idx: InstanceIndex, clen: int) -> list[dict[str, Any]]:
    """Occurrence per consensus match position; failed instances count as source unavailable, never as negatives."""
    obs_by = {(o["seq_name"], o["label_seq_id"]): o for r in records for o in r["observations"]}
    positions = []
    for pos in range(1, clen + 1):
        n_pos = n_neg = n_none = n_unavail = n_nocoord = 0
        reasons: Counter = Counter(); ccd_inst: dict[str, set] = defaultdict(set); ccd_entry: dict[str, set] = defaultdict(set); entries: set = set(); support: list[str] = []
        for r in records:
            if r["status"] != "ok":
                n_unavail += 1; continue
            m = idx.match_residue(r["seq_name"], pos)
            if m is None:
                n_none += 1; continue
            sid = int(m["polymer_seq_position"]); st = r["_status_by_seq"].get(sid, "missing")
            if st == "missing":
                reasons[r["_missing_reason_by_seq"].get(sid, "not_in_mapping")] += 1; continue
            if sid in r["label_seq_ids_without_coordinates"]:
                n_nocoord += 1
            if st == "positive":
                o = obs_by[(r["seq_name"], sid)]
                n_pos += 1; ccd_inst[o["ccd_id"]].add(r["seq_name"]); ccd_entry[o["ccd_id"]].add(r["pdb_id"]); entries.add(r["pdb_id"]); support.append(o["observation_id"])
            else:
                n_neg += 1
        den = n_pos + n_neg
        positions.append({"rfam_consensus_pos": pos, "n_positive": n_pos, "n_negative": n_neg, "n_missing_identity": sum(reasons.values()), "missing_reasons": dict(reasons),
                          "n_no_residue": n_none, "n_source_unavailable": n_unavail, "denominator": den, "frequency": round(n_pos / den, 4) if den else None,
                          "n_supporting_entries": len(entries), "n_denominator_without_coordinates": n_nocoord,
                          "ccd": [{"ccd_id": c, "n_instances": len(ccd_inst[c]), "n_entries": len(ccd_entry[c])} for c in sorted(ccd_inst, key=lambda c: (-len(ccd_inst[c]), c))],
                          "supporting_observations": sorted(support)})
    return positions


# ----------------------------------------------------------------------------- optional sources (cached envelopes: status, time, headers, body)
def _envelope_path(cache: Cache, url: str, key: str) -> Path:
    return cache._path(url, ".envelope.json", "POST" if key else "GET", key or None)


def load_pdbe(cache: Cache, entries: list[str], mode: str, seed: Optional[Path] = None) -> dict[str, Any]:
    """PDBe modified_AA_or_NA for all entries in one cached batch request. mode: 'cache' (fetch when not cached), 'offline' (cache only), 'off'."""
    if mode == "off":
        return {"status": "not_checked", "records": {}}
    key = ",".join(sorted(e.lower() for e in entries)); p = _envelope_path(cache, PDBE_MODS_URL, key)
    if not p.exists() and seed and seed.exists():
        env = json.load(open(seed))
        if sorted(env.get("request_ids", [])) == sorted(key.split(",")):
            env = {"url": PDBE_MODS_URL, "method": "POST", "request_ids": env["request_ids"], "http_status": env["http_status"], "retrieved_utc": env["retrieved_utc"], "body": env["body"], "seeded_from": str(seed)}
            p.write_text(json.dumps(env))
    if not p.exists():
        if mode == "offline":
            return {"status": "unavailable", "reason": "not cached (offline)", "records": {}}
        try:
            import requests
            r = requests.post(PDBE_MODS_URL, data=json.dumps(key), headers={"Content-Type": "application/json", "accept": "*/*"}, timeout=cache.timeout)
            if r.status_code not in (200, 404):
                return {"status": "unavailable", "reason": f"HTTP {r.status_code}", "records": {}}
            p.write_text(json.dumps({"url": PDBE_MODS_URL, "method": "POST", "request_ids": key.split(","), "http_status": r.status_code, "retrieved_utc": _now(),
                                     "body": r.json() if r.status_code == 200 else None}))
        except Exception as e:   # network failure: cross-check unavailable, core output unaffected
            return {"status": "unavailable", "reason": str(e), "records": {}}
    env = json.load(open(p))
    recs = env["body"] if env["http_status"] == 200 and isinstance(env["body"], dict) else {}
    return {"status": "checked", "url": PDBE_MODS_URL, "http_status": env["http_status"], "retrieved_utc": env["retrieved_utc"], "records": recs,
            "n_entries_with_records": len(recs), "n_records": sum(len(v) for v in recs.values()), "cache_file": p.name}


def load_nakb(cache: Cache, mode: str, seed_dir: Optional[Path] = None) -> dict[str, Any]:
    if mode == "off":
        return {"status": "not_used", "components": {}}
    p = _envelope_path(cache, NAKB_URL, "")
    if not p.exists() and seed_dir and (seed_dir / "modified_to_change_data.json").exists():
        hdr = {}
        hf = seed_dir / "json_headers.txt"
        if hf.exists():
            for line in hf.read_text().splitlines():
                if ":" in line:
                    k, v = line.split(":", 1); hdr[k.strip().lower()] = v.strip()
        ret = (seed_dir / "retrieved_utc.txt").read_text().strip() if (seed_dir / "retrieved_utc.txt").exists() else _mtime(seed_dir / "modified_to_change_data.json")
        p.write_text(json.dumps({"url": NAKB_URL, "http_status": 200, "retrieved_utc": ret, "last_modified": hdr.get("last-modified"), "etag": hdr.get("etag"),
                                 "body": json.load(open(seed_dir / "modified_to_change_data.json")), "seeded_from": str(seed_dir)}))
    if not p.exists():
        if mode == "offline":
            return {"status": "unavailable", "reason": "not cached (offline)", "components": {}}
        try:
            import requests
            r = requests.get(NAKB_URL, timeout=cache.timeout); r.raise_for_status()
            p.write_text(json.dumps({"url": NAKB_URL, "http_status": r.status_code, "retrieved_utc": _now(), "last_modified": r.headers.get("Last-Modified"), "etag": r.headers.get("ETag"), "body": r.json()}))
        except Exception as e:
            return {"status": "unavailable", "reason": str(e), "components": {}}
    env = json.load(open(p)); body = env["body"]
    return {"status": "used", "url": NAKB_URL, "reference_page": NAKB_PAGE, "retrieved_utc": env["retrieved_utc"], "last_modified": env.get("last_modified"), "etag": env.get("etag"),
            "sha256": hashlib.sha256(json.dumps(body, sort_keys=True).encode()).hexdigest(), "n_components": len(body), "components": body}


def ccd_parent(cache: Cache, ccd: str, offline: bool) -> dict[str, Any]:
    """mon_nstd_parent_comp_id from the wwPDB CCD definition: PDBe file first, RCSB only when PDBe cannot provide it (cached)."""
    import gemmi
    reasons = []
    for url in (u.format(ccd=ccd) for u in CCD_URLS):
        p = cache._path(url, ".cif", "GET", None)
        if offline and not p.exists():
            reasons.append(f"{url}: offline and not cached"); continue
        try:
            b = gemmi.cif.read_string(cache.fetch_text(url, ".cif")).sole_block()
            if b.find_value("_chem_comp.id") is None:
                raise ValueError("no _chem_comp record")
            v = b.find_value("_chem_comp.mon_nstd_parent_comp_id")
            par = gemmi.cif.as_string(v) if v else None
            return {"status": "used", "parent": None if par in (None, "", "?", ".") else par, "url": url, "retrieved_utc": _mtime(p)}
        except Exception as e:
            reasons.append(f"{url}: {e}")
    return {"status": "unavailable", "parent": None, "reason": "; ".join(reasons)}


def component_info(ccd: str, names: Counter, types: Counter, ccd_rec: dict[str, Any], nakb: dict[str, Any]) -> dict[str, Any]:
    nk = nakb["components"].get(ccd) if nakb.get("components") else None
    locs = sorted({c.get("change_location") for c in (nk or {}).get("changes", []) if c.get("change_location")} | {c.get("change_location_2") for c in (nk or {}).get("changes", []) if c.get("change_location_2")})
    nakb_parent = (nk or {}).get("standard_base") or None
    if ccd_rec.get("parent"):
        parent, src = ccd_rec["parent"], "CCD"
    elif nakb_parent:
        parent, src = "/".join(nakb_parent), "NAKB"
    else:
        parent, src = None, None
    return {"ccd_id": ccd, "chem_comp_name": names.most_common(1)[0][0] if names else None, "chem_comp_type": types.most_common(1)[0][0] if types else None,
            "ccd_parent": ccd_rec.get("parent"), "ccd_parent_status": ccd_rec.get("status"),
            "nakb_in_reference": nk is not None if nakb.get("status") == "used" else None, "nakb_parent": nakb_parent, "nakb_change_locations": locs if nk else None,
            "nakb_name": (nk or {}).get("pdb", {}).get("name"), "nakb_modomics": (nk or {}).get("modomics"),
            "parent": parent, "parent_source": src}


def _stream_block(url: str) -> tuple[Any, dict[str, Any]]:
    """Download a whole-entry mmCIF to a temporary file, parse it and delete the file (not cached). Transient transport errors are retried."""
    import gemmi, hashlib, os, tempfile, time as _time
    import requests
    from .cache import RETRIES, RETRY_WAIT_S
    for attempt in range(RETRIES):
        try:
            r = requests.get(url, timeout=300); break
        except (requests.ConnectionError, requests.exceptions.ChunkedEncodingError, requests.Timeout):
            if attempt == RETRIES - 1:
                raise
            _time.sleep(RETRY_WAIT_S * (attempt + 1))
    r.raise_for_status()
    fd, tmp = tempfile.mkstemp(suffix=".cif")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(r.content)
        block = gemmi.cif.read(tmp).sole_block()
    finally:
        os.remove(tmp)
    return block, {"url": url, "retrieved_utc": _now(), "bytes": len(r.content), "md5": hashlib.md5(r.content).hexdigest(), "kept": False}


# ----------------------------------------------------------------------------- build
def build_modifications(acc: str, cache_dir: Path, out_dir: Path, pdbe: str = "cache", nakb: str = "cache", seed_dir: Optional[Path] = None,
                        offline: bool = False, write: bool = True) -> dict[str, Any]:
    import gemmi
    fam_dir = out_dir / "families" / acc
    idx = InstanceIndex(fam_dir); cache = Cache(cache_dir)
    generated = _now(); records: list[dict[str, Any]] = []; blocks: dict[str, Any] = {}; sources: dict[str, dict[str, Any]] = {}
    # families built from per-chain files (e.g. 70S ribosomes) stream each whole-entry mmCIF instead of keeping it: download to a temporary
    # file, parse, delete; the record keeps URL, retrieval time, size and MD5. Only the current entry's parsed block is held in memory.
    stream = family_uses_chain_files(fam_dir)
    for inst in sorted(idx.instances, key=lambda i: i["pdb_id"]) if stream else idx.instances:
        pdb = inst["pdb_id"]
        try:
            if pdb not in blocks:
                url = f"{PDBE_FILES}/{pdb.lower()}_updated.cif"
                if stream:
                    if offline:
                        raise FileNotFoundError("mmCIF is streamed, not cached (offline)")
                    blocks.clear()
                    blocks[pdb], sources[pdb] = _stream_block(url)
                else:
                    path = cache._path(url, ".cif", "GET", None)
                    if offline and not path.exists():
                        raise FileNotFoundError("mmCIF not cached (offline)")
                    path = cache.fetch_bytes(url, suffix=".cif")
                    blocks[pdb] = gemmi.cif.read(str(path)).sole_block(); sources[pdb] = {"url": url, "cache_file": path.name, "retrieved_utc": _mtime(path)}
            rec = extract_instance(blocks[pdb], inst, idx.mapping(inst["seq_name"]))
            rec["mmcif_source"] = sources[pdb]
        except Exception as e:
            log.error("%s failed: %s", inst["seq_name"], e)
            rec = {"seq_name": inst["seq_name"], "pdb_id": pdb, "entity_id": int(inst["entity_id"]), "label_asym_id": inst["label_asym_id"], "auth_asym_id": inst["auth_asym_id"],
                   "status": "source_unavailable", "error": str(e), "observations": [], "missing_identity": []}
        records.append(rec)

    # optional cross-validation and enrichment
    entries = sorted({i["pdb_id"] for i in idx.instances})
    pdbe_res = load_pdbe(cache, entries, "offline" if offline and pdbe == "cache" else pdbe, seed_dir / "pdbe" / "batch_all_entries.json" if seed_dir else None)
    nakb_res = load_nakb(cache, "offline" if offline and nakb == "cache" else nakb, seed_dir / "nakb" if seed_dir else None)
    pdbe_idx: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for e, recs in pdbe_res["records"].items():
        for x in recs:
            pdbe_idx[(e.lower(), x["struct_asym_id"], int(x["residue_number"]))].append(x)
    family_chains = {(i["pdb_id"].lower(), i["label_asym_id"]) for i in idx.instances}
    names: dict[str, Counter] = defaultdict(Counter); types: dict[str, Counter] = defaultdict(Counter)
    matched = set()
    for r in records:
        for o in r["observations"]:
            names[o["ccd_id"]][o["chem_comp_name"]] += 1; types[o["ccd_id"]][o["chem_comp_type"]] += 1
            k = (o["pdb_id"].lower(), o["label_asym_id"], o["label_seq_id"])
            if pdbe_res["status"] != "checked":
                o["pdbe"] = {"status": "not_checked"}
            elif k in pdbe_idx:
                x = pdbe_idx[k]; matched.add(k)
                o["pdbe"] = {"status": "reported", "chem_comp_id": x[0]["chem_comp_id"], "agrees": x[0]["chem_comp_id"] == o["ccd_id"], "n_records": len(x)}
            else:
                o["pdbe"] = {"status": "not_reported"}
    pdbe_unmatched = [{"pdb_id": k[0], "struct_asym_id": k[1], "residue_number": k[2], "chem_comp_id": v[0]["chem_comp_id"]} for k, v in pdbe_idx.items() if (k[0], k[1]) in family_chains and k not in matched]
    comps = {c: component_info(c, names[c], types[c], ccd_parent(cache, c, offline), nakb_res) for c in sorted(names)}
    for r in records:
        for o in r["observations"]:
            ci = comps[o["ccd_id"]]
            o.update({"ccd_parent": ci["ccd_parent"], "nakb_parent": ci["nakb_parent"], "nakb_change_locations": ci["nakb_change_locations"], "parent": ci["parent"], "parent_source": ci["parent_source"]})

    positions = aggregate(records, idx, idx.clen)
    obs = [o for r in records for o in r["observations"]]
    for c, ci in comps.items():
        os_ = [o for o in obs if o["ccd_id"] == c]
        ci.update({"n_observations": len(os_), "n_instances": len({o["seq_name"] for o in os_}), "n_entries": len({o["pdb_id"] for o in os_}),
                   "n_on_consensus_axis": sum(o["on_consensus_axis"] for o in os_), "consensus_positions": sorted({o["rfam_consensus_pos"] for o in os_ if o["on_consensus_axis"]}),
                   "pdbe_reported": sum(o["pdbe"]["status"] == "reported" for o in os_) if pdbe_res["status"] == "checked" else None})
    off = [o for o in obs if not o["on_consensus_axis"]]
    ok = [r for r in records if r["status"] == "ok"]
    counts = {"n_instances": len(records), "n_instances_ok": len(ok), "n_instances_source_unavailable": len(records) - len(ok), "n_entries": len({r["pdb_id"] for r in records}),
              "n_polymer_residues": sum(r["n_polymer_residues"] for r in ok), "n_standard": sum(r["n_standard"] for r in ok), "n_nonstandard": len(obs),
              "n_missing_identity": sum(r["n_missing_identity"] for r in ok), "missing_identity_reasons": dict(Counter(m["reason"] for r in ok for m in r["missing_identity"])),
              "n_residues_without_coordinates": sum(r["n_without_coordinates"] for r in ok),
              "n_instances_with_nonstandard": len({o["seq_name"] for o in obs}), "n_entries_with_nonstandard": len({o["pdb_id"] for o in obs}), "n_distinct_ccd": len(comps),
              "n_on_consensus_axis": len(obs) - len(off), "n_off_axis": len(off), "off_axis_by_state": dict(Counter(o["alignment_state"] or "not_in_mapping" for o in off)),
              "consensus_positions_with_nonstandard": sorted({o["rfam_consensus_pos"] for o in obs if o["on_consensus_axis"]})}
    mp = idx.mapping(idx.instances[0]["seq_name"])["provenance"]
    pdbe_prov = {k: v for k, v in pdbe_res.items() if k != "records"} | {"note": "cross-validation only; a missing PDBe record is not evidence of a standard nucleotide",
                                                                         "n_family_chain_records_without_observation": len(pdbe_unmatched), "family_chain_records_without_observation": pdbe_unmatched}
    nakb_prov = {k: v for k, v in nakb_res.items() if k != "components"} | {"note": "chemistry reference only (parent nucleotide, change locations); not evidence of occurrence in a chain"}
    prov = {"rfam_acc": acc, "rfam_release": mp.get("rfam_release"), "cm_checksum": mp.get("cm_checksum"), "cm_clen": idx.clen, "infernal_version": mp.get("infernal_version"),
            "mapping_pipeline_version": mp.get("pipeline_version"), "mapping_generated": mp.get("generated"), "pipeline_version": __version__, "generated": generated,
            "mmcif_source": f"{PDBE_FILES}/<pdb>_updated.cif (PDBe updated mmCIF, pipeline cache; per-instance retrieval in instances[].mmcif_source)",
            "pdbe_modified_AA_or_NA": pdbe_prov, "nakb": nakb_prov, "ccd": {"urls_in_order": list(CCD_URLS), "rule": "PDBe first; RCSB only when PDBe cannot provide the definition", "field": "_chem_comp.mon_nstd_parent_comp_id"},
            "investigation": "experiments/modifications_RF00167/REPORT.md"}
    for r in records:
        r.pop("_status_by_seq", None); r.pop("_missing_reason_by_seq", None)
    detailed = {"schema_version": SCHEMA_VERSION, "rfam_acc": acc, "clen": idx.clen, "method": METHOD, "counts": counts, "components": comps,
                "observations": obs, "instances": records, "provenance": prov}
    tracks = {"schema_version": SCHEMA_VERSION, "rfam_acc": acc, "clen": idx.clen, "track_id": "nonstandard_occurrence", "category": "Nonstandard nucleotides",
              "measurement": {"name": "Occurrence", "definition": "RNA instances with a nonstandard deposited component at the consensus position / RNA instances with a determined component there",
                              "unit": "fraction of RNA instances"},
              "aggregation_rules": {"positive": "nonstandard component (not A/C/G/U) at the consensus match position", "negative": "A, C, G or U at the position (coordinates not required)",
                                    "missing_identity": "component cannot be established (" + ", ".join(MISSING_REASONS) + "); not in the denominator",
                                    "no_residue": "no mapped residue at the position (deletion or outside the instance's mapped span); not in the denominator",
                                    "source_unavailable": "mmCIF unavailable for the instance; not in the denominator", "counting": "one count per RNA instance per position"},
              "summary": counts | {"off_axis_observations": [{"observation_id": o["observation_id"], "ccd_id": o["ccd_id"], "alignment_state": o["alignment_state"],
                                                                "is_5prime_terminal": o["is_5prime_terminal"], "is_3prime_terminal": o["is_3prime_terminal"], "has_coordinates": o["has_coordinates"]} for o in off]},
              "components": {c: {k: ci[k] for k in ("chem_comp_name", "chem_comp_type", "parent", "parent_source", "ccd_parent", "nakb_parent", "nakb_change_locations", "n_observations", "n_instances", "n_entries", "n_on_consensus_axis", "consensus_positions", "pdbe_reported")} for c, ci in comps.items()},
              "positions": positions, "method": METHOD, "provenance": prov}
    if write:
        validate_json(detailed, "modifications"); validate_json(tracks, "modification_tracks")
        with open(fam_dir / "modifications.json", "w") as f: json.dump(detailed, f, indent=1)
        with open(fam_dir / "modification_tracks.json", "w") as f: json.dump(tracks, f, indent=1)
    return {"acc": acc, "counts": counts, "pdbe": pdbe_res["status"], "nakb": nakb_res["status"], "modifications": detailed, "tracks": tracks}
