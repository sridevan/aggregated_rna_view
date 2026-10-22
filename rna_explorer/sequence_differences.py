"""Reported sequence differences of the family RNA chains, read per chain from the PDBe updated mmCIF and mapped to Rfam consensus positions.

Why mmCIF: the PDBe API mutated_AA_or_NA currently covers protein chains only ("at present it does not provide information about mutated
nucleotides in RNA or DNA chains"; experiments/mutations_investigation/REPORT.md). The same information is deposited in the mmCIF category
_struct_ref_seq_dif, which lists each residue where a polymer's deposited sequence differs from its reference sequence (_struct_ref).
The full method, written so it can be rebuilt elsewhere (e.g. at PDBe), is docs/SEQUENCE_DIFFERENCES_METHOD.md.

Per chain (extract_chain_differences, independent of Rfam):
  * a row belongs to an entity through _struct_ref_seq_dif.align_id -> _struct_ref_seq.ref_id -> _struct_ref.entity_id, and to a chain through
    its author strand (pdbx_pdb_strand_id, else _struct_ref_seq.pdbx_strand_id) -> label asym id via _pdbx_poly_seq_scheme (entity + strand);
  * deposited residue = mon_id at seq_num (= label_seq_id), reference residue = db_mon_id; change kind: substitution (both present),
    added_residue (no reference residue), absent_residue (no deposited residue / no seq_num);
  * the classification is the deposited 'details' text, kept verbatim (not normalised, not reinterpreted);
  * the deposited residue is checked against _pdbx_poly_seq_scheme.mon_id at that label_seq_id (identity check) and the author number against
    pdbx_auth_seq_num / pdb_ins_code there.
Family layer: each row is placed through the chain's own validated mapping (identity and author numbering checked again); match-state residues
get the consensus position; inserts, flanks, local ends, absent residues and rows that fail a check keep a null position and stay listed.
Counting: one observation per row; per position and classification, distinct RNA chain instances and PDB entries. No frequency is computed:
most entities have no external reference sequence, so the absence of a row is not evidence of an unchanged nucleotide.
"""
from __future__ import annotations
import datetime as _dt
import json
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional
from . import __version__
from .annotation_mapping import InstanceIndex
from .cache import Cache
from .export import validate_json
from .modifications import PDBE_FILES, _mtime, _stream_block
from .pdbe import family_uses_chain_files

log = logging.getLogger("rna_explorer.sequence_differences")
SCHEMA_VERSION = "1.0"
MISSING = ("", "?", ".")
CHANGE_KINDS = ("substitution", "added_residue", "absent_residue")
SEQ_DIF_COLS = ["align_id", "pdbx_pdb_strand_id", "seq_num", "pdbx_auth_seq_num", "pdbx_pdb_ins_code", "mon_id", "db_mon_id", "pdbx_seq_db_seq_num", "details", "pdbx_ordinal"]
METHOD = {
    "source": "PDBe updated mmCIF, category _struct_ref_seq_dif (deposited vs reference sequence of each polymer entity)",
    "why_not_api": "PDBe API mutated_AA_or_NA covers protein chains only at present (endpoint description); no RNA records for any family",
    "entity_join": "_struct_ref_seq_dif.align_id -> _struct_ref_seq.ref_id -> _struct_ref.entity_id (reference db_name, db_code, pdbx_db_accession kept)",
    "chain_join": "author strand (pdbx_pdb_strand_id, else _struct_ref_seq.pdbx_strand_id) + entity -> label asym id via _pdbx_poly_seq_scheme",
    "residue": "seq_num = label_seq_id; mon_id = deposited residue; db_mon_id = reference residue",
    "change_kinds": {"substitution": "deposited and reference residues present", "added_residue": "deposited residue without a reference residue",
                     "absent_residue": "reference residue without a deposited residue (no seq_num)"},
    "classification": "the deposited _struct_ref_seq_dif.details text, verbatim (case and wording preserved)",
    "checks": ["deposited residue = _pdbx_poly_seq_scheme.mon_id at label_seq_id", "author number + insertion code = _pdbx_poly_seq_scheme.pdb_seq_num/pdb_ins_code",
               "deposited residue = residue_name of the validated family mapping", "copies of one entity report the same differences"],
    "consensus_positions": "validated family mapping only (match state); insert, flank, local-end, absent and conflicting rows keep a null position",
    "counting": "one observation per row; per position and classification: distinct RNA chain instances and PDB entries",
    "absence": "not a negative: most RNA entities reference their own PDB sequence and cannot record differences",
}


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def _val(v: Optional[str]) -> Optional[str]:
    return None if v is None or v.strip() in MISSING else v.strip()


def _rows(block, cat: str, cols: list[str]) -> list[list[Optional[str]]]:
    """rows of an mmCIF category; columns prefixed with '?' are optional and read as None when the file does not have them"""
    import gemmi
    if not block.find_mmcif_category(cat):
        return []
    t = block.find(cat, cols)
    return [[gemmi.cif.as_string(r[i]) if r.has(i) else None for i in range(len(cols))] for r in t]


def change_kind(deposited: Optional[str], reference: Optional[str]) -> str:
    return "substitution" if deposited and reference else "added_residue" if deposited else "absent_residue"


def extract_chain_differences(block, pdb_id: str) -> dict[str, Any]:
    """Every _struct_ref_seq_dif row of an entry resolved to (entity, label asym id) and checked against the deposited sequence scheme.
    Independent of Rfam: returns rows for all polymer entities, with the entity's polymer type."""
    ref = {r[0]: {"entity_id": r[1], "db_name": _val(r[2]), "db_code": _val(r[3]), "accession": _val(r[4])}
           for r in _rows(block, "_struct_ref.", ["id", "entity_id", "?db_name", "?db_code", "?pdbx_db_accession"])}
    align = {r[0]: {"ref_id": r[1], "strand": _val(r[2])} for r in _rows(block, "_struct_ref_seq.", ["align_id", "ref_id", "?pdbx_strand_id"])}
    poly_type = {r[0]: r[1] for r in _rows(block, "_entity_poly.", ["entity_id", "type"])}
    scheme_asym: dict[tuple[str, str], str] = {}
    scheme: dict[tuple[str, int], dict[str, str]] = {}
    for asym, ent, seq, mon, strand, auth, ins in _rows(block, "_pdbx_poly_seq_scheme.", ["asym_id", "entity_id", "seq_id", "mon_id", "pdb_strand_id", "?pdb_seq_num", "?pdb_ins_code"]):
        scheme_asym[(ent, strand)] = asym
        scheme.setdefault((asym, int(seq)), {"mon_id": mon, "auth": _val(auth), "ins": _val(ins) or ""})   # first of microheterogeneous rows
    rows: list[dict[str, Any]] = []
    for v in (dict(zip(SEQ_DIF_COLS, r)) for r in _rows(block, "_struct_ref_seq_dif.", ["align_id"] + ["?" + c for c in SEQ_DIF_COLS[1:]])):
        a = align.get(v["align_id"], {}); rf = ref.get(a.get("ref_id"), {})
        ent = rf.get("entity_id"); strand = _val(v["pdbx_pdb_strand_id"]) or a.get("strand")
        asym = scheme_asym.get((ent, strand)) if ent and strand else None
        seq = _val(v["seq_num"]); dep = _val(v["mon_id"]); refres = _val(v["db_mon_id"])
        row = {"pdb_id": pdb_id.lower(), "entity_id": ent, "polymer_type": poly_type.get(ent) if ent else None, "auth_asym_id": strand, "label_asym_id": asym,
               "label_seq_id": int(seq) if seq else None, "author_residue_number": _val(v["pdbx_auth_seq_num"]), "author_insertion_code": _val(v["pdbx_pdb_ins_code"]) or "",
               "deposited": dep, "reference": refres, "reference_seq_num": _val(v["pdbx_seq_db_seq_num"]), "classification": _val(v["details"]),
               "change_kind": change_kind(dep, refres), "reference_sequence": {k: rf.get(k) for k in ("db_name", "db_code", "accession")}, "ordinal": _val(v["pdbx_ordinal"])}
        if not ent:
            row["chain_status"] = "entity_unresolved"
        elif not asym:
            row["chain_status"] = "chain_unresolved"
        else:
            row["chain_status"] = "resolved"
            s = scheme.get((asym, row["label_seq_id"])) if row["label_seq_id"] is not None else None
            row["scheme_check"] = ("no_deposited_residue" if row["label_seq_id"] is None else "not_in_scheme" if s is None
                                   else "agrees" if s["mon_id"] == dep and (row["author_residue_number"] is None or (s["auth"] == row["author_residue_number"] and s["ins"] == row["author_insertion_code"]))
                                   else "identity_conflict" if s["mon_id"] != dep else "author_number_conflict")
        rows.append(row)
    return {"rows": rows, "n_rows": len(rows), "reference_sequences": list(ref.values())}


def place_on_family(row: dict[str, Any], mapping_res: Optional[dict[str, Any]]) -> dict[str, Any]:
    """Consensus placement of one resolved family-chain row through the chain's validated mapping (null position unless every check passes)."""
    out = {"alignment_state": None, "rfam_consensus_pos": None, "has_coordinates": None}
    if row["label_seq_id"] is None:
        return out | {"mapping_status": "no_deposited_residue"}
    if row.get("scheme_check") not in ("agrees",):
        return out | {"mapping_status": row.get("scheme_check") or "not_checked"}
    if mapping_res is None:
        return out | {"mapping_status": "not_in_mapping"}
    out |= {"alignment_state": mapping_res["alignment_state"], "has_coordinates": bool(mapping_res["has_coordinates"])}
    if mapping_res["residue_name"] != row["deposited"]:
        return out | {"mapping_status": "mapping_identity_conflict"}
    if row["author_residue_number"] is not None and (str(mapping_res["author_residue_number"]) != str(row["author_residue_number"]) or (mapping_res["author_insertion_code"] or "") != row["author_insertion_code"]):
        return out | {"mapping_status": "mapping_author_conflict"}
    if mapping_res["alignment_state"] != "match":
        return out | {"mapping_status": f"no_consensus_position:{mapping_res['alignment_state']}"}
    return out | {"rfam_consensus_pos": mapping_res["rfam_consensus_pos"], "mapping_status": "mapped"}


def copy_consistency(obs: list[dict[str, Any]], idx: InstanceIndex) -> list[dict[str, Any]]:
    """Copies (chains) of one entity in one entry share a sequence, so they should report the same differences. Lists entities where they do not."""
    chains: dict[tuple[str, str], list[str]] = defaultdict(list)
    for i in idx.instances:
        chains[(i["pdb_id"].lower(), str(i["entity_id"]))].append(i["label_asym_id"])
    sig: dict[tuple[str, str, str], set] = defaultdict(set)
    for o in obs:
        sig[(o["pdb_id"], str(o["entity_id"]), o["label_asym_id"])].add((o["label_seq_id"], o["reference_seq_num"], o["deposited"], o["reference"], o["classification"]))
    out = []
    for (pdb, ent), labs in sorted(chains.items()):
        if len(labs) < 2:
            continue
        sets = {lab: sig.get((pdb, ent, lab), set()) for lab in sorted(labs)}
        if len({frozenset(s) for s in sets.values()}) > 1:
            out.append({"pdb_id": pdb, "entity_id": ent, "rows_per_chain": {lab: len(s) for lab, s in sets.items()}})
    return out


def aggregate(obs: list[dict[str, Any]], clen: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    on = [o for o in obs if o["rfam_consensus_pos"] is not None]
    by_pos: dict[int, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for o in on:
        by_pos[o["rfam_consensus_pos"]][o["classification"] or "(no classification)"].append(o)
    positions = []
    for pos in range(1, clen + 1):
        classes = {}
        for cls, os_ in sorted(by_pos.get(pos, {}).items()):
            classes[cls] = {"n_observations": len(os_), "n_instances": len({o["seq_name"] for o in os_}), "n_entries": len({o["pdb_id"] for o in os_}),
                            "changes": dict(Counter(f"{o['reference'] or '-'}>{o['deposited'] or '-'}" for o in os_)),
                            "change_kinds": dict(Counter(o["change_kind"] for o in os_)), "observation_ids": [o["observation_id"] for o in os_]}
        positions.append({"rfam_consensus_pos": pos, "classifications": classes})
    summary = {}
    for cls in sorted({o["classification"] or "(no classification)" for o in obs}):
        os_ = [o for o in obs if (o["classification"] or "(no classification)") == cls]
        summary[cls] = {"n_observations": len(os_), "n_on_consensus_axis": sum(o["rfam_consensus_pos"] is not None for o in os_), "n_instances": len({o["seq_name"] for o in os_}),
                        "n_entries": len({o["pdb_id"] for o in os_}), "change_kinds": dict(Counter(o["change_kind"] for o in os_)),
                        "consensus_positions": sorted({o["rfam_consensus_pos"] for o in os_ if o["rfam_consensus_pos"] is not None})}
    return positions, summary


def build_sequence_differences(acc: str, cache_dir: Path, out_dir: Path, offline: bool = False, write: bool = True) -> dict[str, Any]:
    import gemmi
    fam_dir = out_dir / "families" / acc
    idx = InstanceIndex(fam_dir); cache = Cache(cache_dir)
    generated = _now(); stream = family_uses_chain_files(fam_dir)
    fam_chain = {(i["pdb_id"].lower(), i["label_asym_id"]): i for i in idx.instances}
    fam_entities = {(i["pdb_id"].lower(), str(i["entity_id"])) for i in idx.instances}
    entries = sorted({i["pdb_id"].lower() for i in idx.instances})
    by_seq: dict[str, dict[int, dict[str, Any]]] = {}
    obs: list[dict[str, Any]] = []; unresolved: list[dict[str, Any]] = []; sources: dict[str, Any] = {}; other_rows = Counter()
    for pdb in entries:
        url = f"{PDBE_FILES}/{pdb}_updated.cif"
        try:
            if stream:
                if offline:
                    raise FileNotFoundError("mmCIF is streamed, not cached (offline)")
                block, src = _stream_block(url)
            else:
                path = cache._path(url, ".cif", "GET", None)
                if offline and not path.exists():
                    raise FileNotFoundError("mmCIF not cached (offline)")
                path = cache.fetch_bytes(url, suffix=".cif")
                block = gemmi.cif.read(str(path)).sole_block(); src = {"url": url, "cache_file": path.name, "retrieved_utc": _mtime(path)}
            res = extract_chain_differences(block, pdb)
        except Exception as e:
            log.error("%s failed: %s", pdb, e)
            sources[pdb] = {"status": "source_unavailable", "error": str(e)}
            continue
        sources[pdb] = {"status": "evaluated", "n_rows_all_entities": res["n_rows"]} | src
        for k, row in enumerate(res["rows"]):
            if (pdb, str(row["entity_id"])) not in fam_entities:
                other_rows[row["polymer_type"] or "unknown"] += 1      # other polymers of the entry (proteins, other RNAs): counted, not shown
                continue
            inst = fam_chain.get((pdb, row["label_asym_id"])) if row["chain_status"] == "resolved" else None
            row["observation_id"] = f"{pdb}|{row['label_asym_id'] or row['auth_asym_id'] or '?'}|{row['label_seq_id'] if row['label_seq_id'] is not None else 'ref' + str(row['reference_seq_num'])}|{row['ordinal'] or k}"
            if inst is None:
                unresolved.append(row | {"mapping_status": row["chain_status"] if row["chain_status"] != "resolved" else "chain_not_in_family"})
                continue
            if inst["seq_name"] not in by_seq:
                by_seq[inst["seq_name"]] = {int(r["polymer_seq_position"]): r for r in idx.residues(inst["seq_name"])}
            mres = by_seq[inst["seq_name"]].get(row["label_seq_id"]) if row["label_seq_id"] is not None else None
            obs.append(row | {"seq_name": inst["seq_name"]} | place_on_family(row, mres))
    positions, summary = aggregate(obs, idx.clen)
    inconsistent = copy_consistency(obs, idx)
    off = [o for o in obs if o["rfam_consensus_pos"] is None]
    evaluated = [e for e, s in sources.items() if s["status"] == "evaluated"]
    counts = {"n_entries": len(entries), "n_entries_evaluated": len(evaluated), "n_entries_source_unavailable": len(entries) - len(evaluated),
              "n_instances": len(idx.instances), "n_instances_evaluated": sum(1 for i in idx.instances if i["pdb_id"].lower() in evaluated),
              "n_observations": len(obs), "n_on_consensus_axis": len(obs) - len(off), "n_off_axis": len(off), "n_unresolved_rows": len(unresolved),
              "n_instances_with_differences": len({o["seq_name"] for o in obs}), "n_entries_with_differences": len({o["pdb_id"] for o in obs}),
              "mapping_status": dict(Counter(o["mapping_status"] for o in obs)), "change_kinds": dict(Counter(o["change_kind"] for o in obs)),
              "classifications": dict(Counter(o["classification"] or "(no classification)" for o in obs)),
              "reference_databases": dict(Counter(o["reference_sequence"]["db_name"] or "?" for o in obs)),
              "rows_on_other_entities_by_polymer_type": dict(other_rows), "n_entities_with_inconsistent_copies": len(inconsistent)}
    mp = idx.mapping(idx.instances[0]["seq_name"])["provenance"]
    prov = {"rfam_acc": acc, "rfam_release": mp.get("rfam_release"), "cm_checksum": mp.get("cm_checksum"), "cm_clen": idx.clen, "infernal_version": mp.get("infernal_version"),
            "mapping_pipeline_version": mp.get("pipeline_version"), "mapping_generated": mp.get("generated"), "pipeline_version": __version__, "generated": generated,
            "mmcif_source": f"{PDBE_FILES}/<pdb>_updated.cif (PDBe updated mmCIF; per-entry retrieval in entry_sources)",
            "method_document": "docs/SEQUENCE_DIFFERENCES_METHOD.md", "investigation": "experiments/mutations_investigation/REPORT.md"}
    detailed = {"schema_version": SCHEMA_VERSION, "rfam_acc": acc, "clen": idx.clen, "method": METHOD, "counts": counts, "classifications": summary,
                "observations": obs, "unresolved_rows": unresolved, "inconsistent_copies": inconsistent, "entry_sources": sources, "provenance": prov}
    tracks = {"schema_version": SCHEMA_VERSION, "rfam_acc": acc, "clen": idx.clen, "category": "Reported sequence differences",
              "measurement": {"name": "Reported sequence difference", "definition": "a _struct_ref_seq_dif row of a family RNA chain at the consensus position, per deposited classification",
                              "unit": "RNA chain instances and PDB entries (no frequency: absence of a row is not evidence of an unchanged nucleotide)"},
              "classifications": summary, "positions": positions,
              "off_axis": [{k: o[k] for k in ("observation_id", "pdb_id", "auth_asym_id", "label_asym_id", "seq_name", "label_seq_id", "author_residue_number", "author_insertion_code",
                                               "reference", "deposited", "change_kind", "classification", "alignment_state", "mapping_status")} for o in off],
              "observations": {o["observation_id"]: {k: o[k] for k in ("pdb_id", "auth_asym_id", "label_asym_id", "seq_name", "label_seq_id", "author_residue_number", "author_insertion_code",
                                                                       "reference", "deposited", "change_kind", "classification", "rfam_consensus_pos")} for o in obs if o["rfam_consensus_pos"] is not None},
              "counts": counts, "method": METHOD, "provenance": prov}
    if write:
        validate_json(detailed, "sequence_differences"); validate_json(tracks, "sequence_difference_tracks")
        with open(fam_dir / "sequence_differences.json", "w") as f: json.dump(detailed, f, indent=1)
        with open(fam_dir / "sequence_difference_tracks.json", "w") as f: json.dump(tracks, f, indent=1)
    return {"acc": acc, "counts": counts, "detailed": detailed, "tracks": tracks}
