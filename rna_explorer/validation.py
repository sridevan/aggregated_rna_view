"""Mapping validation (from poc_rf00177/scripts/04_validate.py): mmCIF identifier agreement, coordinate cross-check, continuity, insertion bookkeeping."""
from __future__ import annotations
from collections import Counter
from typing import Any, Optional
from .infernal import ScoreRow
from .mapping import ChainMapping


def _runs(xs: list[int]) -> list[list[int]]:
    out: list[list[int]] = []
    for x in xs:
        if out and x == out[-1][1] + 1:
            out[-1][1] = x
        else:
            out.append([x, x])
    return out


def validate_chain(m: ChainMapping, scheme: dict[int, dict[str, str]], clen: int, rfam_hit: Optional[dict[str, Any]] = None,
                   score: Optional[ScoreRow] = None) -> dict[str, Any]:
    """Return a validation record for one chain mapping. `scheme` = _pdbx_poly_seq_scheme rows keyed by seq_id."""
    info = m.info
    mism: list[dict[str, Any]] = []
    obs_conf: list[dict[str, Any]] = []
    for r in m.rows:
        s = scheme.get(r.polymer_seq_position)
        if s is None:
            mism.append({"pos": r.polymer_seq_position, "problem": "absent from _pdbx_poly_seq_scheme"})
            continue
        chk = {"auth_asym_id": (s["pdb_strand_id"], r.auth_asym_id), "label_asym_id": (s["asym_id"], r.label_asym_id), "entity_id": (s["entity_id"], r.entity_id),
               "author_residue_number": (s["pdb_seq_num"], r.author_residue_number),
               "author_insertion_code": ("" if s["pdb_ins_code"] in (".", "?") else s["pdb_ins_code"], r.author_insertion_code), "residue_name": (s["mon_id"], r.residue_name)}
        bad = {k: v for k, v in chk.items() if str(v[0]) != str(v[1])}
        if bad:
            mism.append({"pos": r.polymer_seq_position, "mismatch": bad})
        cif_modelled = s["auth_seq_num"] not in ("?", ".")
        api_obs = (r.observed_ratio_api or 0) > 0
        if cif_modelled != bool(r.has_coordinates) or api_obs != bool(r.has_coordinates):
            obs_conf.append({"pos": r.polymer_seq_position, "scheme_modelled": cif_modelled, "atom_site": bool(r.has_coordinates), "api_observed": api_obs})
    match = [(r.polymer_seq_position, r.rfam_consensus_pos) for r in m.rows if r.alignment_state == "match"]
    mono = all(b[1] > a[1] for a, b in zip(match, match[1:]))   # type: ignore[operator]
    ins_ok = anchors_ok = True
    prev_cons = 0
    last_anchor: Optional[int] = None
    last_off = 0
    anchor_after_deletion = 0
    for r in m.rows:
        if r.alignment_state == "match":
            prev_cons = r.rfam_consensus_pos or 0
            last_anchor, last_off = None, 0
        elif r.alignment_state in ("insert", "local_end"):
            a, o = r.insertion_anchor or 0, r.insertion_offset or 0
            if a < prev_cons:
                anchors_ok = False
            if a > prev_cons and last_anchor != a:
                anchor_after_deletion += 1
            if (last_anchor == a and o != last_off + 1) or (last_anchor != a and o != 1) or r.rfam_consensus_pos is not None:
                ins_ok = False
            last_anchor, last_off = a, o
    cpos = Counter(c.rfam_consensus_pos for c in m.cons_rows)
    cons_complete = sorted(cpos) == list(range(1, clen + 1)) and max(cpos.values()) == 1
    unobs = [r.polymer_seq_position for r in m.rows if r.in_pdbe_mapped_region and not r.has_coordinates]
    rec: dict[str, Any] = {
        "seq_name": info.seq_name, "pdb_id": info.pdb_id, "entity_id": info.entity_id, "label_asym_id": info.label_asym_id, "auth_asym_id": info.auth_asym_id,
        "chain_length": len(info.sequence), "scheme_rows": len(scheme), "n_mmcif_mismatches": len(mism), "mmcif_mismatches": mism[:20],
        "n_observed_conflicts": len(obs_conf), "observed_conflicts": obs_conf[:20], "n_modelled_residues": len(info.observed),
        "mapped_residues_without_coordinates": len(unobs), "unmodelled_runs_in_mapped_region": _runs(unobs),
        "consensus_monotonic": mono, "insertion_offsets_consistent": ins_ok, "insertion_anchors_not_before_last_match": anchors_ok,
        "insertion_sites_whose_anchor_follows_deleted_consensus": anchor_after_deletion, "consensus_columns_complete_and_unique": cons_complete,
        "states": dict(Counter(r.alignment_state for r in m.rows)),
        "author_numbering_continuous": all(_int(m.rows[i + 1].author_residue_number) - _int(m.rows[i].author_residue_number) == 1 for i in range(len(m.rows) - 1)),
        "label_vs_auth_chain_differ": info.label_asym_id != info.auth_asym_id,
        "modified_in_mapped_region": {str(k): v.get("three_letter_code") for k, v in info.modified_in_mapped_region.items()},
        "species": {"api": info.api_organism, "mmcif": info.cif_organism, "mmcif_tax_id": info.cif_tax_id}}
    if score:
        rec["cmalign"] = score.to_dict()
    if rfam_hit:
        rec["rfam_structures_endpoint"] = {k: rfam_hit.get(k) for k in ("cm_start", "cm_end", "bit_score", "evalue_score", "pdb_start", "pdb_end")}
        rec["pdbe_region_agrees_with_rfam"] = (info.map_start, info.map_end) == (rfam_hit.get("pdb_start"), rfam_hit.get("pdb_end"))
        if score:
            rec["cm_span_agrees_with_rfam"] = (score.cm_from, score.cm_to) == (rfam_hit.get("cm_start"), rfam_hit.get("cm_end"))
    rec["ok"] = not mism and not obs_conf and mono and ins_ok and anchors_ok and cons_complete
    return rec


def _int(x: Any) -> int:
    try:
        return int(x)
    except (TypeError, ValueError):
        return -10 ** 9
