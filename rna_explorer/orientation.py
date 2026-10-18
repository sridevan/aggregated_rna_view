"""Common reference orientation for the structures of a family: rigid-body (Kabsch) superposition of each chain instance onto the
display representative, using corresponding consensus positions from the validated residue mapping.

Convention (documented for the frontend): x_ref = R @ x_target + t for column vectors, R row-major 3x3 with det(R) = +1 (no reflection,
no scaling). In Mol* this is Mat4.ofRows([[R00,R01,R02,t0],[R10,R11,R12,t1],[R20,R21,R22,t2],[0,0,0,1]]) applied with
TransformStructureConformation (transpose: false), i.e. the target structure is moved into the reference coordinate frame.
"""
from __future__ import annotations
import datetime as _dt
import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional
import numpy as np
from . import __version__

ATOM_PRIORITY: tuple[str, ...] = ("C4'", "P", "C1'")
PP_STRICT = ("*",)                 # >= 0.95
PP_RELAXED = ("*", "9")            # >= 0.85
MIN_PAIRS = 6
RELAX_BELOW = 10
SCHEMA_VERSION = "1.0"


@dataclass
class KabschResult:
    rotation: list[list[float]]
    translation: list[float]
    rmsd: float
    n_pairs: int
    determinant: float


def kabsch(target: np.ndarray, reference: np.ndarray) -> KabschResult:
    """Optimal rotation R and translation t minimising |R p + t - q| over paired points (p = target, q = reference). Reflections are
    rejected (det(R) = +1) so the molecule's handedness is preserved."""
    if target.shape != reference.shape or target.ndim != 2 or target.shape[1] != 3:
        raise ValueError("target and reference must be N x 3 arrays of equal size")
    n = target.shape[0]
    if n < 3:
        raise ValueError("at least 3 paired points are required")
    cp, cq = target.mean(axis=0), reference.mean(axis=0)
    p, q = target - cp, reference - cq
    h = p.T @ q
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    if d == 0:
        d = 1.0
    corr = np.diag([1.0, 1.0, d])
    r = vt.T @ corr @ u.T
    t = cq - r @ cp
    moved = (r @ target.T).T + t
    rmsd = float(np.sqrt(((moved - reference) ** 2).sum(axis=1).mean()))
    return KabschResult(r.tolist(), t.tolist(), rmsd, n, float(np.linalg.det(r)))


@dataclass
class PairSelection:
    pairs: list[tuple[int, int, int, str]]          # (consensus_pos, ref_label_seq_id, target_label_seq_id, atom)
    pp_filter: str
    n_shared_consensus: int
    n_shared_with_coords: int
    excluded: dict[str, int] = field(default_factory=dict)
    atoms_used: dict[str, int] = field(default_factory=dict)


def _match_rows(mapping: dict[str, Any]) -> dict[int, dict[str, Any]]:
    return {r["rfam_consensus_pos"]: r for r in mapping["residues"] if r["alignment_state"] == "match" and r["rfam_consensus_pos"] is not None}


def select_pairs(ref_map: dict[str, Any], tgt_map: dict[str, Any], ref_xyz: dict[int, dict[str, tuple[float, float, float]]],
                 tgt_xyz: dict[int, dict[str, tuple[float, float, float]]]) -> PairSelection:
    """Corresponding nucleotides: same consensus position, match state in both, coordinates in both, reliable PP in both, a common atom."""
    rm, tm = _match_rows(ref_map), _match_rows(tgt_map)
    shared = sorted(set(rm) & set(tm))
    excluded = {"no_coords_reference": 0, "no_coords_target": 0, "low_pp": 0, "missing_atom": 0}
    with_coords = []
    for c in shared:
        r, t = rm[c], tm[c]
        if not r["has_coordinates"]:
            excluded["no_coords_reference"] += 1; continue
        if not t["has_coordinates"]:
            excluded["no_coords_target"] += 1; continue
        with_coords.append(c)

    def build(pp_ok: tuple[str, ...]) -> tuple[list[tuple[int, int, int, str]], dict[str, int], int, int]:
        pairs, atoms_used, low_pp, missing = [], {}, 0, 0
        for c in with_coords:
            r, t = rm[c], tm[c]
            if r["posterior_prob_char"] not in pp_ok or t["posterior_prob_char"] not in pp_ok:
                low_pp += 1; continue
            ra, ta = ref_xyz.get(r["polymer_seq_position"], {}), tgt_xyz.get(t["polymer_seq_position"], {})
            atom = next((a for a in ATOM_PRIORITY if a in ra and a in ta), None)
            if atom is None:
                missing += 1; continue
            pairs.append((c, r["polymer_seq_position"], t["polymer_seq_position"], atom))
            atoms_used[atom] = atoms_used.get(atom, 0) + 1
        return pairs, atoms_used, low_pp, missing

    pairs, atoms_used, low_pp, missing = build(PP_STRICT)
    pp_filter = "*"
    if len(pairs) < RELAX_BELOW:
        pairs, atoms_used, low_pp, missing = build(PP_RELAXED)
        pp_filter = "9"
    excluded["low_pp"], excluded["missing_atom"] = low_pp, missing
    return PairSelection(pairs, pp_filter, len(shared), len(with_coords), excluded, atoms_used)


def fit_quality(rmsd: Optional[float], n_pairs: int) -> str:
    if n_pairs < MIN_PAIRS or rmsd is None:
        return "insufficient"
    if rmsd <= 2.0 and n_pairs >= 20:
        return "good"
    if rmsd <= 4.0:
        return "moderate"
    return "poor"


def orient_instance(ref_map: dict[str, Any], tgt_map: dict[str, Any], ref_xyz: dict[int, dict[str, tuple[float, float, float]]],
                    tgt_xyz: dict[int, dict[str, tuple[float, float, float]]], ref_n_match_coords: int, tgt_n_match_coords: int) -> dict[str, Any]:
    sel = select_pairs(ref_map, tgt_map, ref_xyz, tgt_xyz)
    rec: dict[str, Any] = {
        "target": {k: tgt_map[k] for k in ("pdb_id", "entity_id", "label_asym_id", "auth_asym_id")},
        "n_pairs": len(sel.pairs), "n_shared_consensus": sel.n_shared_consensus, "n_shared_consensus_with_coords": sel.n_shared_with_coords,
        "reference_coverage": round(len(sel.pairs) / ref_n_match_coords, 4) if ref_n_match_coords else None,
        "target_coverage": round(len(sel.pairs) / tgt_n_match_coords, 4) if tgt_n_match_coords else None,
        "pp_filter_used": sel.pp_filter, "atoms_used": sel.atoms_used, "excluded": sel.excluded,
        "paired_consensus_positions": [p[0] for p in sel.pairs],
    }
    if len(sel.pairs) < MIN_PAIRS:
        rec.update({"rotation": None, "translation": None, "rmsd": None, "determinant": None, "fit_quality": "insufficient", "is_identity": False, "check": None})
        return rec
    t_xyz = np.array([tgt_xyz[p[2]][p[3]] for p in sel.pairs])
    r_xyz = np.array([ref_xyz[p[1]][p[3]] for p in sel.pairs])
    k = kabsch(t_xyz, r_xyz)
    is_identity = bool(np.allclose(np.array(k.rotation), np.eye(3), atol=1e-6) and np.allclose(np.array(k.translation), 0, atol=1e-6))
    c0, _, tgt_seq, atom = sel.pairs[0]
    orig = tgt_xyz[tgt_seq][atom]
    moved = (np.array(k.rotation) @ np.array(orig) + np.array(k.translation)).tolist()
    rec.update({"rotation": k.rotation, "translation": k.translation, "rmsd": round(k.rmsd, 4), "determinant": round(k.determinant, 6),
                "fit_quality": fit_quality(k.rmsd, k.n_pairs), "is_identity": is_identity,
                "check": {"consensus_pos": c0, "label_seq_id": tgt_seq, "atom": atom, "original": list(orig), "transformed": moved}})
    return rec


def build_orientation(fam_dir: Path, cif_for: Any, out_path: Optional[Path] = None) -> dict[str, Any]:
    """Compute transforms for every instance of a built family onto its display representative. `cif_for(instance row) -> Path`
    (whole-entry mmCIF, or the instance's own chain file for selection builds)."""
    from .mmcif import atom_site_coords
    family = json.load(open(fam_dir / "family.json"))
    instances = json.load(open(fam_dir / "instances.json"))
    ref = family.get("display_representative") or family.get("representative")
    if not ref:
        raise RuntimeError("no representative to use as reference")
    ref_map = json.load(open(fam_dir / ref["mapping_file"]))
    ref_xyz = atom_site_coords(cif_for(ref), ref["label_asym_id"], ref["entity_id"])
    transforms: dict[str, Any] = {}
    xyz_cache: dict[tuple[str, str, int], dict[int, dict[str, tuple[float, float, float]]]] = {}
    for inst in instances:
        key = (inst["pdb_id"], inst["label_asym_id"], inst["entity_id"])
        if key not in xyz_cache:
            xyz_cache[key] = atom_site_coords(cif_for(inst), inst["label_asym_id"], inst["entity_id"])
        tgt_map = json.load(open(fam_dir / inst["mapping_file"]))
        transforms[inst["seq_name"]] = orient_instance(ref_map, tgt_map, ref_xyz, xyz_cache[key], ref["n_match_with_coordinates"], inst["n_match_with_coordinates"])
    out = {
        "schema_version": SCHEMA_VERSION, "rfam_acc": family["rfam_acc"],
        "reference": {k: ref[k] for k in ("seq_name", "pdb_id", "entity_id", "label_asym_id", "auth_asym_id")} | {"selection_rule": ref.get("selection_rule", "coverage representative")},
        "method": {"algorithm": "Kabsch (SVD) rigid-body superposition, rotation + translation, no scaling, reflections rejected (det R = +1)",
                   "atom_priority": list(ATOM_PRIORITY), "atom_fallback": "per pair: first atom present in both residues in the order C4', P, C1'",
                   "coordinates": "first model (pdbx_PDB_model_num) and first alternate location of the PDBe updated mmCIF",
                   "pair_selection": "same Rfam consensus position, match state in both chains, coordinates in both, PP '*' (>= 0.95) in both; relaxed to PP >= '9' (>= 0.85) when fewer than %d pairs; transform omitted below %d pairs" % (RELAX_BELOW, MIN_PAIRS),
                   "quality_bands": "good: rmsd <= 2.0 A and >= 20 pairs; moderate: rmsd <= 4.0 A; poor: otherwise; insufficient: < %d pairs" % MIN_PAIRS,
                   "convention": "x_ref = R @ x_target + t (column vectors); rotation row-major 3x3; Mol*: Mat4.ofRows([[R00,R01,R02,t0],[R10,R11,R12,t1],[R20,R21,R22,t2],[0,0,0,1]]) with TransformStructureConformation (transpose false); the target structure is moved into the reference frame"},
        "transforms": transforms,
        "provenance": {"pipeline_version": __version__, "generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                       "mapping_provenance": ref_map.get("provenance"), "numpy": np.__version__},
    }
    if out_path:
        with open(out_path, "w") as f:
            json.dump(out, f, indent=1)
    return out
