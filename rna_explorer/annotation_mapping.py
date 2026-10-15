"""Resolve external nucleotide identifiers (BGSU unit ids, PDBe polymer sequence indices) against the validated family mappings.

Canonical residue identifier: (pdb_id, label_asym_id, label_seq_id). Every resolution returns a ResidueRef carrying the
Rfam consensus position (when the residue is a match state), the alignment state, PP confidence and coordinate availability.
Nothing is forced onto a consensus position: insertions, local ends and residues outside the mapped region keep
``rfam_consensus_pos = None`` and a ``mapping_status`` explaining why.
"""
from __future__ import annotations
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Optional
from .mapping import pp_category, pp_value

MAPPING_STATUSES = ("mapped", "mapped_no_coordinates", "insert", "local_end", "outside_mapped_region", "unresolved_identifier")


@dataclass
class UnitId:
    """BGSU unit id ``PDB|model|chain|comp|number[|atom|alt|inscode|sym]``."""
    raw: str
    pdb_id: str
    model: int
    auth_asym_id: str
    comp_id: str
    author_residue_number: str
    insertion_code: str

    @classmethod
    def parse(cls, unit: str) -> "UnitId":
        f = unit.strip().split("|")
        if len(f) < 5:
            raise ValueError(f"malformed unit id: {unit!r}")
        return cls(raw=unit.strip(), pdb_id=f[0].lower(), model=int(f[1]) if f[1].isdigit() else 1, auth_asym_id=f[2], comp_id=f[3],
                   author_residue_number=f[4], insertion_code=(f[7] if len(f) > 7 else "") or "")


@dataclass
class ResidueRef:
    seq_name: Optional[str]
    pdb_id: str
    entity_id: Optional[int]
    label_asym_id: Optional[str]
    auth_asym_id: Optional[str]
    label_seq_id: Optional[int]
    author_residue_number: Optional[int]
    insertion_code: str
    residue_name: Optional[str]
    alignment_state: Optional[str]
    rfam_consensus_pos: Optional[int]
    pp_char: str
    pp_value: Optional[float]
    pp_category: Optional[str]
    has_coordinates: Optional[bool]
    mapping_status: str
    source_identifier: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _status(res: dict[str, Any]) -> str:
    st = res["alignment_state"]
    if st == "match":
        return "mapped" if res["has_coordinates"] else "mapped_no_coordinates"
    if st in ("insert", "local_end"):
        return st
    return "outside_mapped_region"


class InstanceIndex:
    """Instances of one family with their residue mappings, indexed for identifier resolution."""

    def __init__(self, fam_dir: Path) -> None:
        self.fam_dir = Path(fam_dir)
        self.family: dict[str, Any] = json.load(open(self.fam_dir / "family.json"))
        self.instances: list[dict[str, Any]] = json.load(open(self.fam_dir / "instances.json"))
        self.clen: int = int(self.family["clen"])
        self.by_seq: dict[str, dict[str, Any]] = {r["seq_name"]: r for r in self.instances}
        self.by_auth: dict[tuple[str, str], dict[str, Any]] = {(r["pdb_id"], r["auth_asym_id"]): r for r in self.instances}
        self.by_label: dict[tuple[str, str], dict[str, Any]] = {(r["pdb_id"], r["label_asym_id"]): r for r in self.instances}
        self.by_entity: dict[tuple[str, int], list[dict[str, Any]]] = {}
        self.by_pdb: dict[str, list[dict[str, Any]]] = {}
        for r in self.instances:
            self.by_entity.setdefault((r["pdb_id"], int(r["entity_id"])), []).append(r)
            self.by_pdb.setdefault(r["pdb_id"], []).append(r)
        self._maps: dict[str, dict[str, Any]] = {}
        self._auth_idx: dict[str, dict[tuple[str, str], dict[str, Any]]] = {}
        self._seq_idx: dict[str, dict[int, dict[str, Any]]] = {}
        self._cons_idx: dict[str, dict[int, dict[str, Any]]] = {}

    # ---- mapping access -------------------------------------------------
    def mapping(self, seq_name: str) -> dict[str, Any]:
        if seq_name not in self._maps:
            m = json.load(open(self.fam_dir / self.by_seq[seq_name]["mapping_file"]))
            self._maps[seq_name] = m
            self._auth_idx[seq_name] = {(str(r["author_residue_number"]), r["author_insertion_code"] or ""): r for r in m["residues"]}
            self._seq_idx[seq_name] = {int(r["polymer_seq_position"]): r for r in m["residues"]}
            self._cons_idx[seq_name] = {int(r["rfam_consensus_pos"]): r for r in m["residues"] if r["alignment_state"] == "match"}
        return self._maps[seq_name]

    def residues(self, seq_name: str) -> list[dict[str, Any]]:
        return self.mapping(seq_name)["residues"]

    def match_residue(self, seq_name: str, consensus_pos: int) -> Optional[dict[str, Any]]:
        self.mapping(seq_name)
        return self._cons_idx[seq_name].get(consensus_pos)

    def observable(self, seq_name: str, consensus_pos: int) -> bool:
        """True when the instance has an experimentally observed (modelled) nucleotide aligned to this consensus position."""
        r = self.match_residue(seq_name, consensus_pos)
        return bool(r and r["has_coordinates"])

    def polymer_length(self, seq_name: str) -> int:
        return len(self.residues(seq_name))

    def polymer_residue_names(self, seq_name: str) -> set[str]:
        return {r["residue_name"].upper() for r in self.residues(seq_name)}

    # ---- resolution -----------------------------------------------------
    def _ref(self, inst: Optional[dict[str, Any]], res: Optional[dict[str, Any]], source_identifier: str, pdb_id: str,
             auth_asym_id: Optional[str], ins: str = "") -> ResidueRef:
        if inst is None or res is None:
            return ResidueRef(seq_name=inst["seq_name"] if inst else None, pdb_id=pdb_id, entity_id=inst["entity_id"] if inst else None,
                              label_asym_id=inst["label_asym_id"] if inst else None, auth_asym_id=auth_asym_id, label_seq_id=None,
                              author_residue_number=None, insertion_code=ins, residue_name=None, alignment_state=None, rfam_consensus_pos=None,
                              pp_char="", pp_value=None, pp_category=None, has_coordinates=None, mapping_status="unresolved_identifier",
                              source_identifier=source_identifier)
        ch = res["posterior_prob_char"] or ""
        return ResidueRef(seq_name=inst["seq_name"], pdb_id=inst["pdb_id"], entity_id=int(inst["entity_id"]), label_asym_id=inst["label_asym_id"],
                          auth_asym_id=inst["auth_asym_id"], label_seq_id=int(res["polymer_seq_position"]), author_residue_number=res["author_residue_number"],
                          insertion_code=res["author_insertion_code"] or "", residue_name=res["residue_name"], alignment_state=res["alignment_state"],
                          rfam_consensus_pos=res["rfam_consensus_pos"], pp_char=ch, pp_value=pp_value(ch) if ch else None,
                          pp_category=pp_category(ch) if ch else None, has_coordinates=bool(res["has_coordinates"]), mapping_status=_status(res),
                          source_identifier=source_identifier)

    def resolve_unit(self, unit: str | UnitId) -> ResidueRef:
        """BGSU unit id (author chain + author number + insertion code) -> canonical residue."""
        u = UnitId.parse(unit) if isinstance(unit, str) else unit
        inst = self.by_auth.get((u.pdb_id, u.auth_asym_id))
        res = None
        if inst is not None:
            self.mapping(inst["seq_name"])
            res = self._auth_idx[inst["seq_name"]].get((u.author_residue_number, u.insertion_code))
        return self._ref(inst, res, u.raw, u.pdb_id, u.auth_asym_id, u.insertion_code)

    def resolve_seq_index(self, seq_name: str, label_seq_id: int, source_identifier: str = "") -> ResidueRef:
        """PDBe polymer sequence index (label_seq_id) on a known instance -> canonical residue."""
        inst = self.by_seq[seq_name]
        self.mapping(seq_name)
        res = self._seq_idx[seq_name].get(int(label_seq_id))
        return self._ref(inst, res, source_identifier or f"{inst['pdb_id']}|{inst['label_asym_id']}|{label_seq_id}", inst["pdb_id"], inst["auth_asym_id"])
