"""Chain instances, PDBe-mapped region extraction, parent-base normalisation and deduplication."""
from __future__ import annotations
import hashlib
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional


@dataclass
class ChainInstance:
    """One PDBe nucleic-mapping segment (one chain) with the entity sequence and per-residue listing."""
    pdb_id: str
    entity_id: int
    label_asym_id: str
    auth_asym_id: str
    map_start: int                      # polymer (label_seq_id) position, 1-based inclusive
    map_end: int
    map_auth_start: Any
    map_auth_end: Any
    sequence: str                       # PDBe one-letter sequence, modified residues as parent base
    pdb_sequence: str
    residues: list[dict[str, Any]]      # PDBe residue listing rows sorted by residue_number
    modified_positions: dict[str, Any] = field(default_factory=dict)
    observed: set[int] = field(default_factory=set)   # label_seq_id with coordinates (atom_site)
    molecule_name: str = ""
    api_organism: Optional[str] = None
    api_tax_id: Optional[int] = None
    cif_organism: Optional[str] = None
    cif_tax_id: Optional[str] = None
    cif_src_category: Optional[str] = None

    @property
    def mapped_sequence(self) -> str:
        return self.sequence[self.map_start - 1:self.map_end]

    @property
    def sequence_key(self) -> str:
        return sequence_key(self.mapped_sequence)

    @property
    def seq_name(self) -> str:
        return f"{self.pdb_id}_{self.entity_id}_{self.auth_asym_id}_{self.label_asym_id}"

    @property
    def modified_in_mapped_region(self) -> dict[int, Any]:
        return {int(k): v for k, v in self.modified_positions.items() if self.map_start <= int(k) <= self.map_end}

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["observed"] = sorted(self.observed)
        return d


def sequence_key(mapped_sequence: str) -> str:
    return hashlib.md5(mapped_sequence.encode()).hexdigest()


def chain_instance_from_pdbe(pdb_id: str, segment: dict[str, Any], entity: dict[str, Any], residue_listing: dict[str, Any],
                             observed: set[int], source: Optional[dict[str, Optional[str]]] = None) -> ChainInstance:
    """Build a ChainInstance from a nucleic_mappings segment + molecules entity + residue_listing (+ atom_site observed set)."""
    ent = int(segment["entity_id"])
    chain, asym = segment["chain_id"], segment["struct_asym_id"]
    listing = [c for m in residue_listing.get("molecules", []) if int(m["entity_id"]) == ent
               for c in m["chains"] if c["chain_id"] == chain and c["struct_asym_id"] == asym]
    if len(listing) != 1:
        raise ValueError(f"{pdb_id} entity {ent} chain {chain}/{asym}: {len(listing)} residue listings")
    residues = sorted(listing[0]["residues"], key=lambda r: r["residue_number"])
    seq = entity["sequence"]
    if not (len(seq) == len(residues) == int(entity["length"])):
        raise ValueError(f"{pdb_id}/{chain}: sequence length {len(seq)} != residues {len(residues)} != entity length {entity['length']}")
    src = (entity.get("source") or [{}])[0]
    return ChainInstance(pdb_id=pdb_id.lower(), entity_id=ent, label_asym_id=asym, auth_asym_id=chain,
                         map_start=int(segment["start"]["residue_number"]), map_end=int(segment["end"]["residue_number"]),
                         map_auth_start=segment["start"].get("author_residue_number"), map_auth_end=segment["end"].get("author_residue_number"),
                         sequence=seq, pdb_sequence=entity.get("pdb_sequence", seq), residues=residues,
                         modified_positions=entity.get("pdb_sequence_indices_with_multiple_residues") or {}, observed=set(observed),
                         molecule_name="; ".join(entity.get("molecule_name") or []), api_organism=src.get("organism_scientific_name"), api_tax_id=src.get("tax_id"),
                         cif_organism=(source or {}).get("cif_organism"), cif_tax_id=(source or {}).get("cif_tax_id"), cif_src_category=(source or {}).get("cif_src_category"))


@dataclass
class UniqueSequence:
    key: str
    mapped_sequence: str
    members: list[str]      # ChainInstance.seq_name of every chain with this mapped sequence

    @property
    def fasta_name(self) -> str:
        return f"seq_{self.key}"


def unique_sequences(instances: list[ChainInstance]) -> dict[str, UniqueSequence]:
    out: dict[str, UniqueSequence] = {}
    for inst in instances:
        k = inst.sequence_key
        if k not in out:
            out[k] = UniqueSequence(k, inst.mapped_sequence, [])
        out[k].members.append(inst.seq_name)
    return out


def non_acgu_letters(seq: str) -> list[str]:
    return sorted(set(seq) - set("ACGU"))


def write_fasta(uniques: dict[str, UniqueSequence], path: Path | str) -> None:
    with open(path, "w") as f:
        for u in uniques.values():
            f.write(f">{u.fasta_name}\n{u.mapped_sequence}\n")
