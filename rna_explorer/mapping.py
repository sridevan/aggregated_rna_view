"""Consensus residue mapping: aligned-sequence states (match / insert / local_end) and per-chain mapping rows.
Logic validated in poc_rf00177 (03_parse_stockholm.py) and against BGSU RNA 3D Hub."""
from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict, field
from typing import Any, Optional
from .sequences import ChainInstance
from .stockholm import ColumnMap, GAPS

# Infernal PP annotation bins: '*' >= 0.95, '9' = [0.85,0.95), ..., '0' = [0,0.05). Value = bin midpoint.
PP_BINS: dict[str, float] = {"*": 0.975, "0": 0.025, **{str(i): i / 10 for i in range(1, 10)}}
PP_THRESHOLDS = {"high": 0.95, "medium": 0.80}   # provisional display thresholds on bin midpoints


def pp_value(ch: str) -> Optional[float]:
    return PP_BINS.get(ch)


def pp_category(ch: str) -> Optional[str]:
    v = pp_value(ch)
    if v is None:
        return None
    return "high" if v >= PP_THRESHOLDS["high"] else "medium" if v >= PP_THRESHOLDS["medium"] else "low"


@dataclass
class AlignedResidue:
    index: int                     # 1-based position within the aligned (mapped-region) sequence
    aligned_char: str
    state: str                     # match | insert | local_end
    consensus_pos: Optional[int]
    anchor: Optional[int]          # consensus position preceding the insert/local-end column (0 before the first)
    offset: Optional[int]          # 1..n within the run after that anchor
    pp_char: str
    column: int                    # 1-based alignment column


@dataclass
class AlignedConsensus:
    consensus_pos: int
    state: str                     # match | deletion (inside aligned span) | uncovered (outside span)
    index: Optional[int]
    column: int
    pp_char: str
    rf: str
    ss_cons: str


@dataclass
class AlignedSequence:
    name: str
    length: int
    residues: list[AlignedResidue]
    consensus: list[AlignedConsensus]
    cons_min: int
    cons_max: int

    @property
    def match_by_consensus(self) -> dict[int, AlignedResidue]:
        return {r.consensus_pos: r for r in self.residues if r.state == "match" and r.consensus_pos is not None}


def map_aligned_sequence(name: str, aln_row: str, pp_row: str, cm: ColumnMap) -> AlignedSequence:
    """Assign a state to every residue of an aligned row and a state to every consensus position."""
    residues: list[AlignedResidue] = []
    cons: list[AlignedConsensus] = []
    idx = 0
    offset = 0
    for col, a in enumerate(aln_row):
        ctype = cm.col_type[col]
        pp = pp_row[col] if pp_row else ""
        if a in GAPS:
            if ctype == "match":
                offset = 0
                cons.append(AlignedConsensus(cm.col_cons[col] or 0, "deletion", None, col + 1, "", cm.rf[col], cm.ss_cons[col]))
            continue
        idx += 1
        if ctype == "match":
            offset = 0
            residues.append(AlignedResidue(idx, a, "match", cm.col_cons[col], None, None, pp, col + 1))
            cons.append(AlignedConsensus(cm.col_cons[col] or 0, "match", idx, col + 1, pp, cm.rf[col], cm.ss_cons[col]))
        else:
            offset += 1
            residues.append(AlignedResidue(idx, a, ctype, None, cm.prev_cons[col], offset, pp, col + 1))
    matched = [r.consensus_pos for r in residues if r.state == "match" and r.consensus_pos]
    if not matched:
        raise ValueError(f"{name}: no residue aligned to a consensus column")
    lo, hi = min(matched), max(matched)
    for c in cons:
        if c.state == "deletion" and not (lo <= c.consensus_pos <= hi):
            c.state = "uncovered"
    cons.sort(key=lambda c: c.consensus_pos)
    return AlignedSequence(name, idx, residues, cons, lo, hi)


@dataclass
class ResidueRow:
    pdb_id: str
    entity_id: int
    label_asym_id: str
    auth_asym_id: str
    polymer_seq_position: int
    author_residue_number: Any
    author_insertion_code: str
    residue_name: str
    nucleotide_parent: str
    is_modified: int
    alignment_state: str
    rfam_consensus_pos: Optional[int]
    insertion_anchor: Optional[int]
    insertion_offset: Optional[int]
    aligned_char: str
    aln_column: Optional[int]
    rf: str
    ss_cons: str
    posterior_prob: str
    posterior_prob_value: Optional[float]
    has_coordinates: int
    observed_ratio_api: Any
    in_pdbe_mapped_region: int
    rfam_acc: str
    seq_name: str
    sequence_key: str
    organism_api: Optional[str]
    organism_mmcif: Optional[str]
    tax_id_mmcif: Optional[str]


@dataclass
class ConsensusRow:
    pdb_id: str
    auth_asym_id: str
    seq_name: str
    rfam_consensus_pos: int
    aln_column: Optional[int]
    state: str
    polymer_seq_position: Optional[int]
    nucleotide_parent: str
    residue_name: str
    has_coordinates: Optional[int]
    posterior_prob: str
    ss_cons: str
    rf: str


@dataclass
class ChainMapping:
    info: ChainInstance
    aligned: AlignedSequence
    rows: list[ResidueRow]
    cons_rows: list[ConsensusRow]
    summary: dict[str, Any] = field(default_factory=dict)

    @property
    def key(self) -> str:
        return self.info.seq_name


def build_chain_mapping(info: ChainInstance, aligned: AlignedSequence, cm: ColumnMap, rfam_acc: str) -> ChainMapping:
    """Join the unique-sequence alignment with chain-specific residue identifiers and coordinate availability."""
    if aligned.length != info.map_end - info.map_start + 1:
        raise ValueError(f"{info.seq_name}: aligned length {aligned.length} != mapped region length")
    res_by_num = {int(r["residue_number"]): r for r in info.residues}
    common = dict(pdb_id=info.pdb_id, entity_id=info.entity_id, label_asym_id=info.label_asym_id, auth_asym_id=info.auth_asym_id, rfam_acc=rfam_acc,
                  seq_name=info.seq_name, sequence_key=info.sequence_key, organism_api=info.api_organism, organism_mmcif=info.cif_organism, tax_id_mmcif=info.cif_tax_id)

    def base(p: int) -> dict[str, Any]:
        r = res_by_num[p]
        return dict(common, polymer_seq_position=p, author_residue_number=r["author_residue_number"], author_insertion_code=r.get("author_insertion_code") or "",
                    residue_name=r["residue_name"], nucleotide_parent=info.sequence[p - 1], is_modified=int(r["residue_name"] not in ("A", "C", "G", "U")),
                    has_coordinates=int(p in info.observed), observed_ratio_api=r.get("observed_ratio"), in_pdbe_mapped_region=int(info.map_start <= p <= info.map_end))

    rows: list[ResidueRow] = []
    for ar in aligned.residues:
        p = info.map_start - 1 + ar.index
        rows.append(ResidueRow(**base(p), alignment_state=ar.state, rfam_consensus_pos=ar.consensus_pos, insertion_anchor=ar.anchor, insertion_offset=ar.offset,
                               aligned_char=ar.aligned_char, aln_column=ar.column, rf=cm.rf[ar.column - 1], ss_cons=cm.ss_cons[ar.column - 1],
                               posterior_prob=ar.pp_char, posterior_prob_value=pp_value(ar.pp_char)))
    for p in range(1, len(info.sequence) + 1):
        if p < info.map_start or p > info.map_end:
            rows.append(ResidueRow(**base(p), alignment_state="not_submitted_flank", rfam_consensus_pos=None, insertion_anchor=None, insertion_offset=None,
                                   aligned_char="", aln_column=None, rf="", ss_cons="", posterior_prob="", posterior_prob_value=None))
    rows.sort(key=lambda r: r.polymer_seq_position)
    by_pos = {r.polymer_seq_position: r for r in rows}
    cons_rows: list[ConsensusRow] = []
    for c in aligned.consensus:
        if c.state == "match" and c.index is not None:
            p = info.map_start - 1 + c.index
            r = by_pos[p]
            cons_rows.append(ConsensusRow(info.pdb_id, info.auth_asym_id, info.seq_name, c.consensus_pos, c.column, "match", p, r.nucleotide_parent, r.residue_name, r.has_coordinates, c.pp_char, c.ss_cons, c.rf))
        else:
            cons_rows.append(ConsensusRow(info.pdb_id, info.auth_asym_id, info.seq_name, c.consensus_pos, c.column, c.state, None, "", "", None, "", c.ss_cons, c.rf))
    mapped = [r for r in rows if r.alignment_state != "not_submitted_flank"]
    ins_after: Counter[int] = Counter(r.insertion_anchor for r in mapped if r.alignment_state in ("insert", "local_end") and r.insertion_anchor is not None)
    ppv = [r.posterior_prob for r in mapped if r.posterior_prob]
    n_match = sum(r.alignment_state == "match" for r in mapped)
    summary: dict[str, Any] = {
        "mapped_region_len": len(mapped), "n_match": n_match, "n_insert": sum(r.alignment_state == "insert" for r in mapped),
        "n_local_end": sum(r.alignment_state == "local_end" for r in mapped), "cons_span": [aligned.cons_min, aligned.cons_max],
        "n_deleted_in_span": sum(c.state == "deletion" for c in cons_rows), "n_uncovered_outside_span": sum(c.state == "uncovered" for c in cons_rows),
        "fraction_CLEN_matched": round(n_match / cm.clen, 4), "n_insertion_sites": len(ins_after),
        "longest_insertion": max(ins_after.values()) if ins_after else 0, "longest_insertion_anchor": max(ins_after, key=ins_after.get) if ins_after else None,
        "pp_star_fraction": round(sum(p == "*" for p in ppv) / len(ppv), 4) if ppv else None, "pp_below_0.7_count": sum(p not in ("*", "9", "8", "7") for p in ppv),
        "n_mapped_with_coordinates": sum(r.has_coordinates for r in mapped), "n_match_with_coordinates": sum(r.has_coordinates for r in mapped if r.alignment_state == "match"),
        "n_modified_in_mapped": sum(r.is_modified for r in mapped), "pp_category_counts": dict(Counter(pp_category(p) for p in ppv))}
    return ChainMapping(info, aligned, rows, cons_rows, summary)


def rows_to_dicts(rows: list[Any]) -> list[dict[str, Any]]:
    return [asdict(r) for r in rows]
