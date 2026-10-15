"""Frontend export: static JSON (validated against JSON Schemas) plus CSV twins for regression."""
from __future__ import annotations
import csv, json
from dataclasses import asdict
from pathlib import Path
from typing import Any
from .mapping import ChainMapping

SCHEMA_DIR = Path(__file__).parent / "schemas"


def _dump(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)


def validate_json(obj: Any, schema_name: str) -> None:
    try:
        import jsonschema
    except ImportError:   # pragma: no cover
        return
    schema = json.load(open(SCHEMA_DIR / f"{schema_name}.schema.json"))
    jsonschema.validate(obj, schema)


def mapping_json(m: ChainMapping, provenance: dict[str, Any], alignment: dict[str, Any] | None) -> dict[str, Any]:
    info = m.info
    return {"pdb_id": info.pdb_id, "entity_id": info.entity_id, "label_asym_id": info.label_asym_id, "auth_asym_id": info.auth_asym_id, "rfam_acc": provenance["rfam_acc"],
            "sequence_key": info.sequence_key, "mapped_region": {"start": info.map_start, "end": info.map_end}, "alignment": alignment, "summary": m.summary,
            "residues": [{"polymer_seq_position": r.polymer_seq_position, "author_residue_number": r.author_residue_number, "author_insertion_code": r.author_insertion_code,
                          "residue_name": r.residue_name, "nucleotide_parent": r.nucleotide_parent, "is_modified": bool(r.is_modified), "alignment_state": r.alignment_state,
                          "rfam_consensus_pos": r.rfam_consensus_pos, "insertion_anchor": r.insertion_anchor, "insertion_offset": r.insertion_offset,
                          "posterior_prob_char": r.posterior_prob, "posterior_prob": r.posterior_prob_value, "has_coordinates": bool(r.has_coordinates)} for r in m.rows],
            "consensus_columns": [{"rfam_consensus_pos": c.rfam_consensus_pos, "state": c.state, "polymer_seq_position": c.polymer_seq_position} for c in m.cons_rows],
            "provenance": provenance}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})


def write_family(out_dir: Path, family: dict[str, Any], instances: list[dict[str, Any]], chains: list[ChainMapping], stats: dict[str, Any],
                 validation: list[dict[str, Any]], provenance: dict[str, Any], alignments: dict[str, dict[str, Any]]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    validate_json(family, "family"); validate_json(instances, "instances"); validate_json(stats, "consensus_stats")
    _dump(out_dir / "family.json", family)
    _dump(out_dir / "instances.json", instances)
    _dump(out_dir / "consensus_stats.json", stats)
    _dump(out_dir / "validation.json", validation)
    res_rows: list[dict[str, Any]] = []
    cons_rows: list[dict[str, Any]] = []
    for m in chains:
        mj = mapping_json(m, provenance, alignments.get(m.info.sequence_key))
        validate_json(mj, "mapping")
        _dump(out_dir / "mappings" / f"{m.info.seq_name}.json", mj)
        res_rows.extend(asdict(r) for r in m.rows)
        cons_rows.extend(asdict(c) for c in m.cons_rows)
    write_csv(out_dir / "residue_mapping.csv", res_rows)
    write_csv(out_dir / "consensus_columns.csv", cons_rows)
    return {"family": out_dir / "family.json", "instances": out_dir / "instances.json", "stats": out_dir / "consensus_stats.json"}
