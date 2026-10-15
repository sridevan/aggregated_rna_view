"""Built family JSON files validate against the schemas and are internally consistent (runs only when a build exists)."""
import json
from pathlib import Path
import pytest
from rna_explorer.export import validate_json

ROOT = Path(__file__).resolve().parent.parent
FAM = ROOT / "data/families/RF00167"
pytestmark = pytest.mark.skipif(not (FAM / "family.json").exists(), reason="RF00167 build not present")


def test_schemas_and_consistency() -> None:
    fam = json.load(open(FAM / "family.json")); inst = json.load(open(FAM / "instances.json")); st = json.load(open(FAM / "consensus_stats.json"))
    validate_json(fam, "family"); validate_json(inst, "instances"); validate_json(st, "consensus_stats")
    assert fam["counts"]["instances"] == len(inst) == st["structure_level"]["n_instances"]
    assert len(st["sequence_level"]["coverage"]) == fam["clen"] == len(fam["consensus_sequence"])
    for row in inst[:10]:
        m = json.load(open(FAM / row["mapping_file"]))
        validate_json(m, "mapping")
        assert m["pdb_id"] == row["pdb_id"] and m["auth_asym_id"] == row["auth_asym_id"]
        matches = [r for r in m["residues"] if r["alignment_state"] == "match"]
        assert len(matches) == row["n_match"]
        assert all(r["rfam_consensus_pos"] is None for r in m["residues"] if r["alignment_state"] in ("insert", "local_end"))
        assert sorted(c["rfam_consensus_pos"] for c in m["consensus_columns"]) == list(range(1, fam["clen"] + 1))
