"""End-to-end check of a built RF00167 data set against the POC (runs only when data/families/RF00167 exists)."""
import csv
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parent.parent
BUILT = ROOT / "data/families/RF00167/residue_mapping.csv"
POC = ROOT / "poc_rf00167/results/mapped_residue_mapping.csv"
pytestmark = pytest.mark.skipif(not (BUILT.exists() and POC.exists()), reason="RF00167 build or POC not present")


def test_built_rf00167_matches_poc() -> None:
    built = {(r["pdb_id"], r["auth_asym_id"], int(r["polymer_seq_position"])): r for r in csv.DictReader(open(BUILT))}
    n = 0
    for p in csv.DictReader(open(POC)):
        b = built[(p["pdb_id"], p["auth_chain_id"], int(p["polymer_seq_position"]))]
        assert b["alignment_state"] == p["status"]
        assert b["rfam_consensus_pos"] == p["rfam_consensus_pos"]
        assert b["posterior_prob"] == p["posterior_prob"]
        assert b["author_residue_number"] == p["author_residue_number"] and b["residue_name"] == p["residue_name"]
        n += 1
    assert n == 434
