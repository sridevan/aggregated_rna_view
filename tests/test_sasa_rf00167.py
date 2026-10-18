"""Regression of the generated RF00167 SASA files (runs only when present)."""
import json, math
from pathlib import Path
import pytest
from rna_explorer.export import validate_json

FAM = Path(__file__).resolve().parent.parent / "data/families/RF00167"
pytestmark = pytest.mark.skipif(not (FAM / "sasa_tracks.json").exists(), reason="RF00167 SASA not built")


def test_sasa_dataset_invariants():
    d = json.load(open(FAM / "sasa.json")); t = json.load(open(FAM / "sasa_tracks.json"))
    validate_json(d, "sasa"); validate_json(t, "sasa_tracks")
    assert d["counts"]["n_instances_measured"] == 84 and d["failures"] == []
    assert d["method"]["probe_radius_A"] == 1.4 and d["method"]["algorithm"].startswith("Shrake-Rupley") and d["method"]["software"]["freesasa"]
    assert d["counts"]["nucleotide_status_counts"] == {"complete": 6133, "incomplete": 46, "no_coordinates": 422}
    for m in d["instances"]:
        vals = [n["sasa_A2"] for n in m["nucleotides"] if n["sasa_A2"] is not None]
        assert all(math.isfinite(v) and v >= 0 for v in vals)
        assert abs(m["chain_total_sasa_A2"] - m["sum_nucleotide_sasa_A2"]) < 0.05           # per-nucleotide values sum to the chain total
        assert all(n["sasa_A2"] is None for n in m["nucleotides"] if n["status"] == "no_coordinates")   # missing residues never become zero
        assert all(n["rfam_consensus_pos"] is None for n in m["nucleotides"] if n["alignment_state"] in ("insert", "local_end"))
        assert all(n["missing_atoms"] for n in m["nucleotides"] if n["status"] == "incomplete")
    p53 = t["positions"][52]
    assert p53["n_contributing_instances"] == 80 and p53["n_mapped_instances_with_coordinates"] == 80 and p53["median_sasa_A2"] == 134.85
    assert t["summary"]["n_positions_with_value"] == 102 and t["colour_domain_suggestion_A2"]["max"] == 350.0
    for p in t["positions"]:
        assert p["n_contributing_instances"] + p["n_excluded_incomplete"] + p["n_excluded_unsupported"] == p["n_mapped_instances_with_coordinates"]
        assert p["n_contributing_instances"] <= 84
        if p["median_sasa_A2"] is not None:
            assert p["min_sasa_A2"] <= p["q1_sasa_A2"] <= p["median_sasa_A2"] <= p["q3_sasa_A2"] <= p["max_sasa_A2"]
