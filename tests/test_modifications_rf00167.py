"""Regression of the generated RF00167 nonstandard-nucleotide files against the validated investigation (experiments/modifications_RF00167)."""
import json
from collections import Counter
from pathlib import Path
import pytest
from rna_explorer.export import validate_json

ROOT = Path(__file__).resolve().parent.parent
FAM = ROOT / "data/families/RF00167"
pytestmark = pytest.mark.skipif(not (FAM / "modification_tracks.json").exists(), reason="RF00167 modifications not built")


def _load():
    return json.load(open(FAM / "modifications.json")), json.load(open(FAM / "modification_tracks.json"))


def test_counts_and_inventory():
    d, t = _load(); validate_json(d, "modifications"); validate_json(t, "modification_tracks")
    c = d["counts"]
    assert (c["n_instances"], c["n_instances_ok"], c["n_entries"], c["n_polymer_residues"], c["n_missing_identity"]) == (84, 84, 66, 6601, 0)
    assert (c["n_nonstandard"], c["n_distinct_ccd"], c["n_instances_with_nonstandard"], c["n_entries_with_nonstandard"]) == (13, 7, 9, 8)
    assert (c["n_on_consensus_axis"], c["n_off_axis"], c["off_axis_by_state"], c["consensus_positions_with_nonstandard"]) == (6, 7, {"not_submitted_flank": 7}, [22, 32, 35, 84, 87])
    assert set(d["components"]) == {"GTP", "GDP", "CCC", "23G", "IU", "50N", "50L"}
    assert {k: v["n_observations"] for k, v in d["components"].items()} == {"GTP": 4, "GDP": 1, "CCC": 2, "23G": 2, "IU": 2, "50N": 1, "50L": 1}
    # parents: CCD first, NAKB for components without a CCD parent; the mapping's 'X' is never used
    assert {k: (v["parent"], v["parent_source"]) for k, v in d["components"].items()} == {"GTP": ("G", "CCD"), "GDP": ("G", "CCD"), "CCC": ("C", "CCD"), "23G": ("G", "CCD"), "IU": ("U", "CCD"), "50N": ("C", "NAKB"), "50L": ("G", "NAKB")}


def test_23g_and_pdbe_crosscheck():
    d, _ = _load(); o = [x for x in d["observations"] if x["ccd_id"] == "23G"]
    assert {(x["seq_name"], x["label_seq_id"], x["auth_seq_id"]) for x in o} == {("3ivn_1_A_A", 70, "74"), ("3ivn_1_B_B", 70, "74")}
    assert all(x["rfam_consensus_pos"] is None and not x["on_consensus_axis"] and x["outside_mapped_region"] and x["is_3prime_terminal"] and not x["has_coordinates"] for x in o)
    assert all(x["pdbe"]["status"] == "reported" and x["pdbe"]["agrees"] for x in o)
    not_rep = sorted((x["ccd_id"], x["seq_name"], x["label_seq_id"]) for x in d["observations"] if x["pdbe"]["status"] != "reported")
    assert not_rep == [("50L", "5c7w_1_C_A", 64), ("50N", "5c7w_1_C_A", 4)]                      # included although absent from PDBe
    assert d["provenance"]["pdbe_modified_AA_or_NA"]["n_records"] == 11 and d["provenance"]["pdbe_modified_AA_or_NA"]["family_chain_records_without_observation"] == []
    assert d["provenance"]["nakb"]["status"] == "used" and d["provenance"]["nakb"]["last_modified"]


def test_occurrence_regression_and_reconstruction():
    d, t = _load(); P = {p["rfam_consensus_pos"]: p for p in t["positions"]}
    expected = {22: ("50N", 1, 83), 32: ("IU", 1, 84), 35: ("IU", 1, 84), 84: ("50L", 1, 83), 87: ("CCC", 2, 77)}
    for pos, (ccd, k, n) in expected.items():
        assert (P[pos]["n_positive"], P[pos]["denominator"], [x["ccd_id"] for x in P[pos]["ccd"]]) == (k, n, [ccd])
        assert P[pos]["frequency"] == round(k / n, 4)
    zero = [p for p in t["positions"] if p["rfam_consensus_pos"] not in expected]
    assert len(zero) == 97 and all(p["n_positive"] == 0 and p["denominator"] > 0 and p["frequency"] == 0.0 for p in zero)
    assert P[87]["n_no_residue"] == 7 and P[87]["n_denominator_without_coordinates"] == 14
    # independent reconstruction from the deposited components of each instance + the validated mapping files
    inst = {i["seq_name"]: i for i in d["instances"]}
    for p in t["positions"]:
        pos_, neg, none = Counter(), 0, 0
        for sn, r in inst.items():
            m = json.load(open(FAM / "mappings" / f"{sn}.json"))
            res = next((x for x in m["residues"] if x["alignment_state"] == "match" and x["rfam_consensus_pos"] == p["rfam_consensus_pos"]), None)
            if res is None: none += 1; continue
            comp = r["deposited_components"][res["polymer_seq_position"] - 1]
            if comp in ("A", "C", "G", "U"): neg += 1
            else: pos_[comp] += 1
        assert (sum(pos_.values()), neg, none) == (p["n_positive"], p["n_negative"], p["n_no_residue"])
        assert {x["ccd_id"]: x["n_instances"] for x in p["ccd"]} == dict(pos_)
