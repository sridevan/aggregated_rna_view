"""Regression of the generated RF00167 pocket files (runs only when present)."""
import json
from pathlib import Path
import pytest
from rna_explorer.export import validate_json
from rna_explorer.pockets import qualifying

FAM = Path(__file__).resolve().parent.parent / "data/families/RF00167"
pytestmark = pytest.mark.skipif(not (FAM / "pocket_tracks.json").exists(), reason="RF00167 pockets not built")


def test_pockets_dataset_invariants():
    P = json.load(open(FAM / "pockets.json")); T = json.load(open(FAM / "pocket_tracks.json"))
    validate_json(P, "pockets"); validate_json(T, "pocket_tracks")
    assert P["counts"]["n_instances"] == 84 and P["counts"]["status_counts"] == {"pockets": 65, "no_pockets": 19} and P["counts"]["n_pockets_total"] == 70
    assert P["method"]["fpocket_parameters"] == {"m": 3.0, "M": 5.7, "i": 42, "D": 1.65, "A": 3, "p": 0.0}
    for r in P["instances"]:
        assert r["input"]["roundtrip_ok"] and not r["warnings"]
        for p in r["pockets"]:
            assert p["nt_match_csv"] and p["map_consistent"] and p["unresolved_nucleotides"] == 0 and p["score"] >= 0
            assert len({n["label_seq_id"] for n in p["nucleotides"]}) == p["n_nucleotides"]            # one membership per nucleotide
            assert all(n["rfam_consensus_pos"] is None for n in p["nucleotides"] if n["alignment_state"] in ("insert", "local_end", "not_submitted_flank"))
    # representative mapping checks (pilot values)
    by = {r["seq_name"]: r for r in P["instances"]}
    for sn in ("4fe5_1_B_A", "1y26_1_X_A"):
        p = by[sn]["pockets"][0]
        assert p["consensus_positions"] == [25, 26, 51, 54, 55, 56, 79, 80, 81] and sorted(p["pocket_nt_author_csv"]) == [21, 22, 47, 50, 51, 52, 73, 74, 75]
    assert by["4fe5_1_B_A"]["pockets"][0]["score"] == 0.55 and by["1y26_1_X_A"]["pockets"][0]["score"] == 0.53
    assert by["3ivn_1_A_A"]["status"] == "no_pockets" and by["9ndd_1_A_A"]["status"] == "no_pockets"
    # tracks: denominators, filters, known site
    for flt, tr in T["tracks"].items():
        for p in tr:
            assert p["n_supporting_instances"] <= p["n_evaluable_instances"] <= 84 and (p["frequency"] is None) == (p["n_evaluable_instances"] == 0)
            assert p["n_supporting_entries"] <= p["n_supporting_instances"] and p["n_distinct_pockets"] >= p["n_supporting_instances"]
    a = T["tracks"]["all_score_ge_0"]
    assert a[50]["n_supporting_instances"] == 59 and a[50]["n_evaluable_instances"] == 75 and a[50]["frequency"] == 0.7867   # consensus 51
    assert a[55]["n_supporting_instances"] == 61 and a[55]["n_evaluable_instances"] == 80                                   # consensus 56
    assert T["tracks"]["top1"][50]["n_supporting_instances"] == 59 and T["tracks"]["score_ge_0.4"][50]["n_supporting_instances"] == 52
    # independent reconstruction of one position from the pocket records
    sup = {r["seq_name"] for r in P["instances"] if r["status"] != "failed" and any(56 in p["consensus_positions"] for p in qualifying(r["pockets"], "all_score_ge_0"))}
    assert len(sup) == a[55]["n_supporting_instances"]
