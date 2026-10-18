"""Structural clustering for the Structures tab: the ported UPGMA/cut/medoid code (synthetic) and the generated RF00167 file (regression)."""
import json
from pathlib import Path
import numpy as np
import pytest
from rna_explorer.export import validate_json
from rna_explorer.structure_clustering import average_linkage, check_matrix, cluster_metric, cophenetic_correlation, cut_tree, medoid, read_matrix

ROOT = Path(__file__).resolve().parent.parent
FAM = ROOT / "data/families/RF00167"


def _rand(n=12, seed=3):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, 3)); X[: n // 3] += 6
    D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
    return D


def test_upgma_matches_scipy_and_is_deterministic():
    sp = pytest.importorskip("scipy.cluster.hierarchy")
    from scipy.spatial.distance import squareform
    D = _rand()
    Z, order, leaves = average_linkage(D)
    Zs = sp.linkage(squareform(D, checks=False), "average")
    assert np.allclose(sorted(z[2] for z in Z), sorted(Zs[:, 2]))                                  # identical merge heights
    for k in (2, 3, 4):
        parts, *_ = cut_tree(Z, len(D), k)
        mine = sorted(sorted(leaves[p]) for p in parts)
        theirs = sp.fcluster(Zs, k, "maxclust")
        assert mine == sorted(sorted(np.where(theirs == c)[0].tolist()) for c in set(theirs))   # identical partitions
    assert abs(cophenetic_correlation(Z, leaves, D) - np.corrcoef(squareform(D, checks=False), sp.cophenet(Zs))[0, 1]) < 1e-9
    assert average_linkage(D) == (Z, order, leaves) and sorted(order) == list(range(len(D)))


def test_ties_cut_and_medoid():
    D = np.array([[0, 1, 1, 4], [1, 0, 1, 4], [1, 1, 0, 4], [4, 4, 4, 0]], float)   # three-way tie at 1
    Z, order, leaves = average_linkage(D)
    assert Z[0][:2] == [0, 1]                       # first pair met wins the tie
    assert order[:3] == [0, 1, 2] and order[-1] == 3  # larger child drawn first
    parts, cut, low, high = cut_tree(Z, 4, 2)
    assert sorted(len(leaves[p]) for p in parts) == [1, 3] and cut == (low + high) / 2 and low == 4.0 and high == 1.0
    assert medoid([0, 1, 2], D) == (0, 1.0)         # all equal: lowest index
    res = cluster_metric(["a", "b", "c", "d"], {"a": "1x", "b": "1x", "c": "2y", "d": "3z"}, D, 2)
    g1, g2 = res["groups"]
    assert (g1["label"], g1["n_instances"], g1["n_entries"], g1["representative"]) == ("Structural group 1", 3, 2, "a")
    assert g2["members"] == ["d"] and g2["within"]["mean"] is None and g2["small"] and not g2["instance_level_view"]
    assert g1["within"] == {"mean": 1.0, "max": 1.0, "n_pairs": 3} and g1["across"] == {"mean": 4.0, "n_pairs": 3} and g1["across_to_within_ratio"] == 4.0
    assert res["assignments"] == {"a": 1, "b": 1, "c": 1, "d": 2} and res["cut"]["is_largest_gap_of_tree"]


def test_matrix_checks(tmp_path):
    p = tmp_path / "m.csv"; p.write_text(",a,b\na,0,1\nb,1,0\n")
    assert read_matrix(p, ["a", "b"]).tolist() == [[0, 1], [1, 0]]
    with pytest.raises(ValueError):
        read_matrix(p, ["b", "a"])
    with pytest.raises(ValueError):
        check_matrix(np.array([[0, 1], [2, 0]], float), 0.0, "x")
    with pytest.raises(ValueError):
        check_matrix(np.array([[1, 0.5], [0.5, 0.9]], float), 1.0, "tm")
    with pytest.raises(ValueError):
        check_matrix(np.array([[0, np.nan], [np.nan, 0]]), 0.0, "x")


@pytest.mark.skipif(not (FAM / "structure_clustering.json").exists(), reason="RF00167 clustering not built")
def test_rf00167_regression():
    d = json.load(open(FAM / "structure_clustering.json")); validate_json(d, "structure_clustering")
    ds = d["dataset"]; inst = json.load(open(FAM / "instances.json"))
    assert (ds["n_mapped_instances"], ds["n_mapped_entries"], ds["n_clustered_instances"], ds["n_clustered_entries"], ds["n_excluded_instances"]) == (84, 66, 76, 61, 8)
    assert sorted(e["seq_name"] for e in ds["excluded"]) == ["2xo1_1_A_A", "3fo4_1_A_A", "6p2h_1_A_A", "9nbc_1_A_A", "9nbh_1_A_A", "9nbh_1_B_B", "9nbh_1_C_C", "9nbh_1_D_D"]
    assert set(ds["instances"]) | {e["seq_name"] for e in ds["excluded"]} == {i["seq_name"] for i in inst} and len(set(ds["instances"])) == 76
    assert ds["n_shared_positions"] == 57 and ds["coverage_threshold"] == 64 and d["default_metric"] == "discrepancy"
    n = 76
    exp = {"discrepancy": ([72, 4], ["6uc8_1_B_A", "5e54_1_B_B"], 2, "stable"), "tm": ([51, 17, 4, 4], ["6uc8_1_B_A", "9v4v_1_A_A", "5e54_1_B_B", "9ndd_1_A_A"], 4, "provisional")}
    for key, m in d["metrics"].items():
        c = m["clustering"]
        assert len(m["values_upper"]) == n * (n - 1) // 2 and all(np.isfinite(m["values_upper"]))
        assert sorted(c["leaf_order"]) == list(range(n)) and len(c["linkage"]) == n - 1
        assert [g["n_instances"] for g in c["groups"]] == exp[key][0] and [g["representative"] for g in c["groups"]] == exp[key][1]
        assert c["cut"]["k"] == exp[key][2] and c["validation"] == exp[key][3] and c["cut"]["is_largest_gap_of_tree"]
        assert sorted(c["assignments"]) == sorted(ds["instances"]) and sum(g["n_instances"] for g in c["groups"]) == n
        assert all(g["label"] == f"Structural group {g['group_id']}" for g in c["groups"]) and [g["group_id"] for g in c["groups"]] == list(range(1, len(c["groups"]) + 1))
        # rebuild the matrix and check the tree, the medoids and the contiguity of each group in leaf order
        M = np.zeros((n, n)); iu = np.triu_indices(n, 1); M[iu] = m["values_upper"]; M = M + M.T
        np.fill_diagonal(M, 1.0 if key == "tm" else 0.0)
        dist = 1 - M if key == "tm" else M
        np.fill_diagonal(dist, 0.0)
        assert (m["kind"], m["distance"]) == (("similarity", "1 - TM-score") if key == "tm" else ("dissimilarity", "D"))
        Z, order, leaves = average_linkage(dist)
        assert order == c["leaf_order"] and [[int(a), int(b), round(h, 6), int(s)] for a, b, h, s in Z] == c["linkage"]
        ids = ds["instances"]
        for g in c["groups"]:
            idx = [ids.index(s) for s in g["members"]]
            assert ids[medoid(idx, dist)[0]] == g["representative"] and g["representative"] in g["members"]
            assert [ids[i] for i in c["leaf_order"][g["leaf_range"][0]: g["leaf_range"][1] + 1]] == g["members"]
    dg = d["metrics"]["discrepancy"]["clustering"]["groups"]
    assert sorted(dg[1]["members"]) == ["5e54_1_A_A", "5e54_1_B_B", "5swd_1_A_A", "5swd_1_B_B"] and dg[1]["n_entries"] == 2
    assert abs(d["metrics"]["discrepancy"]["clustering"]["cut"]["height"] - 0.495389) < 1e-6 and d["metrics"]["discrepancy"]["clustering"]["cophenetic_correlation"] == 0.9845
    assert d["metrics"]["tm"]["clustering"]["cophenetic_correlation"] == 0.9562                    # equal to the research values
    assert d["metrics"]["tm"]["value_range"]["max"] <= 1.0 and d["metrics"]["tm"]["label"] == "TM-score"


# --- per-run research summary layout (RF00162 and later) and the RF00162 groups (2026-10-17) ---
from rna_explorer.structure_clustering import research_evidence  # noqa: E402

FAM162 = Path(__file__).resolve().parents[1] / "data/families/RF00162"


def test_research_evidence_per_run_layout(tmp_path):
    summary = {"primary_run": "T91", "runs": {"T91": {"metrics": {"disc_base": {"noise_floor": {"n": 2, "p95": 0.0602}, "cophenetic_correlation": 0.97, "partitions": {"3": {"silhouette": 0.56}}}}}},
               "stability": {"position_bootstrap": {"disc_base": {"n_bootstrap": 50, "ARI_k3_p05": 1.0}},
                             "element_jackknife": {"H1": {"n_positions_kept": 70, "disc_base_ARI_k3": 1.0, "rmsd_c4_ARI_k3": 1.0}}},
               "run_sensitivity_ARI": {"disc_base|T91_vs_Tte_only|k3": 0.0, "one_minus_tm_fixed|T91_vs_Tte_only|k4": 0.327}}
    p = tmp_path / "clustering_summary.json"; p.write_text(json.dumps(summary))
    d = research_evidence(p, "discrepancy")
    assert d["layout"] == "per_run" and d["primary"] == "T91" and d["noise_floor_p95"] == 0.0602 and d["noise_floor_n_pairs"] == 2
    assert d["cophenetic_correlation"] == 0.97 and d["silhouette_by_k"] == {"3": 0.56}
    assert d["position_bootstrap"] == {"n_bootstrap": 50, "ARI_k3_p05": 1.0} and d["element_jackknife"] == {"H1": {"disc_base_ARI_k3": 1.0}}
    assert d["run_sensitivity_ARI"] == {"T91_vs_Tte_only|k3": 0.0}
    t = research_evidence(p, "tm")
    assert "position_bootstrap" not in t and t["run_sensitivity_ARI"] == {"T91_vs_Tte_only|k4": 0.327}


@pytest.mark.skipif(not (FAM162 / "structure_clustering.json").exists(), reason="RF00162 clustering not built")
def test_rf00162_groups():
    d = json.load(open(FAM162 / "structure_clustering.json"))
    validate_json(d, "structure_clustering")
    assert d["dataset"]["n_clustered_instances"] == 29 and d["dataset"]["n_excluded_instances"] == 0 and d["dataset"]["n_shared_positions"] == 86
    disc, tm = d["metrics"]["discrepancy"]["clustering"], d["metrics"]["tm"]["clustering"]
    assert disc["cut"]["k"] == 3 and disc["validation"] == "stable" and disc["cut"]["is_largest_gap_of_tree"]
    assert [g["n_instances"] for g in disc["groups"]] == [26, 2, 1]
    assert sorted(disc["groups"][1]["members"]) == ["3v7e_2_C_C", "3v7e_2_D_D"] and disc["groups"][2]["members"] == ["4kqy_1_A_A"]
    assert tm["cut"]["k"] == 4 and tm["validation"] == "provisional" and [g["n_instances"] for g in tm["groups"]] == [22, 4, 2, 1]
    assert sorted(tm["groups"][1]["members"]) == ["2ygh_1_A_A", "4b5r_1_A_A", "5fjc_1_A_A", "5fk4_1_A_A"]
    for c in (disc, tm):   # RF00162 floor omitted by configuration (only two duplicate pairs); the research value is kept as evidence
        assert c["floor"]["value"] is None and "Not shown" in c["floor"]["definition"] and "only two duplicate pairs" in c["floor"]["definition"]
        assert c["research_evidence"]["layout"] == "per_run" and c["research_evidence"]["noise_floor_n_pairs"] == 2
    assert disc["research_evidence"]["position_bootstrap"]["ARI_k3_p05"] == 1.0


FAM059 = Path(__file__).resolve().parents[1] / "data/families/RF00059"


@pytest.mark.skipif(not (FAM059 / "structure_clustering.json").exists(), reason="RF00059 clustering not built")
def test_rf00059_apo_vs_bound_groups():
    d = json.load(open(FAM059 / "structure_clustering.json"))
    validate_json(d, "structure_clustering")
    assert d["dataset"]["n_clustered_instances"] == 37 and [e["seq_name"] for e in d["dataset"]["excluded"]] == ["9c6j_1_A_A"]
    assert d["dataset"]["n_shared_positions"] == 42
    for key, validation, floor in (("discrepancy", "stable", 0.1268), ("tm", "provisional", 0.2078)):
        c = d["metrics"][key]["clustering"]
        assert c["cut"]["k"] == 2 and c["validation"] == validation and c["cut"]["is_largest_gap_of_tree"]
        assert [g["n_instances"] for g in c["groups"]] == [35, 2] and sorted(c["groups"][1]["members"]) == ["8f4o_1_A_A", "8f4o_1_B_B"]
        assert c["floor"]["value"] == floor and c["floor"]["n_pairs"] == 15          # the per-run summary's duplicate floor is read


FAM177 = Path(__file__).resolve().parents[1] / "data/families/RF00177"


@pytest.mark.skipif(not (FAM177 / "structure_clustering.json").exists(), reason="RF00177 clustering not built")
def test_rf00177_head_rotation_groups():
    d = json.load(open(FAM177 / "structure_clustering.json"))
    validate_json(d, "structure_clustering")
    assert d["dataset"]["n_clustered_instances"] == 60 and d["dataset"]["n_excluded_instances"] == 0 and d["dataset"]["n_shared_positions"] == 1469
    groups = {}
    for key in ("discrepancy", "tm"):
        c = d["metrics"][key]["clustering"]
        assert c["cut"]["k"] == 2 and c["validation"] == "provisional" and c["cut"]["is_largest_gap_of_tree"]
        assert [g["n_instances"] for g in c["groups"]] == [50, 10] and c["floor"]["value"] is None and "Not shown" in c["floor"]["definition"]
        groups[key] = sorted(c["groups"][1]["members"])
    assert groups["discrepancy"] == groups["tm"]                                        # both metrics: the same 10 most rotated chains


from rna_explorer.structure_clustering import internal_gaps  # noqa: E402


def test_internal_gaps_and_gap_free_medoid():
    res = lambda p, num, obs: {"polymer_seq_position": p, "author_residue_number": num, "author_insertion_code": "", "has_coordinates": obs}
    # unmodelled at the ends are not breaks; inside they are
    assert internal_gaps({"residues": [res(1, 10, False), res(2, 11, True), res(3, 12, False), res(4, 13, True), res(5, 14, False)]}) == ["12"]
    assert internal_gaps({"residues": [res(1, 10, False), res(2, 11, True), res(3, 12, True)]}) == []
    D = np.array([[0, .1, .2, .3], [.1, 0, .2, .3], [.2, .2, 0, .4], [.3, .3, .4, 0]], dtype=float)
    assert medoid([0, 1, 2, 3], D)[0] == 0                       # plain medoid
    assert medoid([0, 1, 2, 3], D, candidates=[2, 3])[0] == 2     # restricted choice, mean still over all members
    ids, entries = ["a", "b", "c", "d"], {"a": "1", "b": "1", "c": "2", "d": "3"}
    g = cluster_metric(ids, entries, D, 1, gaps={"a": ["5"], "b": ["7"], "c": [], "d": []})["groups"][0]
    assert g["representative"] == "c" and g["medoid"] == "a" and g["representative_rule"] == "medoid among members without internal coordinate gaps" and g["representative_internal_gaps"] == []
    g = cluster_metric(ids, entries, D, 1, gaps={k: ["1"] for k in ids})["groups"][0]
    assert g["representative"] == "a" and g["representative_rule"] == "medoid (every member has internal coordinate gaps)" and g["representative_internal_gaps"] == ["1"]
    assert cluster_metric(ids, entries, D, 1)["groups"][0]["representative_rule"] == "medoid"   # no gap information: plain medoid


FAM005 = Path(__file__).resolve().parents[1] / "data/families/RF00005"


@pytest.mark.skipif(not (FAM005 / "structure_clustering.json").exists(), reason="RF00005 clustering not built")
def test_rf00005_groups_and_gap_free_representatives():
    d = json.load(open(FAM005 / "structure_clustering.json"))
    validate_json(d, "structure_clustering")
    for key, sizes in (("discrepancy", [56, 2, 1, 1]), ("tm", [55, 2, 1, 1, 1])):
        c = d["metrics"][key]["clustering"]
        assert c["validation"] == "provisional" and [g["n_instances"] for g in c["groups"]] == sizes and c["floor"]["value"] is None
        for g in c["groups"]:                                                         # every representative has a continuous backbone
            assert internal_gaps(json.load(open(FAM005 / "mappings" / f"{g['representative']}.json"))) == []
    g1 = d["metrics"]["discrepancy"]["clustering"]["groups"][0]
    assert g1["medoid"] == "9tnl_6_Y_F" and g1["representative"] == "6ufg_1_B_A"     # plain medoid lacks 4 internal nucleotides
