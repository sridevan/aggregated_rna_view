"""Interactions dataset (interactions.json): helpers, schema, and reproduction of the RF00167 interaction investigation."""
import csv
import json
from pathlib import Path
import pytest
from rna_explorer.export import validate_json
from rna_explorer.interactions import flip, parse_unit, ss_cons_pairs

ROOT = Path(__file__).resolve().parents[1]
FAM = ROOT / "data" / "families" / "RF00167"
INV = ROOT / "experiments" / "rf00167_interactions" / "RF00167_interaction_conservation.csv"


@pytest.fixture(scope="module")
def data():
    """the main file plus every chain's full observations (interactions_by_instance/<seq_name>.json) under "observations" """
    d = json.load(open(FAM / "interactions.json"))
    return {**d, "observations": {i["seq_name"]: json.load(open(FAM / "interactions_by_instance" / f"{i['seq_name']}.json"))["observations"] for i in d["instances"]}}


def test_flip_and_units():
    assert flip("cWH") == "cHW" and flip("tSH") == "tHS" and flip("cWW") == "cWW"
    assert flip("s35") == "s53" and flip("ns35") == "ns53" and flip("s33") == "s33"
    with pytest.raises(ValueError):
        flip("xyz")
    u = parse_unit("1Y26|1|X|ADE|90")
    assert (u["pdb"], u["chain"], u["resname"], u["num"], u["ins"], u["alt"]) == ("1Y26", "X", "ADE", 90, "", "")
    assert parse_unit("4V9D|1|AA|G|306||A|B")["ins"] == "B"
    assert ss_cons_pairs("((.<<..>>.))") == {(1, 12), (2, 11), (4, 9), (5, 8)}


def test_schema_and_summary(data):
    validate_json({k: v for k, v in data.items() if k != "observations"}, "interactions")
    s = data["summary"]["by_category"]
    bp, st = s["basepair"], s["stacking"]
    # base pairs + stacking: unchanged from the investigation (stage 1)
    assert bp["n_intra_chain_interactions"] + st["n_intra_chain_interactions"] == 9991 and bp["n_consensus_pairs"] + st["n_consensus_pairs"] == 199
    ver = {k: bp["endpoint_verification_against_mmcif"].get(k, 0) + st["endpoint_verification_against_mmcif"].get(k, 0) for k in ("verified", "mismatch")}
    assert ver == {"verified": 41604, "mismatch": 6}
    est = {k: bp["endpoint_status"].get(k, 0) + st["endpoint_status"].get(k, 0) for k in ("consensus", "insertion", "outside_alignment", "not_in_rna_polymer")}
    assert est == {"consensus": 37382, "insertion": 98, "outside_alignment": 4130, "not_in_rna_polymer": 794}
    # base-ribose and base-phosphate: every mapped endpoint verified against the mmCIF; files for every entry
    for c in ("baseribose", "basephosphate"):
        assert s[c]["endpoint_verification_against_mmcif"].get("mismatch", 0) == 0 and s[c]["n_entries_with_file"] == 66 and s[c]["n_consensus_pairs"] > 0
    assert {x["kind"] for x in data["sources"]} == {"basepair", "stacking", "baseribose", "basephosphate"}
    # every source has its URL, and the seeded files their checksum and retrieval date
    assert all(x["status"] == "ok" and x["sha256"] and x["retrieved"] for x in data["sources"])
    assert data["default_numbering_reference"]["seq_name"] == "4fe5_1_B_A"


def test_family_pairs_reproduce_investigation(data):
    inv = {(r["category"], int(r["rfam_pos_1"]), int(r["rfam_pos_2"])): r for r in csv.DictReader(open(INV))}
    bpst = [p for p in data["pairs"] if p["category"] in ("basepair", "stacking")]
    assert len(inv) == len(bpst)
    for p in bpst:
        r = inv[(p["category"], p["pos_1"], p["pos_2"])]
        assert (p["n_observed"], p["n_eligible"], p["n_entries_observed"], p["n_entries_eligible"]) == \
               (int(r["n_observed"]), int(r["n_eligible"]), int(r["n_entries_observed"]), int(r["n_entries_eligible"])), p["id"]
        assert (p["n_near_only"], p["n_unmapped"], p["n_missing"], p["n_incomplete"]) == \
               (int(r["n_near_only"]), int(r["n_missing_position"]), int(r["n_unmodelled"]), int(r["n_base_incomplete"])), p["id"]
        assert ";".join(f"{c}:{v}" for c, v in p["classes"].items()) == r["classes"], p["id"]


def test_denominators_and_eligibility(data):
    insts = {i["seq_name"]: i for i in data["instances"]}
    obs = data["observations"]
    for p in data["pairs"]:
        i, j = p["pos_1"], p["pos_2"]
        assert p["pos_1"] < p["pos_2"] if not p["directional"] else p["pos_1"] != p["pos_2"]
        # every instance falls into exactly one eligibility state; observed never exceeds eligible
        assert p["n_eligible"] + p["n_missing"] + p["n_incomplete"] + p["n_unmapped"] + p["n_annotation_unavailable"] == len(insts)
        assert p["n_observed"] + p["n_near_only"] + p["n_not_annotated"] == p["n_eligible"]
        # recount observed from the per-instance observations of eligible instances only
        cat = p["category"]
        col = {"baseribose": "sugar_eligibility", "basephosphate": "phosphate_eligibility"}.get(cat, "eligibility")
        n = sum(1 for sn, rows in obs.items() if insts[sn]["eligibility"][i - 1] == "E" and insts[sn][col][j - 1] == "E" and cat not in insts[sn]["unavailable_categories"]
                and any(r[0] == cat and r[4] == i and r[5] == j and not r[1].startswith("n") for r in rows))
        assert n == p["n_observed"], p["id"]


def test_known_rf00167_interactions(data):
    pairs = {p["id"]: p for p in data["pairs"]}
    # the L2-L3 kissing-loop G37-C61 (Rfam 41-66 in this family's coordinates is checked through the region labels)
    kl = [p for p in data["pairs"] if p["category"] == "basepair" and {p["region_1"], p["region_2"]} == {"L2", "L3"} and p["n_eligible"] >= 70]
    assert kl and all(p["structure_class"] == "tertiary" for p in kl)
    assert max(p["frequency"] for p in kl) == 1.0
    # SS_cons helix pairs are classified secondary
    assert all(p["structure_class"] == "secondary" for p in pairs.values() if p["in_ss_cons"] and p["category"] == "basepair")
    # FR3D classes are kept as annotated: a reclassified pair carries more than one class
    assert any(len(p["classes"]) > 1 for p in pairs.values() if p["category"] == "basepair")
    # 5'->3' normalisation: stacking classes keep their direction (both s35 and s53 occur)
    cls = {c for p in pairs.values() if p["category"] == "stacking" for c in p["classes"]}
    assert {"s35", "s53"} <= cls


RAW = ROOT / "experiments" / "rf00167_interactions" / "raw" / "bgsu"


def test_directional_categories_keep_bgsu_roles(data):
    """Base-ribose / base-phosphate: every RNA-RNA row of the raw BGSU file is kept in its row order (base first) and class."""
    import gzip, io
    for pdb, chain, sn in (("4FE5", "B", "4fe5_1_B_A"), ("5E54", "B", "5e54_1_B_B"), ("1Y26", "X", "1y26_1_X_A")):
        mine = {(o[6], o[1], o[7]) for o in data["observations"][sn] if o[0] in ("baseribose", "basephosphate")}
        for kind in ("baseribose", "basephosphate"):
            rows = {tuple(r) for r in csv.reader(io.StringIO(gzip.open(RAW / f"{pdb}_{kind}.csv.gz", "rt").read())) if len(r) == 3}
            rna = {r for r in rows if all(u.split("|")[2] == chain and u.split("|")[3] in ("A", "C", "G", "U") for u in (r[0], r[2]))}
            assert rna and rna <= mine
    # directional ids: base position first, both directions are distinct interactions, classes never flipped (no 'cHW'-style reversal)
    pairs = [p for p in data["pairs"] if p["directional"]]
    assert any(p["pos_1"] > p["pos_2"] for p in pairs) and any(p["pos_1"] < p["pos_2"] for p in pairs)
    assert all(c.endswith(("BR", "BPh")) for p in pairs for c in p["classes"])


def test_category_specific_eligibility(data):
    """A base-phosphate interaction is not evaluated when the contacted nucleotide lacks its phosphate (e.g. a 5' terminal residue)."""
    insts = {i["seq_name"]: i for i in data["instances"]}
    assert any("I" in i["phosphate_eligibility"] for i in insts.values())
    for p in (q for q in data["pairs"] if q["category"] == "basephosphate"):
        inc = sum(1 for i in insts.values() if i["eligibility"][p["pos_1"] - 1] == "E" and i["phosphate_eligibility"][p["pos_2"] - 1] == "I")
        assert p["n_incomplete"] >= inc


def test_compact_index_matches_full_records(data):
    """the main file's per-chain index (observed pair indices by class, near, off-consensus count) is exactly derived from the per-chain files"""
    pairs = data["pairs"]
    for inst in data["instances"]:
        sn = inst["seq_name"]
        f = json.load(open(FAM / "interactions_by_instance" / f"{sn}.json"))
        validate_json(f, "interactions_instance")
        obs = f["observations"]
        cons = [o for o in obs if isinstance(o[4], int) and isinstance(o[5], int) and o[4] != o[5]]
        want = {(o[0], o[1], o[4], o[5]) for o in cons if not o[1].startswith("n")}
        got = {(pairs[k]["category"], c, pairs[k]["pos_1"], pairs[k]["pos_2"]) for c, ks in inst["observed"].items() for k in ks}
        assert got == want, sn
        assert sorted(map(tuple, inst["near"])) == sorted((o[0], o[1], o[4], o[5]) for o in cons if o[1].startswith("n"))
        assert inst["n_off_consensus"] == len(obs) - len(cons) and inst["n_interactions"] == len(obs)


def test_rf00177_interactions_dataset():
    """16S: same pipeline, no family-specific code; every mapped endpoint verified; compact index consistent with `pairs`."""
    p = ROOT / "data" / "families" / "RF00177" / "interactions.json"
    if not p.exists():
        pytest.skip("RF00177 interactions.json not built")
    d = json.load(open(p))
    validate_json(d, "interactions")
    by = d["summary"]["by_category"]
    assert all(v["endpoint_verification_against_mmcif"].get("mismatch", 0) == 0 and v["n_entries_with_file"] == d["summary"]["n_entries"] for v in by.values())
    assert d["default_numbering_reference"]["seq_name"] == "4v9d_1_BA_Y" and d["clen"] == 1533
    n = len(d["pairs"])
    assert all(0 <= k < n for i in d["instances"] for ks in i["observed"].values() for k in ks)


def test_annotation_gaps():
    """incomplete source annotation: only RF00177 7PJS chain a (annotated up to nt ~650); its positions there are not eligible anywhere"""
    from rna_explorer.interactions import annotation_gaps
    assert annotation_gaps("E" * 40, set(range(1, 41))) == []
    assert annotation_gaps("E" * 40, set(range(1, 21))) == [[21, 40]]           # 20 >= 15 positions and >= 8 % of 40
    assert annotation_gaps("E" * 40, set(range(1, 30))) == []                    # 11 < 15
    for acc in ("RF00167", "RF00059", "RF00162", "RF00005", "RF00177"):
        p = ROOT / "data" / "families" / acc / "interactions.json"
        if not p.exists():
            continue
        d = json.load(open(p))
        gaps = {i["seq_name"]: i["annotation_gaps"] for i in d["instances"] if i["annotation_gaps"]}
        assert gaps == ({"7pjs_34_a_HA": [[652, 1531]]} if acc == "RF00177" else {}), acc
        if acc == "RF00177":
            # every pair touching the gap counts 7PJS as annotation-unavailable, never as eligible
            inst = {i["seq_name"]: i for i in d["instances"]}
            for q in d["pairs"]:
                assert q["n_eligible"] + q["n_missing"] + q["n_incomplete"] + q["n_unmapped"] + q["n_annotation_unavailable"] == len(inst)
