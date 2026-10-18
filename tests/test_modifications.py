"""Unit tests for rna_explorer/modifications.py on a synthetic mmCIF (no network, no family dataset)."""
import json
from pathlib import Path
import gemmi
import pytest
from rna_explorer.cache import Cache
from rna_explorer.export import validate_json
from rna_explorer.modifications import PDBE_FILES, aggregate, build_modifications, classify_position, extract_instance, read_polymer

# Entity 1 (RNA, chains A and B) = GTP A PSU C 5MC [G|A] N U ; entity 2 = non-polymer GTP ligand (asym C).
#   seq 1 GTP: 5' terminal, outside the mapped region          seq 5 5MC: insert state (no consensus position)
#   seq 2 A  : standard, cons 1                                 seq 6 G/A: microheterogeneity (missing identity), cons 4
#   seq 3 PSU: nonstandard, cons 2                              seq 7 N  : unknown component (missing identity), cons 5
#   seq 4 C  : standard, no coordinates in chain A, cons 3      seq 8 U  : standard, 3' terminal, cons 6 (chain B: deleted)
SEQ = ["GTP", "A", "PSU", "C", "5MC", "G", "N", "U"]


def _cif() -> str:
    eps = [(i + 1, m, "n") for i, m in enumerate(SEQ)]; eps[5] = (6, "G", "y"); eps.insert(6, (6, "A", "y"))
    scheme = []
    for asym in "AB":
        for num, mon, het in eps:
            scheme.append(f"{asym} 1 {num} {mon} {num + 10} . {het}")
    atoms, k = [], 1
    for asym in "AB":
        for i, mon in enumerate(SEQ, start=1):
            if asym == "A" and i == 4:      # chain A residue 4 unresolved
                continue
            atoms.append(f"HETATM {k} {mon if i != 6 else 'G'} {asym} {i} 1"); k += 1
    atoms.append(f"HETATM {k} GTP C . 1")                                   # non-polymer GTP ligand
    atoms.append(f"HETATM {k + 1} PSU A 3 2")                               # model 2 is ignored
    return "\n".join([
        "data_1TST", "loop_", "_entity.id", "_entity.type", "_entity.pdbx_description", "1 polymer RNA", "2 non-polymer GTP",
        "loop_", "_entity_poly.entity_id", "_entity_poly.type", "1 polyribonucleotide",
        "loop_", "_entity_poly_seq.entity_id", "_entity_poly_seq.num", "_entity_poly_seq.mon_id", "_entity_poly_seq.hetero", *[f"1 {n} {m} {h}" for n, m, h in eps],
        "loop_", "_pdbx_poly_seq_scheme.asym_id", "_pdbx_poly_seq_scheme.entity_id", "_pdbx_poly_seq_scheme.seq_id", "_pdbx_poly_seq_scheme.mon_id",
        "_pdbx_poly_seq_scheme.pdb_seq_num", "_pdbx_poly_seq_scheme.pdb_ins_code", "_pdbx_poly_seq_scheme.hetero", *scheme,
        "loop_", "_pdbx_nonpoly_scheme.asym_id", "_pdbx_nonpoly_scheme.entity_id", "_pdbx_nonpoly_scheme.mon_id", "C 2 GTP",
        "loop_", "_chem_comp.id", "_chem_comp.type", "_chem_comp.name",
        *[f"{c} '{t}' '{n}'" for c, t, n in [("A", "RNA linking", "ADENOSINE-5'-MONOPHOSPHATE"), ("C", "RNA linking", "CYTIDINE"), ("G", "RNA linking", "GUANOSINE"), ("U", "RNA linking", "URIDINE"),
                                             ("N", "RNA linking", "ANY 5'-MONOPHOSPHATE NUCLEOTIDE"), ("PSU", "RNA linking", "PSEUDOURIDINE-5'-MONOPHOSPHATE"),
                                             ("5MC", "RNA linking", "5-METHYLCYTIDINE-5'-MONOPHOSPHATE"), ("GTP", "non-polymer", "GUANOSINE-5'-TRIPHOSPHATE")]],
        "loop_", "_atom_site.group_PDB", "_atom_site.id", "_atom_site.label_comp_id", "_atom_site.label_asym_id", "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_model_num", *atoms, ""])


def _block():
    return gemmi.cif.read_string(_cif()).sole_block()


STATES = {1: ("not_submitted_flank", None), 2: ("match", 1), 3: ("match", 2), 4: ("match", 3), 5: ("insert", None), 6: ("match", 4), 7: ("match", 5), 8: ("match", 6)}


def _mapping(chain_b: bool = False) -> dict:
    res = []
    for i, mon in enumerate(SEQ, start=1):
        st, pos = STATES[i]
        if chain_b and i == 8:
            st, pos = "not_submitted_flank", None   # chain B: no residue at consensus 6
        res.append({"polymer_seq_position": i, "residue_name": mon if i != 6 else "G", "nucleotide_parent": "X" if mon == "N" else mon[0], "alignment_state": st,
                    "rfam_consensus_pos": pos, "insertion_anchor": 2 if st == "insert" else None, "author_residue_number": i + 10, "author_insertion_code": ""})
    return {"residues": res, "provenance": {"rfam_release": "15.1", "cm_checksum": "1", "pipeline_version": "0.1.0", "generated": "t", "infernal_version": "x"}}


INST = [{"seq_name": "1tst_1_A_A", "pdb_id": "1tst", "entity_id": 1, "label_asym_id": "A", "auth_asym_id": "A", "mapping_file": "mappings/1tst_1_A_A.json"},
        {"seq_name": "1tst_1_B_B", "pdb_id": "1tst", "entity_id": 1, "label_asym_id": "B", "auth_asym_id": "B", "mapping_file": "mappings/1tst_1_B_B.json"}]


class FakeIdx:
    def __init__(self):
        self.maps = {"1tst_1_A_A": _mapping(), "1tst_1_B_B": _mapping(True)}
    def match_residue(self, sn, pos):
        return next((r for r in self.maps[sn]["residues"] if r["alignment_state"] == "match" and r["rfam_consensus_pos"] == pos), None)


def test_polymer_membership_ignores_ligand_and_other_models():
    poly = read_polymer(_block(), "A", "1")
    assert sorted(poly["scheme"]) == list(range(1, 9)) and poly["entity_poly_type"] == "polyribonucleotide"
    assert poly["coords"][1] == {"GTP"} and 4 not in poly["coords"] and poly["coords"][3] == {"PSU"}     # model 2 and the C-asym GTP ligand never read
    assert poly["chem_comp"]["GTP"]["type"] == "non-polymer"        # CCD type says non-polymer, yet it is a polymer residue by the sequence scheme
    with pytest.raises(ValueError):
        read_polymer(_block(), "A", "2")                             # asym A is not entity 2


def test_extract_statuses_and_observations():
    r = extract_instance(_block(), INST[0], _mapping())
    st = r["_status_by_seq"]
    assert [st[i] for i in range(1, 9)] == ["positive", "negative", "positive", "negative", "positive", "missing", "missing", "negative"]
    assert r["_missing_reason_by_seq"] == {6: "microheterogeneity", 7: "unknown_component"}
    assert (r["n_standard"], r["n_nonstandard"], r["n_missing_identity"], r["n_without_coordinates"]) == (3, 3, 2, 1)
    assert r["deposited_components"] == ["GTP", "A", "PSU", "C", "5MC", "A|G", "N", "U"] and r["label_seq_ids_without_coordinates"] == [4]
    o = {x["label_seq_id"]: x for x in r["observations"]}
    assert set(o) == {1, 3, 5}                                       # one observation per nonstandard polymer residue (no ligand, no duplicates)
    assert o[1]["ccd_id"] == "GTP" and o[1]["is_5prime_terminal"] and o[1]["outside_mapped_region"] and o[1]["rfam_consensus_pos"] is None and not o[1]["on_consensus_axis"]
    assert o[3]["rfam_consensus_pos"] == 2 and o[3]["on_consensus_axis"] and o[3]["auth_seq_id"] == "13" and o[3]["chem_comp_name"].startswith("PSEUDOURIDINE")
    assert o[5]["alignment_state"] == "insert" and o[5]["rfam_consensus_pos"] is None and o[5]["insertion_anchor"] == 2   # insert never forced onto the axis
    assert {m["label_seq_id"]: m["reason"] for m in r["missing_identity"]} == {6: "microheterogeneity", 7: "unknown_component"}


def test_classify_conflicts():
    poly = {"scheme": {1: [{"mon": "A", "auth_seq_id": "1", "ins": "", "hetero": "n"}]}, "entity_poly_seq": {1: [("A", "n")]}, "coords": {1: {"A"}}, "chem_comp": {}}
    assert classify_position(1, poly, {"residue_name": "A"})["status"] == "negative"
    assert classify_position(1, poly, None)["missing_reason"] == "not_in_mapping"
    assert classify_position(1, poly, {"residue_name": "G"})["missing_reason"] == "mapping_identity_conflict"
    assert classify_position(1, poly | {"coords": {1: {"G"}}}, {"residue_name": "A"})["missing_reason"] == "coordinate_identity_conflict"
    assert classify_position(1, poly | {"entity_poly_seq": {1: [("G", "n")]}}, {"residue_name": "A"})["missing_reason"] == "sequence_record_conflict"
    assert classify_position(1, poly | {"coords": {}}, {"residue_name": "A"}) | {} == {"label_seq_id": 1, "component": "A", "components": ["A"], "has_coordinates": False, "status": "negative", "missing_reason": None}


def test_aggregation_counts_and_denominators():
    recs = [extract_instance(_block(), INST[0], _mapping()), extract_instance(_block(), INST[1], _mapping(True)),
            {"seq_name": "2bad_1_A_A", "pdb_id": "2bad", "status": "source_unavailable", "observations": []}]
    a = {p["rfam_consensus_pos"]: p for p in aggregate(recs, FakeIdx(), 7)}
    p2 = a[2]   # PSU in both chains sharing the entity: two instances, one entry
    assert (p2["n_positive"], p2["n_negative"], p2["denominator"], p2["frequency"], p2["n_supporting_entries"]) == (2, 0, 2, 1.0, 1)
    assert p2["ccd"] == [{"ccd_id": "PSU", "n_instances": 2, "n_entries": 1}] and p2["supporting_observations"] == ["1tst_1_A_A:3", "1tst_1_B_B:3"] and p2["n_source_unavailable"] == 1
    p3 = a[3]   # chain A unresolved but identity known: negative, counted in the denominator
    assert (p3["n_negative"], p3["denominator"], p3["frequency"], p3["n_denominator_without_coordinates"]) == (2, 2, 0.0, 1)
    p4 = a[4]   # microheterogeneity: missing in both, not in the denominator
    assert (p4["n_missing_identity"], p4["missing_reasons"], p4["denominator"], p4["frequency"]) == (2, {"microheterogeneity": 2}, 0, None)
    assert a[5]["missing_reasons"] == {"unknown_component": 2}
    p6 = a[6]   # chain B has no residue there: excluded, not a negative
    assert (p6["n_negative"], p6["n_no_residue"], p6["denominator"]) == (1, 1, 1)
    assert a[7]["n_no_residue"] == 2 and a[7]["frequency"] is None
    for p in a.values():
        assert p["n_positive"] + p["n_negative"] + p["n_missing_identity"] + p["n_no_residue"] + p["n_source_unavailable"] == 3   # every instance accounted once


def _family(tmp_path: Path) -> tuple[Path, Path]:
    out = tmp_path / "data"; fam = out / "families" / "RFTEST"; (fam / "mappings").mkdir(parents=True)
    json.dump({"clen": 7}, open(fam / "family.json", "w")); json.dump(INST, open(fam / "instances.json", "w"))
    json.dump(_mapping(), open(fam / "mappings/1tst_1_A_A.json", "w")); json.dump(_mapping(True), open(fam / "mappings/1tst_1_B_B.json", "w"))
    cache = tmp_path / "cache"; c = Cache(cache); url = f"{PDBE_FILES}/1tst_updated.cif"
    c._path(url, ".cif", "GET", None).write_text(_cif())
    return out, cache


def test_build_without_enrichment_writes_valid_json(tmp_path):
    out, cache = _family(tmp_path)
    res = build_modifications("RFTEST", cache, out, pdbe="cache", nakb="off", offline=True)
    fam = out / "families" / "RFTEST"
    d = json.load(open(fam / "modifications.json")); t = json.load(open(fam / "modification_tracks.json"))
    validate_json(d, "modifications"); validate_json(t, "modification_tracks")
    assert res["pdbe"] == "unavailable" and res["nakb"] == "not_used"                      # optional sources missing: core output still produced
    assert d["counts"]["n_nonstandard"] == 6 and d["counts"]["n_on_consensus_axis"] == 2 and d["counts"]["n_off_axis"] == 4 and d["counts"]["off_axis_by_state"] == {"not_submitted_flank": 2, "insert": 2}
    assert all(o["pdbe"] == {"status": "not_checked"} for o in d["observations"]) and all(o["parent"] is None and o["parent_source"] is None for o in d["observations"])
    assert d["components"]["PSU"]["nakb_in_reference"] is None and d["provenance"]["pdbe_modified_AA_or_NA"]["status"] == "unavailable"
    assert [p["n_positive"] for p in t["positions"]] == [0, 2, 0, 0, 0, 0, 0] and t["summary"]["n_off_axis"] == 4
    # reconstruction: occurrence counts follow from the detailed records (deposited components + validated mapping)
    for p in t["positions"]:
        assert len(p["supporting_observations"]) == p["n_positive"]


def test_build_with_seeded_pdbe_and_nakb(tmp_path):
    out, cache = _family(tmp_path)
    seed = tmp_path / "seed"; (seed / "pdbe").mkdir(parents=True); (seed / "nakb").mkdir()
    json.dump({"request_ids": ["1tst"], "http_status": 200, "retrieved_utc": "2026-01-01T00:00:00Z",
               "body": {"1tst": [{"struct_asym_id": "A", "residue_number": 3, "chem_comp_id": "PSU"}, {"struct_asym_id": "A", "residue_number": 3, "chem_comp_id": "PSU"}]}},
              open(seed / "pdbe/batch_all_entries.json", "w"))   # duplicated API record: still one observation
    json.dump({"PSU": {"standard_base": ["U"], "changes": [{"change_location": "base"}], "pdb": {"name": "PSEUDOURIDINE"}}}, open(seed / "nakb/modified_to_change_data.json", "w"))
    (seed / "nakb/json_headers.txt").write_text("Last-Modified: Fri, 17 Apr 2026 22:31:29 GMT\nETag: \"abc\"\n")
    res = build_modifications("RFTEST", cache, out, seed_dir=seed, offline=True)
    d = res["modifications"]; o = {x["observation_id"]: x for x in d["observations"]}
    assert res["pdbe"] == "checked" and res["nakb"] == "used"
    assert o["1tst_1_A_A:3"]["pdbe"] == {"status": "reported", "chem_comp_id": "PSU", "agrees": True, "n_records": 2} and o["1tst_1_B_B:3"]["pdbe"]["status"] == "not_reported"
    assert o["1tst_1_A_A:3"]["parent"] == "U" and o["1tst_1_A_A:3"]["parent_source"] == "NAKB" and o["1tst_1_A_A:3"]["nakb_change_locations"] == ["base"]   # CCD unavailable offline
    assert d["components"]["GTP"]["nakb_in_reference"] is False and d["provenance"]["nakb"]["etag"] == '"abc"'
    assert len(d["observations"]) == 6
