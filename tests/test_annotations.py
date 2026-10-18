"""Unit tests for the Milestone 3B annotation backend on a small synthetic family with a fake fetcher.

Synthetic family RFTEST (clen 10):
  1abc chain A (entity 1, single chain): residues 1..6 (auth 11..16); consensus 3,4,-,5,6,7 (residue 3 is an insertion).
  2xyz chains A and B (entity 1, two chains): residues 1..5 -> consensus 3..7; chain B residue 5 has no coordinates.
  3nnn chain A (entity 1): every external source fails (HTTP error) -> 'source unavailable'.
"""
import json
from pathlib import Path
import pytest
from rna_explorer.annotation_mapping import InstanceIndex, UnitId
from rna_explorer.annotation_sources import SourceMeta, parse_motifs_html, parse_loops_csv, parse_reference_annotations
from rna_explorer.annotations import AnnotationBuilder, load_config
from rna_explorer.export import validate_json

ACC = "RFTEST"


# ----------------------------------------------------------------------------- fixtures
def _residue(i, auth, name, state, cons, coords=True, pp="*"):
    return {"polymer_seq_position": i, "author_residue_number": auth, "author_insertion_code": "", "residue_name": name, "nucleotide_parent": name,
            "is_modified": False, "alignment_state": state, "rfam_consensus_pos": cons, "insertion_anchor": 4 if state == "insert" else None,
            "insertion_offset": 1 if state == "insert" else None, "posterior_prob_char": pp if state == "match" else "", "posterior_prob": 0.975 if state == "match" else None,
            "has_coordinates": coords}


def _mapping(pdb, entity, label, auth, residues, key):
    cons = {r["rfam_consensus_pos"]: r["polymer_seq_position"] for r in residues if r["alignment_state"] == "match"}
    return {"pdb_id": pdb, "entity_id": entity, "label_asym_id": label, "auth_asym_id": auth, "rfam_acc": ACC, "sequence_key": key,
            "mapped_region": {"start": 1, "end": len(residues)}, "alignment": None, "summary": {}, "residues": residues,
            "consensus_columns": [{"rfam_consensus_pos": p, "state": "match" if p in cons else "uncovered", "polymer_seq_position": cons.get(p)} for p in range(1, 11)],
            "provenance": {"rfam_acc": ACC, "rfam_release": "15.1", "cm_name": "T", "cm_checksum": "1", "cm_clen": 10, "pipeline_version": "0.1.0", "generated": "2026-10-15T00:00:00+00:00"}}


@pytest.fixture
def fam_dir(tmp_path: Path) -> Path:
    d = tmp_path / "families" / ACC
    (d / "mappings").mkdir(parents=True)
    json.dump({"rfam_acc": ACC, "clen": 10, "description": "test"}, open(d / "family.json", "w"))
    r1 = [_residue(1, 11, "G", "match", 3), _residue(2, 12, "G", "match", 4), _residue(3, 13, "A", "insert", None), _residue(4, 14, "U", "match", 5),
          _residue(5, 15, "C", "match", 6, pp="8"), _residue(6, 16, "C", "match", 7)]
    rA = [_residue(i, 10 + i, "GGUCC"[i - 1], "match", i + 2) for i in range(1, 6)]
    rB = [_residue(i, 10 + i, "GGUCC"[i - 1], "match", i + 2, coords=(i != 5)) for i in range(1, 6)]
    r3 = [_residue(i, i, "GGUCC"[i - 1], "match", i + 2) for i in range(1, 6)]
    maps = {"1abc_1_A_A": _mapping("1abc", 1, "A", "A", r1, "k1"), "2xyz_1_A_A": _mapping("2xyz", 1, "A", "A", rA, "k2"),
            "2xyz_1_B_B": _mapping("2xyz", 1, "B", "B", rB, "k2"), "3nnn_1_A_A": _mapping("3nnn", 1, "A", "A", r3, "k3")}
    inst = []
    for sn, m in maps.items():
        json.dump(m, open(d / "mappings" / f"{sn}.json", "w"))
        inst.append({"pdb_id": m["pdb_id"], "entity_id": m["entity_id"], "label_asym_id": m["label_asym_id"], "auth_asym_id": m["auth_asym_id"], "seq_name": sn,
                     "sequence_key": m["sequence_key"], "mapping_file": f"mappings/{sn}.json"})
    json.dump(inst, open(d / "instances.json", "w"))
    return d


def _row(loop_id, location, ann, group=None, mapped=None):
    link = f'<a href="https://rna.bgsu.edu/rna3dhub/motif/view/{group}">{group}</a>' if group else ""
    ann_cell = f"{ann}<br>{link}" if mapped is None else ann
    map_cell = "" if mapped is None else f"{mapped[0]}<br>{mapped[1]} {link}"
    return (f"<tr><td class=\"loop\"><input id='{loop_id}'></td><td><a href=\"https://rna.bgsu.edu/rna3dhub/loops/view/{loop_id}\">{loop_id}</a></td>"
            f"<td>{location}</td><td>{ann_cell}</td><td>{map_cell}</td></tr>")


def _prow(loop_id, problem, ann):
    return (f"<tr><td class=\"loop\"><input id='{loop_id}'></td><td><a href=\"https://rna.bgsu.edu/rna3dhub/loops/view/{loop_id}\">{loop_id}</a></td>"
            f"<td><label>{problem}</label></td><td>{ann}</td></tr>")


def _page(valid_rows, prob_rows=()):
    return ("<html><h3>Valid loops</h3><table><tbody>" + "".join(valid_rows) + "</tbody></table>"
            "<h3>Problematic loops</h3><table><tbody>" + "".join(prob_rows) + "</tbody></table></html>")


REF_FILE = "# following https://example.org/paper\n1ABC|1|A|G|11, Helix P1\n1ABC|1|A|G|12, Helix P1\n1ABC|1|A|A|13, L1\n1ABC|1|A|U|14, L1\n1ABC|1|A|C|15, L1\n1ABC|1|A|C|16, Helix P1\n"
CHAIN_2XYZ_A = "2XYZ|1|A|G|11\tHelix P1\n2XYZ|1|A|G|12\tHelix P1\n2XYZ|1|A|U|13\tL1\n2XYZ|1|A|C|14\tL1\n2XYZ|1|A|C|15\tHelix P1\n"
LOOPS_1ABC = ('"HL_1ABC_001","1ABC|1|A|G|12,1ABC|1|A|A|13,1ABC|1|A|U|14"\n'
              '"J3_1ABC_001","1ABC|1|A|G|11,1ABC|1|A|G|12,1ABC|1|A|C|15,1ABC|1|A|C|16"\n'
              '"J4_1ABC_001","1ABC|1|A|G|11,1ABC|1|A|C|16"\n')
LOOPS_2XYZ = ('"IL_2XYZ_001","2XYZ|1|A|G|12,2XYZ|1|A|U|13,2XYZ|1|B|G|12"\n'
              '"HL_2XYZ_001","2XYZ|1|B|U|13,2XYZ|1|B|C|14,2XYZ|1|B|C|15"\n'
              '"HL_2XYZ_002","2XYZ|1|Z|U|13,2XYZ|1|Z|C|14"\n')
MOTIFS_1ABC = _page([_row("HL_1ABC_001", "1/A/12:14", "Kink", "HL_1.1"), _row("J3_1ABC_001", "1/A/11:12, 1/A/15:16", "No text annotation", "J3_1.1"),
                     _row("J4_1ABC_001", "1/A/11, 1/A/16", "No text annotation")])
MOTIFS_2XYZ = _page([_row("HL_2XYZ_001", "1/B/13:15", "No text annotation", "HL_1.1", ("homologous", "HL_1ABC_001")),
                     _row("HL_2XYZ_002", "1/Z/13:14", "No text annotation", "HL_9.9", ("geometric", "HL_1ABC_001"))],
                    [_prow("IL_2XYZ_001", "poor quality numbers", "No text annotation")])
SITES_1ABC = {"1abc": {"length": 6, "dataType": "LIGAND BINDING SITES", "data": [
    {"name": "ADENINE", "accession": "ADE", "additionalData": {"bestChainId": "A", "entityId": 1, "scaffoldId": "x", "residueCount": 3}, "dataType": "Ligand",
     "residues": [{"startIndex": 1, "endIndex": 2, "startCode": "G", "endCode": "G", "indexType": "PDB", "additionalData": {"boundMoleculeId": "bm1"}},
                  {"startIndex": 2, "endIndex": 4, "startCode": "G", "endCode": "U", "indexType": "PDB", "additionalData": {"boundMoleculeId": "bm2"}}]},
    {"name": "MAGNESIUM ION", "accession": "MG", "additionalData": {"bestChainId": "A", "entityId": 1}, "dataType": "Ligand",
     "residues": [{"startIndex": 5, "endIndex": 5, "startCode": "C", "endCode": "C", "indexType": "PDB", "additionalData": {"boundMoleculeId": "bm3"}}]},
    {"name": "GTP", "accession": "GTP", "additionalData": {"bestChainId": "A", "entityId": 1}, "dataType": "Ligand",
     "residues": [{"startIndex": 1, "endIndex": 1, "startCode": "G", "endCode": "G", "indexType": "PDB", "additionalData": {"boundMoleculeId": "bm4"}}]}]}}
SITES_2XYZ = {"2xyz": {"length": 5, "data": [
    {"name": "ADENINE", "accession": "ADE", "additionalData": {"bestChainId": "B", "entityId": 1}, "dataType": "Ligand",
     "residues": [{"startIndex": 2, "endIndex": 3, "startCode": "G", "endCode": "U", "indexType": "PDB", "additionalData": {"boundMoleculeId": "bm1"}}]}]}}


class FakeFetcher:
    """Serves canned responses by URL; anything mentioning 3NNN/3nnn fails; unknown URLs are empty."""
    def __init__(self):
        self.calls = []

    def text(self, url, suffix, seed_name=None):
        self.calls.append(url)
        if "3NNN" in url or "3nnn" in url:
            return None, SourceMeta(url, "http_error", None, "404")
        table = {f"nucleotide_annotations/{ACC}_annotations.txt": REF_FILE, "chain=1ABC%7C1%7CA": REF_FILE.replace(", ", "\t"), "chain=2XYZ%7C1%7CA": CHAIN_2XYZ_A,
                 "chain=2XYZ%7C1%7CB": "", "loops/download/1ABC": LOOPS_1ABC, "loops/download/2XYZ": LOOPS_2XYZ, "pdb/1ABC/motifs": MOTIFS_1ABC, "pdb/2XYZ/motifs": MOTIFS_2XYZ,
                 "binding_sites/1abc/1": json.dumps(SITES_1ABC), "binding_sites/2xyz/1": json.dumps(SITES_2XYZ), "release/hl/current": "<h1>Release 4.15</h1>",
                 "release/il/current": "Release 4.15", "release/j3/current": "Release 4.15", "release/j4/current": "Release 4.15"}
        for k, v in table.items():
            if k in url:
                return v, SourceMeta(url, "ok" if v.strip() else "empty", "2026-10-15")
        return "", SourceMeta(url, "empty", "2026-10-15")


@pytest.fixture
def built(fam_dir):
    b = AnnotationBuilder(fam_dir, FakeFetcher(), load_config())
    summary = b.run()
    return b, fam_dir, summary


# ----------------------------------------------------------------------------- parsers and identifier mapping
def test_unit_id_and_resolution(fam_dir):
    u = UnitId.parse("1ABC|1|A|G|12")
    assert (u.pdb_id, u.auth_asym_id, u.author_residue_number, u.insertion_code) == ("1abc", "A", "12", "")
    assert UnitId.parse("1ABC|1|A|G|12|||B").insertion_code == "B"
    idx = InstanceIndex(fam_dir)
    r = idx.resolve_unit("1ABC|1|A|G|12")
    assert (r.seq_name, r.label_asym_id, r.label_seq_id, r.rfam_consensus_pos, r.mapping_status) == ("1abc_1_A_A", "A", 2, 4, "mapped")
    assert idx.resolve_unit("1ABC|1|A|A|13").mapping_status == "insert"
    assert idx.resolve_unit("1ABC|1|A|A|99").mapping_status == "unresolved_identifier"
    assert idx.resolve_unit("1ABC|1|Q|A|11").seq_name is None
    assert idx.resolve_seq_index("2xyz_1_B_B", 5).mapping_status == "mapped_no_coordinates"
    assert idx.resolve_seq_index("2xyz_1_B_B", 5).rfam_consensus_pos == 7


def test_parsers():
    rows, comments = parse_reference_annotations(REF_FILE)
    assert len(rows) == 6 and rows[0].labels == ["Helix P1"] and comments[0].startswith("following")
    assert parse_reference_annotations("X|1|A|G|1\tHelix P1; L1")[0][0].labels == ["Helix P1", "L1"]
    loops = parse_loops_csv(LOOPS_1ABC)
    assert [L.loop_id for L in loops] == ["HL_1ABC_001", "J3_1ABC_001", "J4_1ABC_001"] and loops[1].loop_type == "J3" and len(loops[1].units) == 4
    m = {x.loop_id: x for x in parse_motifs_html(MOTIFS_2XYZ)}
    assert m["HL_2XYZ_001"].motif_group == "HL_1.1" and m["HL_2XYZ_001"].assignment_type == "homologous" and m["HL_2XYZ_001"].mapped_via_loop_id == "HL_1ABC_001"
    assert m["HL_2XYZ_002"].assignment_type == "geometric"
    assert m["IL_2XYZ_001"].atlas_status == "problematic" and m["IL_2XYZ_001"].problem == "poor quality numbers" and m["IL_2XYZ_001"].motif_group is None
    m1 = {x.loop_id: x for x in parse_motifs_html(MOTIFS_1ABC)}
    assert m1["HL_1ABC_001"].assignment_type == "member" and m1["HL_1ABC_001"].annotation_text == "Kink"
    assert m1["J4_1ABC_001"].motif_group is None and m1["J4_1ABC_001"].atlas_status == "valid"


# ----------------------------------------------------------------------------- reference elements
def test_reference_transfer(built):
    b, d, _ = built
    ref = json.load(open(d / "reference_structural_elements.json"))
    assert ref["track"] == [None, None, "Helix P1", "Helix P1", "L1", "L1", "Helix P1", None, None, None]
    assert ref["reference"]["pdb_id"] == "4feo" or ref["reference"]["pdb_id"] == "1abc"
    assert ref["reference"]["seq_name"] == "1abc_1_A_A" and ref["n_consensus_positions_annotated"] == 5
    assert [x["unit_id"] for x in ref["reference_nucleotides_not_on_consensus"]] == ["1ABC|1|A|A|13"]      # insertion kept, not forced
    assert ref["elements"][0] == {"element": "Helix P1", "spans": [[3, 4], [7, 7]], "n_positions": 3}
    it = ref["instance_transfer"]
    assert it["n_instances_evaluated"] == 2 and sorted(it["unavailable_instances"]) == ["2xyz_1_B_B", "3nnn_1_A_A"]
    assert it["n_labels_disagree"] == 0 and it["n_labels_agree_with_family_track"] == 10
    inst = json.load(open(d / "annotations" / "2xyz_1_A_A.json"))
    assert inst["sources"]["reference_structural_elements"]["status"] == "evaluated"
    assert [r["reference_element"] for r in inst["residues"]] == ["Helix P1", "Helix P1", "L1", "L1", "Helix P1"]
    assert json.load(open(d / "annotations" / "2xyz_1_B_B.json"))["sources"]["reference_structural_elements"]["status"] == "unavailable"


# ----------------------------------------------------------------------------- loops
def test_loop_membership_strands_and_insertions(built):
    b, d, _ = built
    loops = {L["loop_id"]: L for L in json.load(open(d / "loops.json"))["loops"]}
    j3 = loops["J3_1ABC_001"]
    assert len(j3["strands"]) == 2 and [s["n_nucleotides"] for s in j3["strands"]] == [2, 2]
    assert j3["consensus_positions"] == [3, 4, 6, 7] and j3["motif_group"] == "J3_1.1" and j3["motif_group_assignment"] == "member"
    hl = loops["HL_1ABC_001"]
    assert hl["mapping_status_counts"] == {"mapped": 2, "insert": 1} and hl["consensus_positions"] == [4, 5]
    assert [n["mapping_status"] for n in hl["nucleotides"]] == ["mapped", "insert", "mapped"]
    assert hl["nucleotides"][1]["rfam_consensus_pos"] is None and hl["nucleotides"][1]["label_seq_id"] == 3
    assert loops["J4_1ABC_001"]["motif_group"] is None and loops["J4_1ABC_001"]["atlas_status"] == "valid"
    assert loops["IL_2XYZ_001"]["atlas_status"] == "problematic" and loops["IL_2XYZ_001"]["problem"] == "poor quality numbers"
    # loop with no nucleotide in a family chain: retained as a compact record (no nucleotide-level data), not counted as unresolved
    assert "HL_2XYZ_002" not in loops
    out = {L["loop_id"]: L for L in b.loops_outside}
    assert out["HL_2XYZ_002"] == {"loop_id": "HL_2XYZ_002", "pdb_id": "2xyz", "loop_type": "HL", "auth_asym_ids": ["Z"], "n_nucleotides": 2, "motif_group": "HL_9.9",
                                  "atlas_status": out["HL_2XYZ_002"]["atlas_status"], "loop_url": out["HL_2XYZ_002"]["loop_url"], "reason": "no nucleotide in a family chain"}
    assert not any(u["identifier"] == "2XYZ|1|Z|U|13" for u in b.unresolved) and b.loop_source["2xyz"]["n_loops_outside_family_chains"] == 1
    # an inter-chain loop touching a family chain is unchanged (full record)
    assert "IL_2XYZ_001" in loops


def test_multi_chain_loop_is_one_observation(built):
    b, d, _ = built
    loops = {L["loop_id"]: L for L in json.load(open(d / "loops.json"))["loops"]}
    il = loops["IL_2XYZ_001"]
    assert il["multi_chain"] and [c["auth_asym_id"] for c in il["chains"]] == ["A", "B"] and sorted(il["instances"]) == ["2xyz_1_A_A", "2xyz_1_B_B"]
    t = json.load(open(d / "tracks.json"))["positions"][3]      # consensus 4: A|12 and B|12
    il_t = t["internal_loops"]
    assert il_t["n_loop_observations"] == 1 and il_t["n_instances"] == 2 and il_t["n_entries"] == 1
    assert il_t["by_motif_group"] == {"unassigned": {"n_loops": 1, "n_instances": 2, "n_entries": 1}} and il_t["n_loops_without_motif_group"] == 1


def test_j3_j4_overlap_union_not_sum(built):
    b, d, _ = built
    p = json.load(open(d / "tracks.json"))["positions"][2]   # consensus 3: residue 11 in J3 and J4 of 1abc
    assert p["junction_loops_J3"]["n_instances"] == 1 and p["junction_loops_J4"]["n_instances"] == 1
    assert p["junction_loops_union"]["n_instances"] == 1 and p["junction_loops_union"]["n_loop_observations"] == 2
    inst = json.load(open(d / "annotations" / "1abc_1_A_A.json"))
    assert sorted(x["loop_type"] for x in inst["residues"][0]["loops"]) == ["J3", "J4"]     # both kept independently


def test_motif_group_counts_and_homologous_assignment(built):
    b, d, _ = built
    p = json.load(open(d / "tracks.json"))["positions"][4]   # consensus 5: HL_1ABC_001 (U|14), HL_2XYZ_001 (B|13)
    hl = p["hairpin_loops"]
    assert hl["n_loop_observations"] == 2 and hl["n_distinct_motif_groups"] == 1 and hl["by_motif_group"]["HL_1.1"]["n_loops"] == 2 and hl["by_motif_group"]["HL_1.1"]["n_entries"] == 2
    assert hl["confidence"] == {"high": 2, "medium": 0, "low": 0}


# ----------------------------------------------------------------------------- ligands
def test_entity_level_ligand_attribution_and_dedup(built):
    b, d, _ = built
    lig = json.load(open(d / "ligands.json"))
    recs = {(r["pdb_id"], r["ccd_id"]): r for r in lig["bound_molecules"]}
    ade = recs[("1abc", "ADE")]
    assert ade["chain_attribution"] == "resolved" and ade["n_entity_chains"] == 1
    assert [c["label_seq_id"] for c in ade["contacts"]] == [1, 2, 3, 4] and ade["contacts"][1]["bound_molecule_ids"] == ["bm1", "bm2"]   # deduplicated, ids merged
    assert ade["contacts"][2]["mapping_status"] == "insert" and ade["contacts"][2]["rfam_consensus_pos"] is None
    assert ade["consensus_positions"] == [3, 4, 5] and ade["n_bound_molecule_instances"] == 2
    assert recs[("1abc", "MG")]["ligand_class"] == "metal_ion" and not recs[("1abc", "MG")]["in_default_track"] and recs[("1abc", "MG")]["exclusion_reason"] == "class:metal_ion"
    gtp = recs[("1abc", "GTP")]
    assert gtp["flags"] == ["contacts_only_polymer_termini"] and gtp["ligand_class"] == "ligand"
    x = recs[("2xyz", "ADE")]
    assert x["chain_attribution"] == "unresolved" and x["best_chain_label_asym_id"] == "B" and sorted(x["entity_chains"]) == ["2xyz_1_A_A", "2xyz_1_B_B"]
    a = json.load(open(d / "annotations" / "2xyz_1_A_A.json"))
    bm = a["bound_molecules"][0]
    assert bm["attribution"] == "entity_unresolved" and bm["is_best_chain"] is False
    assert a["residues"][1]["bound_molecules"][0]["attribution"] == "entity_unresolved"
    assert lig["counts"]["n_distinct_ccd_ids"] == 3 and lig["counts"]["by_class"] == {"ligand": 2, "metal_ion": 1}
    assert lig["entry_entity_sources"]["3nnn:1"]["status"] == "unavailable"


def test_ligand_track_entry_level_counts_and_denominators(built):
    b, d, _ = built
    cfg = load_config()
    excl_termini = next(r for r in cfg["artefact_rules"] if r["id"] == "contacts_only_polymer_termini").get("exclude_from_default_track", False)
    p = json.load(open(d / "tracks.json"))["positions"][3]   # consensus 4: 1abc residue 2 (ADE), 2xyz residue 2 (ADE, entity-level)
    lg = p["ligand_binding_sites"]
    assert lg["n_entry_entities"] == 2 and lg["by_ccd_n_entry_entities"] == {"ADE": 2}
    den = lg["denominators"]
    # observable entities at pos 4: 1abc, 2xyz, 3nnn; evaluated: 1abc, 2xyz; 3nnn unavailable
    assert den == {"counting_unit": "pdb_entry_entity", "n_observable_entry_entities": 3, "n_source_evaluated_entry_entities": 2, "n_observable_and_evaluated": 2, "n_positive": 2,
                   "n_positive_observable_and_evaluated": 2, "n_not_annotated_by_source": 0, "n_source_unavailable": 1, "n_positive_not_observable": 0}
    assert lg["frequency"] == 1.0
    p3 = json.load(open(d / "tracks.json"))["positions"][2]  # consensus 3: ADE at 1abc residue 1 (+GTP flagged)
    assert p3["ligand_binding_sites"]["by_ccd_n_entry_entities"] == ({"ADE": 1} if excl_termini else {"ADE": 1, "GTP": 1})
    assert p3["ligand_binding_sites"]["by_class_n_entry_entities"] == {"ligand": 1}


# ----------------------------------------------------------------------------- denominators / missing sources
def test_loop_denominators_and_missing_source(built):
    b, d, _ = built
    tr = json.load(open(d / "tracks.json"))
    assert tr["source_availability"]["motif_atlas_loops"] == {"n_instances_evaluated": 3, "n_instances_unavailable": 1}
    p = tr["positions"][6]    # consensus 7: observable 1abc(6), 2xyz A(5); 2xyz B(5) has no coordinates; 3nnn(5) unavailable source
    assert p["observable"]["n_instances"] == 3
    hl = p["hairpin_loops"]["denominators"]
    # HL_2XYZ_001 contains B|15 (consensus 7) -> positive on a non-observable instance (no coordinates) is reported, not hidden
    assert hl["n_observable_instances"] == 3 and hl["n_source_evaluated_instances"] == 3 and hl["n_observable_and_evaluated"] == 2
    assert hl["n_positive"] == 1 and hl["n_positive_not_observable"] == 1 and hl["n_positive_observable_and_evaluated"] == 0
    assert hl["n_not_annotated_by_source"] == 2 and hl["n_source_unavailable"] == 1 and p["hairpin_loops"]["frequency"] == 0.0
    # consensus 1 and 2: nobody observable -> frequency null, not 0
    assert tr["positions"][0]["hairpin_loops"]["frequency"] is None and tr["positions"][0]["ligand_binding_sites"]["frequency"] is None
    inst = json.load(open(d / "annotations" / "3nnn_1_A_A.json"))
    assert {k: v["status"] for k, v in inst["sources"].items()} == {"reference_structural_elements": "unavailable", "motif_atlas_loops": "unavailable", "pdbe_binding_sites": "unavailable"}


# ----------------------------------------------------------------------------- schemas and traceability
def test_schemas_and_traceability(built):
    b, d, summary = built
    for name, schema in [("reference_structural_elements.json", "reference_structural_elements"), ("loops.json", "loops"), ("ligands.json", "ligands"), ("tracks.json", "tracks")]:
        validate_json(json.load(open(d / name)), schema)
    for f in (d / "annotations").glob("*.json"):
        validate_json(json.load(open(f)), "annotation_instance")
    loops = {L["loop_id"]: L for L in json.load(open(d / "loops.json"))["loops"]}
    tr = json.load(open(d / "tracks.json"))
    for p in tr["positions"]:
        for track in ("hairpin_loops", "internal_loops", "junction_loops_J3", "junction_loops_J4"):
            ids = p[track]["loop_ids"]
            assert all(p["rfam_consensus_pos"] in loops[i]["consensus_positions"] for i in ids)
            assert p[track]["n_loop_observations"] == len(ids)
    assert summary["n_unresolved_identifiers"] == 0     # the only unresolved units belonged to a loop with no family-chain nucleotide (kept as a compact record)


# ----------------------------------------------------------------------------- overlapping loop classes (regression for RF00167 consensus 53 semantics)
LOOPS_OVERLAP = ('"HL_1ABC_002","1ABC|1|A|G|12,1ABC|1|A|A|13,1ABC|1|A|U|14"\n'          # hairpin inside the junction strand (consensus 4, insert, 5)
                 '"J3_1ABC_002","1ABC|1|A|G|11,1ABC|1|A|G|12,1ABC|1|A|A|13,1ABC|1|A|U|14,1ABC|1|A|C|15"\n'   # J3 strand covering the hairpin
                 '"J3_1ABC_003","1ABC|1|A|G|11,1ABC|1|A|G|12,1ABC|1|A|A|13,1ABC|1|A|U|14"\n'                 # second J3 of the same chain (alternative boundaries)
                 '"J4_1ABC_002","1ABC|1|A|G|11,1ABC|1|A|C|15"\n')                                              # J4 excluding the hairpin strand
MOTIFS_OVERLAP = _page([_row("HL_1ABC_002", "1/A/12:14", "Purine riboswitch", "HL_9.1"), _row("J3_1ABC_002", "1/A/11:15", "No text annotation", "J3_9.1"),
                        _row("J3_1ABC_003", "1/A/11:14", "No text annotation", "J3_8.1"), _row("J4_1ABC_002", "1/A/11, 1/A/15", "No text annotation", "J4_9.1")])


class OverlapFetcher(FakeFetcher):
    def text(self, url, suffix, seed_name=None):
        if "loops/download/1ABC" in url:
            return LOOPS_OVERLAP, SourceMeta(url, "ok", "2026-10-16")
        if "pdb/1ABC/motifs" in url:
            return MOTIFS_OVERLAP, SourceMeta(url, "ok", "2026-10-16")
        return super().text(url, suffix, seed_name)


def test_overlapping_hl_and_j3_are_both_kept_and_counted_per_instance(fam_dir):
    b = AnnotationBuilder(fam_dir, OverlapFetcher(), load_config())
    b.run()
    tr = json.load(open(fam_dir / "tracks.json"))
    p4 = tr["positions"][3]      # consensus 4 = residue 12: in HL_1ABC_002, J3_1ABC_002, J3_1ABC_003 (2 J3 loops, one chain)
    assert p4["hairpin_loops"]["n_instances"] == 1 and p4["junction_loops_J3"]["n_instances"] == 1 and p4["junction_loops_J4"]["n_instances"] == 0
    assert p4["junction_loops_J3"]["n_loop_observations"] == 2            # two J3 loops, still one instance
    assert p4["junction_loops_union"]["n_instances"] == 1 and p4["junction_loops_union"]["n_loop_observations"] == 2
    p3 = tr["positions"][2]      # consensus 3 = residue 11: J3 ×2 and J4 -> union counts the instance once, loops 3
    assert p3["junction_loops_J3"]["n_instances"] == 1 and p3["junction_loops_J4"]["n_instances"] == 1 and p3["junction_loops_union"]["n_instances"] == 1
    assert p3["junction_loops_union"]["n_loop_observations"] == 3 and p3["hairpin_loops"]["n_instances"] == 0
    # overlapping classes are preserved: both loops remain in loops.json with their original types and the insertion kept in the hairpin membership
    loops = {L["loop_id"]: L for L in json.load(open(fam_dir / "loops.json"))["loops"]}
    assert loops["HL_1ABC_002"]["loop_type"] == "HL" and loops["J3_1ABC_002"]["loop_type"] == "J3"
    assert set(loops["HL_1ABC_002"]["consensus_positions"]) <= set(loops["J3_1ABC_002"]["consensus_positions"])
    assert [n["mapping_status"] for n in loops["HL_1ABC_002"]["nucleotides"]] == ["mapped", "insert", "mapped"]
    inst = json.load(open(fam_dir / "annotations" / "1abc_1_A_A.json"))
    assert sorted(x["loop_type"] for x in inst["residues"][1]["loops"]) == ["HL", "J3", "J3"]


# --- BGSU per-chain annotations for a family without a BGSU annotation file (RF00162, 2026-10-17) ---
from rna_explorer.annotation_sources import fetch_chain_annotations  # noqa: E402


class _MessageFetcher:
    """BGSU replies HTTP 200 with a plain-text message (no unit_id, label rows)."""
    def __init__(self, body):
        self.body = body

    def text(self, url, suffix, seed_name=None):
        return self.body, SourceMeta(url, "ok", "2026-10-17")


def test_reference_parser_skips_lines_without_separator():
    rows, comments = parse_reference_annotations("# c\nCould not find annotation file RF00162_annotations.txt\nX|1|A|G|1\tHelix P1\n")
    assert [r.unit_id for r in rows] == ["X|1|A|G|1"] and rows[0].labels == ["Helix P1"] and comments == ["c"]


def test_chain_annotations_bgsu_message_is_source_not_available():
    rows, meta = fetch_chain_annotations(_MessageFetcher("Could not find annotation file RF00162_annotations.txt"), "5fjc", "A")
    assert rows == [] and meta.status == "not_available" and meta.detail.startswith("Could not find annotation file")
    assert "chain=5FJC%7C1%7CA" in meta.url and meta.retrieved == "2026-10-17"


def test_chain_annotations_unparseable_reply_is_parse_error_and_valid_rows_still_parse():
    rows, meta = fetch_chain_annotations(_MessageFetcher("<html>maintenance</html>"), "5fjc", "A")
    assert rows == [] and meta.status == "parse_error"
    rows, meta = fetch_chain_annotations(_MessageFetcher("5FJC|1|A|G|1\tHelix P1\n"), "5fjc", "A")
    assert len(rows) == 1 and meta.status == "ok"


def test_loops_by_instance_files(built):
    b, fam_dir, summary = built
    loops = json.load(open(fam_dir / "loops.json"))["loops"]
    for sn in ["1abc_1_A_A", "2xyz_1_A_A", "2xyz_1_B_B"]:
        f = fam_dir / "loops_by_instance" / f"{sn}.json"
        if not f.exists():
            continue
        d = json.load(open(f)); validate_json(d, "loops_instance")
        assert d["seq_name"] == sn and sorted(L["loop_id"] for L in d["loops"]) == sorted(L["loop_id"] for L in loops if sn in L["instances"])
    files = sorted(p.stem for p in (fam_dir / "loops_by_instance").glob("*.json"))
    assert files == sorted(i["seq_name"] for i in b.idx.instances)          # one file per chain instance
    # the inter-chain internal loop appears in both chains' files
    a, bb = (json.load(open(fam_dir / "loops_by_instance" / f"2xyz_1_{c}_{c}.json"))["loops"] for c in "AB")
    assert "IL_2XYZ_001" in {L["loop_id"] for L in a} and "IL_2XYZ_001" in {L["loop_id"] for L in bb}
    assert summary["file_sizes_bytes"]["loops_by_instance/ (total)"] > 0


# ----------------------------------------------------------------------------- higher-order junctions (J5) are part of the junction union
LOOPS_J5 = LOOPS_OVERLAP + '"J5_1ABC_001","1ABC|1|A|G|11,1ABC|1|A|C|15,1ABC|1|A|G|16"\n'
MOTIFS_J5 = _page([_row("HL_1ABC_002", "1/A/12:14", "Purine riboswitch", "HL_9.1"), _row("J3_1ABC_002", "1/A/11:15", "No text annotation", "J3_9.1"),
                   _row("J3_1ABC_003", "1/A/11:14", "No text annotation", "J3_8.1"), _row("J4_1ABC_002", "1/A/11, 1/A/15", "No text annotation", "J4_9.1"),
                   _row("J5_1ABC_001", "1/A/11, 1/A/15:16", "No text annotation", "J5_9.1")])


class J5Fetcher(FakeFetcher):
    def text(self, url, suffix, seed_name=None):
        if "loops/download/1ABC" in url:
            return LOOPS_J5, SourceMeta(url, "ok", "2026-10-18")
        if "pdb/1ABC/motifs" in url:
            return MOTIFS_J5, SourceMeta(url, "ok", "2026-10-18")
        return super().text(url, suffix, seed_name)


def test_j5_junction_loops_are_in_the_junction_union(fam_dir):
    AnnotationBuilder(fam_dir, J5Fetcher(), load_config()).run()
    tr = json.load(open(fam_dir / "tracks.json"))
    assert tr["tracks"]["junction_loops_union"]["loop_types"] == ["J3", "J4", "J5"]
    p3 = tr["positions"][2]      # consensus 3 = residue 11: J3 x2, J4 and J5, one chain
    u = p3["junction_loops_union"]
    assert u["n_loop_observations"] == 4 and u["n_instances"] == 1            # the J5 loop counted; the instance once (union, not a sum)
    assert u["by_type_n_instances"] == {"J3": 1, "J4": 1, "J5": 1} and u["loop_types"] == ["J3", "J4", "J5"]
    assert "J5" in u["rule"]
    assert "J5_1ABC_001" in {L["loop_id"] for L in json.load(open(fam_dir / "loops.json"))["loops"]}
