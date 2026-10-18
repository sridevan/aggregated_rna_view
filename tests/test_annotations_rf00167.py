"""Regression of the generated RF00167 annotation files against the Milestone 3A assessment (runs only when the files exist)."""
import json
from collections import Counter
from pathlib import Path
import pytest
from rna_explorer.export import validate_json

ROOT = Path(__file__).resolve().parent.parent
FAM = ROOT / "data/families/RF00167"
pytestmark = pytest.mark.skipif(not (FAM / "tracks.json").exists(), reason="RF00167 annotations not built")


@pytest.fixture(scope="module")
def data():
    return {"ref": json.load(open(FAM / "reference_structural_elements.json")), "loops": json.load(open(FAM / "loops.json")),
            "lig": json.load(open(FAM / "ligands.json")), "tracks": json.load(open(FAM / "tracks.json")), "summary": json.load(open(FAM / "annotation_summary.json")),
            "instances": json.load(open(FAM / "instances.json"))}


def test_schemas(data):
    validate_json(data["ref"], "reference_structural_elements"); validate_json(data["loops"], "loops"); validate_json(data["lig"], "ligands"); validate_json(data["tracks"], "tracks")
    files = sorted((FAM / "annotations").glob("*.json"))
    assert len(files) == len(data["instances"]) == 84
    for f in files[:12]:
        validate_json(json.load(open(f)), "annotation_instance")


def test_milestone_3a_counts(data):
    ref, loops, lig, summ = data["ref"], data["loops"], data["lig"], data["summary"]
    assert ref["reference"]["pdb_id"] == "4feo" and ref["reference"]["auth_asym_id"] == "B" and ref["reference"]["n_annotated_nucleotides"] == 67
    assert ref["n_consensus_positions_annotated"] == 67 and ref["instance_transfer"]["n_instances_evaluated"] == 80
    assert ref["instance_transfer"]["n_labels_disagree"] == 0 and ref["conflicts"] == []
    assert [e["element"] for e in ref["elements"]] == ["Helix P1", "J1/2", "Helix P2", "L2", "J2/3", "Helix P3", "L3", "J3/1"]
    c = loops["counts"]
    assert c["n_loops"] == 496 and c["n_with_motif_group"] == 456 and c["n_distinct_motif_groups"] == 29
    assert {t: v["n_loops"] for t, v in c["by_type"].items()} == {"HL": 234, "IL": 82, "J3": 114, "J4": 66}
    assert c["nucleotide_mapping_status_counts"] == {"mapped": 5419, "outside_mapped_region": 168, "insert": 68}
    assert sum(1 for L in loops["loops"] if L["multi_chain"]) == 4
    assert lig["counts"]["n_distinct_ccd_ids"] == 37 and lig["counts"]["by_class"] == {"ligand": 30, "metal_ion": 5, "additive": 2}
    assert lig["counts"]["n_entries_in_default_track"] == 63
    assert summ["n_unresolved_identifiers"] == 0


def test_tracks_traceable_to_loops(data):
    loops = {L["loop_id"]: L for L in data["loops"]["loops"]}
    for p in data["tracks"]["positions"]:
        pos = p["rfam_consensus_pos"]
        for track, t in (("hairpin_loops", "HL"), ("internal_loops", "IL"), ("junction_loops_J3", "J3"), ("junction_loops_J4", "J4")):
            expected = sorted(i for i, L in loops.items() if L["loop_type"] == t and pos in L["consensus_positions"])
            assert p[track]["loop_ids"] == expected
            inst = {n["seq_name"] for i in expected for n in loops[i]["nucleotides"] if n["rfam_consensus_pos"] == pos and n["mapping_status"] in ("mapped", "mapped_no_coordinates")}
            assert p[track]["n_instances"] == len(inst)
            groups = Counter(loops[i]["motif_group"] or "unassigned" for i in expected)
            assert {g: v["n_loops"] for g, v in p[track]["by_motif_group"].items()} == dict(groups)
        u = p["junction_loops_union"]
        assert max(p["junction_loops_J3"]["n_instances"], p["junction_loops_J4"]["n_instances"]) <= u["n_instances"] <= p["junction_loops_J3"]["n_instances"] + p["junction_loops_J4"]["n_instances"]
    # the J2/3 strand: FR3D hairpin HL_56046.1 and J3 overlap at consensus 53 and are both kept
    p = data["tracks"]["positions"][52]
    assert "HL_56046.1" in p["hairpin_loops"]["by_motif_group"] and p["junction_loops_J3"]["n_instances"] > 60 and p["junction_loops_J4"]["n_instances"] == 0


def test_ligand_track_traceable_and_entry_level(data):
    recs = data["lig"]["bound_molecules"]
    p = data["tracks"]["positions"][24]   # consensus 25
    expected = {f"{r['pdb_id']}:{r['entity_id']}" for r in recs if r["in_default_track"] and 25 in r["consensus_positions"]}
    assert set(p["ligand_binding_sites"]["entry_entities"]) == expected and p["ligand_binding_sites"]["n_entry_entities"] == 62
    assert p["ligand_binding_sites"]["denominators"]["n_observable_and_evaluated"] == 64
    multi = [r for r in recs if r["n_entity_chains"] > 1]
    assert multi and all(r["chain_attribution"] == "unresolved" for r in multi)
    assert all(r["chain_attribution"] == "resolved" for r in recs if r["n_entity_chains"] == 1)
    gtp = [r for r in recs if r["ccd_id"] == "GTP"]
    assert len(gtp) == 4 and all("contacts_only_polymer_termini" in r["flags"] for r in gtp)
    # contacts deduplicated per label_seq_id
    for r in recs:
        ids = [c["label_seq_id"] for c in r["contacts"]]
        assert len(ids) == len(set(ids))


def test_denominators_distinguish_missing_from_negative(data):
    tr = data["tracks"]
    assert tr["source_availability"]["reference_structural_elements"] == {"n_instances_evaluated": 80, "n_instances_unavailable": 4}
    for p in tr["positions"]:
        for track in ("hairpin_loops", "internal_loops", "junction_loops_J3", "junction_loops_J4"):
            d = p[track]["denominators"]
            assert d["n_observable_and_evaluated"] == d["n_positive_observable_and_evaluated"] + d["n_not_annotated_by_source"]
            assert d["n_observable_instances"] == d["n_observable_and_evaluated"] + d["n_source_unavailable"]
            assert d["n_observable_instances"] <= 84
        d = p["ligand_binding_sites"]["denominators"]
        assert d["n_observable_and_evaluated"] == d["n_positive_observable_and_evaluated"] + d["n_not_annotated_by_source"]
    scaffold = tr["positions"][0]      # consensus 1 is only observed in the 9NBH/9NDD scaffold constructs
    assert scaffold["observable"]["n_instances"] == 8 and scaffold["hairpin_loops"]["frequency"] == 0.0


def test_consensus_53_overlapping_loop_classes(data):
    """Validated audit values (experiments/loop_overlap_pos53): HL 66, IL 13, J3 76, J4 0, JL 76, HL∩JL 62, any 80, evaluable 80 of 84."""
    p = data["tracks"]["positions"][52]
    assert p["hairpin_loops"]["n_instances"] == 66 and p["internal_loops"]["n_instances"] == 13
    assert p["junction_loops_J3"]["n_instances"] == 76 and p["junction_loops_J4"]["n_instances"] == 0 and p["junction_loops_union"]["n_instances"] == 76
    assert p["hairpin_loops"]["denominators"]["n_observable_and_evaluated"] == 80 and data["tracks"]["n_instances"] == 84
    by = {"HL": set(), "IL": set(), "J3": set(), "J4": set()}
    for L in data["loops"]["loops"]:
        for n in L["nucleotides"]:
            if n["rfam_consensus_pos"] == 53 and n["mapping_status"] == "mapped" and L["loop_type"] in by:
                by[L["loop_type"]].add(n["seq_name"])
    assert {k: len(v) for k, v in by.items()} == {"HL": 66, "IL": 13, "J3": 76, "J4": 0}
    jl = by["J3"] | by["J4"]
    assert len(by["HL"] & jl) == 62 and len(by["HL"] | by["IL"] | jl) == 80
    assert p["junction_loops_J3"]["n_loop_observations"] == 93    # alternative J3 boundaries in the same chain: observations > instances
