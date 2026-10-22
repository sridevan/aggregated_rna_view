"""Reported sequence differences (mmCIF _struct_ref_seq_dif) per chain, placed on Rfam positions through the validated mapping."""
import json
from pathlib import Path
import gemmi
import pytest
from rna_explorer.sequence_differences import change_kind, copy_consistency, extract_chain_differences, place_on_family
from rna_explorer.export import validate_json

ROOT = Path(__file__).resolve().parent.parent

# entity 1: RNA, two copies (strands A, B -> label asym A, B), referenced to GenBank; entity 2: protein (strand P -> label C)
CIF = """data_TEST
loop_
_entity_poly.entity_id
_entity_poly.type
1 polyribonucleotide
2 'polypeptide(L)'
loop_
_struct_ref.id
_struct_ref.entity_id
_struct_ref.db_name
_struct_ref.db_code
_struct_ref.pdbx_db_accession
1 1 GB REF1 X12345.1
2 2 UNP PROT_X P99999
loop_
_struct_ref_seq.align_id
_struct_ref_seq.ref_id
_struct_ref_seq.pdbx_strand_id
1 1 A
2 1 B
3 2 P
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
_pdbx_poly_seq_scheme.pdb_strand_id
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.pdb_ins_code
A 1 1 G A 10 .
A 1 2 G A 11 .
A 1 3 U A 12 A
A 1 4 C A 13 .
B 1 1 G B 10 .
B 1 2 G B 11 .
B 1 3 U B 12 A
B 1 4 C B 13 .
C 2 1 ALA P 1 .
loop_
_struct_ref_seq_dif.align_id
_struct_ref_seq_dif.pdbx_pdb_strand_id
_struct_ref_seq_dif.seq_num
_struct_ref_seq_dif.pdbx_auth_seq_num
_struct_ref_seq_dif.pdbx_pdb_ins_code
_struct_ref_seq_dif.mon_id
_struct_ref_seq_dif.db_mon_id
_struct_ref_seq_dif.pdbx_seq_db_seq_num
_struct_ref_seq_dif.details
_struct_ref_seq_dif.pdbx_ordinal
1 A 1 10 ? G ? ? 'Engineered mutation' 1
1 A 3 12 A U C 30 'Engineered mutation' 2
2 B 1 10 ? G ? ? 'Engineered mutation' 3
2 B 3 12 A U C 30 'Engineered mutation' 4
1 A ? ? ? ? U 31 deletion 5
1 A 4 13 ? A C 32 Conflict 6
1 Z 2 11 ? G A 29 variant 7
3 P 1 1 ? ALA GLY 1 'Expression tag' 8
"""


@pytest.fixture
def rows():
    return extract_chain_differences(gemmi.cif.read_string(CIF).sole_block(), "1TST")["rows"]


def test_rows_resolve_to_entity_and_label_chain_with_verbatim_classification(rows):
    by = {r["ordinal"]: r for r in rows}
    sub = by["2"]
    assert (sub["entity_id"], sub["polymer_type"], sub["auth_asym_id"], sub["label_asym_id"], sub["label_seq_id"]) == ("1", "polyribonucleotide", "A", "A", 3)
    assert (sub["reference"], sub["deposited"], sub["change_kind"], sub["classification"]) == ("C", "U", "substitution", "Engineered mutation")
    assert (sub["author_residue_number"], sub["author_insertion_code"], sub["scheme_check"]) == ("12", "A", "agrees")
    assert sub["reference_sequence"] == {"db_name": "GB", "db_code": "REF1", "accession": "X12345.1"}
    assert by["4"]["label_asym_id"] == "B"                                                   # the second copy is its own chain
    assert (by["1"]["change_kind"], by["1"]["reference"]) == ("added_residue", None)
    assert (by["5"]["change_kind"], by["5"]["label_seq_id"], by["5"]["scheme_check"]) == ("absent_residue", None, "no_deposited_residue")
    assert by["6"]["scheme_check"] == "identity_conflict"                                    # file says A, the deposited sequence has C
    assert by["7"]["chain_status"] == "chain_unresolved" and by["7"]["label_asym_id"] is None   # strand Z is not in the entity: never guessed
    assert by["8"]["polymer_type"] == "polypeptide(L)" and by["8"]["classification"] == "Expression tag"   # kept with its type; filtered by the family layer
    assert change_kind("U", "C") == "substitution" and change_kind("G", None) == "added_residue" and change_kind(None, "U") == "absent_residue"


def _res(state, pos, name="U", auth="12", ins="A"):
    return {"alignment_state": state, "rfam_consensus_pos": pos, "residue_name": name, "author_residue_number": auth, "author_insertion_code": ins, "has_coordinates": True}


def test_placement_only_when_every_check_passes(rows):
    by = {r["ordinal"]: r for r in rows}
    assert place_on_family(by["2"], _res("match", 80)) | {} == {"alignment_state": "match", "rfam_consensus_pos": 80, "has_coordinates": True, "mapping_status": "mapped"}
    assert place_on_family(by["2"], _res("insert", None))["mapping_status"] == "no_consensus_position:insert"
    assert place_on_family(by["2"], _res("not_submitted_flank", None))["rfam_consensus_pos"] is None
    assert place_on_family(by["2"], None)["mapping_status"] == "not_in_mapping"
    assert place_on_family(by["2"], _res("match", 80, name="C"))["mapping_status"] == "mapping_identity_conflict"
    assert place_on_family(by["2"], _res("match", 80, auth="13"))["mapping_status"] == "mapping_author_conflict"
    assert place_on_family(by["5"], None)["mapping_status"] == "no_deposited_residue"
    assert place_on_family(by["6"], _res("match", 81, name="C", auth="13", ins=""))["mapping_status"] == "identity_conflict"


def test_copies_of_one_entity_must_agree():
    class Idx:
        instances = [{"pdb_id": "1tst", "entity_id": 1, "label_asym_id": "A"}, {"pdb_id": "1tst", "entity_id": 1, "label_asym_id": "B"}]
    o = lambda lab, seq: {"pdb_id": "1tst", "entity_id": "1", "label_asym_id": lab, "label_seq_id": seq, "reference_seq_num": "30", "deposited": "U", "reference": "C", "classification": "Engineered mutation"}
    assert copy_consistency([o("A", 3), o("B", 3)], Idx()) == []
    assert copy_consistency([o("A", 3)], Idx()) == [{"pdb_id": "1tst", "entity_id": "1", "rows_per_chain": {"A": 1, "B": 0}}]


FAM = ROOT / "data/families/RF00167"


@pytest.mark.skipif(not (FAM / "sequence_differences.json").exists(), reason="RF00167 sequence differences not built")
def test_rf00167_regression():
    d = json.load(open(FAM / "sequence_differences.json")); t = json.load(open(FAM / "sequence_difference_tracks.json"))
    validate_json(d, "sequence_differences"); validate_json(t, "sequence_difference_tracks")
    c = d["counts"]
    assert (c["n_entries_evaluated"], c["n_observations"], c["n_on_consensus_axis"], c["n_off_axis"], c["n_unresolved_rows"]) == (66, 85, 77, 8, 0)
    assert c["classifications"] == {"Engineered mutation": 35, "SEE REMARK 999": 4, "Conflict": 6, "variant": 40} and d["inconsistent_copies"] == []
    # 3FO4 C74U: the guanine->adenine specificity switch, at consensus 80
    o = [x for x in d["observations"] if x["pdb_id"] == "3fo4"]
    assert [(x["reference"], x["deposited"], x["author_residue_number"], x["rfam_consensus_pos"], x["classification"]) for x in o] == [("C", "U", "74", 80, "Engineered mutation")]
    # 6VWT: both copies carry the same 10 rows, each on its own chain
    v = [x for x in d["observations"] if x["pdb_id"] == "6vwt"]
    assert {x["label_asym_id"] for x in v} == {"A", "B"} and sum(x["label_asym_id"] == "A" for x in v) == 10
    p80 = t["positions"][79]["classifications"]["Engineered mutation"]
    assert p80["n_entries"] == 3 and p80["changes"] == {"C>U": 2, "U>C": 1}                     # 2B57 and 3FO4 C74U; 4LX6 U74C (reverse direction)
