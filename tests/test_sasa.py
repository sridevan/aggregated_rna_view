"""SASA calculation and aggregation tests (rna_explorer/sasa.py) on synthetic structures; the RF00167 regression lives in test_sasa_rf00167.py."""
import math
from pathlib import Path
import pytest
pytest.importorskip("freesasa"); pytest.importorskip("gemmi"); pytest.importorskip("numpy")
import gemmi
from rna_explorer.sasa import PROBE_RADIUS, RADII, N_POINTS, aggregate, compute_atom_sasa, load_chain, measure_instance, CcdExpectedAtoms, _stats


def _atom(name, el, x, y, z, occ=1.0, alt="\0"):
    a = gemmi.Atom(); a.name = name; a.element = gemmi.Element(el); a.pos = gemmi.Position(x, y, z); a.occ = occ; a.altloc = alt; a.b_iso = 20.0
    return a


def _structure(tmp_path: Path) -> Path:
    """Two chains: RNA chain A (label asym 'A', two residues) with a second altloc copy and a hydrogen; protein-like chain B (label asym 'B') far away; a water."""
    st = gemmi.Structure(); st.name = "1tst"; st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    ch = gemmi.Chain("A")
    r1 = gemmi.Residue(); r1.name = "G"; r1.seqid = gemmi.SeqId(5, " "); r1.label_seq = 1; r1.subchain = "A"; r1.entity_type = gemmi.EntityType.Polymer
    r1.add_atom(_atom("P", "P", 0, 0, 0)); r1.add_atom(_atom("OP1", "O", 1.5, 0, 0)); r1.add_atom(_atom("OP2", "O", -1.5, 0, 0)); r1.add_atom(_atom("H1", "H", 0, 0.9, 0))
    r1.add_atom(_atom("C1'", "C", 0, 0, 1.6, occ=0.7, alt="A")); r1.add_atom(_atom("C1'", "C", 0, 0, 50.0, occ=0.3, alt="B"))   # altloc B far away must be dropped
    r2 = gemmi.Residue(); r2.name = "U"; r2.seqid = gemmi.SeqId(6, " "); r2.label_seq = 2; r2.subchain = "A"; r2.entity_type = gemmi.EntityType.Polymer
    r2.add_atom(_atom("N1", "N", 30, 0, 0))   # isolated atom: SASA = 4 pi (1.55+1.4)^2
    ch.add_residue(r1); ch.add_residue(r2); model.add_chain(ch)
    chb = gemmi.Chain("B")
    rb = gemmi.Residue(); rb.name = "ALA"; rb.seqid = gemmi.SeqId(1, " "); rb.label_seq = 1; rb.subchain = "B"; rb.entity_type = gemmi.EntityType.Polymer
    rb.add_atom(_atom("CA", "C", 0.1, 0.1, 0.1))   # would bury chain A atoms if included
    chb.add_residue(rb); model.add_chain(chb)
    chw = gemmi.Chain("B"); rw = gemmi.Residue(); rw.name = "HOH"; rw.seqid = gemmi.SeqId(101, " "); rw.subchain = "C"; rw.entity_type = gemmi.EntityType.Water
    rw.add_atom(_atom("O", "O", 30.0, 2.0, 0)); chw.add_residue(rw); model.add_chain(chw)
    st.add_model(model)
    for name, et, pt in (("1", gemmi.EntityType.Polymer, gemmi.PolymerType.Rna), ("2", gemmi.EntityType.Polymer, gemmi.PolymerType.PeptideL), ("3", gemmi.EntityType.Water, gemmi.PolymerType.Unknown)):
        e = gemmi.Entity(name); e.entity_type = et; e.polymer_type = pt; e.subchains = [{"1": "A", "2": "B", "3": "C"}[name]]; st.entities.append(e)
    p = tmp_path / "1tst.cif"; st.make_mmcif_document().write_file(str(p))
    return p


def test_chain_isolation_hydrogens_altloc_and_sphere_area(tmp_path):
    chain = load_chain(_structure(tmp_path), "A")
    assert [r["residue_name"] for r in chain["residues"]] == ["G", "U"] and chain["polymer_type"] == "Rna"
    g = chain["residues"][0]
    assert sorted(a["name"] for a in g["atoms"]) == ["C1'", "OP1", "OP2", "P"]        # hydrogen removed, one C1' kept
    c1 = next(a for a in g["atoms"] if a["name"] == "C1'"); assert c1["altloc"] == "A" and c1["z"] == 1.6 and g["n_altloc_copies_dropped"] == 1   # highest occupancy wins
    areas, unsupported = compute_atom_sasa(chain["residues"])
    assert unsupported == [] and len(areas) == 5 and all(math.isfinite(a) and a >= 0 for a in areas)
    n1 = chain["residues"][1]["atoms"][0]
    exact = 4 * math.pi * (RADII["N"] + PROBE_RADIUS) ** 2                                  # isolated atom: protein chain B and water are not in the calculation
    assert abs(n1["sasa"] - exact) / exact < 0.01
    assert N_POINTS >= 100 and PROBE_RADIUS == 1.4


def test_unsupported_element_is_flagged_not_zero(tmp_path):
    p = _structure(tmp_path)
    chain = load_chain(p, "A")
    chain["residues"][1]["atoms"][0]["element"] = "XX"
    areas, unsupported = compute_atom_sasa(chain["residues"])
    assert unsupported == ["XX"] and chain["residues"][1]["atoms"][0]["sasa"] is None and len(areas) == 4


class _Idx:
    def __init__(self, mapping): self._m = mapping
    def mapping(self, _): return self._m


class _Ccd:
    def __init__(self, table): self.t = table; self._memo = {}
    def get(self, ccd): return self.t.get(ccd)


def test_measure_instance_completeness_and_no_coordinates(tmp_path):
    p = _structure(tmp_path)
    mapping = {"residues": [
        {"polymer_seq_position": 1, "author_residue_number": 5, "author_insertion_code": "", "residue_name": "G", "nucleotide_parent": "G", "is_modified": False, "alignment_state": "match", "rfam_consensus_pos": 3, "posterior_prob_char": "*", "has_coordinates": True},
        {"polymer_seq_position": 2, "author_residue_number": 6, "author_insertion_code": "", "residue_name": "U", "nucleotide_parent": "U", "is_modified": False, "alignment_state": "insert", "rfam_consensus_pos": None, "posterior_prob_char": "", "has_coordinates": True},
        {"polymer_seq_position": 3, "author_residue_number": 7, "author_insertion_code": "", "residue_name": "A", "nucleotide_parent": "A", "is_modified": False, "alignment_state": "match", "rfam_consensus_pos": 4, "posterior_prob_char": "*", "has_coordinates": False}]}
    ccd = _Ccd({"G": {"expected_atoms": ["C1'", "OP1", "OP2", "P", "N9"], "leaving_atoms": ["OP3"], "type": "RNA LINKING", "parent": "?"}, "U": {"expected_atoms": ["N1"], "leaving_atoms": ["OP3"], "type": "RNA LINKING", "parent": "?"}})
    import rna_explorer.sasa as S
    orig = S.instance_cif; S.instance_cif = lambda cache, inst, chain_files: p
    try:
        m = measure_instance(_Idx(mapping), None, {"seq_name": "1tst_1_A_A", "pdb_id": "1tst", "entity_id": 1, "label_asym_id": "A", "auth_asym_id": "A"}, ccd)
    finally:
        S.instance_cif = orig
    by = {n["label_seq_id"]: n for n in m["nucleotides"]}
    assert by[1]["status"] == "incomplete" and by[1]["missing_atoms"] == ["N9"] and by[1]["incomplete_reason"] == "missing_heavy_atoms" and by[1]["sasa_A2"] > 0
    assert by[2]["status"] == "complete" and by[2]["alignment_state"] == "insert" and by[2]["rfam_consensus_pos"] is None
    assert by[3]["status"] == "no_coordinates" and by[3]["sasa_A2"] is None and by[3]["n_observed_heavy_atoms"] == 0     # never zero
    assert abs(m["chain_total_sasa_A2"] - m["sum_nucleotide_sasa_A2"]) < 0.01 and m["status_counts"] == {"incomplete": 1, "complete": 1, "no_coordinates": 1}
    # aggregation: incomplete excluded from the primary value but counted; insert never aggregated
    pos = aggregate([m], 5)
    assert pos[2]["n_mapped_instances_with_coordinates"] == 1 and pos[2]["n_excluded_incomplete"] == 1 and pos[2]["median_sasa_A2"] is None and pos[2]["n_contributing_instances"] == 0
    assert all(p["n_contributing_instances"] == 0 for p in pos) and pos[3]["n_mapped_instances_with_coordinates"] == 0


def _meas(seq, values, status="complete"):
    return {"seq_name": seq, "nucleotides": [{"label_seq_id": i + 1, "residue_name": "A", "alignment_state": "match", "rfam_consensus_pos": i + 1, "status": status if v is not None else "no_coordinates", "sasa_A2": v} for i, v in enumerate(values)]}


def test_aggregation_one_per_instance_median_quartiles_and_reproducibility():
    ms = [_meas("i1", [10.0, 50.0, None]), _meas("i2", [20.0, 60.0, 5.0]), _meas("i3", [30.0, 70.0, 7.0]), _meas("i4", [40.0, None, 9.0]), _meas("i5", [100.0, 80.0, 11.0], status="incomplete")]
    pos = aggregate(ms, 3)
    p1 = pos[0]
    assert p1["n_contributing_instances"] == 4 and p1["n_mapped_instances_with_coordinates"] == 5 and p1["n_excluded_incomplete"] == 1
    assert p1["median_sasa_A2"] == 25.0 and p1["mean_sasa_A2"] == 25.0 and p1["q1_sasa_A2"] == 17.5 and p1["q3_sasa_A2"] == 32.5 and p1["min_sasa_A2"] == 10.0 and p1["max_sasa_A2"] == 40.0
    assert pos[1]["n_contributing_instances"] == 3 and pos[1]["median_sasa_A2"] == 60.0 and pos[1]["nucleotide_identity"] == {"A": 3}
    assert pos[2]["n_contributing_instances"] == 3 and pos[2]["n_mapped_instances_with_coordinates"] == 4   # i1 has no coordinates: not counted, not zero
    assert aggregate(ms, 3) == pos                                                                             # deterministic
    # a chain with two residues claiming the same position contributes once
    dup = {"seq_name": "d", "nucleotides": [{"label_seq_id": 1, "residue_name": "G", "alignment_state": "match", "rfam_consensus_pos": 1, "status": "complete", "sasa_A2": 1.0},
                                            {"label_seq_id": 2, "residue_name": "G", "alignment_state": "match", "rfam_consensus_pos": 1, "status": "complete", "sasa_A2": 99.0}]}
    assert aggregate([dup], 1)[0]["n_contributing_instances"] == 1 and aggregate([dup], 1)[0]["median_sasa_A2"] == 1.0
    s = _stats([1.0, 2.0, 3.0, 4.0]); assert s["median"] == 2.5 and s["q1"] == 1.75 and s["q3"] == 3.25 and abs(s["sd"] - 1.29) < 0.01


def test_ccd_expected_atoms_parsing(tmp_path):
    class FakeCache:
        def fetch_text(self, url, suffix):
            return ("data_X\n_chem_comp.type 'RNA LINKING'\n_chem_comp.mon_nstd_parent_comp_id G\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.pdbx_leaving_atom_flag\n"
                    "OP3 O Y\nP P N\nOP1 O N\nC1' C N\nH1 H N\n")
    c = CcdExpectedAtoms(FakeCache())
    e = c.get("X")
    assert e["expected_atoms"] == ["C1'", "OP1", "P"] and e["leaving_atoms"] == ["OP3"] and e["parent"] == "G"
