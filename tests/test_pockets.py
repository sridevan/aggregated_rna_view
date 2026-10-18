"""Unit tests for rna_explorer/pockets.py on synthetic fpocket/fpocketR outputs (no fpocketR execution, no family dataset)."""
import csv
from pathlib import Path
import pytest
from rna_explorer.pockets import FILTERS, aggregate, lining_atoms, map_pockets, parse_run, qualifying, read_info


def _write_run(d: Path, name: str, pockets: list[dict]):
    """pockets: [{rank, score, nts: [(chain, resnum, resname, [atoms])]}]"""
    (d / "pockets").mkdir(parents=True, exist_ok=True)
    cols = ["Parameters", "Name", "PDB", "State", "Type", "Filter", "Pocket", "Score", "Drug_score", "a-sphere", "SASA", "Volume", "Hydrophobic_density", "Apolar_a-sphere_proportion", "Hydrophobicity_score", "Polarity_score",
            "PocketNT", "Pocket_NPR1", "Pocket_NPR2", "Pocket_shape", "Ligand_ID", "Pocket_overlap", "Ligand_overlap", "Center_criteria", "QED_score", "Ligand_NPR1", "Ligand_NPR2", "Ligand_shape"]
    with open(d / f"{name}_out_pocket_characteristics.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for p in pockets:
            nts = sorted({n[1] for n in p["nts"]})
            w.writerow({"Parameters": "-m 3.0 -M 5.7 -i 42 -D 1.65 -A 3 -p 0.0", "Name": name, "PDB": name[:4], "State": "", "Type": "Novel", "Filter": "Pass" if p["score"] > 0 else "Fail", "Pocket": p["rank"], "Score": p["score"], "Drug_score": 0.001,
                        "a-sphere": 45, "SASA": 30.0, "Volume": 140.0, "Hydrophobic_density": 0.0, "Apolar_a-sphere_proportion": 0.0, "Hydrophobicity_score": 1.5, "Polarity_score": 3, "PocketNT": str(nts), "Pocket_NPR1": 0.6, "Pocket_NPR2": 0.8, "Pocket_shape": "Balanced",
                        "Ligand_ID": "", "Pocket_overlap": "", "Ligand_overlap": "", "Center_criteria": "", "QED_score": "", "Ligand_NPR1": "", "Ligand_NPR2": "", "Ligand_shape": ""})
    with open(d / f"{name}_clean_info.txt", "w") as f:
        for p in pockets:
            f.write(f"Pocket {p['rank']} :\n\tScore : \t{p['score']:.3f}\n\tDruggability Score : \t0.001\n\tNumber of Alpha Spheres : \t45\n\tTotal SASA : \t30.000\n\tVolume : \t140.000\n\n")
    for p in pockets:
        with open(d / "pockets" / f"pocket{p['rank']}_atm.pdb", "w") as f:
            f.write("HEADER fake\n"); k = 1
            for ch, num, rn, atoms in p["nts"]:
                for a in atoms:
                    f.write(f"ATOM  {k:5d} {a:<4s} {rn:>3s} {ch}{num:4d}    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           {a[0]}\n"); k += 1


class FakeRef:
    def __init__(self, label_seq_id, cons, state="match"):
        self.label_seq_id = label_seq_id; self.rfam_consensus_pos = cons; self.alignment_state = state; self.mapping_status = "mapped" if label_seq_id is not None else "unresolved_identifier"


class FakeIdx:
    """author number n -> label_seq_id n-10; consensus = label_seq_id + 2 for label 1..6 (match), label 7 insert; label 4 has no coordinates."""
    clen = 10
    def resolve_unit(self, unit):
        num = int(unit.split("|")[4]); lab = num - 10
        if lab < 1 or lab > 7: return FakeRef(None, None, None)
        if lab == 7: return FakeRef(7, None, "insert")
        return FakeRef(lab, lab + 2)
    def match_residue(self, sn, pos): return {"x": 1} if 3 <= pos <= 8 else None
    def observable(self, sn, pos): return 3 <= pos <= 8 and pos != 6   # label 4 (consensus 6) lacks coordinates


RMAP = [{"pdb_chain": "A", "pdb_resnum": 10 + i, "pdb_icode": "", "resname": "G", "label_asym_id": "A", "label_seq_id": i} for i in range(1, 8)]
INST = {"seq_name": "1tst_1_A_A", "pdb_id": "1tst", "label_asym_id": "A", "auth_asym_id": "A"}


def test_parse_lining_dedup_and_mapping(tmp_path):
    d = tmp_path / "1tst_1_A_A_clean_out"
    _write_run(d, "1tst_1_A_A", [{"rank": 1, "score": 0.55, "nts": [("A", 11, "G", ["N1", "C2", "O6"]), ("A", 12, "G", ["N7"]), ("A", 17, "G", ["C8"])]},
                                 {"rank": 2, "score": 0.30, "nts": [("A", 12, "G", ["O6"]), ("A", 13, "G", ["N2", "N3"])]}])
    pk = parse_run(d, "1tst_1_A_A")
    assert [p["rank"] for p in pk] == [1, 2] and pk[0]["score"] == 0.55 and pk[0]["pocket_nt_author_csv"] == [11, 12, 17]
    la = lining_atoms(d / "pockets" / "pocket1_atm.pdb"); assert la[("A", 11, "", "G")] == ["N1", "C2", "O6"]      # three atoms, one residue
    info = read_info(d / "1tst_1_A_A_clean_info.txt"); assert info[1]["Score"] == 0.55 and info[2]["Volume"] == 140.0
    mapped = map_pockets(pk, RMAP, FakeIdx(), INST, incomplete_nts={2}, n_residues=7)
    p1 = mapped[0]
    assert [n["label_seq_id"] for n in p1["nucleotides"]] == [1, 2, 7] and p1["n_nucleotides"] == 3        # one membership per nucleotide despite 3 lining atoms
    assert p1["consensus_positions"] == [3, 4] and p1["n_non_consensus_nucleotides"] == 1                    # insert kept, not forced onto a consensus position
    assert p1["nt_match_csv"] and p1["map_consistent"] and p1["unresolved_nucleotides"] == 0
    assert set(p1["flags"]) == {"includes_incomplete_nucleotide", "touches_chain_terminus", "includes_non_consensus_nucleotide"}
    assert mapped[1]["consensus_positions"] == [4, 5] and mapped[1]["flags"] == ["includes_incomplete_nucleotide", "touches_chain_terminus"]   # overlapping pockets stay distinct; label 2 is incomplete and within the terminal window
    assert p1["pocket_id"] == "1tst_1_A_A:pocket1"


def test_filters_and_aggregation_denominators():
    def rec(sn, status, pockets):
        return {"seq_name": sn, "pdb_id": sn[:4], "status": status, "pockets": pockets}
    P = lambda rank, score, cons: {"pocket_id": f"p{rank}", "rank": rank, "score": score, "consensus_positions": cons}
    recs = [rec("1aaa_1_A_A", "pockets", [P(1, 0.55, [3, 4]), P(2, 0.30, [4, 5])]),      # two pockets overlapping at 4
            rec("2bbb_1_A_A", "pockets", [P(1, 0.20, [5])]),                          # only a low-score pocket
            rec("3ccc_1_A_A", "no_pockets", []),                                      # evaluated, zero pockets
            rec("4ddd_1_A_A", "failed", [])]                                          # failed run: excluded from denominators
    assert [p["rank"] for p in qualifying(recs[0]["pockets"], "all_score_ge_0")] == [1, 2]
    assert [p["rank"] for p in qualifying(recs[0]["pockets"], "top1")] == [1]
    assert [p["rank"] for p in qualifying(recs[0]["pockets"], "score_ge_0.4")] == [1] and qualifying(recs[1]["pockets"], "score_ge_0.4") == []
    tr = aggregate(recs, FakeIdx(), 10)
    assert set(tr) == set(FILTERS)
    a = {p["rfam_consensus_pos"]: p for p in tr["all_score_ge_0"]}
    assert a[4]["n_supporting_instances"] == 1 and a[4]["n_distinct_pockets"] == 2 and a[4]["n_evaluable_instances"] == 3      # one instance counted once, two pockets
    assert a[4]["frequency"] == round(1 / 3, 4) and a[4]["unevaluable_reasons"] == {"failed_run": 1} and a[4]["score_summary"] == {"min": 0.3, "median": 0.55, "max": 0.55}
    assert a[5]["n_supporting_instances"] == 2 and a[5]["n_supporting_entries"] == 2
    assert a[6]["n_evaluable_instances"] == 0 and a[6]["frequency"] is None and a[6]["unevaluable_reasons"] == {"no_coordinates": 3, "failed_run": 1}   # no coordinates: outside the denominator, never zero
    assert a[1]["n_evaluable_instances"] == 0 and a[1]["n_unevaluable_instances"] == 0                                                           # position absent from these chains
    t = {p["rfam_consensus_pos"]: p for p in tr["top1"]}; s = {p["rfam_consensus_pos"]: p for p in tr["score_ge_0.4"]}
    assert t[5]["n_supporting_instances"] == 1 and s[5]["n_supporting_instances"] == 0 and s[5]["n_evaluable_instances"] == 3     # zero pockets after filtering still counts in the denominator
    assert all(p["n_supporting_instances"] <= p["n_evaluable_instances"] for f in tr.values() for p in f)


def test_zero_pocket_run_parses_to_empty(tmp_path):
    d = tmp_path / "x_clean_out"; _write_run(d, "x", [])
    assert parse_run(d, "x") == [] and read_info(d / "x_clean_info.txt") == {}
