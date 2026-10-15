"""Offline regression: rebuild chain mappings from the POC Stockholm/chains files and compare with the POC CSVs (RF00167 and RF00177)."""
import csv, json, re
from pathlib import Path
import pytest
from rna_explorer.mapping import build_chain_mapping, map_aligned_sequence
from rna_explorer.mmcif import atom_site_observed
from rna_explorer.sequences import ChainInstance
from rna_explorer.stockholm import parse_stockholm, type_columns

ROOT = Path(__file__).resolve().parent.parent
CASES = {
    "RF00167": (ROOT / "poc_rf00167/results/rf00167_mapped.sto", ROOT / "poc_rf00167/results/chains.json", ROOT / "poc_rf00167/data", ROOT / "poc_rf00167/results/mapped_residue_mapping.csv", "poc167"),
    "RF00177": (ROOT / "poc_rf00177/intermediate/rf00177_mapped.sto", ROOT / "poc_rf00177/intermediate/chains.json", ROOT / "poc_rf00177/data/cif", ROOT / "poc_rf00177/results/rf00177_residue_mapping.csv", "poc177"),
}


def _instance(c: dict, cif_dir: Path) -> ChainInstance:
    return ChainInstance(pdb_id=c["pdb_id"], entity_id=int(c["entity_id"]), label_asym_id=c["label_asym_id"], auth_asym_id=c["auth_chain_id"], map_start=c["map_start"], map_end=c["map_end"],
                         map_auth_start=c["map_auth_start"], map_auth_end=c["map_auth_end"], sequence=c["sequence"], pdb_sequence=c["pdb_sequence"], residues=c["residues"],
                         modified_positions=c.get("modified_positions", {}), observed=atom_site_observed(cif_dir / f"{c['pdb_id']}_updated.cif", c["label_asym_id"], int(c["entity_id"])),
                         api_organism=c.get("api_organism"), cif_organism=c.get("cif_organism"), cif_tax_id=c.get("cif_tax_id"))


def _poc_rows(path: Path, flavour: str) -> dict:
    out = {}
    for r in csv.DictReader(open(path)):
        key = (r["pdb_id"], int(r["polymer_seq_position"]))
        if flavour == "poc167":
            m = re.match(r"(insert|local_end):after_cons_(\d+):#(\d+)", r["insertion_info"] or "")
            out[key] = {"state": r["status"], "cons": r["rfam_consensus_pos"], "anchor": m.group(2) if m else "", "offset": m.group(3) if m else "", "pp": r["posterior_prob"],
                        "auth": r["author_residue_number"], "ins": r["author_insertion_code"], "resn": r["residue_name"], "nt": r["nucleotide"]}
        else:
            out[key] = {"state": r["alignment_state"], "cons": r["rfam_consensus_pos"], "anchor": r["insertion_anchor"], "offset": r["insertion_offset"], "pp": r["posterior_prob"],
                        "auth": r["author_residue_number"], "ins": r["author_insertion_code"], "resn": r["residue_name"], "nt": r["nucleotide_parent"], "coords": r["has_coordinates"]}
    return out


@pytest.mark.parametrize("acc", list(CASES))
def test_regression_against_poc(acc: str) -> None:
    sto, chains_json, cif_dir, poc_csv, flavour = CASES[acc]
    if not (sto.exists() and chains_json.exists() and poc_csv.exists()):
        pytest.skip("POC files not present")
    aln = parse_stockholm(sto)
    cm = type_columns(aln.gc["RF"], aln.gc.get("SS_cons", ""))
    poc = _poc_rows(poc_csv, flavour)
    n = 0
    for c in json.load(open(chains_json)):
        inst = _instance(c, cif_dir)
        name = f"{inst.seq_name}/{inst.map_start}-{inst.map_end}"
        aligned = map_aligned_sequence(name, aln.seqs[name], aln.gr[name]["PP"], cm)
        m = build_chain_mapping(inst, aligned, cm, acc)
        for r in m.rows:
            p = poc[(r.pdb_id, r.polymer_seq_position)]
            got = {"state": r.alignment_state, "cons": "" if r.rfam_consensus_pos is None else str(r.rfam_consensus_pos), "anchor": "" if r.insertion_anchor is None else str(r.insertion_anchor),
                   "offset": "" if r.insertion_offset is None else str(r.insertion_offset), "pp": r.posterior_prob, "auth": str(r.author_residue_number), "ins": r.author_insertion_code,
                   "resn": r.residue_name, "nt": r.nucleotide_parent}
            if "coords" in p:
                got["coords"] = str(r.has_coordinates)
            assert got == {k: p[k] for k in got}, (acc, r.pdb_id, r.polymer_seq_position)
            n += 1
        assert m.summary["n_match"] == sum(1 for r in m.rows if r.alignment_state == "match")
    assert n == len(poc)
