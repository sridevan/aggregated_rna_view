"""Write isolated single-chain PDB inputs for fpocketR from the cached PDBe updated mmCIF files, with a reversible identifier map.
Chain = all polymer residues of the instance's label_asym_id (observed atoms only, construct extensions kept, modified residues kept as ATOM
records), heavy atoms, altloc = highest occupancy per atom name (ties: first altloc label), model 1. Optional 'holo' variant appends the
ligand residue(s) of the given CCD that belong to the same entry (needed only for fpocketR's own Known/Novel classification via --ligand).
Project python (gemmi) is used read-only; nothing in the project is modified."""
import csv, json, sys
from pathlib import Path
import gemmi
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from rna_explorer.cache import Cache
from rna_explorer.pdbe import get_updated_cif
OUT = Path(__file__).resolve().parents[1] / "inputs"; OUT.mkdir(exist_ok=True)
inst = {r["seq_name"]: r for r in json.load(open(ROOT / "data/families/RF00167/instances.json"))}
cache = Cache(ROOT / "cache")

def pick_atoms(res):
    by = {}
    for a in res:
        if a.element.name.upper() in ("H", "D"): continue
        by.setdefault(a.name, []).append(a)
    return [sorted(v, key=lambda a: (-a.occ, a.altloc))[0] for v in by.values()], sum(len(v) - 1 for v in by.values())

def write(seq_name, ligand_ccd=None):
    r = inst[seq_name]
    st = gemmi.read_structure(str(get_updated_cif(cache, r["pdb_id"]))); st.setup_entities(); model = st[0]
    sub = model.get_subchain(r["label_asym_id"])
    out = gemmi.Structure(); out.name = r["pdb_id"].upper(); m = gemmi.Model("1"); ch = gemmi.Chain(r["auth_asym_id"])
    rows = []; dropped = 0; n_atoms = 0
    for res in sub:
        atoms, d = pick_atoms(res); dropped += d
        nr = gemmi.Residue(); nr.name = res.name; nr.seqid = res.seqid; nr.het_flag = "A"   # polymer residues (incl. modified) as ATOM records
        for a in atoms:
            na = gemmi.Atom(); na.name = a.name; na.element = a.element; na.pos = a.pos; na.occ = a.occ; na.b_iso = a.b_iso; na.altloc = "\0"; nr.add_atom(na); n_atoms += 1
        ch.add_residue(nr)
        rows.append({"pdb_id": r["pdb_id"], "seq_name": seq_name, "pdb_chain": r["auth_asym_id"], "pdb_resnum": res.seqid.num, "pdb_icode": res.seqid.icode.strip(), "resname": res.name,
                     "label_asym_id": r["label_asym_id"], "label_seq_id": res.label_seq, "n_heavy_atoms": len(atoms), "altloc_copies_dropped": d})
    m.add_chain(ch)
    lig_info = None
    if ligand_ccd:
        # only the ligand copy closest to this chain (multi-copy entries have one ligand per RNA copy); its original chain/number are recorded
        chain_pos = [a.pos for res in sub for a in res]
        cands = [(min(a.pos.dist(cp) for a in res for cp in chain_pos), c.name, res) for c in model for res in c if res.name == ligand_ccd and res.het_flag == "H"]
        cands.sort(key=lambda x: x[0])
        lch = gemmi.Chain(r["auth_asym_id"]); n = 0
        if cands:
            dmin, cname, res = cands[0]
            atoms, _ = pick_atoms(res); nr = gemmi.Residue(); nr.name = res.name; nr.seqid = res.seqid; nr.het_flag = "H"
            for a in atoms:
                na = gemmi.Atom(); na.name = a.name; na.element = a.element; na.pos = a.pos; na.occ = a.occ; na.b_iso = a.b_iso; nr.add_atom(na)
            lch.add_residue(nr); n = 1
            lig_info = {"ccd": ligand_ccd, "n_copies_in_entry": len(cands), "copy_used": f"{cname}/{res.seqid.num}", "min_distance_to_chain_A": round(dmin, 2), "n_copies_added": n}
        m.add_chain(lch)
    out.add_model(m); out.setup_entities()
    tag = f"{seq_name}{'_holo_' + ligand_ccd if ligand_ccd else ''}"
    pdb_path = OUT / f"{tag}.pdb"; out.write_pdb(str(pdb_path))
    with open(OUT / f"{tag}_residue_map.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    # validation: re-read the PDB and compare residue numbering / atom counts
    chk = gemmi.read_structure(str(pdb_path)); cch = chk[0][r["auth_asym_id"]]
    back = [(res.seqid.num, res.seqid.icode.strip(), res.name, len(res)) for res in cch if res.het_flag != "H" or ligand_ccd is None]
    back = [b for b in back if b[2] != ligand_ccd]
    ok = [(x["pdb_resnum"], x["pdb_icode"], x["resname"], x["n_heavy_atoms"]) for x in rows] == back
    rep = {"seq_name": seq_name, "file": pdb_path.name, "n_residues": len(rows), "n_heavy_atoms": n_atoms, "altloc_copies_dropped": dropped, "modified_residues": sorted({x["resname"] for x in rows if x["resname"] not in ("A", "C", "G", "U")}),
           "resnum_range": [rows[0]["pdb_resnum"], rows[-1]["pdb_resnum"]], "label_seq_range": [rows[0]["label_seq_id"], rows[-1]["label_seq_id"]], "roundtrip_ok": ok, "ligand": lig_info}
    print(json.dumps(rep)); return rep

if __name__ == "__main__":
    reps = [write("4fe5_1_B_A"), write("1y26_1_X_A"), write("5e54_1_A_A"), write("3ivn_1_A_A"), write("9ndd_1_A_A"), write("1y27_1_X_A"),
            write("4fe5_1_B_A", "HPA"), write("1y26_1_X_A", "ADE"), write("9ndd_1_A_A", "OXG"), write("1y27_1_X_A", "GUN")]
    json.dump(reps, open(OUT / "input_validation.json", "w"), indent=1)
