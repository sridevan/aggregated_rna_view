"""Read-only: every non-polymer copy (model 1, deposited coordinates) in each RF00167 entry, with its minimum heavy-atom distance to each
family RNA chain and the RNA residues (label_seq_id) within 4 and 5 A. Input: cached PDBe updated mmCIF; output: CSV/JSON in this folder."""
import glob, json, sys, csv
from pathlib import Path
import gemmi
ROOT = Path("/Users/sria/Desktop/aggregated_rna_view"); OUT = Path(__file__).parent.parent / "data"
acc = sys.argv[1] if len(sys.argv) > 1 else "RF00167"
inst = json.load(open(ROOT / f"data/families/{acc}/instances.json"))
by_pdb = {}
for i in inst: by_pdb.setdefault(i["pdb_id"], []).append(i)
rows, contacts = [], []
for pdb, chains in sorted(by_pdb.items()):
    f = glob.glob(str(ROOT / f"cache/www_ebi_ac_uk/pdbe_entry-files_download_{pdb}_updated.cif-*.cif"))
    if not f: rows.append({"pdb_id": pdb, "status": "no_cached_mmcif"}); continue
    doc = gemmi.cif.read(f[0]); st = gemmi.make_structure_from_block(doc[0]); st.setup_entities()
    model = st[0]
    labels = {c["label_asym_id"] for c in chains}
    rna_atoms = {}   # label_asym -> list of (pos, label_seq)
    lig = []         # (label_asym, auth_chain, auth_seq, ccd, atoms)
    for ch in model:
        for sub in ch.subchains():
            sid = sub.subchain_id()
            for r in sub:
                if sid in labels and r.entity_type == gemmi.EntityType.Polymer:
                    rna_atoms.setdefault(sid, []).extend((a.pos, r.label_seq) for a in r if a.element.name not in ("H", "D"))
                elif r.entity_type in (gemmi.EntityType.NonPolymer, gemmi.EntityType.Branched) and r.name not in ("HOH", "DOD"):
                    lig.append((sid, ch.name, r.seqid.num, r.name, [a.pos for a in r if a.element.name not in ("H", "D")]))
    for sid, ach, anum, ccd, atoms in lig:
        d = {}
        for lab, ra in rna_atoms.items():
            best = 1e9; near4, near5 = set(), set()
            for p, seq in ra:
                for q in atoms:
                    x = p.dist(q)
                    if x < best: best = x
                    if x <= 5.0: near5.add(seq)
                    if x <= 4.0: near4.add(seq)
            d[lab] = round(best, 2)
            contacts.append({"pdb_id": pdb, "ligand_label_asym": sid, "ccd": ccd, "rna_label_asym": lab, "min_dist": round(best, 2), "res_4A": sorted(near4), "res_5A": sorted(near5)})
        nearest = min(d, key=d.get) if d else None
        rows.append({"pdb_id": pdb, "status": "ok", "ligand_label_asym": sid, "ligand_auth": f"{ach}/{anum}", "ccd": ccd, "n_heavy_atoms": len(atoms),
                     "nearest_family_chain": nearest, "min_dist_nearest": d.get(nearest), "min_dist_by_chain": json.dumps(d)})
json.dump({"ligands": rows, "contacts": contacts}, open(OUT / f"{acc}_inventory.json", "w"))
print(len(rows), "ligand copies;", sum(r.get("status") == "no_cached_mmcif" for r in rows), "entries without mmCIF")
