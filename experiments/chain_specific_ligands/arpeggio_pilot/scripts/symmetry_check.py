"""Document contacts potentially omitted because the pilot uses deposited asymmetric-unit coordinates without assembly or symmetry expansion.

Two independent checks per entry:
 A. Geometric: gemmi NeighborSearch with the crystal cell/space group. For every nonpolymer ligand atom, RNA atoms of a *symmetry image*
    (image_idx != 0) within 5 A are reported as contacts that exist in the crystal lattice but not in the asymmetric unit (X-ray only;
    cryo-EM entries have no cell).
 B. Source comparison: PDBe graph-api bound_molecule_interactions (computed on PDBe's preferred assembly, symmetry mates included) partner
    residues that the local Arpeggio run does not list for the same ligand instance; partner chains absent from the deposited file are
    symmetry/assembly copies. Ligand instances PDBe did not process (no bm id: outside the preferred assembly) are listed as covered only locally.

Run with the Arpeggio venv python:  <venv>/bin/python3 scripts/symmetry_check.py --family RF00167
Output: results/symmetry_check_<family>.json and .csv
"""
from __future__ import annotations
import argparse, csv, json
from collections import defaultdict
from pathlib import Path
import gemmi

PILOT = Path(__file__).resolve().parents[1]
WATER = {"HOH", "DOD"}


def geometric(pdb: str, rec: dict, radius: float = 5.0) -> dict:
    st = gemmi.read_structure(str(PILOT / "cif" / f"{pdb.upper()}.cif")); st.setup_entities()
    if not st.cell.a or not st.spacegroup_hm or rec["coordinate_context"].get("experimental_method", "").upper() != "X-RAY DIFFRACTION":
        return {"applicable": False, "reason": f"not X-ray ({rec['coordinate_context'].get('experimental_method')}); cell {st.cell.a} {st.spacegroup_hm}"}
    model = st[0]
    ns = gemmi.NeighborSearch(model, st.cell, radius).populate(include_h=False)
    targets = {c["auth_asym_id"] for c in rec["target_chains"]}
    ligs = {(l["auth_asym_id"], l["auth_seq_num"], l["ins_code"]): l for l in rec["ligand_instances"]}
    hits: dict[tuple, dict] = {}
    for ch in model:
        for res in ch:
            L = ligs.get((ch.name, res.seqid.num, res.seqid.icode.strip()))
            if not L:
                continue
            for atom in res:
                for m in ns.find_atoms(atom.pos, "\0", radius=radius):
                    if m.image_idx == 0:
                        continue
                    cra = m.to_cra(model)
                    if cra.chain.name not in targets or cra.residue.name in WATER or cra.residue.het_flag == "H":
                        continue
                    d = st.cell.find_nearest_image(atom.pos, cra.atom.pos).dist()
                    k = (L["label_asym_id"], cra.chain.name, cra.residue.seqid.num, cra.residue.seqid.icode.strip())
                    h = hits.setdefault(k, {"ligand_ccd": L["ccd_id"], "ligand_class": L["ligand_class"], "ligand_label_asym_id": L["label_asym_id"], "ligand_auth": f"{L['auth_asym_id']}/{L['auth_seq_num']}",
                                            "rna_auth_asym_id": cra.chain.name, "rna_auth_seq_num": cra.residue.seqid.num, "rna_comp_id": cra.residue.name, "min_distance": d, "n_atom_pairs": 0, "image_idx": m.image_idx})
                    h["n_atom_pairs"] += 1; h["min_distance"] = min(h["min_distance"], d)
    for h in hits.values():
        h["min_distance"] = round(h["min_distance"], 2)
    return {"applicable": True, "spacegroup": st.spacegroup_hm, "radius": radius, "n_symmetry_contacts": len(hits), "contacts": sorted(hits.values(), key=lambda h: (h["ligand_label_asym_id"], h["rna_auth_asym_id"], h["rna_auth_seq_num"]))}


def pdbe_comparison(pdb: str, rec: dict) -> dict:
    P = PILOT / "pdbe"
    try:
        bms = json.load(open(P / f"bound_molecules_{pdb}.json"))[pdb]
        mons = json.load(open(P / f"ligand_monomers_{pdb}.json"))[pdb]
        summ = json.load(open(P / f"summary_{pdb}.json"))[pdb][0]
    except Exception as e:
        return {"available": False, "reason": str(e)}
    preferred = next((a["assembly_id"] for a in summ.get("assemblies", []) if a.get("preferred")), None)
    deposited_chains = {c["auth_asym_id"] for c in rec["target_chains"]}
    lig_keys = {(l["auth_asym_id"], l["auth_seq_num"]) for l in rec["ligand_instances"]}
    local = defaultdict(set)
    for x in rec["residue_contacts"]:
        local[(x["ligand_auth_asym_id"], x["ligand_auth_seq_num"])].add((x["rna_auth_asym_id"], x["rna_auth_seq_num"]))
    out = {"available": True, "preferred_assembly": preferred, "n_ligand_instances_total": len(mons), "n_processed_by_pdbe": sum(1 for m in mons if m["bm_id"]),
           "unprocessed_by_pdbe": [f"{m['chem_comp_id']} {m['chain_id']}/{m['author_residue_number']}" for m in mons if not m["bm_id"]], "per_bm": []}
    for b in bms:
        lig = b["composition"]["ligands"][0]
        try:
            d = json.load(open(P / f"bmi_{pdb}_{b['bm_id']}.json"))[pdb][0]
        except Exception:
            out["per_bm"].append({"bm_id": b["bm_id"], "error": "no interactions file"}); continue
        partners = {(i["end"]["chain_id"], i["end"]["author_residue_number"]) for i in d["interactions"] if i["end"]["chem_comp_id"] not in WATER and i["end"]["chain_id"] in deposited_chains
                    and (i["end"]["chain_id"], i["end"]["author_residue_number"]) not in lig_keys}
        foreign = {(i["end"]["chain_id"], i["end"]["author_residue_number"], i["end"]["chem_comp_id"]) for i in d["interactions"] if i["end"]["chem_comp_id"] not in WATER and i["end"]["chain_id"] not in deposited_chains}
        loc = local.get((lig["chain_id"], lig["author_residue_number"]), set())
        out["per_bm"].append({"bm_id": b["bm_id"], "ligand": f"{lig['chem_comp_id']} {lig['chain_id']}/{lig['author_residue_number']}", "n_pdbe_rna_partners": len(partners), "n_local_rna_partners": len(loc),
                              "pdbe_not_local": sorted(f"{c}/{n}" for c, n in partners - loc), "local_not_pdbe": sorted(f"{c}/{n}" for c, n in loc - partners),
                              "pdbe_partners_on_non_target_chains": sorted(f"{c}/{n} {r}" for c, n, r in foreign)})
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--family", required=True); a = ap.parse_args()
    R = PILOT / "results"
    report, rows = {}, []
    for f in sorted(R.glob("*.json")):
        if f.name.startswith(a.family) or f.name.startswith("symmetry"):
            continue
        rec = json.load(open(f)); pdb = rec["pdb_id"]
        g = geometric(pdb, rec); p = pdbe_comparison(pdb, rec)
        report[pdb] = {"geometric_symmetry_contacts": g, "pdbe_comparison": p}
        for h in g.get("contacts", []):
            rows.append({"pdb_id": pdb, "check": "geometric_symmetry_image", **h})
        for b in p.get("per_bm", []):
            for x in b.get("pdbe_not_local", []):
                rows.append({"pdb_id": pdb, "check": "pdbe_partner_missing_locally", "ligand_auth": b["ligand"], "rna_auth": x})
    n_x = sum(1 for v in report.values() if v["geometric_symmetry_contacts"].get("applicable"))
    summary = {"n_entries": len(report), "n_xray_entries_checked": n_x,
               "n_entries_with_symmetry_contacts": sum(1 for v in report.values() if v["geometric_symmetry_contacts"].get("n_symmetry_contacts")),
               "n_symmetry_contacts": sum(v["geometric_symmetry_contacts"].get("n_symmetry_contacts", 0) for v in report.values()),
               "symmetry_contacts_by_class": dict(sum((__import__("collections").Counter(h["ligand_class"] for h in v["geometric_symmetry_contacts"].get("contacts", [])) for v in report.values()), __import__("collections").Counter())),
               "n_ligand_instances_total": sum(v["pdbe_comparison"].get("n_ligand_instances_total", 0) for v in report.values()),
               "n_ligand_instances_unprocessed_by_pdbe": sum(len(v["pdbe_comparison"].get("unprocessed_by_pdbe", [])) for v in report.values()),
               "n_pdbe_partners_missing_locally": sum(len(b.get("pdbe_not_local", [])) for v in report.values() for b in v["pdbe_comparison"].get("per_bm", [])),
               "n_local_partners_not_in_pdbe": sum(len(b.get("local_not_pdbe", [])) for v in report.values() for b in v["pdbe_comparison"].get("per_bm", []))}
    json.dump({"summary": summary, "entries": report}, open(R / f"symmetry_check_{a.family}.json", "w"), indent=1)
    if rows:
        keys = []; [keys.append(k) for r in rows for k in r if k not in keys]
        with open(R / f"symmetry_check_{a.family}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
