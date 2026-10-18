"""Compare chain-specific RNA-ligand contact sources for RF00167 test entries (1Y26, 3RKF, 9NDD, 9NBH).

Sources (all responses cached in ../cache/, retrieved 2026-10-15):
  v2        PDBe api/v2/pdb/entry/binding_sites/{pdb}/{entity}   (Milestone 3B; entity level, label_seq_id indices, bound-molecule ids)
  bm        PDBe graph-api/pdb/bound_molecule_interactions/{pdb}/{bm_id}   (Arpeggio, residue level, interaction types, author ids)
  lig       PDBe graph-api/pdb/bound_ligand_interactions/{pdb}/{chain}/{resnum}   (Arpeggio, atom level with distances, author ids)
  rcsb      RCSB data API core/polymer_entity_instance/{PDB}/{asym}  rcsb_ligand_neighbors   (distance based <= 5 A, label ids)
  arp       local pdbe-arpeggio 1.4.4 on the PDBe updated mmCIF (asymmetric unit), INTER contacts, author ids

Every record is resolved through the validated family mapping (rna_explorer.annotation_mapping.InstanceIndex, read-only) to
(pdb_id, label_asym_id, label_seq_id) and the Rfam consensus position. Output: ../CONTACT_COMPARISON.csv and a summary printed to stdout.
"""
from __future__ import annotations
import csv, glob, json, sys
from collections import defaultdict
from pathlib import Path
HERE = Path(__file__).resolve().parent
EXP = HERE.parent
ROOT = EXP.parent.parent
sys.path.insert(0, str(ROOT))
from rna_explorer.annotation_mapping import InstanceIndex   # noqa: E402

C = EXP / "cache"
idx = InstanceIndex(ROOT / "data/families/RF00167")
ENTRIES = ["1y26", "3rkf", "9ndd", "9nbh"]
WATER = {"HOH", "DOD"}


def jload(name):
    return json.load(open(C / name))


def mapping_by_auth(seq_name):
    idx.mapping(seq_name)
    return idx._auth_idx[seq_name]


def resolve_auth(pdb, auth_chain, auth_num, ins=""):
    """author chain + number -> ResidueRef (via the BGSU-style unit id resolver)."""
    unit = f"{pdb.upper()}|1|{auth_chain}|X|{auth_num}" + (f"|||{ins}" if ins and ins.strip() else "")
    return idx.resolve_unit(unit)


rows: list[dict] = []
ligand_instances: dict[tuple, dict] = {}     # (pdb, label_asym of ligand) -> info


def add(pdb, source, lig_ccd, lig_auth_chain, lig_auth_num, lig_label_asym, bm_id, ref, detail: dict, attribution="chain"):
    rows.append({"pdb_id": pdb, "source": source, "ligand_ccd": lig_ccd, "ligand_auth_chain": lig_auth_chain, "ligand_auth_resnum": lig_auth_num,
                 "ligand_label_asym_id": lig_label_asym, "bound_molecule_id": bm_id, "attribution": attribution,
                 "rna_seq_name": ref.seq_name, "rna_label_asym_id": ref.label_asym_id, "rna_auth_asym_id": ref.auth_asym_id, "rna_label_seq_id": ref.label_seq_id,
                 "rna_auth_resnum": ref.author_residue_number, "rna_ins_code": ref.insertion_code, "rna_residue": ref.residue_name, "mapping_status": ref.mapping_status,
                 "rfam_consensus_pos": ref.rfam_consensus_pos, "pp_char": ref.pp_char, **detail})


for pdb in ENTRIES:
    P = pdb.upper()
    # ---- ligand instances (PDBe ligand_monomers: struct_asym_id == label_asym_id, bm_id or None)
    mons = jload(f"www.ebi.ac.uk_pdbe_api_pdb_entry_ligand_monomers_{pdb}.json")[pdb]
    by_auth_lig = {(m["chain_id"], m["author_residue_number"]): m for m in mons}
    by_asym_lig = {m["struct_asym_id"]: m for m in mons}
    for m in mons:
        ligand_instances[(pdb, m["struct_asym_id"])] = {"ccd": m["chem_comp_id"], "auth_chain": m["chain_id"], "auth_num": m["author_residue_number"], "bm_id": m["bm_id"],
                                                       "entity_id": m["entity_id"], "altconf": m["alternate_conformers"]}
    # ---- v2 binding sites (entity level)
    v2 = json.load(open(ROOT / "annotation_assessment/cache/pdbe" / f"binding_sites_{pdb}.json"))[pdb]
    chains = idx.by_entity[(pdb, 1)]
    for site in v2["data"]:
        bms = sorted({r["additionalData"]["boundMoleculeId"] for r in site["residues"]})
        for seg in site["residues"]:
            for i in range(seg["startIndex"], seg["endIndex"] + 1):
                ref = idx.resolve_seq_index(chains[0]["seq_name"], i)
                # entity level: chain unknown -> blank chain fields, keep best chain as metadata
                ref2 = ref.to_dict(); ref2.update(seq_name=None, label_asym_id=None, auth_asym_id=None)
                from rna_explorer.annotation_mapping import ResidueRef
                r = ResidueRef(**ref2)
                add(pdb, "pdbe_binding_sites_v2", site["accession"], None, None, None, seg["additionalData"]["boundMoleculeId"], r,
                    {"interaction_types": "", "distance": "", "ligand_atoms": "", "rna_atoms": "", "note": f"entity-level; bestChainId={site['additionalData'].get('bestChainId')}; bms={','.join(bms)}"},
                    attribution="entity")
    # ---- graph-api bound_molecule_interactions (residue level)
    bms = jload(f"www.ebi.ac.uk_pdbe_graph-api_pdb_bound_molecules_{pdb}.json")[pdb]
    for b in bms:
        lig = b["composition"]["ligands"][0]
        d = jload(f"www.ebi.ac.uk_pdbe_graph-api_pdb_bound_molecule_interactions_{pdb}_{b['bm_id']}.json")[pdb][0]
        lm = by_auth_lig.get((lig["chain_id"], lig["author_residue_number"]), {})
        for it in d["interactions"]:
            e = it["end"]
            if e["chem_comp_id"] in WATER:
                continue
            ref = resolve_auth(pdb, e["chain_id"], e["author_residue_number"], e.get("author_insertion_code", ""))
            types = sorted({t for k, v in it["interactions"].items() for t in v})
            add(pdb, "pdbe_graphapi_bm_interactions", lig["chem_comp_id"], lig["chain_id"], lig["author_residue_number"], lm.get("struct_asym_id"), b["bm_id"], ref,
                {"interaction_types": ";".join(types), "distance": "", "ligand_atoms": "", "rna_atoms": "", "note": "residue-level Arpeggio; preferred assembly only"})
    # ---- graph-api bound_ligand_interactions (atom level) for the files we fetched
    for f in glob.glob(str(C / f"www.ebi.ac.uk_pdbe_graph-api_pdb_bound_ligand_interactions_{pdb}_*.json")):
        d = json.load(open(f))
        if pdb not in d:
            continue
        rec = d[pdb][0]
        lg = rec["ligand"]
        lm = by_auth_lig.get((lg["chain_id"], lg["author_residue_number"]), {})
        per_res: dict[tuple, dict] = {}
        for it in rec["interactions"]:
            e = it["end"]
            if e["chem_comp_id"] in WATER:
                continue
            k = (e["chain_id"], e["author_residue_number"], e.get("author_insertion_code", ""))
            x = per_res.setdefault(k, {"types": set(), "dmin": None, "latoms": set(), "ratoms": set(), "n": 0})
            x["types"].update(it["interaction_details"]); x["n"] += 1
            x["dmin"] = it["distance"] if x["dmin"] is None else min(x["dmin"], it["distance"])
            x["latoms"].update(it["ligand_atoms"]); x["ratoms"].update(e["atom_names"] or [])
        for (ch, num, ins), x in per_res.items():
            ref = resolve_auth(pdb, ch, num, ins)
            add(pdb, "pdbe_graphapi_ligand_interactions", lg["chem_comp_id"], lg["chain_id"], lg["author_residue_number"], lm.get("struct_asym_id"), lm.get("bm_id"), ref,
                {"interaction_types": ";".join(sorted(x["types"])), "distance": x["dmin"], "ligand_atoms": ";".join(sorted(x["latoms"])), "rna_atoms": ";".join(sorted(x["ratoms"])),
                 "note": f"atom-level Arpeggio, {x['n']} atom pairs; preferred assembly only"})
    # ---- RCSB ligand neighbors (polymer instance, label ids)
    for inst in idx.by_pdb[pdb]:
        f = C / f"data.rcsb.org_rest_v1_core_polymer_entity_instance_{P}_{inst['label_asym_id']}.json"
        if not f.exists():
            continue
        d = json.load(open(f)).get("rcsb_ligand_neighbors", [])
        per: dict[tuple, dict] = {}
        for n in d:
            if n["ligand_comp_id"] in WATER:
                continue
            k = (n["ligand_asym_id"], n["seq_id"])
            x = per.setdefault(k, {"comp": n["ligand_comp_id"], "dmin": n["distance"], "latoms": set(), "ratoms": set(), "n": 0, "bound": n.get("ligand_is_bound")})
            x["dmin"] = min(x["dmin"], n["distance"]); x["latoms"].add(n["ligand_atom_id"]); x["ratoms"].add(n["atom_id"]); x["n"] += 1
        for (lasym, seq), x in per.items():
            ref = idx.resolve_seq_index(inst["seq_name"], seq)
            lm = by_asym_lig.get(lasym, {})
            add(pdb, "rcsb_ligand_neighbors", x["comp"], lm.get("chain_id"), lm.get("author_residue_number"), lasym, lm.get("bm_id"), ref,
                {"interaction_types": "proximity<=5A" + (";covalent_or_metal" if x["bound"] else ""), "distance": x["dmin"], "ligand_atoms": ";".join(sorted(x["latoms"])),
                 "rna_atoms": ";".join(sorted(x["ratoms"])), "note": f"{x['n']} atom pairs; asymmetric unit"})
    # ---- local Arpeggio (author ids)
    for f in glob.glob(str(C / f"arpeggio_local_{pdb}.json")):
        d = json.load(open(f))
        per: dict[tuple, dict] = {}
        for c in d:
            if c["interacting_entities"] != "INTER":
                continue
            b, e = c["bgn"], c["end"]
            # orient: ligand side = the one whose (chain,num) is a ligand monomer
            if (b["auth_asym_id"], b["auth_seq_id"]) in by_auth_lig:
                lig, rna = b, e
            elif (e["auth_asym_id"], e["auth_seq_id"]) in by_auth_lig:
                lig, rna = e, b
            else:
                continue
            if rna["label_comp_id"] in WATER:
                continue
            k = (lig["auth_asym_id"], lig["auth_seq_id"], rna["auth_asym_id"], rna["auth_seq_id"], rna["pdbx_PDB_ins_code"].strip())
            x = per.setdefault(k, {"ccd": lig["label_comp_id"], "types": set(), "dmin": None, "latoms": set(), "ratoms": set(), "n": 0})
            x["types"].update(c["contact"]); x["n"] += 1
            x["dmin"] = c["distance"] if x["dmin"] is None else min(x["dmin"], c["distance"])
            x["latoms"].add(lig["auth_atom_id"]); x["ratoms"].add(rna["auth_atom_id"])
        for (lch, lnum, rch, rnum, ins), x in per.items():
            ref = resolve_auth(pdb, rch, rnum, ins)
            lm = by_auth_lig.get((lch, lnum), {})
            specific = sorted(t for t in x["types"] if t != "proximal")
            add(pdb, "arpeggio_local", x["ccd"], lch, lnum, lm.get("struct_asym_id"), lm.get("bm_id"), ref,
                {"interaction_types": ";".join(sorted(x["types"])), "distance": x["dmin"], "ligand_atoms": ";".join(sorted(x["latoms"])), "rna_atoms": ";".join(sorted(x["ratoms"])),
                 "note": f"{x['n']} atom pairs; asymmetric unit; specific={'yes' if specific else 'proximal-only'}"})

keys = list(rows[0].keys())
with open(EXP / "CONTACT_COMPARISON.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
    for r in sorted(rows, key=lambda r: (r["pdb_id"], r["ligand_ccd"], str(r["ligand_auth_chain"]), str(r["ligand_auth_resnum"]), r["source"], str(r["rna_label_asym_id"]), r["rna_label_seq_id"] or 0)):
        w.writerow({k: ("" if v is None else v) for k, v in r.items()})

# ---- summary: per ligand instance, residue sets per source
print("rows:", len(rows))
summary = defaultdict(lambda: defaultdict(set))
for r in rows:
    if r["ligand_ccd"] in ("HOH",):
        continue
    key = (r["pdb_id"], r["ligand_ccd"], r["ligand_auth_chain"], r["ligand_auth_resnum"])
    summary[key][r["source"]].add((r["rna_label_asym_id"], r["rna_label_seq_id"], r["rfam_consensus_pos"]))
for key in sorted(summary, key=lambda k: (k[0], k[1], str(k[2]), str(k[3]))):
    print(key)
    for src, s in sorted(summary[key].items()):
        print(f"   {src:34s} n={len(s):3d}  consensus={sorted({c for _, _, c in s if c})}")
