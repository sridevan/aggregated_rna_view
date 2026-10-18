"""Map fpocketR pockets to canonical identifiers and Rfam consensus positions, compare with chain-level Arpeggio ligand contacts, and run the filter sensitivity analysis."""
import csv, json, re, sys
from collections import defaultdict
from pathlib import Path
P = Path(__file__).resolve().parents[1]; ROOT = P.parents[1]
sys.path.insert(0, str(ROOT))
from rna_explorer.annotation_mapping import InstanceIndex
idx = InstanceIndex(ROOT / "data/families/RF00167")
RUNS = ["4fe5_1_B_A", "1y26_1_X_A", "5e54_1_A_A", "3ivn_1_A_A", "1y27_1_X_A", "9ndd_1_A_A"]
HOLO = {"4fe5_1_B_A": "HPA", "1y26_1_X_A": "ADE", "1y27_1_X_A": "GUN", "9ndd_1_A_A": "OXG"}
LIG_CONTACT_CCD = {"4fe5_1_B_A": "HPA", "1y26_1_X_A": "ADE", "1y27_1_X_A": "GUN", "9ndd_1_A_A": "OXG"}

def read_info(path):
    out = {}; cur = None
    for line in open(path):
        m = re.match(r"Pocket (\d+) :", line)
        if m: cur = int(m.group(1)); out[cur] = {}; continue
        if cur and ":" in line:
            k, v = line.split(":", 1); out[cur][k.strip()] = float(v.strip()) if re.fullmatch(r"-?[\d.]+(e[+-]?\d+)?", v.strip()) else v.strip()
    return out

def lining_atoms(atm_pdb):
    res = defaultdict(list)
    for line in open(atm_pdb):
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21].strip(), int(line[22:26]), line[26].strip(), line[17:20].strip()); res[key].append(line[12:16].strip())
    return res

def analyze(run_name, seq_name, variant):
    d = P / "runs" / f"{run_name}_clean_out"
    csvp = d / f"{run_name}_out_pocket_characteristics.csv"; info = read_info(d / f"{run_name}_clean_info.txt")
    rmap = {(r["pdb_chain"], int(r["pdb_resnum"]), r["pdb_icode"]): r for r in csv.DictReader(open(P / "inputs" / f"{seq_name}_residue_map.csv"))}
    inst = idx.by_seq[seq_name]
    pockets = []
    rows = list(csv.DictReader(open(csvp))) if csvp.exists() else []
    for r in rows:
        n = int(r["Pocket"]); la = lining_atoms(d / "pockets" / f"pocket{n}_atm.pdb")
        nts = []
        for (ch, num, ic, rn), atoms in sorted(la.items(), key=lambda x: x[0][1]):
            m = rmap.get((ch, num, ic)); ref = idx.resolve_unit(f"{inst['pdb_id'].upper()}|1|{ch}|{rn}|{num}" + (f"|||{ic}" if ic else ""))
            nts.append({"pdb_chain": ch, "pdb_resnum": num, "icode": ic, "resname": rn, "label_asym_id": m["label_asym_id"] if m else None, "label_seq_id": int(m["label_seq_id"]) if m else None,
                        "resolver_label_seq_id": ref.label_seq_id, "rfam_consensus_pos": ref.rfam_consensus_pos, "alignment_state": ref.alignment_state, "mapping_status": ref.mapping_status, "n_lining_atoms": len(atoms), "lining_atoms": atoms})
        pocket_nt_csv = json.loads(r["PocketNT"])
        pockets.append({"pocket_id": f"{seq_name}:{variant}:pocket{n}", "rank": n, "score": float(r["Score"]), "drug_score": float(r["Drug_score"]), "volume_A3": float(r["Volume"]), "n_alpha_spheres": int(r["a-sphere"]), "pocket_sasa_A2": float(r["SASA"]),
                        "hydrophobicity_score": float(r["Hydrophobicity_score"]), "polarity_score": float(r["Polarity_score"]), "npr1": float(r["Pocket_NPR1"]) if r["Pocket_NPR1"] else None, "npr2": float(r["Pocket_NPR2"]) if r["Pocket_NPR2"] else None, "shape": r["Pocket_shape"],
                        "filter_pass": r["Filter"] == "Pass", "type": r["Type"], "ligand_id": r["Ligand_ID"] or None, "pocket_overlap": float(r["Pocket_overlap"]) if r["Pocket_overlap"] else None, "ligand_overlap": float(r["Ligand_overlap"]) if r["Ligand_overlap"] else None, "center_criteria_A": float(r["Center_criteria"]) if r["Center_criteria"] else None,
                        "fpocket_info": info.get(n, {}), "pocket_nt_author": pocket_nt_csv, "pocket_nt_from_atm_pdb": sorted({x["pdb_resnum"] for x in nts}), "nt_match": sorted(pocket_nt_csv) == sorted({x["pdb_resnum"] for x in nts}),
                        "consensus_positions": sorted({x["rfam_consensus_pos"] for x in nts if x["rfam_consensus_pos"] is not None}), "n_insert_or_flank_nt": sum(1 for x in nts if x["rfam_consensus_pos"] is None), "nucleotides": nts,
                        "map_consistent": all(x["label_seq_id"] == x["resolver_label_seq_id"] for x in nts)})
    return {"run": run_name, "seq_name": seq_name, "variant": variant, "n_pockets": len(pockets), "pockets": pockets}

def ligand_contacts(seq_name, ccd):
    f = ROOT / "experiments/chain_specific_ligands/arpeggio_pilot/results" / f"{idx.by_seq[seq_name]['pdb_id']}.json"
    if not f.exists(): return None
    rec = json.load(open(f)); c = [x for x in rec["residue_contacts"] if x["seq_name"] == seq_name and x["ligand_ccd"] == ccd]
    return {"any": sorted({x["rfam_consensus_pos"] for x in c if x["rfam_consensus_pos"]}), "specific": sorted({x["rfam_consensus_pos"] for x in c if x["rfam_consensus_pos"] and x["specific"]}), "n_residues": len(c)}

results = []
for sn in RUNS:
    results.append(analyze(sn, sn, "apo_input"))
    if sn in HOLO: results.append(analyze(f"{sn}_holo_{HOLO[sn]}", sn, f"holo_input_{HOLO[sn]}"))
# comparison with experimental contacts (chain-level Arpeggio on deposited coordinates)
comp = []
for res in results:
    if res["variant"] != "apo_input": continue
    sn = res["seq_name"]; ccd = LIG_CONTACT_CCD.get(sn); lc = ligand_contacts(sn, ccd) if ccd else None
    for pk in res["pockets"]:
        if lc:
            pset = set(pk["consensus_positions"]); s = set(lc["specific"]); a = set(lc["any"])
            comp.append({"seq_name": sn, "pocket": pk["rank"], "score": pk["score"], "pocket_positions": sorted(pset), "contacts_specific": lc["specific"], "contacts_any": lc["any"],
                         "jaccard_specific": round(len(pset & s) / len(pset | s), 3) if pset | s else None, "recall_specific": round(len(pset & s) / len(s), 3) if s else None, "precision_vs_any": round(len(pset & a) / len(pset), 3) if pset else None})
        else:
            comp.append({"seq_name": sn, "pocket": pk["rank"], "score": pk["score"], "pocket_positions": pk["consensus_positions"], "contacts_specific": None, "note": "apo structure (no default-track ligand)"})
    if not res["pockets"]: comp.append({"seq_name": sn, "pocket": None, "note": "no pocket predicted", "contacts_specific": lc["specific"] if lc else None})
# sensitivity analysis over the apo-input runs (one indicator per instance per consensus position)
def positives(res, mode):
    pk = [p for p in res["pockets"] if p["filter_pass"]]
    if mode == "top1": pk = pk[:1]
    if mode == "score>=0.4": pk = [p for p in pk if p["score"] >= 0.4]
    return set().union(*[set(p["consensus_positions"]) for p in pk]) if pk else set(), len(pk)
known = {25, 26, 51, 55, 56, 80, 81}   # family-level pocket from PDBe entity-level contacts (>=54 of 64 entries)
sens = {}
for mode in ("all_score>=0", "top1", "score>=0.4"):
    per = {}; freq = defaultdict(int); evaluable = defaultdict(int)
    for res in results:
        if res["variant"] != "apo_input": continue
        pos, n = positives(res, mode); per[res["seq_name"]] = {"n_pockets_retained": n, "n_consensus_positions": len(pos), "positions": sorted(pos), "overlap_with_known_site": len(pos & known)}
        sn = res["seq_name"]
        for p in range(1, idx.clen + 1):
            if idx.observable(sn, p): evaluable[p] += 1
            if p in pos: freq[p] += 1
    sens[mode] = {"per_instance": per, "positions_with_any_occurrence": sorted(freq), "occurrence_counts": {str(p): f"{freq[p]}/{evaluable[p]}" for p in sorted(freq)}}
out = {"method": {"fpocketR": "1.3.6", "fpocket": "4.0.3 (conda-forge h8a118fb_2, x86_64 under Rosetta 2)", "parameters": "-m 3.0 -M 5.7 -i 42 -D 1.65 -A 3 -p 0.0 (published RNA-optimised set; fpocketR defaults)", "pocket_nucleotides": "residues of fpocket pockets/pocketN_atm.pdb (atoms contacted by the pocket's Voronoi vertices); no extra distance cut-off",
                  "known_site_reference": "chain-level Arpeggio contacts on deposited coordinates (experiments/chain_specific_ligands/arpeggio_pilot) and the family-level PDBe pocket positions 25,26,51,55,56,80,81"},
       "runs": results, "comparison_with_contacts": comp, "sensitivity": sens}
json.dump(out, open(P / "pockets_pilot.json", "w"), indent=1)
for res in results: print(res["run"], "pockets", res["n_pockets"], [(p["rank"], p["score"], p["type"], p["consensus_positions"], "nt_match" if p["nt_match"] else "NT MISMATCH", "map_ok" if p["map_consistent"] else "MAP MISMATCH", p["n_insert_or_flank_nt"]) for p in res["pockets"]])
print(json.dumps(comp, indent=1)); print(json.dumps({k: {"per_instance": {s: (v["n_pockets_retained"], v["n_consensus_positions"], v["overlap_with_known_site"]) for s, v in m["per_instance"].items()}, "occurrence": m["occurrence_counts"]} for k, m in sens.items()}, indent=1))
