"""Milestone 3A: map BGSU helix/loop annotations and PDBe bound molecules onto canonical residue ids and Rfam consensus positions."""
import csv, json, re, collections, glob
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
FAM = ROOT.parent / "data" / "families" / "RF00167"
C = ROOT / "cache"
inst = json.load(open(FAM / "instances.json"))
by_pdb_auth = {(r["pdb_id"], r["auth_asym_id"]): r for r in inst}
by_pdb_entity = collections.defaultdict(list)
for r in inst: by_pdb_entity[(r["pdb_id"], r["entity_id"])].append(r)
maps = {r["seq_name"]: json.load(open(FAM / r["mapping_file"])) for r in inst}
PP_HIGH = {"*"}; PP_MED = {"9", "8"}

def residue_index(m):
    by_auth = {}
    for res in m["residues"]:
        by_auth[(str(res["author_residue_number"]), res["author_insertion_code"] or "")] = res
    return by_auth, {res["polymer_seq_position"]: res for res in m["residues"]}

def classify_residue(res, m):
    """Canonical residue -> consensus status."""
    if res is None: return "unresolved_identifier", None
    st = res["alignment_state"]
    if st == "match":
        pp = res["posterior_prob_char"]
        conf = "high" if pp in PP_HIGH else "medium" if pp in PP_MED else "low"
        return ("mapped" if res["has_coordinates"] else "mapped_no_coordinates"), {"consensus": res["rfam_consensus_pos"], "pp": pp, "confidence": conf}
    if st in ("insert", "local_end"): return st, {"anchor": res["insertion_anchor"], "offset": res["insertion_offset"]}
    return "outside_mapped_region", None

issues, helix_rows, loop_rows, lig_rows = [], [], [], []
cov = {}
# ---------------- 1. BGSU helix / element annotations ----------------
helix = collections.defaultdict(list)   # seq_name -> rows
stat = collections.Counter(); cons_set = set(); entries = set(); chains = set(); unresolved = 0; missing = []
for r in inst:
    fn = C / "bgsu" / f"nucleotide_annotation_{r['pdb_id'].upper()}_1_{r['auth_asym_id']}.txt"
    if not fn.exists() or fn.stat().st_size == 0: missing.append(r["seq_name"]); continue
    m = maps[r["seq_name"]]; by_auth, _ = residue_index(m)
    n = 0
    for line in open(fn):
        if not line.strip() or line.startswith("#"): continue
        unit, label = line.rstrip("\n").split("\t", 1)
        f = unit.split("|"); auth_num, ins = f[4], (f[7] if len(f) > 7 else "")
        labels = [x.strip() for x in re.split(r";|\|\|", label)]   # multiple assignments would be separated; none observed
        res = by_auth.get((auth_num, ins))
        status, info = classify_residue(res, m)
        stat[status] += 1; n += 1
        row = {"source": "bgsu_helix", "pdb_id": r["pdb_id"], "auth_asym_id": r["auth_asym_id"], "label_asym_id": r["label_asym_id"], "unit_id": unit,
               "label_seq_id": res["polymer_seq_position"] if res else "", "author_residue_number": auth_num, "nucleotide": f[3], "element": labels[0], "all_elements": "; ".join(labels),
               "n_assignments": len(labels), "status": status, "rfam_consensus_pos": info.get("consensus", "") if info else "", "pp": info.get("pp", "") if info else "", "confidence": info.get("confidence", "") if info else "",
               "insertion_anchor": info.get("anchor", "") if info else "", "reference_chain": "4FEO|1|B (RF00167_annotations.txt)"}
        helix[r["seq_name"]].append(row)
        if status == "mapped": cons_set.add(info["consensus"]); entries.add(r["pdb_id"]); chains.add(r["seq_name"])
        if status in ("unresolved_identifier",): unresolved += 1
        if status != "mapped": issues.append({**{k: row[k] for k in ("source", "pdb_id", "auth_asym_id", "unit_id", "label_seq_id", "status")}, "detail": json.dumps(info) if info else ""})
    if n == 0: missing.append(r["seq_name"])
cov["bgsu_helix"] = {"pdb_entries_with_annotation": len(entries), "instances_with_annotation": len(chains), "annotated_nucleotides": sum(stat.values()), "unique_consensus_positions": len(cons_set),
                     "unresolved_identifiers": unresolved, "missing_source_records(instances)": len(missing), "status_breakdown": dict(stat), "missing_instances": missing}
helix_rows = [row for rows in helix.values() for row in rows]
# ---------------- 2. BGSU loops / Motif Atlas ----------------
loops = list(csv.DictReader(open(C / "bgsu" / "loops_all.csv")))
lstat = collections.Counter(); lcons = collections.defaultdict(set); lentries = collections.defaultdict(set); lchains = collections.defaultdict(set); lunres = 0; lmotifs = collections.defaultdict(set)
loops_in_family_chains = 0; loops_other_chains = 0
for L in loops:
    units = L["units"].split(";")
    chains_in = {u.split("|")[2] for u in units}
    t = L["loop_type"]
    seg_rows = []
    for u in units:
        f = u.split("|"); pdb, chain, auth_num, ins = f[0].lower(), f[2], f[4], (f[7] if len(f) > 7 else "")
        r = by_pdb_auth.get((pdb, chain))
        if not r: seg_rows.append((u, None, "chain_not_in_rf00167_instances", None)); continue
        m = maps[r["seq_name"]]; by_auth, _ = residue_index(m)
        res = by_auth.get((auth_num, ins)); status, info = classify_residue(res, m)
        seg_rows.append((u, r, status, info))
    if all(s[1] is None for s in seg_rows): loops_other_chains += 1; continue
    loops_in_family_chains += 1
    if L["motif_group"]: lmotifs[t].add(L["motif_group"])
    for u, r, status, info in seg_rows:
        lstat[(t, status)] += 1
        if status == "unresolved_identifier": lunres += 1
        if status == "mapped": lcons[t].add(info["consensus"]); lentries[t].add(L["pdb_id"]); lchains[t].add(r["seq_name"])
        loop_rows.append({"source": "bgsu_loop", "loop_id": L["loop_id"], "loop_type": t, "motif_group": L["motif_group"], "annotation_text": L["annotation_text"], "pdb_id": L["pdb_id"],
                          "auth_asym_id": u.split("|")[2], "label_asym_id": r["label_asym_id"] if r else "", "unit_id": u, "label_seq_id": info and None or "", "status": status,
                          "rfam_consensus_pos": info.get("consensus", "") if info else "", "pp": info.get("pp", "") if info else "", "confidence": info.get("confidence", "") if info else "", "n_nucleotides_in_loop": len(units), "multi_chain_loop": len(chains_in) > 1})
        if r and status != "mapped": issues.append({"source": "bgsu_loop", "pdb_id": L["pdb_id"], "auth_asym_id": u.split("|")[2], "unit_id": u, "label_seq_id": "", "status": status, "detail": f"{L['loop_id']} {json.dumps(info) if info else ''}"})
# fix label_seq_id in loop rows
for row in loop_rows:
    r = by_pdb_auth.get((row["pdb_id"], row["auth_asym_id"]))
    if r:
        by_auth, _ = residue_index(maps[r["seq_name"]]); f = row["unit_id"].split("|"); res = by_auth.get((f[4], f[7] if len(f) > 7 else ""))
        row["label_seq_id"] = res["polymer_seq_position"] if res else ""
cov["bgsu_loops"] = {"loops_total_in_entries": len(loops), "loops_on_rf00167_chains": loops_in_family_chains, "loops_on_other_chains": loops_other_chains,
                     "by_type": {t: {"loops": sum(1 for L in loops if L["loop_type"] == t), "with_motif_group": sum(1 for L in loops if L["loop_type"] == t and L["motif_group"]), "without_motif_group": sum(1 for L in loops if L["loop_type"] == t and not L["motif_group"]),
                                    "distinct_motif_groups": len(lmotifs[t]), "pdb_entries": len(lentries[t]), "instances": len(lchains[t]), "unique_consensus_positions": len(lcons[t])} for t in ("HL", "IL", "J3", "J4")},
                     "annotated_nucleotides(loop memberships)": sum(lstat.values()), "unresolved_identifiers": lunres, "status_breakdown": {f"{t}:{s}": n for (t, s), n in sorted(lstat.items())},
                     "pdb_entries_with_loops": len({L["pdb_id"] for L in loops}), "missing_source_records(entries)": 66 - len({L["pdb_id"] for L in loops})}
# ---------------- 3. PDBe bound molecules ----------------
WATER = {"HOH", "DOD"}
METALS = {"MG", "NA", "K", "CA", "MN", "ZN", "CO", "NI", "CD", "SR", "BA", "CS", "RB", "TL", "PB", "IR", "FE", "CU", "LI", "HG", "AU", "PT", "AG", "OS", "SM", "LU", "GD", "EU", "YB", "TB", "3CO", "NCO", "IRI", "6MO"}
ADDITIVES = {"CL", "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "P6G", "MPD", "ACT", "FMT", "DMS", "BME", "TRS", "EPE", "MES", "CAC", "NH4", "ACY", "IPA", "BR", "IOD", "SCN", "NO3", "SPM", "SPD", "SPK", "BU3", "PDO", "MLI", "CIT", "OXL", "TAR", "BTB", "IMD", "GSH", "MRD", "2PE", "HEZ", "PE4", "PEG"}
def lig_class(acc, name):
    if acc in WATER: return "water"
    if acc in METALS: return "metal_ion"
    if acc in ADDITIVES: return "additive"
    return "ligand"
lcount = collections.Counter(); lentries2 = collections.defaultdict(set); linst = collections.defaultdict(set); lcons2 = collections.defaultdict(set); lclass = {}; lnames = {}; lstat2 = collections.Counter(); lunres2 = 0; entries_with_sites = set(); entries_without = []
for pdb in sorted({r["pdb_id"] for r in inst}):
    d = json.load(open(C / "pdbe" / f"binding_sites_{pdb}.json"))
    ent = d.get(pdb)
    if not ent or not ent.get("data"): entries_without.append(pdb); continue
    entries_with_sites.add(pdb)
    rows_for_entity = by_pdb_entity.get((pdb, 1), [])
    for site in ent["data"]:
        acc, name = site["accession"], site["name"]; cls = lig_class(acc, name); lclass[acc] = cls; lnames[acc] = name
        best = site.get("additionalData", {}).get("bestChainId"); entity = site.get("additionalData", {}).get("entityId")
        chain_rows = [r for r in rows_for_entity if r["label_asym_id"] == best] or rows_for_entity[:1]
        for seg in site["residues"]:
            for idx in range(seg["startIndex"], seg["endIndex"] + 1):
                for r in chain_rows:
                    _, by_pos = residue_index(maps[r["seq_name"]]); res = by_pos.get(idx)
                    status, info = classify_residue(res, maps[r["seq_name"]])
                    # sequence-level PDBe index: the residue exists in sequence; coordinates may be absent in this chain
                    lstat2[(cls, status)] += 1
                    if status == "unresolved_identifier": lunres2 += 1
                    if status.startswith("mapped"): lcons2[acc].add(info["consensus"])
                    lcount[acc] += 1; lentries2[acc].add(pdb); linst[acc].add(r["seq_name"])
                    lig_rows.append({"source": "pdbe_binding_sites", "pdb_id": pdb, "entity_id": entity, "best_label_asym_id": best, "auth_asym_id": r["auth_asym_id"], "label_asym_id": r["label_asym_id"],
                                     "ccd_id": acc, "ligand_name": name, "ligand_class": cls, "bound_molecule_id": seg.get("additionalData", {}).get("boundMoleculeId", ""),
                                     "label_seq_id": idx, "index_type": seg["indexType"], "sequence_code": seg["startCode"] if idx == seg["startIndex"] else "", "status": status,
                                     "rfam_consensus_pos": info.get("consensus", "") if info else "", "pp": info.get("pp", "") if info else "", "confidence": info.get("confidence", "") if info else "", "observation_level": "entity (bestChainId reported; residues pooled per ligand)"})
                    if status not in ("mapped",): issues.append({"source": "pdbe_binding_sites", "pdb_id": pdb, "auth_asym_id": r["auth_asym_id"], "unit_id": f"{pdb}|{entity}|{best}|{acc}|seq{idx}", "label_seq_id": idx, "status": status, "detail": json.dumps(info) if info else ""})
cov["pdbe_bound_molecules"] = {"pdb_entries_with_binding_sites": len(entries_with_sites), "pdb_entries_without": entries_without, "distinct_ccd_ids": len(lcount), "by_class": dict(collections.Counter(lclass.values())),
                               "ligands(default display)": {acc: {"name": lnames[acc], "pdb_entries": len(lentries2[acc]), "instances(best chain)": len(linst[acc]), "binding_residue_rows": lcount[acc], "consensus_positions": sorted(lcons2[acc])} for acc in sorted(lcount) if lclass[acc] == "ligand"},
                               "excluded_from_default": {acc: {"class": lclass[acc], "pdb_entries": len(lentries2[acc])} for acc in sorted(lcount) if lclass[acc] != "ligand"},
                               "status_breakdown": {f"{c}:{s}": n for (c, s), n in sorted(lstat2.items())}, "unresolved_identifiers": lunres2,
                               "unique_consensus_positions(ligand class)": len({c for acc in lcons2 if lclass[acc] == "ligand" for c in lcons2[acc]})}
# ---------------- write ----------------
def dump(name, rows):
    if not rows: return
    keys = []; [keys.append(k) for r in rows for k in r if k not in keys]
    with open(ROOT / name, "w", newline="") as f: w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
dump("helix_annotations_sample.csv", [r for r in helix_rows if r["pdb_id"] in ("1y26", "4fe5", "9ndd", "6p2h", "3rkf")])
dump("helix_annotations_all.csv", helix_rows)
dump("loop_annotations_sample.csv", [r for r in loop_rows if r["pdb_id"] in ("1y26", "4fe5", "9ndd", "6p2h", "3rkf")])
dump("loop_annotations_all.csv", loop_rows)
dump("ligand_annotations_sample.csv", [r for r in lig_rows if r["pdb_id"] in ("1y26", "4fe5", "9ndd", "6p2h", "3rkf", "4lx5")])
dump("ligand_annotations_all.csv", lig_rows)
dump("identifier_mapping_issues.csv", issues)
covrows = []
for src, v in cov.items():
    flat = {k: (json.dumps(x) if isinstance(x, (dict, list)) else x) for k, x in v.items()}
    covrows.append({"source": src, **flat})
dump("annotation_coverage.csv", covrows)
json.dump(cov, open(ROOT / "annotation_coverage.json", "w"), indent=1)
print(json.dumps(cov, indent=1)[:6000])
print("issues by (source,status):", collections.Counter((i["source"], i["status"]) for i in issues))
