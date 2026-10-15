"""Add Rfam hit confidence (bits, E-value, model span) and SSU-clan competition to rf00177_instances.csv."""
import csv, json, sys, collections
from pathlib import Path
DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
CLAN = {"RF00177": "SSU_rRNA_bacteria", "RF01959": "SSU_rRNA_archaea", "RF01960": "SSU_rRNA_eukarya", "RF02542": "SSU_rRNA_microsporidia", "RF02545": "SSU_trypano_mito"}
hits = collections.defaultdict(list)   # (pdb, auth_chain) -> [(fam, bits, evalue, seq_start, seq_end, cm_start, cm_end, source)]
# Primary source: Rfam family 'structures' endpoint for every SSU-clan family (current, includes entries newer than the FTP table).
for fam in CLAN:
    for m in json.load(open(DATA/f"api/{fam.lower()}_rfam_structures.json"))["mapping"]:
        hits[(m["pdb_id"].lower(), m["chain"])].append((fam, float(m["bit_score"]), m["evalue_score"], m["pdb_start"], m["pdb_end"], m["cm_start"], m["cm_end"], "rfam_structures_endpoint"))
rows = list(csv.DictReader(open(RES/"rf00177_instances.csv")))
n_not_best, n_missing = 0, 0
for r in rows:
    h = hits.get((r["pdb_id"], r["auth_asym_id"]), [])
    mine = [x for x in h if x[0] == "RF00177"]
    best = max(h, key=lambda x: x[1]) if h else None
    r["rfam_bits"] = mine[0][1] if mine else ""; r["rfam_evalue"] = mine[0][2] if mine else ""
    r["rfam_cm_start"] = mine[0][5] if mine else ""; r["rfam_cm_end"] = mine[0][6] if mine else ""
    r["rfam_seq_start"] = mine[0][3] if mine else ""; r["rfam_seq_end"] = mine[0][4] if mine else ""
    r["ssu_clan_families_hit"] = ";".join(f"{x[0]}:{x[1]:.0f}" for x in sorted(h, key=lambda x: -x[1]))
    r["ssu_clan_best_family"] = best[0] if best else ""
    r["rf00177_is_clan_best"] = int(bool(best) and best[0] == "RF00177") if h else ""
    if not mine: n_missing += 1
    elif best[0] != "RF00177": n_not_best += 1
with open(RES/"rf00177_instances.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("instances:", len(rows), "| RF00177 hit absent from Rfam pdb_full_region:", n_missing, "| RF00177 not clan-best:", n_not_best)
print("--- not clan-best (should be excluded from RF00177 aggregation):")
for r in rows:
    if r["rf00177_is_clan_best"] == 0: print(f"  {r['pdb_id']} {r['auth_asym_id']} {r['organism']} | {r['molecule_name'][:40]} | {r['ssu_clan_families_hit']}")
print("--- clan-best but with competing clan hits:")
for r in rows:
    if r["rf00177_is_clan_best"] == 1 and ";" in r["ssu_clan_families_hit"]: print(f"  {r['pdb_id']} {r['auth_asym_id']} {r['organism']} | {r['ssu_clan_families_hit']}")
bits = sorted((float(r["rfam_bits"]), r["pdb_id"], r["auth_asym_id"], r["organism"], r["mapped_length"]) for r in rows if r["rfam_bits"] != "")
print("--- lowest 12 RF00177 bit scores:"); [print("  ", b) for b in bits[:12]]
print("--- bit score distribution:", {k: sum(1 for b in bits if lo <= b[0] < hi) for k, (lo, hi) in {"<100":(0,100),"100-300":(100,300),"300-600":(300,600),"600-1200":(600,1200),">=1200":(1200,1e9)}.items()})
s = json.load(open(RES/"rf00177_instances_summary.json")); s["rfam_hit_absent"] = n_missing; s["rf00177_not_ssu_clan_best"] = n_not_best
s["rfam_bits_distribution"] = {k: sum(1 for b in bits if lo <= b[0] < hi) for k, (lo, hi) in {"<100":(0,100),"100-300":(100,300),"300-600":(300,600),"600-1200":(600,1200),">=1200":(1200,1e9)}.items()}
json.dump(s, open(RES/"rf00177_instances_summary.json", "w"), indent=1)
