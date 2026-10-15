"""Normalise BGSU align_chains output (unit ids) to PDB polymer positions via mmCIF _pdbx_poly_seq_scheme; cross-check with BGSU SeqtoUnitMapping.
Unit id: PDB|model|chain|comp|number[|atom|alt|ins|sym]."""
import csv, json, sys, re, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts")); from common import cif_category
BG, CIF, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(exist_ok=True, parents=True)
CHAINS = {"7K00|1|A": ("7k00", "A"), "5J7L|1|AA": ("5j7l", "AA"), "1J5E|1|A": ("1j5e", "A"), "4Y4O|1|1a": ("4y4o", "1a"), "6ERI|1|BA": ("6eri", "BA"),
          "7P2E|1|A": ("7p2e", "A"), "5MRC|1|aa": ("5mrc", "aa"), "6YW5|1|aa": ("6yw5", "aa"), "9GYT|1|2": ("9gyt", "2"), "8APN|1|B3": ("8apn", "B3")}
# mmCIF lookup: (auth chain, auth num, ins) -> (seq_id, mon_id, asym) ; keep only polymer RNA asym matching strand id
lookup, stats = {}, {}
for ife, (pdb, chain) in CHAINS.items():
    rows = [r for r in cif_category(CIF/f"{pdb}_updated.cif", "_pdbx_poly_seq_scheme") if r["pdb_strand_id"] == chain]
    asyms = collections.Counter(r["asym_id"] for r in rows)
    d = {}
    for r in rows:
        key = (r["pdb_seq_num"], "" if r["pdb_ins_code"] in (".", "?") else r["pdb_ins_code"])
        d.setdefault(key, []).append((int(r["seq_id"]), r["mon_id"], r["asym_id"], r["auth_seq_num"] not in ("?", ".")))
    lookup[ife] = d; stats[ife] = {"scheme_rows": len(rows), "asym_ids_for_strand": dict(asyms), "ambiguous_keys": sum(1 for v in d.values() if len(v) > 1)}
def parse_unit(u):
    f = u.split("|"); return {"pdb": f[0].lower(), "model": f[1], "chain": f[2], "comp": f[3], "num": f[4], "ins": f[7] if len(f) > 7 else "", "raw": u}
def resolve(u):
    p = parse_unit(u); ife = f"{p['pdb'].upper()}|{p['model']}|{p['chain']}"
    cands = lookup.get(ife, {}).get((p["num"], p["ins"]), [])
    if not cands: return p, None, "unresolved"
    if len(cands) > 1: return p, None, "ambiguous"
    seq_id, mon, asym, modelled = cands[0]
    status = "ok" if mon == p["comp"] else "comp_mismatch"   # BGSU comp is the mmCIF comp_id (3-letter for modified)
    return p, {"polymer_pos": seq_id, "mon_id": mon, "label_asym_id": asym, "modelled": modelled}, status
# --- SeqtoUnitMapping cross-check: sequence position -> unit id (BGSU's own numbering) ---
s2u_check = {}
for ife in CHAINS:
    f = BG/f"seq2unit_{ife.replace('|','_')}.txt"
    if not f.exists(): continue
    txt = f.read_text().replace("</br>", "\n"); n_ok = n_bad = n_null = n_unres = 0; bad = []
    for line in txt.splitlines():
        m = re.match(r"(\S+)\|sequence\|(\S+?)\|(\S+)\|(\d+) observed_as (\S+)", line.strip())
        if not m: continue
        seqpos, unit = int(m.group(4)), m.group(5)
        if unit == "NULL": n_null += 1; continue
        p, res, st = resolve(unit)
        if res is None: n_unres += 1; continue
        if res["polymer_pos"] == seqpos: n_ok += 1
        else: n_bad += 1; bad.append((seqpos, unit, res["polymer_pos"]))
    s2u_check[ife] = {"seqpos_equals_mmcif_seq_id": n_ok, "differs": n_bad, "unobserved_NULL": n_null, "unresolved": n_unres, "examples_differ": bad[:3]}
# --- pairwise/multi alignment files ---
files = sorted(BG.glob("pair_*.txt")) + [BG/"multi_all10_column.txt", BG/"bgsu_column_5J7L_1J5E.txt"]
out_rows, file_stats = [], {}
for f in files:
    lines = [l.rstrip("\n").split("\t") for l in open(f) if l.strip()]
    ncol = len(lines[0]); has_col = lines[0][0].isdigit()
    chains_in_file = [None] * (ncol - (1 if has_col else 0))
    for t in lines:
        for i, u in enumerate(t[1:] if has_col else t):
            if u != "NULL" and chains_in_file[i] is None: chains_in_file[i] = "|".join(u.split("|")[:3])
    st = collections.Counter()
    for t in lines:
        col = int(t[0]) if has_col else ""
        units = t[1:] if has_col else t
        rec = {"file": f.name, "bgsu_column": col}
        for i, u in enumerate(units):
            ife = chains_in_file[i]
            if u == "NULL": rec[f"c{i}_ife"] = ife; rec[f"c{i}_unit"] = ""; rec[f"c{i}_pos"] = ""; rec[f"c{i}_status"] = "NULL"; st["NULL"] += 1; continue
            p, res, status = resolve(u); st[status] += 1
            rec[f"c{i}_ife"] = ife; rec[f"c{i}_unit"] = u; rec[f"c{i}_pos"] = res["polymer_pos"] if res else ""; rec[f"c{i}_status"] = status
            rec[f"c{i}_comp"] = p["comp"]; rec[f"c{i}_auth"] = p["num"] + p["ins"]
        out_rows.append(rec)
    file_stats[f.name] = {"rows": len(lines), "chains": chains_in_file, "unit_resolution": dict(st)}
# write long table of resolved rows (generic columns)
keys = sorted({k for r in out_rows for k in r}, key=lambda k: (k != "file", k != "bgsu_column", k))
with open(OUT/"bgsu_alignments_normalised.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(out_rows)
# pairwise correspondences table (2-chain files only)
pairs = []
for r in out_rows:
    if r["file"].startswith("multi"): continue
    if r.get("c0_status") not in ("NULL", None) and r.get("c1_status") not in ("NULL", None) and r["c0_pos"] != "" and r["c1_pos"] != "":
        pairs.append({"file": r["file"], "bgsu_column": r["bgsu_column"], "ife1": r["c0_ife"], "pdb1": r["c0_ife"].split("|")[0].lower(), "chain1": r["c0_ife"].split("|")[2], "pos1": r["c0_pos"], "unit1": r["c0_unit"], "comp1": r["c0_comp"],
                      "ife2": r["c1_ife"], "pdb2": r["c1_ife"].split("|")[0].lower(), "chain2": r["c1_ife"].split("|")[2], "pos2": r["c1_pos"], "unit2": r["c1_unit"], "comp2": r["c1_comp"]})
with open(OUT/"bgsu_correspondences_normalised.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(pairs[0].keys())); w.writeheader(); w.writerows(pairs)
json.dump({"mmcif_lookup": stats, "seq2unit_crosscheck": s2u_check, "files": file_stats}, open(OUT/"bgsu_normalisation_stats.json", "w"), indent=1)
print(json.dumps({"seq2unit_crosscheck": s2u_check}, indent=1)); print(json.dumps(file_stats, indent=1)); print("pairs:", len(pairs))
