"""(a) BGSU column -> CM consensus mapping via match residues; consistency across pairwise files. (b) Nature of BGSU pairs inside our unnumbered regions."""
import csv, json, sys, collections
from pathlib import Path
RES, MAIN = Path(sys.argv[1]), Path(sys.argv[2])
ours = {}
for f in [MAIN/"rf00177_residue_mapping.csv", RES/"repro/rf00177_residue_mapping.csv"]:
    for r in csv.DictReader(open(f)): ours[(r["pdb_id"], int(r["polymer_seq_position"]))] = r
norm = list(csv.DictReader(open(RES/"bgsu_alignments_normalised.csv")))
# (a) column consistency: for every (pdb,pos) collect the BGSU columns seen across all column-bearing files
cols = collections.defaultdict(set)
for r in norm:
    if r["bgsu_column"] == "": continue
    for i in range(10):
        if r.get(f"c{i}_pos") not in (None, ""): cols[(r[f"c{i}_ife"].split("|")[0].lower(), int(r[f"c{i}_pos"]))].add(int(r["bgsu_column"]))
incons = {k: v for k, v in cols.items() if len(v) > 1}
print("residues seen in >=2 files:", sum(1 for k, v in cols.items()), "| residues with inconsistent BGSU column across files:", len(incons), list(incons.items())[:5])
# column -> consensus via match residues of all chains
col2cons = collections.defaultdict(collections.Counter); col_states = collections.defaultdict(collections.Counter)
for (pdb, pos), cs in cols.items():
    o = ours.get((pdb, pos)); 
    if not o: continue
    for c in cs:
        col_states[c][o["alignment_state"]] += 1
        if o["alignment_state"] == "match": col2cons[c][int(o["rfam_consensus_pos"])] += 1
conflict = {c: dict(v) for c, v in col2cons.items() if len(v) > 1}
match_cols = sorted(col2cons); cons_seen = {list(v)[0] for v in col2cons.values()}
mixed = {c: dict(v) for c, v in col_states.items() if "match" in v and (("insert" in v) or ("local_end" in v))}
print(f"BGSU columns holding match residues: {len(match_cols)} covering {len(cons_seen)} distinct consensus positions; columns mapping to >1 consensus position: {len(conflict)} {list(conflict.items())[:5]}")
print(f"columns holding both match and insert/local-end residues (mixed): {len(mixed)} {list(mixed.items())[:6]}")
print("max BGSU column:", max(cols[k] for k in cols) if cols else None, "| columns used by any residue:", len({c for v in cols.values() for c in v}))
# monotonic? consensus increases with column
seq = sorted((c, list(v)[0]) for c, v in col2cons.items() if len(v) == 1)
print("column->consensus monotonic:", all(b[1] > a[1] for a, b in zip(seq, seq[1:])))
with open(RES/"bgsu_column_to_consensus.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["bgsu_column", "rfam_consensus_pos", "n_match_residues", "n_insert_residues", "n_local_end_residues"])
    for c in sorted(col_states): w.writerow([c, list(col2cons[c])[0] if len(col2cons[c]) == 1 else ("CONFLICT" if col2cons[c] else ""), col_states[c]["match"], col_states[c]["insert"], col_states[c]["local_end"]])
# (b) BGSU pairs inside unnumbered regions: are they offset-identical (same k-th residue after same anchor)?
disc = [d for d in csv.DictReader(open(RES/"discrepancies.csv")) if d["category"] in ("bgsu_pairs_insert_residue", "bgsu_pairs_local_end_residue")]
out = []
for d in disc:
    a, b = ours[(d["pdbA"], int(d["posA"]))], ours[(d["pdbB"], int(d["posB"]))]
    out.append({"file": d["file"], "pdbA": d["pdbA"], "posA": d["posA"], "stateA": a["alignment_state"], "anchorA": a["insertion_anchor"], "offsetA": a["insertion_offset"], "ntA": a["nucleotide_parent"], "ppA": a["posterior_prob"], "coordsA": a["has_coordinates"],
                "pdbB": d["pdbB"], "posB": d["posB"], "stateB": b["alignment_state"], "anchorB": b["insertion_anchor"], "offsetB": b["insertion_offset"], "ntB": b["nucleotide_parent"], "ppB": b["posterior_prob"], "coordsB": b["has_coordinates"],
                "same_anchor": a["insertion_anchor"] == b["insertion_anchor"], "same_offset": a["insertion_offset"] == b["insertion_offset"], "same_nt": a["nucleotide_parent"] == b["nucleotide_parent"],
                "bgsu_column": next((r["bgsu_column"] for r in norm if r["file"] == d["file"] and r.get("c0_pos") == d["posA"] and r.get("c1_pos") == d["posB"]), "")})
with open(RES/"unnumbered_bgsu_pairs_detail.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
summ = collections.defaultdict(collections.Counter)
for o in out: summ[(o["file"], o["stateA"] + "/" + o["stateB"])][f"same_anchor={o['same_anchor']},same_offset={o['same_offset']}"] += 1; summ[(o["file"], o["stateA"] + "/" + o["stateB"])]["same_nt"] += o["same_nt"]; summ[(o["file"], o["stateA"] + "/" + o["stateB"])]["n"] += 1
for k, v in summ.items(): print(k, dict(v))
# insertion run lengths at the shared anchors in each chain (to see whether BGSU pairs the k-th residue of unequal-length runs)
runs = collections.defaultdict(int)
for (pdb, pos), o in ours.items():
    if o["alignment_state"] in ("insert", "local_end"): runs[(pdb, o["alignment_state"], o["insertion_anchor"])] += 1
for key in sorted({(o["pdbA"], o["stateA"], o["anchorA"]) for o in out} | {(o["pdbB"], o["stateB"], o["anchorB"]) for o in out}): print("run length", key, runs[key])
