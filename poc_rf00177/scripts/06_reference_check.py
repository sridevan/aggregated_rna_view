"""Step 7: correspondence checks independent of alignment score.
(a) E. coli (7k00) numbering -> consensus position lookup table (E. coli numbering is the field's reference frame).
(b) Universally conserved 16S motifs located in E. coli; nucleotides each instance places at the same consensus positions.
(c) Identity to E. coli in 50-consensus-position windows; (d) base-pair consistency with SS_cons per window -> flags low-confidence regions."""
import csv, json, sys, collections
from pathlib import Path
RES = Path(sys.argv[1]); REF = "7k00"
rows = list(csv.DictReader(open(RES/"rf00177_residue_mapping.csv")))
summ = json.load(open(RES/"rf00177_summary.json")); CLEN = summ["CLEN"]
cons = collections.defaultdict(dict)   # pdb -> cons_pos -> row
for r in rows:
    if r["alignment_state"] == "match": cons[r["pdb_id"]][int(r["rfam_consensus_pos"])] = r
pp_of = lambda r: 1.0 if r["posterior_prob"] == "*" else (int(r["posterior_prob"]) + 0.5) / 10 if r["posterior_prob"].isdigit() else None
# (a) E. coli table
ref = {int(r["polymer_seq_position"]): r for r in rows if r["pdb_id"] == REF}
with open(RES/"ecoli_7k00_numbering_to_consensus.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["ecoli_residue_number(7k00 auth)", "nucleotide", "residue_name", "alignment_state", "rfam_consensus_pos", "insertion_anchor", "insertion_offset", "posterior_prob"])
    for p in sorted(ref): r = ref[p]; w.writerow([r["author_residue_number"], r["nucleotide_parent"], r["residue_name"], r["alignment_state"], r["rfam_consensus_pos"], r["insertion_anchor"], r["insertion_offset"], r["posterior_prob"]])
# SS_cons pair map from consensus_columns (match columns only)
cc = [c for c in csv.DictReader(open(RES/"rf00177_consensus_columns.csv")) if c["pdb_id"] == REF and c["state"] in ("match", "deletion")]
ss = {int(c["rfam_consensus_pos"]): c["ss_cons"] for c in cc}
stack, pairs = [], {}
for i in sorted(ss):
    ch = ss[i]
    if ch in "<([{": stack.append(i)
    elif ch in ">)]}": j = stack.pop(); pairs[i] = j; pairs[j] = i
canon = {"AU", "UA", "GC", "CG", "GU", "UG"}
# (b) motifs (E. coli 16S, universally conserved); located by sequence search in 7k00, never assumed numbering
ecoli_seq = "".join(ref[p]["nucleotide_parent"] for p in sorted(ref))
motifs = {"530_loop": "GUGCCAGCAGCCGCGGUAA", "790_loop": "GGAUUAGAUACCC", "1400_region": "GUGAAGUCGUAACAAGG", "anti_SD_3prime": "GAUCACCUCCUUA", "h44_decoding_1492_1493": "GUCGUAACAAGGUA"}
motif_report = {}
for mname, mseq in motifs.items():
    idx = ecoli_seq.find(mseq)
    if idx < 0: motif_report[mname] = {"found_in_7k00": False}; continue
    ecoli_range = (idx + 1, idx + len(mseq)); cps = [ref[p]["rfam_consensus_pos"] for p in range(idx + 1, idx + len(mseq) + 1)]
    per = {}
    for pdb in cons:
        got = "".join(cons[pdb][int(c)]["nucleotide_parent"] if c and int(c) in cons[pdb] else "-" for c in cps)
        ident = sum(a == b for a, b in zip(got, mseq)) / len(mseq)
        per[pdb] = {"sequence_at_same_consensus": got, "identity": round(ident, 3), "mean_pp": round(sum(pp_of(cons[pdb][int(c)]) or 0 for c in cps if c and int(c) in cons[pdb]) / max(1, sum(1 for c in cps if c and int(c) in cons[pdb])), 3)}
    motif_report[mname] = {"found_in_7k00": True, "ecoli_range": ecoli_range, "consensus_positions": [cps[0], cps[-1]], "all_match_cols": all(cps), "per_instance": per}
# (c)/(d) windows
W = 50; win_rows = []
for pdb in sorted(cons):
    for lo in range(1, CLEN + 1, W):
        hi = min(CLEN, lo + W - 1); cols = range(lo, hi + 1)
        matched = [c for c in cols if c in cons[pdb]]
        ident = [cons[pdb][c]["nucleotide_parent"] == cons[REF][c]["nucleotide_parent"] for c in matched if c in cons[REF]]
        pp = [pp_of(cons[pdb][c]) for c in matched]; pp = [x for x in pp if x is not None]
        pr = [c for c in matched if c in pairs and pairs[c] in cons[pdb]]
        ok = [cons[pdb][c]["nucleotide_parent"] + cons[pdb][pairs[c]]["nucleotide_parent"] in canon for c in pr]
        obs = [cons[pdb][c]["has_coordinates"] == "1" for c in matched]
        ins = [r for r in rows if r["pdb_id"] == pdb and r["alignment_state"] in ("insert", "local_end") and r["insertion_anchor"] and lo <= int(r["insertion_anchor"]) <= hi]
        flag = []
        if matched and len(matched) / len(cols) < 0.5: flag.append("sparse_coverage")
        if ident and sum(ident) / len(ident) < 0.5: flag.append("low_identity_to_ecoli")
        if pp and sum(pp) / len(pp) < 0.8: flag.append("low_posterior")
        if ok and sum(ok) / len(ok) < 0.7: flag.append("pairing_inconsistent")
        if len(ins) >= 20: flag.append("long_insertion")
        win_rows.append({"pdb_id": pdb, "cons_from": lo, "cons_to": hi, "n_matched": len(matched), "coverage": round(len(matched) / len(cols), 3),
                         "identity_to_ecoli": round(sum(ident) / len(ident), 3) if ident else "", "mean_pp": round(sum(pp) / len(pp), 3) if pp else "",
                         "n_pairs_evaluable": len(pr), "pair_consistency": round(sum(ok) / len(ok), 3) if ok else "", "fraction_with_coordinates": round(sum(obs) / len(obs), 3) if obs else "",
                         "inserted_or_EL_residues_anchored_here": len(ins), "flags": ";".join(flag)})
with open(RES/"rf00177_window_confidence.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(win_rows[0].keys())); w.writeheader(); w.writerows(win_rows)
per_inst = {}
for pdb in sorted(cons):
    wr = [w for w in win_rows if w["pdb_id"] == pdb]
    allpr = [c for c in cons[pdb] if c in pairs and pairs[c] in cons[pdb]]
    per_inst[pdb] = {"identity_to_ecoli_overall": round(sum(cons[pdb][c]["nucleotide_parent"] == cons[REF][c]["nucleotide_parent"] for c in cons[pdb] if c in cons[REF]) / max(1, sum(1 for c in cons[pdb] if c in cons[REF])), 3),
                     "pair_consistency_overall": round(sum(cons[pdb][c]["nucleotide_parent"] + cons[pdb][pairs[c]]["nucleotide_parent"] in canon for c in allpr) / max(1, len(allpr)), 3), "n_pairs_evaluable": len(allpr) // 2,
                     "flagged_windows": [(w["cons_from"], w["cons_to"], w["flags"]) for w in wr if w["flags"]]}
json.dump({"motifs": motif_report, "per_instance": per_inst, "n_consensus_pairs_in_SS_cons": len(pairs) // 2}, open(RES/"reference_correspondence.json", "w"), indent=1)
print("SS_cons pairs:", len(pairs) // 2)
for m, v in motif_report.items():
    print(m, v.get("ecoli_range"), v.get("consensus_positions"), "all_match_cols" , v.get("all_match_cols"))
    for pdb, x in v.get("per_instance", {}).items(): print(f"   {pdb}: {x['sequence_at_same_consensus']} id={x['identity']} pp={x['mean_pp']}")
for pdb, v in per_inst.items(): print(pdb, {k: v[k] for k in ["identity_to_ecoli_overall", "pair_consistency_overall", "n_pairs_evaluable"]}, "flagged:", len(v["flagged_windows"]), v["flagged_windows"][:6])
