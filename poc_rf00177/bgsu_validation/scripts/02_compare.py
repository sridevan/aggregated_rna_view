"""Compare BGSU nucleotide correspondences with our shared-consensus-position pairs."""
import csv, json, sys, collections
from pathlib import Path
RES, MAIN = Path(sys.argv[1]), Path(sys.argv[2])
ours = {}
for f in [MAIN/"rf00177_residue_mapping.csv", RES/"repro/rf00177_residue_mapping.csv"]:
    for r in csv.DictReader(open(f)): ours[(r["pdb_id"], r["auth_asym_id"], int(r["polymer_seq_position"]))] = r
winflags = collections.defaultdict(dict)
for w in csv.DictReader(open(MAIN/"rf00177_window_confidence.csv")):
    if w["flags"]: winflags[w["pdb_id"]][int(w["cons_from"])] = w["flags"]
def flagged(pdb, cons): return winflags.get(pdb, {}).get(((int(cons) - 1) // 50) * 50 + 1, "")
bg = list(csv.DictReader(open(RES/"bgsu_correspondences_normalised.csv")))
files = sorted({r["file"] for r in bg})
metrics, disc, regions, unnumbered = [], [], [], []
for f in files:
    rows = [r for r in bg if r["file"] == f]
    pdbA, chA, pdbB, chB = rows[0]["pdb1"], rows[0]["chain1"], rows[0]["pdb2"], rows[0]["chain2"]
    B = {(int(r["pos1"]), int(r["pos2"])) for r in rows}
    partnerA = {a: b for a, b in B}; partnerB = {b: a for a, b in B}
    A_rows = {k[2]: v for k, v in ours.items() if k[0] == pdbA and k[1] == chA}; B_rows = {k[2]: v for k, v in ours.items() if k[0] == pdbB and k[1] == chB}
    byconsB = {int(v["rfam_consensus_pos"]): p for p, v in B_rows.items() if v["alignment_state"] == "match"}
    O = {(p, byconsB[int(v["rfam_consensus_pos"])]) for p, v in A_rows.items() if v["alignment_state"] == "match" and int(v["rfam_consensus_pos"]) in byconsB}
    obs = lambda pr: A_rows[pr[0]]["has_coordinates"] == "1" and B_rows[pr[1]]["has_coordinates"] == "1"
    O_obs = {pr for pr in O if obs(pr)}
    both = B & O
    hi = lambda pr: A_rows[pr[0]]["posterior_prob"] == "*" and B_rows[pr[1]]["posterior_prob"] == "*"
    cons_of = lambda pr: A_rows[pr[0]]["rfam_consensus_pos"] if A_rows[pr[0]]["alignment_state"] == "match" else B_rows[pr[1]].get("rfam_consensus_pos", "")
    # classify disagreements
    cat = collections.Counter()
    for pr in sorted(B - O):
        a, b = A_rows.get(pr[0]), B_rows.get(pr[1])
        if a is None or b is None: c = "identifier_unresolved"
        else:
            sa, sb = a["alignment_state"], b["alignment_state"]
            if sa == "match" and sb == "match": c = "shift_both_match_different_consensus"
            elif "local_end" in (sa, sb): c = "bgsu_pairs_local_end_residue"
            elif "insert" in (sa, sb): c = "bgsu_pairs_insert_residue"
            elif "not_submitted_flank" in (sa, sb): c = "bgsu_pairs_flank_outside_pdbe_region"
            else: c = "other"
        cat[c] += 1
        disc.append({"file": f, "direction": "bgsu_only", "category": c, "pdbA": pdbA, "posA": pr[0], "stateA": a["alignment_state"] if a else "", "consA": a["rfam_consensus_pos"] if a else "", "anchorA": a["insertion_anchor"] if a else "", "ppA": a["posterior_prob"] if a else "", "coordsA": a["has_coordinates"] if a else "",
                     "pdbB": pdbB, "posB": pr[1], "stateB": b["alignment_state"] if b else "", "consB": b["rfam_consensus_pos"] if b else "", "anchorB": b["insertion_anchor"] if b else "", "ppB": b["posterior_prob"] if b else "", "coordsB": b["has_coordinates"] if b else "",
                     "delta_consensus": (int(a["rfam_consensus_pos"]) - int(b["rfam_consensus_pos"])) if a and b and a["alignment_state"] == "match" and b["alignment_state"] == "match" else "",
                     "window_flags_B": flagged(pdbB, b["rfam_consensus_pos"]) if b and b["alignment_state"] == "match" else ""})
    for pr in sorted(O_obs - B):
        a, b = A_rows[pr[0]], B_rows[pr[1]]
        pa, pb = partnerA.get(pr[0]), partnerB.get(pr[1])
        c = "bgsu_no_partner_for_either" if pa is None and pb is None else ("bgsu_partner_differs(shift)" if (pa is not None and pa != pr[1]) or (pb is not None and pb != pr[0]) else "bgsu_no_partner_for_one")
        cat[c] += 1
        disc.append({"file": f, "direction": "ours_only(observed)", "category": c, "pdbA": pdbA, "posA": pr[0], "stateA": "match", "consA": a["rfam_consensus_pos"], "anchorA": "", "ppA": a["posterior_prob"], "coordsA": a["has_coordinates"],
                     "pdbB": pdbB, "posB": pr[1], "stateB": "match", "consB": b["rfam_consensus_pos"], "anchorB": "", "ppB": b["posterior_prob"], "coordsB": b["has_coordinates"], "delta_consensus": 0, "window_flags_B": flagged(pdbB, b["rfam_consensus_pos"]),
                     "bgsu_partner_of_A": pa if pa is not None else "", "bgsu_partner_of_B": pb if pb is not None else ""})
    n_hi, n_hi_ok = sum(1 for pr in O_obs if hi(pr)), sum(1 for pr in O_obs if hi(pr) and pr in B)
    n_lo, n_lo_ok = sum(1 for pr in O_obs if not hi(pr)), sum(1 for pr in O_obs if not hi(pr) and pr in B)
    n_fl = sum(1 for pr in O_obs if flagged(pdbB, cons_of(pr))); n_fl_ok = sum(1 for pr in O_obs if flagged(pdbB, cons_of(pr)) and pr in B)
    # Step 4: BGSU partners for our unnumbered residues of chain B (and A)
    for side, rowsX, partner, pdbX, other in (("A", A_rows, partnerA, pdbA, pdbB), ("B", B_rows, partnerB, pdbB, pdbA)):
        for st in ("insert", "local_end"):
            xs = [p for p, v in rowsX.items() if v["alignment_state"] == st]
            paired = [p for p in xs if p in partner]
            if xs:
                byanchor = collections.Counter(rowsX[p]["insertion_anchor"] for p in paired)
                unnumbered.append({"file": f, "chain": f"{pdbX}", "state": st, "n_residues": len(xs), "n_with_coordinates": sum(rowsX[p]["has_coordinates"] == "1" for p in xs),
                                   "n_given_partner_by_bgsu": len(paired), "partner_chain": other, "partners_by_anchor": ";".join(f"{a}:{n}" for a, n in sorted(byanchor.items(), key=lambda x: -x[1])[:8])})
    metrics.append({"file": f, "chainA": f"{pdbA}/{chA}", "chainB": f"{pdbB}/{chB}", "bgsu_pairs": len(B), "our_pairs_all": len(O), "our_pairs_observed": len(O_obs), "matching_pairs": len(both),
                    "precision_vs_bgsu(observed)": round(len(both) / len(O_obs), 4) if O_obs else "", "recall_vs_bgsu": round(len(both) / len(B), 4) if B else "",
                    "our_observed_pairs_not_in_bgsu": len(O_obs - B), "bgsu_pairs_not_ours": len(B - O),
                    "high_conf_pairs(PP* both)": n_hi, "high_conf_agreement": round(n_hi_ok / n_hi, 4) if n_hi else "", "low_conf_pairs": n_lo, "low_conf_agreement": round(n_lo_ok / n_lo, 4) if n_lo else "",
                    "pairs_in_flagged_windows": n_fl, "flagged_window_agreement": round(n_fl_ok / n_fl, 4) if n_fl else "", **{f"disc:{k}": v for k, v in cat.items()}})
    # regions (50 consensus positions, by consensus of chain A residue)
    for lo in range(1, 1534, 50):
        inwin = lambda pr: lo <= int(cons_of(pr) or 0) <= lo + 49
        ob = [pr for pr in O_obs if inwin(pr)]; bb = [pr for pr in B if (A_rows.get(pr[0], {}).get("alignment_state") == "match" and inwin(pr))]
        ag = [pr for pr in ob if pr in B]
        if ob or bb:
            regions.append({"file": f, "chainB": f"{pdbB}/{chB}", "cons_from": lo, "cons_to": lo + 49, "our_observed_pairs": len(ob), "bgsu_pairs_with_A_match": len(bb), "agreeing": len(ag),
                            "agreement_of_ours": round(len(ag) / len(ob), 3) if ob else "", "insert_residues_B": sum(1 for p, v in B_rows.items() if v["alignment_state"] == "insert" and v["insertion_anchor"] and lo <= int(v["insertion_anchor"]) <= lo + 49),
                            "local_end_residues_B": sum(1 for p, v in B_rows.items() if v["alignment_state"] == "local_end" and v["insertion_anchor"] and lo <= int(v["insertion_anchor"]) <= lo + 49),
                            "deleted_consensus_B": 50 - sum(1 for c in byconsB if lo <= c <= lo + 49), "window_flags_B": flagged(pdbB, lo)})
def dump(name, rows):
    keys = []; [keys.append(k) for r in rows for k in r if k not in keys]
    with open(RES/name, "w", newline="") as fh: w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
dump("pair_metrics.csv", metrics); dump("discrepancies.csv", disc); dump("region_windows.csv", regions); dump("unnumbered_regions_bgsu.csv", unnumbered)
for m in metrics: print(json.dumps(m))
print("--- unnumbered"); [print(json.dumps(u)) for u in unnumbered]
