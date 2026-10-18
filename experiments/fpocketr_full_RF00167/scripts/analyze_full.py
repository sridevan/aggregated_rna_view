"""Family-wide analysis of pockets.json / pocket_tracks.json for the report: statistics, filter sensitivity, known-site comparison with
chain-level Arpeggio contacts, apo/holo, construct/quality variability, reproducibility vs the repeat run, and an independent reconstruction
of the aggregated counts from the pocket records. Writes analysis.json and prints a summary."""
import json, statistics, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
E = Path(__file__).resolve().parents[1]; ROOT = E.parents[1]
sys.path.insert(0, str(ROOT))
from rna_explorer.annotation_mapping import InstanceIndex
from rna_explorer.pockets import qualifying, FILTERS
F = ROOT / "data/families/RF00167"
idx = InstanceIndex(F)
P = json.load(open(F / "pockets.json")); T = json.load(open(F / "pocket_tracks.json")); inst = {r["seq_name"]: r for r in idx.instances}
KNOWN = [25, 26, 51, 54, 55, 56, 79, 80, 81]            # pilot-derived reference set
CORE = [25, 26, 51, 55, 56, 80, 81]                      # specific-interaction set
out: dict = {}
recs = P["instances"]
ok = [r for r in recs if r["status"] != "failed"]
# ---- A execution
out["execution"] = {"n_instances": len(recs), "status_counts": P["counts"]["status_counts"], "failures": [(r["seq_name"], r.get("failure_reason")) for r in recs if r["status"] == "failed"],
                    "warnings": [(r["seq_name"], r["warnings"]) for r in recs if r["warnings"]], "total_seconds": P["counts"]["total_seconds"],
                    "seconds_per_run": {k: round(v, 1) for k, v in zip(("min", "median", "max"), (lambda s: (min(s), statistics.median(s), max(s)))([r["run"]["seconds"] for r in recs if r.get("run", {}).get("seconds")]))} if any(r.get("run", {}).get("seconds") for r in recs) else None,
                    "input_checks": {"roundtrip_ok": sum(1 for r in recs if r.get("input", {}).get("roundtrip_ok")), "altloc_copies_dropped_total": sum(r.get("input", {}).get("altloc_copies_dropped", 0) for r in recs),
                                     "modified_residues": dict(Counter(m for r in recs for m in r.get("input", {}).get("modified_residues", []))), "insertion_codes": sum(1 for r in recs if r.get("input", {}).get("has_insertion_codes"))}}
# ---- B pocket statistics
allp = [(r, p) for r in ok for p in r["pockets"]]
scores = [p["score"] for _, p in allp]; vols = [p["volume_A3"] for _, p in allp]; nts = [p["n_nucleotides"] for _, p in allp]
out["pockets"] = {"n_pockets": len(allp), "pockets_per_instance": dict(sorted(Counter(r["n_pockets"] for r in ok).items())),
                  "score": {"min": min(scores), "q1": float(np.percentile(scores, 25)), "median": statistics.median(scores), "q3": float(np.percentile(scores, 75)), "max": max(scores)} if scores else None,
                  "volume_A3": {"min": min(vols), "median": statistics.median(vols), "max": max(vols)} if vols else None, "nucleotides_per_pocket": {"min": min(nts), "median": statistics.median(nts), "max": max(nts)} if nts else None,
                  "n_pockets_score_ge_0.4": sum(1 for s in scores if s >= 0.4), "n_pockets_score_lt_0.4": sum(1 for s in scores if s < 0.4),
                  "pockets_with_flags": dict(Counter(f for _, p in allp for f in p["flags"])), "pockets_entirely_non_consensus": sum(1 for _, p in allp if not p["consensus_positions"]),
                  "mapping_checks": {"nt_match_csv_all": all(p["nt_match_csv"] for _, p in allp), "map_consistent_all": all(p["map_consistent"] for _, p in allp), "unresolved_nucleotides": sum(p["unresolved_nucleotides"] for _, p in allp)},
                  "shape": dict(Counter(p["shape"] for _, p in allp))}
pos_count = Counter(c for _, p in allp for c in p["consensus_positions"])
out["pockets"]["most_frequent_positions_by_pocket_count"] = pos_count.most_common(15)
# ---- C / D aggregated occurrence and filters
def profile(flt):
    tr = T["tracks"][flt]
    nz = [p for p in tr if p["n_supporting_instances"]]
    return {"n_positions_nonzero": len(nz), "n_positions_evaluable": sum(1 for p in tr if p["n_evaluable_instances"]),
            "profile": {str(p["rfam_consensus_pos"]): f"{p['n_supporting_instances']}/{p['n_evaluable_instances']}={p['frequency']}" for p in nz},
            "known_site": {str(k): (tr[k - 1]["n_supporting_instances"], tr[k - 1]["n_evaluable_instances"], tr[k - 1]["frequency"]) for k in KNOWN},
            "positions_ge_0.5": [p["rfam_consensus_pos"] for p in tr if (p["frequency"] or 0) >= 0.5], "positions_0.1_0.5": [p["rfam_consensus_pos"] for p in tr if 0.1 <= (p["frequency"] or 0) < 0.5],
            "positions_lt_0.1_nonzero": [p["rfam_consensus_pos"] for p in nz if (p["frequency"] or 0) < 0.1],
            "n_instances_with_qualifying_pocket": sum(1 for r in ok if qualifying(r["pockets"], flt)), "total_qualifying_pockets": sum(len(qualifying(r["pockets"], flt)) for r in ok),
            "pockets_per_instance": dict(sorted(Counter(len(qualifying(r["pockets"], flt)) for r in ok).items())),
            "outside_known_positions_nonzero": [p["rfam_consensus_pos"] for p in nz if p["rfam_consensus_pos"] not in KNOWN]}
out["filters"] = {flt: profile(flt) for flt in FILTERS}
# ---- independent reconstruction of the aggregated counts from the pocket records
recon_ok = True
for flt in FILTERS:
    for pos in range(1, idx.clen + 1):
        sup = {r["seq_name"] for r in ok if idx.observable(r["seq_name"], pos) and any(pos in p["consensus_positions"] for p in qualifying(r["pockets"], flt))}
        ev = {r["seq_name"] for r in ok if idx.match_residue(r["seq_name"], pos) is not None and idx.observable(r["seq_name"], pos)}
        t = T["tracks"][flt][pos - 1]
        if (len(sup), len(ev)) != (t["n_supporting_instances"], t["n_evaluable_instances"]): recon_ok = False
out["reconstruction_matches_tracks"] = recon_ok
# ---- E known-site comparison with chain-level Arpeggio contacts (ligand class, same chain)
ARP = ROOT / "experiments/chain_specific_ligands/arpeggio_pilot/results"
def contacts(sn):
    f = ARP / f"{inst[sn]['pdb_id']}.json"
    if not f.exists(): return None
    c = [x for x in json.load(open(f))["residue_contacts"] if x["seq_name"] == sn and x["ligand_class"] == "ligand" and x["rfam_consensus_pos"]]
    if not c: return {"has_ligand": False}
    return {"has_ligand": True, "ccds": sorted({x["ligand_ccd"] for x in c}), "specific": sorted({x["rfam_consensus_pos"] for x in c if x["specific"]}), "any": sorted({x["rfam_consensus_pos"] for x in c})}
holo, apo, unknown = [], [], []
for r in ok:
    c = contacts(r["seq_name"]); sn = r["seq_name"]; meta = inst[sn]
    row = {"seq_name": sn, "status": r["status"], "n_pockets": r["n_pockets"], "method": meta["experimental_method"], "resolution": meta["resolution"], "chain_length": meta["chain_length"]}
    if c is None: unknown.append(row); continue
    if not c["has_ligand"]:
        apo.append({**row, "pockets": [(p["rank"], p["score"], p["consensus_positions"]) for p in r["pockets"]], "known_overlap_top1": len(set(r["pockets"][0]["consensus_positions"]) & set(KNOWN)) if r["pockets"] else 0}); continue
    best = None
    for p in r["pockets"]:
        ps = set(p["consensus_positions"]); s = set(c["specific"]); a = set(c["any"])
        m = {"rank": p["rank"], "score": p["score"], "volume_A3": p["volume_A3"], "n_alpha_spheres": p["n_alpha_spheres"], "shape": p["shape"], "positions": p["consensus_positions"],
             "recall_specific": round(len(ps & s) / len(s), 3) if s else None, "precision_any": round(len(ps & a) / len(ps), 3) if ps else None, "jaccard_any": round(len(ps & a) / len(ps | a), 3) if ps | a else None}
        if best is None or (m["jaccard_any"] or 0) > (best["jaccard_any"] or 0): best = m
    holo.append({**row, "ccds": c["ccds"], "contacts_specific": c["specific"], "contacts_any": c["any"], "best_pocket": best, "top1_is_best": bool(best and best["rank"] == 1), "n_pockets_overlapping_site": sum(1 for p in r["pockets"] if set(p["consensus_positions"]) & set(c["any"]))})
def summ(vals): vals = [v for v in vals if v is not None]; return {"n": len(vals), "median": statistics.median(vals), "min": min(vals), "max": max(vals)} if vals else None
out["holo"] = {"n": len(holo), "n_with_pocket_overlapping_site": sum(1 for h in holo if h["best_pocket"] and (h["best_pocket"]["jaccard_any"] or 0) > 0), "n_top1_overlaps": sum(1 for h in holo if h["top1_is_best"] and (h["best_pocket"]["jaccard_any"] or 0) > 0),
               "n_without_any_pocket": sum(1 for h in holo if h["n_pockets"] == 0), "recall_specific": summ([h["best_pocket"]["recall_specific"] for h in holo if h["best_pocket"]]), "precision_any": summ([h["best_pocket"]["precision_any"] for h in holo if h["best_pocket"]]),
               "jaccard_any": summ([h["best_pocket"]["jaccard_any"] for h in holo if h["best_pocket"]]), "best_score": summ([h["best_pocket"]["score"] for h in holo if h["best_pocket"]]), "best_volume": summ([h["best_pocket"]["volume_A3"] for h in holo if h["best_pocket"]]),
               "no_pocket_cases": [(h["seq_name"], h["method"], h["resolution"], h["chain_length"], h["ccds"]) for h in holo if h["n_pockets"] == 0], "rows": holo}
out["apo"] = {"n": len(apo), "rows": apo}; out["no_contact_data"] = unknown
# ---- F construct / quality variability
def group(rows, key):
    g = defaultdict(list)
    for r in rows: g[key(r)].append(r)
    return {str(k): {"n": len(v), "with_pockets": sum(1 for x in v if x["n_pockets"]), "mean_pockets": round(sum(x["n_pockets"] for x in v) / len(v), 2)} for k, v in sorted(g.items(), key=lambda kv: str(kv[0]))}
rows = [{"seq_name": r["seq_name"], "n_pockets": r["n_pockets"], "method": inst[r["seq_name"]]["experimental_method"], "resolution": inst[r["seq_name"]]["resolution"], "chain_length": inst[r["seq_name"]]["chain_length"],
         "scaffold": inst[r["seq_name"]]["chain_length"] > 100, "modified": bool(r.get("input", {}).get("modified_residues"))} for r in ok]
out["variability"] = {"by_method": group(rows, lambda r: r["method"]), "by_resolution_band": group(rows, lambda r: "<=2.0" if (r["resolution"] or 9) <= 2.0 else "2.0-2.5" if r["resolution"] <= 2.5 else "2.5-3.0" if r["resolution"] <= 3.0 else ">3.0"),
                      "by_length_band": group(rows, lambda r: "<=72" if r["chain_length"] <= 72 else "73-110" if r["chain_length"] <= 110 else ">110 (scaffold)"), "by_modified": group(rows, lambda r: r["modified"]),
                      "flagged_pockets": dict(Counter(f for _, p in allp for f in p["flags"])), "scaffold_instances": [r["seq_name"] for r in rows if r["scaffold"]],
                      "pockets_on_scaffold_instances": [(r["seq_name"], p["rank"], p["score"], p["consensus_positions"], p["n_non_consensus_nucleotides"]) for r in ok if inst[r["seq_name"]]["chain_length"] > 100 for p in r["pockets"]]}
# sensitivity subset: aptamer-only constructs (<= 72 nt), all filter
sub = [r for r in ok if inst[r["seq_name"]]["chain_length"] <= 72]
subprof = {}
for pos in KNOWN:
    ev = [r["seq_name"] for r in sub if idx.match_residue(r["seq_name"], pos) and idx.observable(r["seq_name"], pos)]
    sup = [s for s in ev if any(pos in p["consensus_positions"] for r in sub if r["seq_name"] == s for p in qualifying(r["pockets"], "all_score_ge_0"))]
    subprof[str(pos)] = f"{len(sup)}/{len(ev)}={round(len(sup)/len(ev),3) if ev else None}"
out["variability"]["aptamer_only_subset_known_site"] = {"n_instances": len(sub), "profile": subprof}
# ---- repeat run comparison
rep = E / "repeat" / "work" / "runs"
cmp = []
if rep.exists():
    from rna_explorer.pockets import parse_run
    for r in ok:
        d = rep / f"{r['seq_name']}_clean_out"
        if not d.exists(): continue
        p2 = parse_run(d, r["seq_name"])
        same_scores = [p["score"] for p in r["pockets"]] == [p["score"] for p in p2]; same_as = [p["n_alpha_spheres"] for p in r["pockets"]] == [p["n_alpha_spheres"] for p in p2]
        same_nt = [sorted(p["pocket_nt_author_csv"]) for p in r["pockets"]] == [sorted(p["pocket_nt_author_csv"]) for p in p2]
        vol = [(p["volume_A3"], q["volume_A3"]) for p, q in zip(r["pockets"], p2)]
        cmp.append({"seq_name": r["seq_name"], "n_pockets": (len(r["pockets"]), len(p2)), "scores_identical": same_scores, "alpha_spheres_identical": same_as, "lining_nucleotides_identical": same_nt,
                    "volume_rel_diff_max": round(max((abs(a - b) / a for a, b in vol if a), default=0.0), 4)})
out["reproducibility"] = cmp
json.dump(out, open(E / "analysis.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ("execution", "pockets", "reconstruction_matches_tracks")}, indent=1)[:6000])
for flt, pr in out["filters"].items(): print(flt, {k: v for k, v in pr.items() if k not in ("profile",)})
print("HOLO", {k: v for k, v in out["holo"].items() if k != "rows"})
print("APO", [(a["seq_name"], a["n_pockets"], a["pockets"], a["known_overlap_top1"]) for a in out["apo"]["rows"]])
print("VAR", json.dumps({k: v for k, v in out["variability"].items() if k != "pockets_on_scaffold_instances"}, indent=0)[:3000])
print("SCAFFOLD POCKETS", out["variability"]["pockets_on_scaffold_instances"][:12])
print("REPRO", out["reproducibility"])
