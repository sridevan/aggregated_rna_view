"""Audit of RNA loop counts at Rfam consensus position POS for RF00167 (investigation only; reads existing outputs and raw BGSU files)."""
import csv, json, sys
from collections import defaultdict
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
F = ROOT / "data/families/RF00167"
POS = int(sys.argv[1]) if len(sys.argv) > 1 else 53
tracks = json.load(open(F / "tracks.json")); loops = json.load(open(F / "loops.json")); inst = json.load(open(F / "instances.json"))
out = {}
# 1. displayed values
t = tracks["positions"][POS - 1]
out["tracks_json"] = {k: {"n_instances": t[k]["n_instances"], "n_loop_observations": t[k]["n_loop_observations"], "n_entries": t[k]["n_entries"],
                          "denominators": t[k]["denominators"]} for k in ("hairpin_loops", "internal_loops", "junction_loops_J3", "junction_loops_J4", "junction_loops_union")}
# 2. recompute from loops.json (per loop record, exact nucleotide membership)
by = defaultdict(set); obs = defaultdict(set); loop_ids = defaultdict(list)
for L in loops["loops"]:
    nts = [n for n in L["nucleotides"] if n["rfam_consensus_pos"] == POS and n["mapping_status"] in ("mapped", "mapped_no_coordinates")]
    if not nts: continue
    for n in nts:
        by[L["loop_type"]].add(n["seq_name"]); loop_ids[L["loop_type"]].append((L["loop_id"], n["seq_name"], n["unit_id"], L["motif_group"]))
    obs[L["loop_type"]].add(L["loop_id"])
JL = by["J3"] | by["J4"]; ANY = by["HL"] | by["IL"] | JL
# evaluable instances at POS: match residue with coordinates (observable) and loops source evaluated (all 66 entries)
observable = set()
for r in inst:
    m = json.load(open(F / r["mapping_file"]))
    if any(x["alignment_state"] == "match" and x["rfam_consensus_pos"] == POS and x["has_coordinates"] for x in m["residues"]): observable.add(r["seq_name"])
out["recomputed_from_loops_json"] = {"HL": len(by["HL"]), "IL": len(by["IL"]), "J3": len(by["J3"]), "J4": len(by["J4"]), "JL_union": len(JL), "J3_plus_J4_sum": len(by["J3"]) + len(by["J4"]),
    "HL_and_JL_overlap": len(by["HL"] & JL), "HL_only": len(by["HL"] - JL), "JL_only": len(JL - by["HL"]), "any_loop": len(ANY), "no_loop_but_observable": len(observable - ANY),
    "observable_instances": len(observable), "total_instances": len(inst), "loop_observations": {k: len(v) for k, v in obs.items()},
    "memberships_per_instance_HL": max(sum(1 for x in loop_ids["HL"] if x[1] == s) for s in by["HL"]) if by["HL"] else 0}
out["instances_without_any_loop_at_pos"] = sorted(observable - ANY)
out["overlap_instances"] = sorted(by["HL"] & JL)
out["positive_not_observable"] = sorted(ANY - observable)
# 3. representative overlapping instances with raw evidence
def raw_loop(pdb, loop_id):
    for row in csv.reader(open(ROOT / "annotation_assessment/cache/bgsu" / f"loops_{pdb.upper()}.csv")):
        if row and row[0] == loop_id: return row[1].split(",")
    return None
samples = []
for sn in ["1y26_1_X_A", "4fe5_1_B_A", "3rkf_1_A_A", "9ndd_1_A_A"]:
    if sn not in (by["HL"] & JL): continue
    r = next(i for i in inst if i["seq_name"] == sn); m = json.load(open(F / r["mapping_file"]))
    res = next(x for x in m["residues"] if x["rfam_consensus_pos"] == POS and x["alignment_state"] == "match")
    s = {"seq_name": sn, "pdb_id": r["pdb_id"], "auth_asym_id": r["auth_asym_id"], "label_asym_id": r["label_asym_id"],
         "residue_at_pos": {"label_seq_id": res["polymer_seq_position"], "author_residue_number": res["author_residue_number"], "nucleotide": res["residue_name"], "pp": res["posterior_prob_char"]}, "loops": []}
    for L in loops["loops"]:
        if sn in L["instances"] and POS in L["consensus_positions"]:
            raw = raw_loop(L["pdb_id"], L["loop_id"])
            s["loops"].append({"loop_id": L["loop_id"], "type": L["loop_type"], "motif_group": L["motif_group"], "annotation_text": L["annotation_text"],
                               "raw_bgsu_units": raw, "n_raw_units": len(raw) if raw else None, "n_nucleotides_in_loops_json": L["n_nucleotides"],
                               "consensus_positions": L["consensus_positions"], "strands": [(st["start_unit_id"], st["end_unit_id"], st["n_nucleotides"]) for st in L["strands"]],
                               "unit_at_pos": [n["unit_id"] for n in L["nucleotides"] if n["rfam_consensus_pos"] == POS]})
    samples.append(s)
out["samples"] = samples
# 4. mapping correctness: duplicates within loops, unit→residue uniqueness, author number vs label_seq_id
dups = []; unit_pos = defaultdict(set)
for L in loops["loops"]:
    units = [n["unit_id"] for n in L["nucleotides"]]
    if len(units) != len(set(units)): dups.append(L["loop_id"])
    for n in L["nucleotides"]:
        if n["rfam_consensus_pos"] is not None: unit_pos[n["unit_id"]].add(n["rfam_consensus_pos"])
out["mapping_checks"] = {"loops_with_duplicated_units": dups, "units_mapped_to_more_than_one_consensus_pos": sum(1 for v in unit_pos.values() if len(v) > 1),
                         "raw_vs_json_membership_mismatch": [L["loop_id"] for L in loops["loops"] if (raw_loop(L["pdb_id"], L["loop_id"]) or []) and len(raw_loop(L["pdb_id"], L["loop_id"])) != L["n_nucleotides"]]}
# HL motif groups at POS and consensus spans of the two overlapping loop classes
out["HL_groups_at_pos"] = {g: n for g, n in t["hairpin_loops"]["by_motif_group"].items()}
out["J3_groups_at_pos"] = {g: n for g, n in t["junction_loops_J3"]["by_motif_group"].items()}
json.dump(out, open(Path(__file__).parent / f"audit_pos{POS}.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("samples",)}, indent=1))
print(json.dumps(out["samples"][0], indent=1) if out["samples"] else "no sample")
