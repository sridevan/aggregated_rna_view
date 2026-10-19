"""RF00167 RNA-RNA interaction investigation, stage 2: conservation with eligible denominators, structural regions, structural groups,
ligand status. Reads stage-1 outputs (RF00167_mapped_interactions.csv, base_completeness.csv) and the family's existing data files
(read-only). Writes RF00167_interaction_conservation.csv, RF00167_cluster_interaction_comparison.csv, RF00167_representative_comparison.csv,
RF00167_ligand_interaction_comparison.csv, instance_metadata.csv and stage2_summary.json in this directory.

Eligibility of an instance for a consensus position pair (i, j): both positions are match-state residues of the instance, both bases are
complete (FR3D standard heavy atoms, model 1), and the instance's chain has BGSU/FR3D annotations. Observed = at least one non-near FR3D
interaction of the category (base pair / stacking) between i and j. Near interactions (FR3D 'n' prefix; stacking only in these files) are
counted separately and never as observed. Insertion and flank endpoints are kept in the mapped table but not aggregated by consensus pair.
"""
import csv, json, collections, math
from pathlib import Path
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
FAM = ROOT / "data/families/RF00167"


def ss_pairs(ss):
    st, pairs = collections.defaultdict(list), set()
    for i, ch in enumerate(ss, 1):
        for o, c in ("()", "<>", "[]", "{}"):
            if ch == o: st[o].append(i)
            elif ch == c: pairs.add((st[o].pop(), i))
    return pairs


def main():
    inst = json.load(open(FAM / "instances.json")); by_sn = {i["seq_name"]: i for i in inst}
    fam = json.load(open(FAM / "family.json")); SS = ss_pairs(fam["ss_cons"])
    ref = json.load(open(FAM / "reference_structural_elements.json")); region = {p: (e or "unassigned") for p, e in enumerate(ref["track"], 1)}
    sc = json.load(open(FAM / "structure_clustering.json")); lig = json.load(open(FAM / "ligands.json"))
    obs = list(csv.DictReader(open(OUT / "RF00167_mapped_interactions.csv")))
    bc = collections.defaultdict(dict); nt = collections.defaultdict(dict)
    for r in csv.DictReader(open(OUT / "base_completeness.csv")):
        if r["alignment_state"] == "match":
            bc[r["seq_name"]][int(r["rfam_pos"])] = (r["has_coordinates"] == "True", r["base_complete"] == "True")
    for sn in by_sn:
        m = json.load(open(FAM / "mappings" / f"{sn}.json"))
        for r in m["residues"]:
            if r["alignment_state"] == "match": nt[sn][r["rfam_consensus_pos"]] = r["nucleotide_parent"]
    annotated = {r["seq_name"] for r in obs}                       # chains with at least one FR3D interaction in the BGSU files
    sns = sorted(by_sn)

    def eligible(sn, i, j):
        if sn not in annotated: return "annotation_unavailable"
        a, b = bc[sn].get(i), bc[sn].get(j)
        if a is None or b is None: return "position_not_in_chain"      # deleted / outside the aligned region
        if not (a[0] and b[0]): return "unmodelled"
        if not (a[1] and b[1]): return "base_incomplete"
        return "eligible"

    # observations at consensus pairs
    cons = collections.defaultdict(lambda: collections.defaultdict(list))   # (cat, i, j) -> sn -> [classes]
    near = collections.defaultdict(set); ins_obs = collections.Counter()
    for r in obs:
        if r["state_1"] != "consensus" or r["state_2"] != "consensus":
            ins_obs[(r["category"], "insertion" if "insertion" in (r["state_1"], r["state_2"]) else "outside_alignment")] += 1
            continue
        k = (r["category"], int(r["rfam_pos_1"]), int(r["rfam_pos_2"]))
        if r["near"] == "True": near[k].add(r["seq_name"])
        else: cons[k][r["seq_name"]].append(r["fr3d_class"])

    # instance metadata: structural groups (both metrics), ligands, sequence
    groups = {m: sc["metrics"][m]["clustering"]["assignments"] for m in ("discrepancy", "tm")}
    reps = {m: {g["group_id"]: g["representative"] for g in sc["metrics"][m]["clustering"]["groups"]} for m in ("discrepancy", "tm")}
    ent_lig = collections.defaultdict(set)
    for b in lig["bound_molecules"]:
        if b["in_default_track"]: ent_lig[b["pdb_id"]].add(b["ccd_id"])
    # ligand-contact positions: contacted by the entry's default-track ligand in >= 50% of ligand-bound entries with that position observed
    contact = collections.Counter(); lig_entries_obs = collections.Counter()
    pos_entries = collections.defaultdict(set)
    for b in lig["bound_molecules"]:
        if b["in_default_track"]:
            for p in b["consensus_positions"]: pos_entries[p].add(b["pdb_id"])
    bound_entries = {e for e, s in ent_lig.items() if s}
    obs_entries = collections.defaultdict(set)
    for sn in sns:
        for p, (hc, _) in bc[sn].items():
            if hc: obs_entries[p].add(by_sn[sn]["pdb_id"])
    pocket = sorted(p for p in range(1, fam["clen"] + 1) if len(bound_entries & obs_entries[p]) and len(pos_entries[p]) / len(bound_entries & obs_entries[p]) >= 0.5)
    meta = []
    for sn in sns:
        i = by_sn[sn]
        meta.append({"seq_name": sn, "pdb_id": i["pdb_id"].upper(), "organism": i.get("organism_species") or i.get("organism"), "sequence_key": i["sequence_key"],
                     "group_discrepancy": groups["discrepancy"].get(sn, ""), "group_tm": groups["tm"].get(sn, ""),
                     "ligands": ";".join(sorted(ent_lig.get(i["pdb_id"], ()))) or "none", "annotated": sn in annotated,
                     "seq_at_pocket": "".join(nt[sn].get(p, "-") for p in pocket), "method": i.get("experimental_method"), "resolution": i.get("resolution")})
    with open(OUT / "instance_metadata.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(meta[0].keys())); w.writeheader(); w.writerows(meta)
    md = {m["seq_name"]: m for m in meta}

    def context(cat, i, j):
        ri, rj = region[i], region[j]
        if (i, j) in SS: return "SS_cons helix pair"
        if cat == "stacking" and j - i == 1: return "adjacent stack"
        if {ri, rj} == {"L2", "L3"}: return "L2-L3 (loop-loop tertiary)"
        if ri == rj: return f"within {ri}"
        return f"{ri} x {rj}"

    def tally(cat, i, j, members):
        st = collections.Counter(); classes = collections.Counter(); ent_e, ent_o = set(), set(); obs_inel = 0
        for sn in members:
            e = eligible(sn, i, j); st[e] += 1
            cl = cons[(cat, i, j)].get(sn)
            if e == "eligible":
                ent_e.add(by_sn[sn]["pdb_id"])
                if cl:
                    st["observed"] += 1; ent_o.add(by_sn[sn]["pdb_id"])
                    for c in sorted(set(cl)): classes[c] += 1
                elif sn in near[(cat, i, j)]: st["near_only"] += 1
                else: st["not_annotated"] += 1
            elif cl: obs_inel += 1
        n, k = st["eligible"], st["observed"]
        return {"n_eligible": n, "n_observed": k, "frequency": round(k / n, 4) if n else None, "n_near_only": st["near_only"], "n_not_annotated": st["not_annotated"],
                "n_entries_eligible": len(ent_e), "n_entries_observed": len(ent_o), "classes": ";".join(f"{c}:{v}" for c, v in classes.most_common()),
                "n_missing_position": st["position_not_in_chain"], "n_unmodelled": st["unmodelled"], "n_base_incomplete": st["base_incomplete"],
                "n_annotation_unavailable": st["annotation_unavailable"], "n_observed_but_ineligible": obs_inel}

    rows = []
    for (cat, i, j), per in sorted(cons.items()):
        t = tally(cat, i, j, sns)
        ex = sorted({by_sn[sn]["pdb_id"].upper() for sn in per})[:6]
        rows.append({"category": cat, "rfam_pos_1": i, "rfam_pos_2": j, "consensus_nt": f"{fam['consensus_sequence'][i-1]}{fam['consensus_sequence'][j-1]}",
                     "region_1": region[i], "region_2": region[j], "context": context(cat, i, j), "in_ss_cons": (i, j) in SS,
                     "ligand_contact_position": (i in pocket) or (j in pocket), **t, "example_entries": " ".join(ex)})
    with open(OUT / "RF00167_interaction_conservation.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    # structural groups: per metric, each group vs the rest (eligible denominators); Fisher exact as a descriptive aid only
    crow = []
    for m in ("discrepancy", "tm"):
        gids = sorted(set(groups[m].values()))
        for (cat, i, j), _ in sorted(cons.items()):
            per = {g: tally(cat, i, j, [sn for sn in sns if groups[m].get(sn) == g]) for g in gids}
            for g in gids:
                rest = tally(cat, i, j, [sn for sn in sns if sn in groups[m] and groups[m][sn] != g])
                a, n1, b, n2 = per[g]["n_observed"], per[g]["n_eligible"], rest["n_observed"], rest["n_eligible"]
                if n1 == 0 or n2 == 0: continue
                f1, f2 = a / n1, b / n2
                if abs(f1 - f2) < 0.2 or (a + b) == 0: continue
                p = fisher_exact([[a, n1 - a], [b, n2 - b]])[1]
                crow.append({"metric": m, "group": g, "group_size": sum(1 for v in groups[m].values() if v == g), "category": cat, "rfam_pos_1": i, "rfam_pos_2": j,
                             "context": context(cat, i, j), "group_observed": a, "group_eligible": n1, "group_frequency": round(f1, 3),
                             "rest_observed": b, "rest_eligible": n2, "rest_frequency": round(f2, 3), "difference": round(f1 - f2, 3), "fisher_p_descriptive": round(p, 6),
                             "group_classes": per[g]["classes"], "rest_classes": rest["classes"]})
    crow.sort(key=lambda r: (r["metric"], r["group"], r["fisher_p_descriptive"], -abs(r["difference"])))
    with open(OUT / "RF00167_cluster_interaction_comparison.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(crow[0].keys())); w.writeheader(); w.writerows(crow)

    # representatives: pairwise comparison of their consensus-mapped interactions
    def iset(sn):
        return {(c, i, j): set(per[sn]) for (c, i, j), per in cons.items() if sn in per}
    rrows = []
    for m in ("discrepancy", "tm"):
        rl = sorted(reps[m].items())
        for x in range(len(rl)):
            for y in range(x + 1, len(rl)):
                (ga, a), (gb, b) = rl[x], rl[y]
                A, B = iset(a), iset(b)
                for k in sorted(set(A) | set(B)):
                    cat, i, j = k
                    if k in A and k in B: s = "conserved" if A[k] == B[k] else "reclassified"
                    elif k in A: s = f"only_{ga}" if eligible(b, i, j) == "eligible" else f"only_{ga} (other not eligible: {eligible(b, i, j)})"
                    else: s = f"only_{gb}" if eligible(a, i, j) == "eligible" else f"only_{gb} (other not eligible: {eligible(a, i, j)})"
                    rrows.append({"metric": m, "group_a": ga, "rep_a": a, "group_b": gb, "rep_b": b, "category": cat, "rfam_pos_1": i, "rfam_pos_2": j, "context": context(cat, i, j),
                                  "status": s, "classes_a": ";".join(sorted(A.get(k, ()))), "classes_b": ";".join(sorted(B.get(k, ())))})
    with open(OUT / "RF00167_representative_comparison.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rrows[0].keys())); w.writeheader(); w.writerows(rrows)

    summary = {"pocket_positions": pocket, "n_consensus_pairs": len(cons), "insertion_or_flank_interactions": {f"{a}|{b}": n for (a, b), n in ins_obs.items()},
               "annotated_instances": len(annotated), "instances": len(sns), "ss_cons_pairs": sorted(SS)}
    json.dump(summary, open(OUT / "stage2_summary.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "ss_cons_pairs"}, indent=1))


if __name__ == "__main__":
    main()
