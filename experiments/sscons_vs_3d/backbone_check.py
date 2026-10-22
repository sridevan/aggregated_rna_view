"""Do base-ribose / base-phosphate (and stacking) contacts between i and j cover SS_cons pairs without a base pair? (investigation only)

Input: followup/near/ss_cons_near.csv (every evaluable SS_cons observation without a true BGSU base pair, with the local near status)
and the production BGSU observations in data/families/<ACC>/interactions_by_instance/ (basepair, baseribose, basephosphate, stacking;
BR and BPh files include near classes). Only contacts between exactly i and j (either direction) count.

Outputs in followup/backbone/: obs_contacts.csv (one row per observation), backbone_summary.json.
Usage: python3 experiments/sscons_vs_3d/backbone_check.py
"""
import csv, json, sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).parent; OUT = HERE / "followup" / "backbone"
sys.path.insert(0, str(HERE))
from followup import load, evaluable  # noqa: E402
from analyse import ss_pairs  # noqa: E402

ACCS = ["RF00167", "RF00162", "RF00059", "RF00005", "RF00011", "RF00177"]


def contacts(acc):
    """seq_name -> {(i,j) unordered consensus pair -> {category: set(classes)}} for non-basepair categories"""
    d = ROOT / f"data/families/{acc}"
    out = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    for f in (d / "interactions_by_instance").glob("*.json"):
        for o in json.load(open(f))["observations"]:
            cat, cls, _, _, r1, r2 = o[:6]
            if cat == "basepair" or not (isinstance(r1, int) and isinstance(r2, int)): continue
            out[f.stem][tuple(sorted((r1, r2)))][cat].add(cls)
    return out


def tier(c, near_status):
    """strongest evidence that i and j are in contact, given no true base pair"""
    br = c.get("baseribose", set()) | c.get("basephosphate", set())
    if near_status == "local_true_pair": return "local_true_pair"
    if near_status == "near_pair": return "near_pair"
    if any(not x.startswith("n") for x in br): return "true_BR_or_BPh"
    if br: return "near_BR_or_BPh_only"
    if any(not x.startswith("n") for x in c.get("stacking", set())): return "stacking_only"
    return "nothing"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    near = {(r["acc"], r["seq_name"], int(r["pos_1"]), int(r["pos_2"])): r for r in csv.DictReader(open(HERE / "followup/near/ss_cons_near.csv"))}
    rows, summary = [], {"inputs": ["followup/near/ss_cons_near.csv", "data/families/<ACC>/interactions_by_instance (BGSU FR3D)"], "families": {}}
    for acc in ACCS:
        fam, it, chains = load(acc); con = contacts(acc)
        tiers, total, pp = Counter(), 0, defaultdict(Counter)
        for i, j in ss_pairs(fam["ss_cons"]):
            for sn, ch in chains.items():
                if not evaluable(ch, i, j): continue
                total += 1
                if ch["bp"].get((i, j)):
                    pp[(i, j)]["true_bp"] += 1; continue
                n = near.get((acc, sn, i, j))
                if n is None:  # chain not in the local run (should not happen)
                    pp[(i, j)]["unknown"] += 1; continue
                c = con[sn].get((i, j), {})
                t = tier(c, n["status"]); tiers[t] += 1; pp[(i, j)][t] += 1
                rows.append({"acc": acc, "pos_1": i, "pos_2": j, "seq_name": sn, "nts": n["nts"], "near_status": n["status"], "near_classes": n["near_classes"],
                             "baseribose": ",".join(sorted(c.get("baseribose", ()))), "basephosphate": ",".join(sorted(c.get("basephosphate", ()))),
                             "stacking": ",".join(sorted(c.get("stacking", ()))), "tier": t, "geometry": n["geometry"]})
        true_bp = sum(v["true_bp"] for v in pp.values())
        def sup(keys):
            return round(sum(sum(v[k] for k in keys) for v in pp.values()) / total, 3)
        summary["families"][acc] = {"evaluable_observations": total, "true_bp": true_bp, "without_true_bp": dict(tiers),
                                    "support_true_bp": sup(["true_bp"]),
                                    "support_plus_near_bp": sup(["true_bp", "near_pair", "local_true_pair"]),
                                    "support_plus_BR_BPh_true": sup(["true_bp", "near_pair", "local_true_pair", "true_BR_or_BPh"]),
                                    "support_plus_BR_BPh_near": sup(["true_bp", "near_pair", "local_true_pair", "true_BR_or_BPh", "near_BR_or_BPh_only"]),
                                    "pairs_gaining_from_BR_BPh": sorted(f"{i}-{j}:{v['true_BR_or_BPh'] + v['near_BR_or_BPh_only']}/{sum(v.values())}"
                                                                        for (i, j), v in pp.items() if v["true_BR_or_BPh"] + v["near_BR_or_BPh_only"])}
        print(acc, json.dumps(summary["families"][acc]), file=sys.stderr)
    with open(OUT / "obs_contacts.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    json.dump(summary, open(OUT / "backbone_summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
