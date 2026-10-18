"""RF00177 (60 chains) structural-group cut candidates from the research T1481 matrices (read-only), using the project's own clustering
functions. Each chain is labelled with its computed head swivel (research, calibrated to RAD: r = 0.999) and curated state category."""
import csv, json, sys
from collections import Counter
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from rna_explorer.structure_clustering import average_linkage, partition, read_matrix, silhouette

R = Path.home() / "Desktop/aggregated_rna_clustering/experiments/rna_structural_clustering/results/rf00177_expanded"
SEL = Path(__file__).resolve().parents[3] / "data/selections"
ids = open(R / "T1481/matrices/disc_base.csv").readline().strip().split(",")[1:]
rot = {r["seq_name"]: r for r in csv.DictReader(open(R / "validation/per_structure_head_rotation_60.csv"))}
state = {(r["pdb_id"].lower(), r["ssu_label_asym_id"]): r["state_category"] for r in csv.DictReader(open(SEL / "RF00177_60_states.csv"))}
st = lambda sn: state.get((sn.split("_")[0], sn.split("_")[3]), "?")
out = {}
for metric, file, tm in (("discrepancy", "disc_base.csv", False), ("tm", "ufix_tm_user.csv", True)):
    M = read_matrix(R / "T1481/matrices" / file, ids); D = 1 - M if tm else M; np.fill_diagonal(D, 0)
    Z, order, leaves = average_linkage(D); h = sorted(z[2] for z in Z)
    res = {"top_merge_heights": [round(x, 4) for x in h[-8:]], "k": {}}
    for k in range(2, 8):
        lab = np.array(partition(Z, len(ids), k, leaves))
        groups = sorted([[ids[j] for j in range(len(ids)) if lab[j] == g] for g in set(lab)], key=len, reverse=True)
        res["k"][k] = {"gap": round(h[-(k - 1)] - h[-k], 4), "silhouette": round(silhouette(D, lab), 4),
                       "groups": [{"n": len(g), "swivel_range": [min(float(rot[s]["computed_swivel_deg"]) for s in g), max(float(rot[s]["computed_swivel_deg"]) for s in g)],
                                   "states": dict(Counter(st(s) for s in g)), "members": g if len(g) <= 12 else None} for g in groups]}
    out[metric] = res
json.dump(out, open(Path(__file__).with_name("candidates.json"), "w"), indent=1)
for m, res in out.items():
    print(f"== {m}: top merge heights {res['top_merge_heights']}")
    for k, v in res["k"].items():
        print(f"  k={k} gap={v['gap']} sil={v['silhouette']} sizes={[g['n'] for g in v['groups']]} swivel={[tuple(round(x,1) for x in g['swivel_range']) for g in v['groups']]}")
