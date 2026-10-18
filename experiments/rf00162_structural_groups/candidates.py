"""Structural-group cut candidates (family and research run as arguments) (read-only on the research worktree). Uses the project's own clustering functions
(rna_explorer.structure_clustering): average linkage, k-cut, silhouette. TM values read from the research pairs_long.csv
(ufix_tm_norm_n), which equals the run's matrices/ufix_tm_user.csv for every RF00162 pair (max difference 0.0, checked 2026-10-17)."""
import csv, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rna_explorer.structure_clustering import average_linkage, cut_tree, partition, silhouette, read_matrix, medoid

ACC, RUN = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("RF00162", "T91")
R = Path.home() / f"Desktop/aggregated_rna_clustering/experiments/rna_structural_clustering/results/{ACC.lower()}"
fam = Path(__file__).resolve().parents[2] / f"data/families/{ACC}"
ids = open(R / f"{RUN}/matrices/disc_base.csv").readline().strip().split(",")[1:]
D = read_matrix(R / f"{RUN}/matrices/disc_base.csv", ids)
pos = {n: i for i, n in enumerate(ids)}; T = np.eye(len(ids))
for r in csv.DictReader(open(R / f"{RUN}/pairs_long.csv")):
    a, b = pos[r["seq_name_1"]], pos[r["seq_name_2"]]; T[a, b] = T[b, a] = float(r["ufix_tm_norm_n"])
inst = {i["seq_name"]: i for i in json.load(open(fam / "instances.json"))}
lig = {}
for o in json.load(open(fam / "ligands.json"))["bound_molecules"]:
    if o["ligand_class"] == "ligand": lig.setdefault(o["pdb_id"], set()).add(o["ccd_id"])
def tag(sn):
    i = inst[sn]; sp = (i.get("organism") or "?").split()[0]
    return f"{sn.split('_')[0].upper()} {sn.split('_')[2]} [{sp}; {'+'.join(sorted(lig.get(i['pdb_id'], ['apo'])))}]"
out = {}
for name, M in (("discrepancy", D), ("tm", 1 - T)):
    Z, order, leaves = average_linkage(M)
    heights = sorted(z[2] for z in Z)
    res = {"top_merge_heights": [round(h, 4) for h in heights[-8:]], "k": {}}
    for k in range(2, 7):
        lab = partition(Z, len(ids), k, leaves)
        groups = sorted([[ids[j] for j in range(len(ids)) if lab[j] == g] for g in set(lab)], key=len, reverse=True)
        res["k"][k] = {"gap": round(heights[-(k - 1)] - heights[-k], 4), "silhouette": None if (s := silhouette(M, np.array(lab))) is None else round(s, 4),
                       "sizes": [len(g) for g in groups], "small_groups": [[tag(s) for s in g] for g in groups[1:]]}
    out[name] = res
json.dump(out, open(Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).with_name("candidates.json"), "w"), indent=1)
for name, res in out.items():
    print(f"== {name}: top merge heights {res['top_merge_heights']}")
    for k, v in res["k"].items(): print(f"  k={k} gap={v['gap']} sil={v['silhouette']} sizes={v['sizes']}\n     " + "\n     ".join(" | ".join(g) for g in v["small_groups"]))
