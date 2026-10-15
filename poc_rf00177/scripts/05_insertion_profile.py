"""Step 8: insertion profile per consensus position: instance-level (8 structures) vs family-level (Rfam seed, 99 sequences); tables + PNG."""
import csv, json, sys, collections
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
rows = list(csv.DictReader(open(RES/"rf00177_residue_mapping.csv")))
summ = json.load(open(RES/"rf00177_summary.json")); CLEN = summ["CLEN"]
chains = json.load(open(Path(sys.argv[3])/"chains.json")); order = [c["pdb_id"] for c in chains]; label = {c["pdb_id"]: f"{c['pdb_id']}/{c['auth_chain_id']} {c['cif_organism'].split(' ')[0][:12]}" for c in chains}
# ---- family-level from seed ----
seqs, gc = collections.defaultdict(str), collections.defaultdict(str)
for line in open(DATA/"cm/RF00177.seed.sto"):
    line = line.rstrip("\n")
    if line.startswith("#=GC"): _, t, x = line.split(None, 2); gc[t] += x
    elif line and not line.startswith("#") and line != "//": n, x = line.split(None, 1); seqs[n] += x.strip()
RF = gc["RF"]; assert all(len(s) == len(RF) for s in seqs.values())
col_cons = []; n = 0
for ch in RF:
    if ch in ".~": col_cons.append(None)
    else: n += 1; col_cons.append(n)
assert n == CLEN, (n, CLEN)
prev = []; last = 0
for c in col_cons:
    if c: last = c
    prev.append(last)
seed_ins = collections.defaultdict(list)   # anchor -> list of insertion lengths (one per seed seq, >0)
for s in seqs.values():
    cnt = collections.Counter()
    for col, a in enumerate(s):
        if col_cons[col] is None and a not in "-._~": cnt[prev[col]] += 1
    for k, v in cnt.items(): seed_ins[k].append(v)
NSEED = len(seqs)
seed_stats = {k: {"seed_fraction_with_insertion": round(len(v)/NSEED, 3), "seed_max_len": max(v), "seed_mean_len": round(sum(v)/len(v), 1)} for k, v in seed_ins.items()}
# ---- instance-level ----
inst = collections.defaultdict(lambda: collections.defaultdict(lambda: {"n_insert": 0, "n_local_end": 0, "pos": [], "auth": [], "coords": 0}))
for r in rows:
    if r["alignment_state"] in ("insert", "local_end"):
        d = inst[r["pdb_id"]][int(r["insertion_anchor"])]
        d["n_insert" if r["alignment_state"] == "insert" else "n_local_end"] += 1
        d["pos"].append(int(r["polymer_seq_position"])); d["auth"].append(r["author_residue_number"] + (r["author_insertion_code"] or "")); d["coords"] += int(r["has_coordinates"])
long_rows, wide = [], collections.defaultdict(dict)
for pdb in order:
    ch = [c for c in chains if c["pdb_id"] == pdb][0]
    for anchor in sorted(inst[pdb]):
        d = inst[pdb][anchor]; tot = d["n_insert"] + d["n_local_end"]; ss = seed_stats.get(anchor, {})
        long_rows.append({"pdb_id": pdb, "auth_asym_id": ch["auth_chain_id"], "label_asym_id": ch["label_asym_id"], "insertion_anchor": anchor, "n_insert": d["n_insert"], "n_local_end": d["n_local_end"], "total_inserted_residues": tot,
                          "polymer_pos_from": min(d["pos"]), "polymer_pos_to": max(d["pos"]), "auth_from": d["auth"][0], "auth_to": d["auth"][-1], "n_with_coordinates": d["coords"],
                          "seed_fraction_with_insertion": ss.get("seed_fraction_with_insertion", 0), "seed_max_len": ss.get("seed_max_len", 0), "seed_mean_len": ss.get("seed_mean_len", 0)})
        wide[anchor][pdb] = tot
with open(RES/"rf00177_insertion_profile.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(long_rows[0].keys())); w.writeheader(); w.writerows(long_rows)
with open(RES/"rf00177_insertion_profile_wide.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["insertion_anchor"] + order + ["n_instances_with_insertion", "seed_fraction_with_insertion", "seed_max_len"])
    for a in sorted(set(wide) | set(seed_stats)):
        w.writerow([a] + [wide[a].get(p, 0) for p in order] + [sum(1 for p in order if wide[a].get(p, 0)), seed_stats.get(a, {}).get("seed_fraction_with_insertion", 0), seed_stats.get(a, {}).get("seed_max_len", 0)])
site_summary = []
for a in sorted(wide):
    ps = [p for p in order if wide[a].get(p)]
    site_summary.append({"insertion_anchor": a, "n_instances": len(ps), "instances": ";".join(f"{p}:{wide[a][p]}" for p in ps), "bacterial_or_chloroplast_have_it": int(any(p in ("7k00", "4y4o", "6eri") for p in ps)),
                         "seed_fraction_with_insertion": seed_stats.get(a, {}).get("seed_fraction_with_insertion", 0), "seed_max_len": seed_stats.get(a, {}).get("seed_max_len", 0)})
with open(RES/"rf00177_insertion_sites_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(site_summary[0].keys())); w.writeheader(); w.writerows(site_summary)
longest = sorted(long_rows, key=lambda r: -r["total_inserted_residues"])[:15]
with open(RES/"rf00177_longest_insertions.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(longest[0].keys())); w.writeheader(); w.writerows(longest)
stats = {"seed_sequences": NSEED, "seed_anchors_with_any_insertion": len(seed_stats), "seed_anchors_insertion_in_>=50%_seqs": sum(1 for v in seed_stats.values() if v["seed_fraction_with_insertion"] >= 0.5),
         "instance_anchors_total": len(wide), "instance_anchors_shared_by_>=2": sum(1 for a in wide if len(wide[a]) >= 2), "instance_anchors_absent_from_seed": sum(1 for a in wide if a not in seed_stats),
         "instance_anchors_with_>=20_residues": sum(1 for r in long_rows if r["total_inserted_residues"] >= 20)}
json.dump(stats, open(RES/"rf00177_insertion_profile_stats.json", "w"), indent=1); print(json.dumps(stats, indent=1))
for r in longest[:10]: print("  ", {k: r[k] for k in ["pdb_id", "auth_asym_id", "insertion_anchor", "n_insert", "n_local_end", "polymer_pos_from", "polymer_pos_to", "auth_from", "auth_to", "n_with_coordinates", "seed_fraction_with_insertion", "seed_max_len"]})
# ---- plot: small multiples, one panel per instance + family panel; insert (blue) vs local-end (orange) ----
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"; C_INS, C_EL, C_SEED = "#2a78d6", "#eb6834", "#1baf7a"
fig, axes = plt.subplots(len(order) + 1, 1, figsize=(12, 1.25 * (len(order) + 1) + 1.2), sharex=True, facecolor=SURF)
ax = axes[0]; ax.set_facecolor(SURF)
xs = sorted(seed_stats); ax.bar(xs, [seed_stats[a]["seed_fraction_with_insertion"] for a in xs], width=3, color=C_SEED, linewidth=0)
ax.set_ylim(0, 1); ax.set_ylabel("seed frac.", color=INK2, fontsize=8); ax.text(0.005, 0.85, f"Family level: fraction of {NSEED} Rfam seed sequences with an insertion after each consensus position", transform=ax.transAxes, fontsize=8.5, color=INK)
ymax = max(r["total_inserted_residues"] for r in long_rows)
for ax, pdb in zip(axes[1:], order):
    ax.set_facecolor(SURF); d = inst[pdb]; xs = sorted(d)
    ax.bar(xs, [d[a]["n_insert"] for a in xs], width=3, color=C_INS, linewidth=0, label="insert (.)")
    ax.bar(xs, [d[a]["n_local_end"] for a in xs], bottom=[d[a]["n_insert"] for a in xs], width=3, color=C_EL, linewidth=0, label="local end (~)")
    ax.set_ylim(0, ymax * 1.15); ax.set_ylabel("nt", color=INK2, fontsize=8)
    top = max(xs, key=lambda a: d[a]["n_insert"] + d[a]["n_local_end"]) if xs else None
    sm = summ["per_seq"][[k for k in summ["per_seq"] if k.startswith(pdb)][0]]
    ax.text(0.005, 0.8, f"{label[pdb]}   match {sm['n_match']}/{CLEN}  insert {sm['n_insert']}  local-end {sm['n_local_end']}  deleted {sm['n_deleted_in_span']}  uncovered {sm['n_uncovered_outside_span']}", transform=ax.transAxes, fontsize=8.5, color=INK)
    if top: ax.annotate(f"{d[top]['n_insert'] + d[top]['n_local_end']} nt after {top}", (top, d[top]["n_insert"] + d[top]["n_local_end"]), xytext=(4, -2), textcoords="offset points", fontsize=7.5, color=INK2)
for ax in axes:
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8); ax.grid(axis="y", color=GRID, linewidth=0.6); ax.set_axisbelow(True)
axes[-1].set_xlabel("RF00177 consensus position (insertion anchor = preceding consensus position)", color=INK2, fontsize=9); axes[-1].set_xlim(0, CLEN + 1)
h, l = axes[1].get_legend_handles_labels(); fig.legend(h, l, loc="upper right", frameon=False, fontsize=8, ncol=2, bbox_to_anchor=(0.99, 0.995))
fig.suptitle("Insertion length per consensus position: RF00177 seed (family) vs 8 PDB instances (structural observations)", x=0.01, ha="left", fontsize=10.5, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.975)); fig.savefig(RES/"rf00177_insertion_profile.png", dpi=150, facecolor=SURF); print("png written")
