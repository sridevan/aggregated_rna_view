"""Per-chain presence of the add riboswitch junction interactions discussed by Stagno et al. 2017 (Nature 541:242), for the 13 V. vulnificus
chains of one sequence (plus 1Y26, a different construct of the same riboswitch, as a reference). Writes RF00167_vvulnificus_junction_states.csv. 1 = FR3D interaction observed (class in brackets), 0 = both bases
complete and annotated but no interaction, '.' = not eligible."""
import csv, collections
meta = {r["seq_name"]: r for r in csv.DictReader(open("instance_metadata.csv"))}
obs = collections.defaultdict(dict)
for r in csv.DictReader(open("RF00167_mapped_interactions.csv")):
    if r["state_1"] == r["state_2"] == "consensus" and r["near"] == "False":
        obs[r["seq_name"]][(r["category"], int(r["rfam_pos_1"]), int(r["rfam_pos_2"]))] = r["fr3d_class"]
bc = collections.defaultdict(dict)
for r in csv.DictReader(open("base_completeness.csv")):
    if r["alignment_state"] == "match": bc[r["seq_name"]][int(r["rfam_pos"])] = r["base_complete"] == "True"
# (category, consensus i, j, label in add numbering, state in which the paper describes it)
COLS = [("basepair", 26, 56, "U22-A52", "bound triple A73.(A52-U22)"), ("basepair", 53, 82, "U49-A76", "bound triple U49.(U20-A76)"),
        ("basepair", 54, 81, "C50-U75", "bound triple C50.(U75-A21)"), ("basepair", 27, 50, "A23-G46", "bound junction"),
        ("basepair", 51, 56, "U47-A52", "apo WC pair"), ("basepair", 53, 79, "U49-A73", "apo WC pair"), ("basepair", 50, 57, "G46-C53", "apo WC pair"),
        ("stacking", 54, 56, "C50/A52", "apo: C50 stacked on U47-A52"), ("stacking", 51, 54, "U47/C50", "apo: C50 stacked on U47-A52"),
        ("basepair", 27, 80, "A23-U74", "apo1: A23 Hoogsteen with U74"), ("basepair", 52, 81, "U48-U75", "apo1: U48-U75 wobble"), ("basepair", 52, 80, "U48-U74", "apo2: U48-U74 wobble")]
key = meta["5e54_1_A_A"]["sequence_key"]                               # the 13 chains of one V. vulnificus add sequence
sns = sorted(s for s, m in meta.items() if m["sequence_key"] == key) + ["1y26_1_X_A"]   # + 1Y26 (same riboswitch, different construct sequence) for reference
with open("RF00167_vvulnificus_junction_states.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["seq_name", "ligands", "group_tm"] + [f"{c[3]} ({c[0]}; {c[4]})" for c in COLS])
    for s in sns:
        cells = []
        for cat, i, j, *_ in COLS:
            if not (bc[s].get(i) and bc[s].get(j)): cells.append(".")
            else:
                c = obs[s].get((cat, i, j)); cells.append(f"1 [{c}]" if c else "0")
        w.writerow([s, meta[s]["ligands"], meta[s]["group_tm"]] + cells)
        print(f"{s:<12} {meta[s]['ligands']:<5} tm{meta[s]['group_tm']:<2}", " ".join(f"{x:>8}" for x in cells))
print(" " * 22, " ".join(f"{c[3]:>8}" for c in COLS))
