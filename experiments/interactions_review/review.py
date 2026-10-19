"""Scientific review of the Interactions datasets for RF00059 (TPP riboswitch) and RF00177 (bacterial 16S rRNA).

Reads data/families/<ACC>/ (interactions.json, interactions_by_instance/, instances.json, ligands.json, structure_clustering.json,
mappings/) read-only and writes the CSV tables used by REVIEW.md into this directory. Usage (from anywhere):
    python3 experiments/interactions_review/review.py
"""
from __future__ import annotations
import collections, csv, json, re, statistics, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Family  # noqa: E402

OUT = Path(__file__).resolve().parent


def write(name: str, rows: list[dict]) -> None:
    if not rows:
        return
    with open(OUT / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def coverage_table(F: Family, ref: str) -> list[dict]:
    rows = []
    for sn in sorted(F.ix):
        i, r = F.ix[sn], F.inst[sn]
        k, n = F.coverage(sn)
        rows.append({"seq_name": sn, "pdb_id": r["pdb_id"].upper(), "method": r["experimental_method"], "resolution": r["resolution"],
                     "group_discrepancy": F.groups.get(sn, ""), "helix_pairs_annotated": k, "helix_pairs_evaluable": n,
                     "coverage": round(k / n, 3) if n else "", "eligible_positions": i["eligibility"].count("E"),
                     "n_interactions": i["n_interactions"], "title": r["title"]})
    return rows


def truncation_table(F: Family, min_run: int = 15) -> list[dict]:
    """chains whose BGSU annotation stops early: longest run of eligible positions with no annotation of any category"""
    rows = []
    for sn, i in F.ix.items():
        obs = json.load(open(F.dir / "interactions_by_instance" / f"{sn}.json"))["observations"]
        ann = {p for x in obs for p in (x[4], x[5]) if isinstance(p, int)}
        el = [k + 1 for k, c in enumerate(i["eligibility"]) if c == "E"]
        best = run = 0; end = None
        for p in el:
            run = 0 if p in ann else run + 1
            if run > best:
                best, end = run, p
        if el and best >= max(min_run, 0.08 * len(el)):
            rows.append({"seq_name": sn, "eligible_positions": len(el), "annotated_span_rfam": f"{min(ann)}-{max(ann)}", "longest_unannotated_eligible_run": best, "run_ends_at_rfam": end})
    return rows


# ====================================================================================================== TPP
def tpp() -> None:
    F = Family("RF00059")
    N2 = F.numbering("2gdi_1_X_A"); N8 = F.numbering("8f4o_1_A_A"); Nt = F.numbering("2cky_1_A_A")
    lab = lambda n, p: n.get(p, f"Rfam{p}")
    write("tpp_chain_coverage.csv", coverage_table(F, "2gdi_1_X_A"))
    # consensus stems (S1.. in 5' order) with 2GDI numbering and frequency
    rows = []
    for i, j in sorted((i, j) for i, j in F.ss.items() if i < j):
        p = F.by_id.get(f"basepair:{i}-{j}")
        k, n = F.freq(p, F.ix) if p else (0, sum(1 for s in F.ix if F.ix[s]["eligibility"][i - 1] == "E" and F.ix[s]["eligibility"][j - 1] == "E"))
        rows.append({"stem": F.helix[i], "rfam": f"{i}-{j}", "2GDI": f"{lab(N2, i)}-{lab(N2, j)}", "observed": k, "eligible": n,
                     "classes": json.dumps(p["classes"]) if p else "", "note": "" if p else "never annotated"})
    write("tpp_consensus_pairs.csv", rows)
    # apo (8F4O) vs holo chains of the identical sequence: every interaction with its status
    key = F.inst["8f4o_1_A_A"]["sequence_key"]
    same = [s for s in F.ix if F.inst[s]["sequence_key"] == key]
    apo = ["8f4o_1_A_A", "8f4o_1_B_B"]; holo = [s for s in same if s not in apo]
    rows = []
    for p in F.pairs:
        k, n = F.freq(p, holo)
        sa = [F.status(s, p) for s in apo]
        if n and (k / n >= 0.5 or "present" in sa):
            rows.append({"category": p["category"], "rfam": p["id"].split(":")[1], "8F4O_numbering": f"{lab(N8, p['pos_1'])}-{lab(N8, p['pos_2'])}",
                         "holo_same_sequence": f"{k}/{n}", "apo_A": sa[0], "apo_B": sa[1], "classes": json.dumps(p["classes"])})
    write("tpp_apo_vs_holo.csv", rows)
    # E. coli vs A. thaliana (all holo chains of each species)
    eco = [s for s in F.ix if (F.inst[s].get("organism_species") or "").startswith("Escherichia") and not s.startswith(("8f4o", "9c6"))]
    ath = [s for s in F.ix if (F.inst[s].get("organism_species") or "").startswith("Arabidopsis")]
    na, nb = F.ix[eco[0]]["nucleotides"], F.ix[ath[0]]["nucleotides"]
    rows = []
    for p in F.pairs:
        a, b = F.freq(p, eco), F.freq(p, ath)
        if a[1] >= 4 and b[1] >= 4 and abs(a[0] / a[1] - b[0] / b[1]) >= 0.5:
            seqdiff = na[p["pos_1"] - 1] != nb[p["pos_1"] - 1] or na[p["pos_2"] - 1] != nb[p["pos_2"] - 1]
            near = any(na[q - 1] != nb[q - 1] for q in range(max(1, p["pos_1"] - 1), min(F.ds["clen"], p["pos_1"] + 1) + 1)) or \
                any(na[q - 1] != nb[q - 1] for q in range(max(1, p["pos_2"] - 1), min(F.ds["clen"], p["pos_2"] + 1) + 1))
            rows.append({"category": p["category"], "rfam": p["id"].split(":")[1], "E.coli(2GDI)": f"{lab(N2, p['pos_1'])}-{lab(N2, p['pos_2'])}",
                         "A.thaliana(2CKY)": f"{lab(Nt, p['pos_1'])}-{lab(Nt, p['pos_2'])}", "E.coli": f"{a[0]}/{a[1]}", "A.thaliana": f"{b[0]}/{b[1]}",
                         "endpoint_nt_differs": seqdiff, "sequence_differs_within_1nt": near})
    write("tpp_species_differences.csv", rows)
    # same-sequence holo: variable base pairs with ligand and resolution of the chains lacking them
    lig = json.load(open(F.dir / "ligands.json")); ent = collections.defaultdict(set)
    for b in lig["bound_molecules"]:
        if b["in_default_track"]:
            ent[b["pdb_id"]].add(b["ccd_id"])
    rows = []
    for p in F.pairs:
        if p["category"] != "basepair":
            continue
        st = {s: F.status(s, p) for s in holo}; ev = [s for s in holo if st[s] != "nc"]; k = sum(st[s] == "present" for s in ev)
        if ev and 2 <= k <= len(ev) - 2:
            rows.append({"rfam": p["id"].split(":")[1], "8F4O_numbering": f"{lab(N8, p['pos_1'])}-{lab(N8, p['pos_2'])}", "present": f"{k}/{len(ev)}",
                         "absent_in": "; ".join(f"{F.inst[s]['pdb_id'].upper()}({','.join(sorted(ent.get(F.inst[s]['pdb_id'], [])))}, {F.inst[s]['resolution']} A, cov {F.coverage(s)[0]}/{F.coverage(s)[1]})" for s in ev if st[s] == "absent")})
    write("tpp_holo_variable.csv", rows)
    write("tpp_truncation.csv", truncation_table(F))


# ====================================================================================================== 16S
def rrna16s() -> None:
    F = Family("RF00177")
    N = F.numbering("4v9d_1_BA_Y"); R = F.rfam_of("4v9d_1_BA_Y")
    lab = lambda p: N.get(p, f"Rfam{p}")
    write("16s_chain_coverage.csv", coverage_table(F, "4v9d_1_BA_Y"))
    trunc = truncation_table(F); write("16s_truncation.csv", trunc)
    bad = {r["seq_name"] for r in trunc}
    ok = [s for s in F.ix if s not in bad]

    def pair(a, b, cat="basepair"):
        i, j = R[str(a)], R[str(b)]
        return F.by_id.get(f"{cat}:{i}-{j}") or F.by_id.get(f"{cat}:{j}-{i}"), i, j

    # literature tertiary interactions (E. coli numbering)
    known = [("central pseudoknot (h2)", (17, 918)), ("central pseudoknot (h2)", (18, 917)), ("central pseudoknot (h2)", (19, 916)),
             ("530-loop pseudoknot (h18)", (505, 526)), ("530-loop pseudoknot (h18)", (506, 525)), ("530-loop pseudoknot (h18)", (507, 524)),
             ("570/866 long-range pairs", (570, 866)), ("570/866 long-range pairs", (571, 865))]
    rows = []
    for name, (a, b) in known:
        p, i, j = pair(a, b)
        k, n = F.freq(p, ok) if p else (0, 0)
        rows.append({"element": name, "E.coli": f"{lab(i)}-{lab(j)}", "rfam": f"{i}-{j}", "observed": k, "eligible": n,
                     "classes": json.dumps(p["classes"]) if p else "", "in_rfam_ss_cons": p["in_ss_cons"] if p else ""})
    write("16s_known_tertiary.csv", rows)
    # conserved long-range non-SS_cons base pairs
    rows = []
    for p in F.pairs:
        if p["category"] == "basepair" and not p["in_ss_cons"] and abs(p["pos_1"] - p["pos_2"]) >= 50:
            k, n = F.freq(p, ok)
            if n and k / n >= 0.9:
                rows.append({"E.coli": f"{lab(p['pos_1'])}-{lab(p['pos_2'])}", "rfam": p["id"].split(":")[1], "observed": k, "eligible": n, "classes": json.dumps(p["classes"])})
    write("16s_conserved_long_range.csv", rows)
    # decoding centre per chain
    keys = {"A1492/A1493 stack": pair(1492, 1493, "stacking")[0], "A1408-A1492 cWW": pair(1408, 1492)[0], "A1492->C518 ribose": pair(1492, 518, "baseribose")[0],
            "A1493-G1494 cSH": pair(1493, 1494)[0], "G1491/A1492 stack": pair(1491, 1492, "stacking")[0]}
    lig = json.load(open(F.dir / "ligands.json")); ent = collections.defaultdict(set)
    for b in lig["bound_molecules"]:
        ent[b["pdb_id"]].add(b["ccd_id"])
    rows = []
    for s in sorted(ok, key=lambda s: F.inst[s]["pdb_id"]):
        t = F.inst[s]["title"] or ""
        a_site = "A-site tRNA" if re.search(r"A-site|A- and P|A-, ", t) else ("A/P hybrid" if "A/P" in t else "")
        rows.append({"pdb_id": F.inst[s]["pdb_id"].upper(), "seq_name": s, "title_A_site": a_site, "ligands": ",".join(sorted(ent.get(F.inst[s]["pdb_id"], []))),
                     **{k: F.status(s, v) if v else "never annotated" for k, v in keys.items()}, "title": t})
    write("16s_decoding_centre.csv", rows)
    # structural groups (head rotation) with coverage control
    cov = {s: (lambda k, n: k / n if n else 0)(*F.coverage(s)) for s in ok}
    g1 = [s for s in ok if F.groups.get(s) == 1]; g2 = [s for s in ok if F.groups.get(s) == 2]
    h1 = [s for s in g1 if cov[s] >= 0.9]; h2 = [s for s in g2 if cov[s] >= 0.9]
    rows = []
    for p in F.pairs:
        if p["category"] != "basepair":
            continue
        a, b = F.freq(p, g1), F.freq(p, g2)
        if a[1] >= 20 and b[1] >= 6 and abs(b[0] / b[1] - a[0] / a[1]) >= 0.4:
            c, d = F.freq(p, h1), F.freq(p, h2)
            rows.append({"E.coli": f"{lab(p['pos_1'])}-{lab(p['pos_2'])}", "rfam": p["id"].split(":")[1], "classes": json.dumps(p["classes"]),
                         "group1": f"{a[0]}/{a[1]}", "group2": f"{b[0]}/{b[1]}", "group1_cov>=0.9": f"{c[0]}/{c[1]}", "group2_cov>=0.9": f"{d[0]}/{d[1]}"})
    write("16s_group_differences.csv", rows)
    # chains with incomplete source annotation: after the fix (annotation_gaps, schema 1.3) the interactions touching the gap are not
    # comparable for that chain; report how many interactions that concerns and the stored gap
    rows = []
    for s in bad:
        gaps = F.ix[s].get("annotation_gaps", [])
        touched = [p for p in F.pairs if any(a <= q <= b for a, b in gaps for q in (p["pos_1"], p["pos_2"]))]
        rows.append({"seq_name": s, "stored_annotation_gaps": json.dumps(gaps), "interactions_not_comparable_because_of_gap": len(touched),
                     "status_now": "; ".join(sorted({F.status(s, p) for p in touched}))})
    write("16s_truncation_impact.csv", rows)


if __name__ == "__main__":
    tpp(); rrna16s()
    print("written:", sorted(p.name for p in OUT.glob("*.csv")))
