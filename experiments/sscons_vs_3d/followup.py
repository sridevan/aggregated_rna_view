"""Follow-up validation of the SS_cons vs 3D comparison (investigation only; see FOLLOWUP_REPORT.md).

Reads the same unchanged inputs as analyse.py plus instances.json, the validated mappings (mappings/<chain>.json) and the cmalign output
(work/aligned.sto). Writes into followup/:
  support_three_measures.csv   per SS_cons pair: any FR3D base pair / any cWW / canonical cWW over the same eligible chains, instance-,
                               sequence- and species-weighted
  focus_pair_classes.csv       class distributions of the geometry examples (SAM-I 46-69, RNase P G.A tandem, ...)
  alt_partner_chains.csv       per chain for each alternative-partner case: residues (label_seq_id, author id, nt, PP), classes, network
  alt_partner_alignment.txt    local cmalign rows (unique sequences) with SS_cons / RF around each case
  unpaired_causes.csv          every SS_cons pair observation without any base pair, with completeness diagnostics
  unannotated_runs.csv         lengths of runs of eligible positions without any BGSU annotation (all chains)
  redundancy.csv               chains / unique sequences / near-identical sequence clusters / species / entries per family
  outlier_helices.csv          per-pair detail for the outlier chains
  hairpin_16s.csv              the 16S 75-86 hairpin per chain (method, resolution, status, eligibility)
  extra_pair_attributes.csv    recurrent 3D-only pairs with independent boolean attributes (overlap check of the kinds)
  followup_summary.json
"""
import csv, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).parent / "followup"; OUT.mkdir(exist_ok=True)
FAMS = ["RF00167", "RF00162", "RF00059", "RF00005", "RF00011", "RF00177"]
CANON = {"AU", "UA", "GC", "CG", "GU", "UG"}
sys.path.insert(0, str(Path(__file__).parent))
from analyse import ss_pairs, crosses, helices  # noqa: E402


def load(acc):
    d = ROOT / f"data/families/{acc}"
    fam = json.load(open(d / "family.json")); it = json.load(open(d / "interactions.json"))
    meta = {x["seq_name"]: x for x in json.load(open(d / "instances.json"))}
    insts = [i for i in it["instances"] if i["annotation"] == "available"]
    chains = {}
    for inst in insts:
        sn = inst["seq_name"]
        m = json.load(open(d / meta[sn]["mapping_file"]))
        by_pos = {r["rfam_consensus_pos"]: r for r in m["residues"] if r["alignment_state"] == "match"}
        by_label = {r["polymer_seq_position"]: r for r in m["residues"]}
        obs = json.load(open(d / f"interactions_by_instance/{sn}.json"))["observations"]
        bp = defaultdict(set)        # unordered (rfam_or_insert, rfam_or_insert) -> non-near base-pair classes (insert endpoints as 'a^o')
        partners = defaultdict(set)  # consensus pos -> {(partner, class)} for non-near base pairs (partner may be an insertion 'a^o' or None)
        annotated = set()            # consensus positions with any BGSU annotation (any category, near included)
        for o in obs:
            cat, cls, l1, l2, r1, r2 = o[:6]
            for r in (r1, r2):
                if isinstance(r, int): annotated.add(r)
            if cat != "basepair" or cls.startswith("n"): continue
            if isinstance(r1, int) and isinstance(r2, int): bp[tuple(sorted((r1, r2)))].add(cls)
            if isinstance(r1, int): partners[r1].add((r2 if r2 is not None else f"label{l2}", cls))
            if isinstance(r2, int): partners[r2].add((r1 if r1 is not None else f"label{l1}", cls))
        chains[sn] = {"inst": inst, "meta": meta[sn], "by_pos": by_pos, "by_label": by_label, "bp": bp, "partners": partners, "annotated": annotated,
                      "gap": {q for a, b in inst.get("annotation_gaps") or [] for q in range(a, b + 1)}}
    return fam, it, chains


def nt(ch, p):
    return ch["inst"]["nucleotides"][p - 1].upper()


def resid(ch, p):
    r = ch["by_pos"].get(p)
    if not r: return "-"
    return f"{r['residue_name']}{r['author_residue_number']}{r['author_insertion_code'] or ''}"


def evaluable(ch, i, j):
    e = ch["inst"]["eligibility"]
    return e[i - 1] == "E" and e[j - 1] == "E" and i not in ch["gap"] and j not in ch["gap"]


# ------------------------------------------------------------------ near-identical sequence clusters (match-state strings)
def seq_clusters(chains, clen, thr=0.95):
    keys = {}
    for sn, ch in chains.items():
        keys.setdefault(ch["meta"]["sequence_key"], (sn, ch["inst"]["nucleotides"]))
    reps, assign = [], {}
    for k, (sn, s) in sorted(keys.items(), key=lambda kv: -sum(c != "-" for c in kv[1][1])):
        for rk, rs in reps:
            shared = [(a, b) for a, b in zip(s, rs) if a not in "-." and b not in "-."]
            if shared and sum(a.upper() == b.upper() for a, b in shared) / len(shared) >= thr and len(shared) >= 0.8 * sum(c not in "-." for c in s):
                assign[k] = rk; break
        else:
            reps.append((k, s)); assign[k] = k
    return {sn: assign[ch["meta"]["sequence_key"]] for sn, ch in chains.items()}


def weighted(vals, groups):
    """mean of per-group means (each group weight 1)"""
    g = defaultdict(list)
    for v, k in zip(vals, groups): g[k].append(v)
    return sum(sum(x) / len(x) for x in g.values()) / len(g) if g else None


# ------------------------------------------------------------------ 1 + 4: three support measures, instance / sequence / species weighted
def part1(acc, fam, it, chains):
    sp = ss_pairs(fam["ss_cons"]); rows = []
    clus = seq_clusters(chains, len(fam["ss_cons"]))
    for i, j in sp:
        anyb, cww, can, seqk, clk, spc, cls = [], [], [], [], [], [], Counter()
        for sn, ch in chains.items():
            if not evaluable(ch, i, j): continue
            c = ch["bp"].get((i, j), set())
            anyb.append(1 if c else 0); cww.append(1 if "cWW" in c else 0); can.append(1 if "cWW" in c and nt(ch, i) + nt(ch, j) in CANON else 0)
            seqk.append(ch["meta"]["sequence_key"]); clk.append(clus[sn]); spc.append(ch["meta"]["organism_species"])
            cls.update(c or {"none"})
        n = len(anyb)
        r = {"acc": acc, "pos_1": i, "pos_2": j, "n_eligible": n, "n_unique_seq": len(set(seqk)), "n_seq_clusters": len(set(clk)), "n_species": len(set(spc)),
             "any_bp": sum(anyb), "cww": sum(cww), "canonical_cww": sum(can),
             "f_any_bp": round(sum(anyb) / n, 3) if n else None, "f_cww": round(sum(cww) / n, 3) if n else None, "f_canonical": round(sum(can) / n, 3) if n else None,
             "seqw_cww": round(weighted(cww, seqk), 3) if n else None, "clusw_cww": round(weighted(cww, clk), 3) if n else None,
             "speciesw_cww": round(weighted(cww, spc), 3) if n else None, "seqw_any_bp": round(weighted(anyb, seqk), 3) if n else None,
             "classes": "; ".join(f"{k}:{v}" for k, v in cls.most_common())}
        rows.append(r)
    return rows, clus


# ------------------------------------------------------------------ 2: alternative partners
ALT = [("RF00162", 50, 64, 66), ("RF00059", 64, 83, 82), ("RF00177", 372, 385, 386), ("RF00177", 648, 737, None), ("RF00059", 21, 27, None), ("RF00059", 20, 28, None)]


def part2(cache):
    rows, txt = [], []
    for acc, i, j, alt in ALT:
        fam, it, chains = cache[acc]
        if alt is None:   # the most frequent cWW partner of i or j other than each other
            c = Counter()
            for ch in chains.values():
                for p, cls in ch["partners"].get(i, ()):
                    if cls == "cWW" and p != j: c[("i", p)] += 1
                for p, cls in ch["partners"].get(j, ()):
                    if cls == "cWW" and p != i: c[("j", p)] += 1
            (side, alt), _ = c.most_common(1)[0]
        else:
            side = "i"
        a, b = (i, alt) if side == "i" else (alt, j)
        for sn, ch in sorted(chains.items()):
            e = ch["inst"]["eligibility"]
            def cl(x, y):
                return ",".join(sorted(ch["bp"].get(tuple(sorted((x, y))), set()))) or "-" if isinstance(y, int) else "-"
            net = {p: ";".join(f"{q}:{k}" for q, k in sorted(ch["partners"].get(p, ()), key=str)) for p in (i, j, alt)}
            # insertions between j and alt (or i and alt) in this chain: residues with alignment_state insert anchored inside the span
            lo, hi = sorted((j, alt) if side == "i" else (i, alt))
            ins = [r for r in ch["by_label"].values() if r["alignment_state"] in ("insert", "local_end") and r["insertion_anchor"] is not None and lo - 1 <= r["insertion_anchor"] <= hi]
            pp = lambda p: (ch["by_pos"].get(p) or {}).get("posterior_prob_char", "")
            rows.append({"acc": acc, "case": f"{i}-{j} vs {a}-{b}", "seq_name": sn, "pdb": ch["meta"]["pdb_id"].upper(), "chain": ch["meta"]["auth_asym_id"],
                         "species": ch["meta"]["organism_species"], "method": ch["meta"]["experimental_method"], "resolution": ch["meta"]["resolution"],
                         "elig_i_j_alt": e[i - 1] + e[j - 1] + e[alt - 1], "res_i": resid(ch, i), "res_j": resid(ch, j), "res_alt": resid(ch, alt),
                         "pp_i_j_alt": pp(i) + pp(j) + pp(alt), "class_ss_pair": cl(i, j), "class_alt_pair": cl(a, b),
                         "insertions_in_span": " ".join(f"{r['residue_name']}{r['author_residue_number']}@{r['insertion_anchor']}^{r['insertion_offset']}" for r in ins),
                         "partners_i": net[i], "partners_j": net[j], "partners_alt": net[alt]})
        # local alignment (unique sequences) from work/aligned.sto
        txt.append(f"=== {acc} SS_cons {i}-{j} vs observed {a}-{b} ===")
        txt += sto_window(acc, sorted({i, j, alt}))
    return rows, txt


def read_sto(acc):
    seqs, gr, gc = {}, {}, {}
    for line in open(ROOT / f"data/families/{acc}/work/aligned.sto"):
        if line.startswith("#=GC"):
            _, k, s = line.split(); gc[k] = gc.get(k, "") + s
        elif line.startswith("#=GR"):
            _, n, k, s = line.split(); gr[n] = gr.get(n, "") + s
        elif line.strip() and not line.startswith("#") and not line.startswith("//"):
            n, s = line.split(); seqs[n] = seqs.get(n, "") + s
    return seqs, gr, gc


def sto_window(acc, poss, flank=4):
    seqs, gr, gc = read_sto(acc)
    rf = gc["RF"]; cols = [k for k, c in enumerate(rf) if c not in ".~"]  # consensus column k -> alignment column (match columns only)
    out = []
    for p in poss:
        a0, a1 = cols[max(p - 1 - flank, 0)], cols[min(p - 1 + flank, len(cols) - 1)]
        out.append(f"-- consensus {max(p - flank, 1)}..{min(p + flank, len(cols))} (centre {p}); columns incl. inserts")
        out.append(f"{'SS_cons':>12} {gc['SS_cons'][a0:a1 + 1]}"); out.append(f"{'RF':>12} {rf[a0:a1 + 1]}")
        seen = Counter()
        for n, s in seqs.items():
            seen[(s[a0:a1 + 1], gr.get(n, '')[a0:a1 + 1])] += 1
        for (w, q), c in seen.most_common():
            out.append(f"{c:>3} seq(s)   {w}   PP {q}")
    return out


# ------------------------------------------------------------------ 3: annotation completeness
def run_len(ch, p, clen):
    if p in ch["annotated"]: return 0
    e = ch["inst"]["eligibility"]; n = 1
    for step in (-1, 1):
        q = p + step
        while 1 <= q <= clen and (e[q - 1] != "E" or q not in ch["annotated"]):
            if e[q - 1] == "E": n += 1
            q += step
    return n


def part3(acc, fam, it, chains):
    sp = ss_pairs(fam["ss_cons"]); clen = len(fam["ss_cons"]); rows = []; elig_counts = Counter()
    for i, j in sp:
        for sn, ch in chains.items():
            e = ch["inst"]["eligibility"]; elig_counts[e[i - 1] + e[j - 1] if "E" != e[i - 1] or "E" != e[j - 1] else "EE"] += 1
            if not evaluable(ch, i, j) or ch["bp"].get((i, j)): continue
            mod = [p for p in (i, j) if (ch["by_pos"].get(p) or {}).get("is_modified")]
            ri, rj = run_len(ch, i, clen), run_len(ch, j, clen)
            pi = sorted({str(q) for q, _ in ch["partners"].get(i, ())}); pj = sorted({str(q) for q, _ in ch["partners"].get(j, ())})
            cause = ("residue without any BGSU annotation" if ri or rj else
                     "both paired elsewhere" if pi and pj else "one paired elsewhere" if pi or pj else "both annotated, neither base-paired")
            rows.append({"acc": acc, "pos_1": i, "pos_2": j, "seq_name": sn, "res_i": resid(ch, i), "res_j": resid(ch, j), "modified": ",".join(map(str, mod)),
                         "unannotated_run_i": ri, "unannotated_run_j": rj, "partners_i": " ".join(pi), "partners_j": " ".join(pj), "cause": cause,
                         "method": ch["meta"]["experimental_method"], "resolution": ch["meta"]["resolution"]})
    runs = Counter()
    for ch in chains.values():
        e = ch["inst"]["eligibility"]; cur = 0
        for p in range(1, clen + 2):
            if p <= clen and e[p - 1] == "E" and p not in ch["annotated"]: cur += 1
            elif p <= clen and e[p - 1] != "E": continue
            else:
                if cur: runs[cur] += 1
                cur = 0
    # modified residues on SS_cons pairs: cWW support with / without a modification at either end
    modstat = Counter()
    for i, j in sp:
        for ch in chains.values():
            if not evaluable(ch, i, j): continue
            m = any((ch["by_pos"].get(p) or {}).get("is_modified") for p in (i, j))
            modstat[(m, "cWW" in ch["bp"].get((i, j), set()))] += 1
    return rows, runs, elig_counts, modstat


# ------------------------------------------------------------------ 5: outliers
OUTLIERS = [("RF00005", "3qsy"), ("RF00005", "6ahu"), ("RF00005", "1j2b"), ("RF00167", "9v50"), ("RF00167", "9nbc")]


def part5(cache):
    rows = []
    for acc, pdb in OUTLIERS:
        fam, it, chains = cache[acc]; sp = ss_pairs(fam["ss_cons"])
        for sn, ch in chains.items():
            if ch["meta"]["pdb_id"].lower() != pdb: continue
            for h in helices(sp):
                for i, j in h:
                    e = ch["inst"]["eligibility"]
                    c = ch["bp"].get((i, j), set())
                    rows.append({"acc": acc, "seq_name": sn, "helix": f"{h[0][0]}-{h[0][1]}", "pos_1": i, "pos_2": j, "elig": e[i - 1] + e[j - 1],
                                 "res_i": resid(ch, i), "res_j": resid(ch, j), "nts": (nt(ch, i) + nt(ch, j)) if e[i - 1] != "-" and e[j - 1] != "-" else "",
                                 "pp": (ch["by_pos"].get(i) or {}).get("posterior_prob_char", "") + (ch["by_pos"].get(j) or {}).get("posterior_prob_char", ""),
                                 "class": ",".join(sorted(c)) or "-", "partners_i": " ".join(f"{q}:{k}" for q, k in sorted(ch["partners"].get(i, ()), key=str)),
                                 "partners_j": " ".join(f"{q}:{k}" for q, k in sorted(ch["partners"].get(j, ()), key=str)),
                                 "annotated_i_j": f"{int(i in ch['annotated'])}{int(j in ch['annotated'])}",
                                 "method": ch["meta"]["experimental_method"], "resolution": ch["meta"]["resolution"], "title": ch["meta"].get("title", "")[:90]})
    # 16S hairpin 75-86 .. 78-83
    fam, it, chains = cache["RF00177"]; hp = [(75, 86), (76, 85), (77, 84), (78, 83)]; hrows = []
    for sn, ch in sorted(chains.items()):
        e = ch["inst"]["eligibility"]
        st = []
        for i, j in hp:
            if not evaluable(ch, i, j): st.append(e[i - 1] + e[j - 1] if not (i in ch["gap"] or j in ch["gap"]) else "gap"); continue
            c = ch["bp"].get((i, j), set()); st.append("cWW" if "cWW" in c else ",".join(sorted(c)) or "none")
        hrows.append({"seq_name": sn, "pdb": ch["meta"]["pdb_id"].upper(), "species": ch["meta"]["organism_species"], "method": ch["meta"]["experimental_method"],
                      "resolution": ch["meta"]["resolution"], "res_75_78": " ".join(resid(ch, p) for p in (75, 76, 77, 78)),
                      "res_83_86": " ".join(resid(ch, p) for p in (83, 84, 85, 86)), "status_75-86_76-85_77-84_78-83": " ".join(st),
                      "n_cww": st.count("cWW"), "annotated_75_86": "".join(str(int(p in ch["annotated"])) for p in range(75, 87)),
                      "elig_70_90": ch["inst"]["eligibility"][69:90], "title": ch["meta"].get("title", "")[:80]})
    return rows, hrows


# ------------------------------------------------------------------ 6: attributes of recurrent 3D-only pairs
def cooccurrence(chains, extra, ssp):
    """chains (both pairs evaluable) where the 3D-only pair and the SS_cons pair sharing a position are formed together (base triple or
    bifurcation) vs only one of them (competing alternatives)"""
    both = only_x = only_s = neither = 0
    for ch in chains.values():
        if not (evaluable(ch, *extra) and evaluable(ch, *ssp)): continue
        x, s = bool(ch["bp"].get(extra)), "cWW" in ch["bp"].get(ssp, set())
        both += x and s; only_x += x and not s; only_s += s and not x; neither += not x and not s
    return both, only_x, only_s, neither


def part6(acc, fam, it, p1rows, chains):
    sp = set(ss_pairs(fam["ss_cons"])); paired = {p for q in sp for p in q}
    low = {(r["pos_1"], r["pos_2"]) for r in p1rows if r["acc"] == acc and r["f_cww"] is not None and r["f_cww"] < 0.5}
    rows = []
    for p in it["pairs"]:
        if p["category"] != "basepair" or p["in_ss_cons"] or p["n_eligible"] < 3 or (p["frequency"] or 0) < 0.5: continue
        i, j = p["pos_1"], p["pos_2"]; cls = p["classes"]; top = max(cls, key=cls.get)
        a = {"crosses_ss_cons": any(crosses((i, j), q) for q in sp),
             "shares_position_with_ss_pair": any(q[0] in (i, j) or q[1] in (i, j) for q in sp),
             "shares_position_with_low_support_ss_pair": any(q[0] in (i, j) or q[1] in (i, j) for q in low),
             "stacks_on_ss_pair": (i - 1, j + 1) in sp or (i + 1, j - 1) in sp,
             "both_positions_unpaired_in_ss_cons": i not in paired and j not in paired,
             "main_class_cww": top == "cWW", "multiple_classes": len(cls) > 1}
        old = ("shift" if any(q in low and (q[0] in (i, j) or q[1] in (i, j)) for q in sp)
               else "helix_extension" if top == "cWW" and a["stacks_on_ss_pair"]
               else "pseudoknot (crosses SS_cons)" if a["crosses_ss_cons"] else "tertiary (nested)")
        # proposed: independent attributes, one primary label by a fixed precedence that is reported with the attributes
        shared = [q for q in sorted(sp) if q[0] in (i, j) or q[1] in (i, j)]
        co = [cooccurrence(chains, (i, j), q) for q in shared]
        a["coexists_with_sharing_ss_pair"] = any(b >= 0.5 * max(1, b + ox) for b, ox, os_, ne in co)      # formed together in >= half the chains with the 3D pair
        a["replaces_sharing_ss_pair"] = any(ox >= 0.5 * max(1, b + ox) for b, ox, os_, ne in co)
        prim = ("alternative partner (replaces an SS_cons pair)" if a["replaces_sharing_ss_pair"]
                else "base triple / bifurcation with an SS_cons pair" if a["coexists_with_sharing_ss_pair"]
                else "helix extension" if a["stacks_on_ss_pair"] and top == "cWW"
                else "crossing (pseudoknot-like)" if a["crosses_ss_cons"] else "nested, outside SS_cons helices")
        rows.append({"acc": acc, "pos_1": i, "pos_2": j, "main_class": top, "classes": json.dumps(cls), "frequency": p["frequency"], "n_eligible": p["n_eligible"],
                     **{k: int(v) for k, v in a.items()}, "n_true_of_first4": sum(a[k] for k in ("crosses_ss_cons", "shares_position_with_ss_pair", "stacks_on_ss_pair", "both_positions_unpaired_in_ss_cons")),
                     "sharing_ss_pairs_both_onlyX_onlySS_neither": "; ".join(f"{q[0]}-{q[1]}:{'/'.join(map(str, c))}" for q, c in zip(shared, co)),
                     "old_kind": old, "proposed_primary": prim})
    return rows


def write(name, rows):
    if not rows: return
    with open(OUT / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    cache = {acc: load(acc) for acc in FAMS}
    p1, p3, p6, red, summ = [], [], [], [], {}
    runs_all = {}
    for acc in FAMS:
        fam, it, chains = cache[acc]
        rows, clus = part1(acc, fam, it, chains); p1 += rows
        n = [r for r in rows if r["n_eligible"]]
        summ[acc] = {"n_pairs": len(rows), "n_evaluable_pairs": len(n),
                     **{f"{m}_ge_0.9": sum(r[m] >= 0.9 for r in n) for m in ("f_any_bp", "f_cww", "f_canonical", "seqw_cww", "clusw_cww", "speciesw_cww")},
                     **{f"{m}_lt_0.5": sum(r[m] < 0.5 for r in n) for m in ("f_any_bp", "f_cww", "f_canonical", "seqw_cww", "clusw_cww", "speciesw_cww")},
                     "pooled_any_bp": round(sum(r["any_bp"] for r in n) / sum(r["n_eligible"] for r in n), 3),
                     "pooled_cww": round(sum(r["cww"] for r in n) / sum(r["n_eligible"] for r in n), 3),
                     "pooled_canonical": round(sum(r["canonical_cww"] for r in n) / sum(r["n_eligible"] for r in n), 3),
                     "pairs_changing_band_seq_weighted": sorted(f"{r['pos_1']}-{r['pos_2']} {r['f_cww']}->{r['clusw_cww']}" for r in n
                                                                if (r["f_cww"] >= 0.9) != (r["clusw_cww"] >= 0.9) or (r["f_cww"] < 0.5) != (r["clusw_cww"] < 0.5))}
        meta = [ch["meta"] for ch in chains.values()]
        red.append({"acc": acc, "n_chains": len(chains), "n_entries": len({m["pdb_id"] for m in meta}), "n_unique_sequences": len({m["sequence_key"] for m in meta}),
                    "n_seq_clusters_95": len(set(clus.values())), "n_species": len({m["organism_species"] for m in meta}),
                    "largest_cluster_chains": max(Counter(clus.values()).values()),
                    "species": "; ".join(f"{k}:{v}" for k, v in Counter(m["organism_species"] for m in meta).most_common()),
                    "methods": "; ".join(f"{k}:{v}" for k, v in Counter(m["experimental_method"] for m in meta).most_common())})
        r3, runs, elig, modstat = part3(acc, fam, it, chains); p3 += r3; runs_all[acc] = runs
        summ[acc]["unpaired_causes"] = dict(Counter(r["cause"] for r in r3))
        summ[acc]["unpaired_with_short_unannotated_run_3_14"] = sum(1 for r in r3 if 3 <= max(r["unannotated_run_i"], r["unannotated_run_j"]) < 15)
        summ[acc]["unpaired_with_modified_end"] = sum(1 for r in r3 if r["modified"])
        summ[acc]["ss_pair_chain_eligibility"] = dict(elig)
        summ[acc]["modified_end_cww"] = {f"modified={k[0]} cWW={k[1]}": v for k, v in modstat.items()}
        p6 += part6(acc, fam, it, rows, chains)
    write("support_three_measures.csv", p1); write("unpaired_causes.csv", p3); write("redundancy.csv", red); write("extra_pair_attributes.csv", p6)
    write("unannotated_runs.csv", [{"acc": a, "run_length": k, "n_runs": v} for a, r in runs_all.items() for k, v in sorted(r.items())])
    focus = [("RF00162", 46, 69), ("RF00011", 324, 342), ("RF00011", 325, 341), ("RF00011", 326, 340), ("RF00011", 129, 162), ("RF00011", 7, 361),
             ("RF00011", 290, 305), ("RF00177", 139, 171), ("RF00177", 712, 723), ("RF00177", 680, 687)]
    write("focus_pair_classes.csv", [r for r in p1 if (r["acc"], r["pos_1"], r["pos_2"]) in focus])
    r2, txt = part2(cache); write("alt_partner_chains.csv", r2); (OUT / "alt_partner_alignment.txt").write_text("\n".join(txt) + "\n")
    r5, h5 = part5(cache); write("outlier_helices.csv", r5); write("hairpin_16s.csv", h5)
    summ["extra_pair_overlap"] = {acc: {"n": sum(1 for r in p6 if r["acc"] == acc), "with_more_than_one_of_4_attributes": sum(1 for r in p6 if r["acc"] == acc and r["n_true_of_first4"] > 1),
                                        "old_vs_proposed": dict(Counter(f"{r['old_kind']} -> {r['proposed_primary']}" for r in p6 if r["acc"] == acc))} for acc in FAMS}
    json.dump(summ, open(OUT / "followup_summary.json", "w"), indent=1)
    print(json.dumps(summ, indent=1))
