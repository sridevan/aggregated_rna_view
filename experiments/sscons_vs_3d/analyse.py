"""Rfam SS_cons base pairs vs FR3D base pairs observed in the 3D structures, on Rfam consensus positions (investigation only).

Inputs (existing outputs, unchanged): data/families/<ACC>/family.json (ss_cons), interactions.json (instances with per-position eligibility
and nucleotides; pairs with frequencies) and interactions_by_instance/<seq_name>.json (every FR3D observation of a chain on Rfam positions).

For every SS_cons pair (i, j) and every chain where both positions are eligible (match residue, coordinates, complete base, FR3D available)
and neither lies in one of the chain's annotation gaps (>= 15 consecutive eligible positions without any BGSU annotation, i.e. incomplete source
annotation; those pairs are not evaluable and are counted separately):
  cww_canonical      FR3D cWW and the nucleotides form AU/UA/GC/CG/GU/UG
  cww_noncanonical   FR3D cWW with another nucleotide combination (e.g. C·A, U·U)
  other_class        another Leontis-Westhof family between i and j (e.g. tWH)
  near_only          only an FR3D near pair (n-prefixed class, from interactions.json instances[].near) between i and j
  shifted            no i-j pair, but i or j forms a cWW with the neighbour of the other partner (j+-1 / i+-1): register shift candidate
  unpaired           none of the above
Support of the SS_cons pair = (cww_canonical + cww_noncanonical) / eligible chains, i.e. the pair forms with Watson-Crick edges; the canonical
share is reported separately (support_canonical), so consistently non-canonical consensus pairs (e.g. 16S G·A, U·U) are not called unsupported.
Recurrent 3D pairs not in SS_cons: pairs list of interactions.json, base pairs with in_ss_cons false, frequency >= 0.5 over >= 3 eligible
chains; classified as crossing the SS_cons (pseudoknot), helix extension (stacks on an SS_cons pair), shift (adjacent to a poorly supported
SS_cons pair sharing one position) or other tertiary.
"""
import csv, json, sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).parent
FAMS = sys.argv[1:] or ["RF00167", "RF00162", "RF00059", "RF00005", "RF00011", "RF00177"]
CANON = {"AU", "UA", "GC", "CG", "GU", "UG"}
STATUSES = ["cww_canonical", "cww_noncanonical", "other_class", "near_only", "shifted", "unpaired"]


def ss_pairs(ss: str) -> list[tuple[int, int]]:
    st = defaultdict(list); out = []
    for i, ch in enumerate(ss, 1):
        for o, c in ("()", "<>", "[]", "{}"):
            if ch == o: st[o].append(i)
            elif ch == c and st[o]: out.append((st[o].pop(), i))
        if ch.isalpha():   # WUSS pseudoknot letters (none in these six families; handled for completeness)
            if ch.isupper(): st["pk" + ch].append(i)
            elif st["pk" + ch.upper()]: out.append((st["pk" + ch.upper()].pop(), i))
    return sorted(out)


def crosses(p, q):
    (a, b), (c, d) = p, q
    return a < c < b < d or c < a < d < b


def helices(pairs: list[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    s = set(pairs); seen = set(); out = []
    for p in sorted(pairs):
        if p in seen: continue
        h = [p]; seen.add(p); i, j = p
        while (i + 1, j - 1) in s: i, j = i + 1, j - 1; h.append((i, j)); seen.add((i, j))
        out.append(h)
    return out


def analyse(acc: str) -> dict:
    fam = json.load(open(ROOT / f"data/families/{acc}/family.json")); it = json.load(open(ROOT / f"data/families/{acc}/interactions.json"))
    ss = fam["ss_cons"]; sp = ss_pairs(ss); sset = set(sp); clen = len(ss)
    insts = [i for i in it["instances"] if i["annotation"] == "available"]
    # per chain: observed base pairs on consensus positions (class by unordered position pair)
    obs = {}
    for inst in insts:
        f = ROOT / f"data/families/{acc}/interactions_by_instance/{inst['seq_name']}.json"
        d = defaultdict(set)
        for o in json.load(open(f))["observations"]:
            if o[0] != "basepair" or not isinstance(o[4], int) or not isinstance(o[5], int): continue
            a, b = sorted((o[4], o[5])); d[(a, b)].add(o[1])
        for n in inst.get("near", []):
            if n[0] == "basepair" and isinstance(n[2], int) and isinstance(n[3], int):
                a, b = sorted((n[2], n[3])); d[(a, b)].add(n[1] if n[1].startswith("n") else "n" + n[1])
        obs[inst["seq_name"]] = d
    def status(inst, i, j):
        d = obs[inst["seq_name"]]; cls = d.get((i, j), set())
        nt = (inst["nucleotides"][i - 1] + inst["nucleotides"][j - 1]).upper()
        if "cWW" in cls: return "cww_canonical" if nt in CANON else "cww_noncanonical"
        if any(not c.startswith("n") for c in cls): return "other_class"
        if cls: return "near_only"
        for (a, b) in ((i, j - 1), (i, j + 1), (i - 1, j), (i + 1, j)):
            if "cWW" in d.get(tuple(sorted((a, b))), set()): return "shifted"
        return "unpaired"
    in_gap = lambda inst, k: any(a <= k <= b for a, b in inst.get("annotation_gaps") or [])
    rows = []; per_chain = defaultdict(Counter)
    for (i, j) in sp:
        c = Counter(); nts = Counter(); gap = 0
        for inst in insts:
            e = inst["eligibility"]
            if e[i - 1] != "E" or e[j - 1] != "E": continue
            if in_gap(inst, i) or in_gap(inst, j): gap += 1; continue
            s = status(inst, i, j); c[s] += 1; per_chain[inst["seq_name"]][s] += 1
            nts[(inst["nucleotides"][i - 1] + inst["nucleotides"][j - 1]).upper()] += 1
        n = sum(c.values())
        rows.append({"acc": acc, "pos_1": i, "pos_2": j, "bracket": ss[i - 1] + ss[j - 1], "region_1": (it["position_regions"][i - 1] if it.get("position_regions") else None),
                     "n_eligible": n, "n_in_annotation_gap": gap, **{k: c[k] for k in STATUSES}, "support": round((c["cww_canonical"] + c["cww_noncanonical"]) / n, 3) if n else None,
                     "support_canonical": round(c["cww_canonical"] / n, 3) if n else None,
                     "nt_pairs": "; ".join(f"{k}:{v}" for k, v in nts.most_common())})
    # helices
    hel = []
    for h in helices(sp):
        hr = [r for r in rows if (r["pos_1"], r["pos_2"]) in set(h)]; elig = sum(r["n_eligible"] for r in hr); ok = sum(r["cww_canonical"] + r["cww_noncanonical"] for r in hr)
        hel.append({"acc": acc, "helix": f"{h[0][0]}-{h[0][1]}", "n_pairs": len(h), "eligible_pair_observations": elig, "support": round(ok / elig, 3) if elig else None,
                    "weakest_pair": min(((r["support"], f"{r['pos_1']}-{r['pos_2']}") for r in hr if r["support"] is not None), default=(None, None))[1]})
    # recurrent 3D pairs not in SS_cons
    extra = []
    low = {(r["pos_1"], r["pos_2"]) for r in rows if r["support"] is not None and r["support"] < 0.5}
    for p in it["pairs"]:
        if p["category"] != "basepair" or p["in_ss_cons"] or p["n_eligible"] < 3 or (p["frequency"] or 0) < 0.5: continue
        i, j = p["pos_1"], p["pos_2"]; cls = p["classes"]; top = max(cls, key=cls.get)
        kind = ("shift" if any(q in low and (q[0] in (i, j) or q[1] in (i, j)) for q in sset)
                else "helix_extension" if top == "cWW" and ((i - 1, j + 1) in sset or (i + 1, j - 1) in sset)
                else "pseudoknot (crosses SS_cons)" if any(crosses((i, j), q) for q in sset) else "tertiary (nested)")
        extra.append({"acc": acc, "pos_1": i, "pos_2": j, "regions": f"{p['region_1']} / {p['region_2']}", "frequency": p["frequency"], "n_observed": p["n_observed"],
                      "n_eligible": p["n_eligible"], "main_class": top, "classes": json.dumps(cls), "nt_pairs": json.dumps(p["nt_pairs"]), "kind": kind,
                      "span": j - i})
    chains = []
    for sn, c in per_chain.items():
        n = sum(c.values())
        chains.append({"acc": acc, "seq_name": sn, "n_eligible_ss_pairs": n, **{k: c[k] for k in STATUSES}, "agreement": round((c["cww_canonical"] + c["cww_noncanonical"]) / n, 3) if n else None})
    tot = Counter(); [tot.update({k: r[k] for k in STATUSES}) for r in rows]
    n_obs = sum(tot.values())
    summary = {"acc": acc, "clen": clen, "n_ss_pairs": len(sp), "n_helices": len(hel), "n_chains": len(insts), "pseudoknot_letters_in_ss_cons": any(ch.isalpha() for ch in ss),
               "pair_observations": n_obs, "status_fractions": {k: round(tot[k] / n_obs, 3) for k in STATUSES} if n_obs else {},
               "ss_pairs_support_ge_0.9": sum(1 for r in rows if r["support"] is not None and r["support"] >= 0.9),
               "ss_pairs_support_lt_0.5": sorted(f"{r['pos_1']}-{r['pos_2']} ({r['support']}, n={r['n_eligible']})" for r in rows if r["support"] is not None and r["support"] < 0.5),
               "ss_pairs_never_eligible": sum(1 for r in rows if not r["n_eligible"]),
               "recurrent_3d_pairs_not_in_ss_cons": dict(Counter(e["kind"] for e in extra)),
               "chain_agreement_median": sorted(c["agreement"] for c in chains if c["agreement"] is not None)[len(chains) // 2] if chains else None,
               "chains_agreement_lt_0.8": sorted(f"{c['seq_name']} ({c['agreement']})" for c in chains if c["agreement"] is not None and c["agreement"] < 0.8)}
    return {"summary": summary, "pairs": rows, "helices": hel, "extra": extra, "chains": chains}


def write(name, rows):
    if not rows: return
    with open(OUT / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    allr = {k: [] for k in ("pairs", "helices", "extra", "chains")}; summ = []
    for acc in FAMS:
        r = analyse(acc); summ.append(r["summary"])
        for k in allr: allr[k] += r[k]
    write("ss_cons_pair_support.csv", allr["pairs"]); write("helix_support.csv", allr["helices"]); write("recurrent_3d_pairs_not_in_ss_cons.csv", allr["extra"]); write("chain_agreement.csv", allr["chains"])
    json.dump(summ, open(OUT / "summary.json", "w"), indent=1)
    for s in summ: print(json.dumps(s))
