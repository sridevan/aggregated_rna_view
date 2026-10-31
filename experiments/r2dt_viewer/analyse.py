"""Summarise R2DT layout coverage per family (reads data/families/<ACC>/r2dt/index.json and the accepted api.json files).

    python3 experiments/r2dt_viewer/analyse.py RF00167 RF00162 ... > experiments/r2dt_viewer/coverage.json

Per family: chain counts by status and layout, templates used, rejection reasons, and layout consistency: for every template shared by
two or more accepted chains, the distance between the drawn positions of the same Rfam consensus position (match states only) in each
chain pair, after removing the best global translation (the canvases differ in size and offset). Also how many chains PDBe has no 2D for
that R2DT covers.
"""
from __future__ import annotations

import itertools
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def consensus_xy(fam: Path, seq_name: str, mapping_file: str) -> dict[int, tuple[float, float]]:
    api = json.loads((fam / "r2dt" / seq_name / "api.json").read_text())
    res = {r["polymer_seq_position"]: r for r in json.loads((fam / mapping_file).read_text())["residues"]}
    out = {}
    for k, lab in enumerate(api["label_seq_ids"][1:-1]):
        p = api["svg_paths"][k + 2].replace("M", "").split(",")   # previous -> current nucleotide; the current one is the end point
        r = res.get(lab)
        if r and r["alignment_state"] == "match" and r["rfam_consensus_pos"]:
            out[r["rfam_consensus_pos"]] = (float(p[2]), float(p[3]))
    return out


def spacing(xy: dict[int, tuple[float, float]]) -> float:
    ks = sorted(xy)
    d = [math.dist(xy[a], xy[b]) for a, b in zip(ks, ks[1:]) if b == a + 1]
    return statistics.median(d) if d else 0.0


def pair_residual(a: dict, b: dict) -> tuple[int, float, float] | None:
    common = sorted(set(a) & set(b))
    if len(common) < 10:
        return None
    tx = statistics.fmean(a[c][0] - b[c][0] for c in common)
    ty = statistics.fmean(a[c][1] - b[c][1] for c in common)
    d = sorted(math.dist(a[c], (b[c][0] + tx, b[c][1] + ty)) for c in common)
    return len(common), statistics.median(d), d[-1]


def family(acc: str) -> dict:
    fam = ROOT / "data" / "families" / acc
    idx = json.loads((fam / "r2dt" / "index.json").read_text())
    inst = {i["seq_name"]: i for i in json.loads((fam / "instances.json").read_text())}
    rows = idx["chains"]
    acc_rows = [r for r in rows if r["status"] == "accepted"]
    reasons = Counter()
    for r in rows:
        if r["status"] == "rejected":
            v = r["validation"]
            if r.get("r2dt_warning"):
                reasons["R2DT warned 'could not map to full sequence'; positions shifted"] += 1
            elif v["n_drawn"] > v["n_chain"]:
                reasons["more positions drawn than the chain has, no R2DT warning"] += 1
            elif v["n_drawn"] < v["n_chain"]:
                reasons["fewer positions drawn and shifted, no R2DT warning"] += 1
            else:
                reasons["base or author number differs: " + "; ".join(v["problems"][:2])] += 1
        elif r["status"] == "failed":
            reasons["failed: " + r.get("reason", "")[:120]] += 1
    by_template = defaultdict(list)
    for r in acc_rows:
        if r["layout"] == "templated":
            by_template[r["template"]].append(r["seq_name"])
    consistency = {}
    for t, names in sorted(by_template.items(), key=lambda kv: -len(kv[1])):
        if len(names) < 2:
            continue
        xy = {n: consensus_xy(fam, n, inst[n]["mapping_file"]) for n in names}
        pairs = {(a, b): x for a, b in itertools.combinations(names, 2) if (x := pair_residual(xy[a], xy[b]))}
        res = list(pairs.values())
        sp = statistics.median(spacing(v) for v in xy.values())
        # outlier chain: its median residual against the other chains exceeds one nucleotide spacing
        per_chain = {n: statistics.median([x[1] for k, x in pairs.items() if n in k] or [0]) for n in names}
        outliers = {n: round(m, 1) for n, m in per_chain.items() if m > sp}
        if res:
            med = [x[1] for x in res]
            rest = [x[1] for k, x in pairs.items() if not set(k) & set(outliers)]
            consistency[t] = {"chains": len(names), "pairs": len(res), "nt_spacing": round(sp, 2),
                              "median_of_pair_medians": round(statistics.median(med), 2),
                              "pairs_median_below_1px": sum(m < 1 for m in med),
                              "pairs_median_above_1_spacing": sum(m > sp for m in med),
                              "max_pair_median": round(max(med), 2),
                              "outlier_chains": outliers,
                              "without_outliers": {"pairs": len(rest), "median_of_pair_medians": round(statistics.median(rest), 2) if rest else None,
                                                   "max_pair_median": round(max(rest), 2) if rest else None,
                                                   "pairs_median_below_half_spacing": sum(m < sp / 2 for m in rest)}}
    no_pdbe = [n for n, i in inst.items() if i.get("pdbe_2d_available") is False]
    return {
        "family": acc, "chains": len(rows),
        "status": dict(Counter(r["status"] for r in rows)),
        "accepted_layout": dict(Counter(r["layout"] for r in acc_rows)),
        "templates": dict(Counter(r["template"] for r in acc_rows if r["layout"] == "templated").most_common()),
        "accepted_on_family_template": sum(r["template"] == acc for r in acc_rows),
        "rejection_reasons": dict(reasons),
        "modified_drawn_as_N_chains": sum(bool(r["validation"]["modified_drawn_as_N"]) for r in acc_rows),
        "accepted_with_matched_positions_not_drawn": {r["seq_name"]: r["validation"]["not_drawn_consensus_matched"]
                                                       for r in acc_rows if r["validation"]["not_drawn_consensus_matched"]},
        "no_pdbe_2d": len(no_pdbe),
        "no_pdbe_2d_covered_by_r2dt_template": sum(any(r["seq_name"] == n and r["status"] == "accepted" and r["layout"] == "templated"
                                                       for r in rows) for n in no_pdbe),
        "layout_consistency": consistency,
        "provenance": idx["provenance"],
    }


if __name__ == "__main__":
    print(json.dumps([family(a) for a in sys.argv[1:]], indent=1))
