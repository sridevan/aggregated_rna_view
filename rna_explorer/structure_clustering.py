"""Structural groups for the Structures tab, from precomputed pairwise matrices (no structural comparison is recomputed here).

Inputs (read-only): a research run directory of the RNA structural clustering benchmark
  <run>/manifest.json                 instance order, shared consensus positions, coverage threshold, exclusions
  <run>/matrices/disc_base.csv        geometric discrepancy D (FR3D base centres + base orientations, fixed correspondence)
  <run>/matrices/ufix_tm_user.csv     US-align TM-score with fixed correspondence (-I), normalised by the n shared positions
and the run's clustering_summary.json for the stability evidence that the research computed.

Method (ported from antibody_antigen_interface/app-react/src/groups/clustering.js, so the page behaves like the reference demo):
  * average linkage (UPGMA) on the metric's dissimilarity (D, or 1 - TM), deterministic: ties go to the first pair met in instance
    order; the larger child is drawn first, then the one holding the lowest index;
  * fixed cut: undo the k - 1 highest merges (k from structure_clustering_config.json); the stored height is midway between the
    lowest merge undone and the highest merge kept;
  * groups numbered by size (largest first), named "Structural group N"; no state labels;
  * representative = medoid (smallest mean dissimilarity to the other members), ties to the earliest instance in matrix order, taken among
    the members whose chain has no internal coordinate gap (no unmodelled nucleotide between its first and last modelled one), so the
    displayed backbone is continuous; the plain medoid is used when every member has a gap, and is always recorded as `medoid`.
"""
from __future__ import annotations
import csv
import datetime as _dt
import hashlib
import json
import logging
from pathlib import Path
from typing import Any, Optional
import numpy as np
from . import __version__
from .annotation_mapping import InstanceIndex
from .config import PROJECT_ROOT
from .export import validate_json
from .seriation import seriation_record

log = logging.getLogger("rna_explorer.structure_clustering")
SCHEMA_VERSION = "1.0"
CONFIG_PATH = Path(__file__).with_name("structure_clustering_config.json")
DETAIL_MAX = 50          # instance-level heatmap/dendrogram limit used by the frontend (visualisation only)
SMALL_GROUP = 5          # groups smaller than this are flagged "small" (as in the reference demo)
METRICS = {
    "discrepancy": {"file": "disc_base.csv", "label": "Geometric discrepancy", "short": "D", "kind": "dissimilarity", "stored_quantity": "D",
                    "distance": "D", "higher_is_more_similar": False,
                    "unit": "dimensionless (base-centre deviation in A and base rotation in rad, weight 1, divided by n)",
                    "definition": "D = sqrt(L + O) / n after Kabsch superposition of the n base centres: L = sum of squared base-centre deviations, "
                                  "O = sum of squared base-rotation angles (FR3D standard base frames); fixed correspondence on the n shared Rfam consensus positions"},
    "tm": {"file": "ufix_tm_user.csv", "label": "TM-score", "short": "TM", "kind": "similarity", "stored_quantity": "TM-score",
           "distance": "1 - TM-score", "higher_is_more_similar": True, "unit": "dimensionless (0-1)",
           "definition": "US-align TM-score with fixed residue correspondence (-mol RNA -I <pair.fasta> -u n), normalised by the n shared positions "
                         "(symmetric); clustering and the dendrogram use 1 - TM-score"},
}


# ----------------------------------------------------------------------------- inputs
def read_matrix(path: Path, ids: list[str]) -> np.ndarray:
    rows = list(csv.reader(open(path)))
    head = rows[0][1:]
    if head != ids or [r[0] for r in rows[1:]] != ids:
        raise ValueError(f"{path.name}: instance order differs from the manifest")
    return np.array([[float(x) for x in r[1:]] for r in rows[1:]], dtype=float)


def check_matrix(M: np.ndarray, diag: float, name: str) -> None:
    if not np.all(np.isfinite(M)):
        raise ValueError(f"{name}: non-finite values")
    if np.max(np.abs(M - M.T)) > 1e-9:
        raise ValueError(f"{name}: not symmetric")
    if np.max(np.abs(np.diag(M) - diag)) > 1e-9:
        raise ValueError(f"{name}: diagonal is not {diag}")


def _sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ----------------------------------------------------------------------------- tree (port of clustering.js)
def average_linkage(D: np.ndarray) -> tuple[list[list[float]], list[int], dict[int, list[int]]]:
    """UPGMA. Returns the scipy-style linkage [[left, right, height, size], ...] (node ids n..2n-2), the leaf drawing order and the
    leaves of every node in drawing order."""
    n = len(D)
    N = 2 * n - 1
    dist = np.zeros((N, N))
    dist[:n, :n] = D
    leaves: dict[int, list[int]] = {i: [i] for i in range(n)}
    active = list(range(n))
    Z: list[list[float]] = []
    nxt = n
    while len(active) > 1:
        ba = bb = -1; bd = np.inf
        for a in range(len(active)):
            ia = active[a]
            for b in range(a + 1, len(active)):
                d = dist[ia, active[b]]
                if d < bd:
                    ba, bb, bd = a, b, d
        A, B = active[ba], active[bb]
        la, lb = leaves[A], leaves[B]
        a_first = len(la) > len(lb) if len(la) != len(lb) else min(la) < min(lb)
        L, R = (A, B) if a_first else (B, A)
        node = nxt; nxt += 1
        leaves[node] = leaves[L] + leaves[R]
        for C in active:
            if C in (A, B):
                continue
            v = (dist[A, C] * len(la) + dist[B, C] * len(lb)) / (len(la) + len(lb))
            dist[node, C] = dist[C, node] = v
        Z.append([L, R, float(bd), len(leaves[node])])
        active = [x for x in active if x not in (A, B)] + [node]
    return Z, leaves[active[0]], leaves


def cut_tree(Z: list[list[float]], n: int, k: int) -> tuple[list[int], float, float, float]:
    """Undo the k - 1 highest merges. Returns the group root nodes in drawing order, the cut height (midway), the lowest merge undone
    and the highest merge kept."""
    height = {n + i: z[2] for i, z in enumerate(Z)}
    child = {n + i: (int(z[0]), int(z[1])) for i, z in enumerate(Z)}
    parts = [2 * n - 2]
    lowest_split = None
    while len(parts) < k:
        at = -1
        for i, p in enumerate(parts):
            if p in child and (at < 0 or height[p] > height[parts[at]]):
                at = i
        if at < 0:
            break
        p = parts[at]
        lowest_split = height[p] if lowest_split is None else min(lowest_split, height[p])
        parts = parts[:at] + list(child[p]) + parts[at + 1:]
    highest_kept = max([0.0] + [height.get(p, 0.0) for p in parts])
    cut = height[2 * n - 2] if lowest_split is None else (lowest_split + highest_kept) / 2
    return parts, cut, lowest_split if lowest_split is not None else height[2 * n - 2], highest_kept


def medoid(members: list[int], D: np.ndarray, candidates: Optional[list[int]] = None) -> tuple[int, float]:
    """Member with the smallest mean dissimilarity to the other members; ties to the lowest matrix index. With `candidates`, only those
    members can be chosen (their mean is still taken over all other members)."""
    if len(members) == 1:
        return members[0], 0.0
    best = None
    for i in sorted(candidates if candidates else members):
        mean = float(np.mean([D[i, j] for j in members if j != i]))
        if best is None or mean < best[1] - 1e-15:
            best = (i, mean)
    return best


def cophenetic_correlation(Z: list[list[float]], leaves: dict[int, list[int]], D: np.ndarray) -> float:
    n = len(D); C = np.zeros_like(D)
    for i, z in enumerate(Z):
        for a in leaves[int(z[0])]:
            for b in leaves[int(z[1])]:
                C[a, b] = C[b, a] = z[2]
    iu = np.triu_indices(n, 1)
    return float(np.corrcoef(D[iu], C[iu])[0, 1])


def silhouette(D: np.ndarray, labels: np.ndarray) -> Optional[float]:
    """Mean silhouette width with a precomputed dissimilarity (singletons score 0)."""
    ks = np.unique(labels)
    if len(ks) < 2:
        return None
    s = []
    for i in range(len(D)):
        own = labels == labels[i]
        if own.sum() == 1:
            s.append(0.0); continue
        a = D[i, own & (np.arange(len(D)) != i)].mean()
        b = min(D[i, labels == k].mean() for k in ks if k != labels[i])
        s.append((b - a) / max(a, b) if max(a, b) > 0 else 0.0)
    return float(np.mean(s))


def partition(Z: list[list[float]], n: int, k: int, leaves: dict[int, list[int]]) -> np.ndarray:
    parts, *_ = cut_tree(Z, n, k)
    lab = np.zeros(n, dtype=int)
    for g, p in enumerate(parts):
        lab[leaves[p]] = g
    return lab


# ----------------------------------------------------------------------------- one metric
def internal_gaps(mapping: dict[str, Any]) -> list[Any]:
    """Author numbers of the chain's unmodelled nucleotides that lie between its first and last modelled nucleotide (breaks in the
    displayed backbone); terminal unmodelled nucleotides are not breaks."""
    res = sorted(mapping["residues"], key=lambda r: r["polymer_seq_position"])
    obs = [i for i, r in enumerate(res) if r["has_coordinates"]]
    if not obs:
        return []
    return [f"{r['author_residue_number']}{r['author_insertion_code']}" for r in res[obs[0]:obs[-1] + 1] if not r["has_coordinates"]]


def cluster_metric(ids: list[str], entries: dict[str, str], D: np.ndarray, k: int, gaps: Optional[dict[str, list[Any]]] = None) -> dict[str, Any]:
    n = len(ids)
    Z, order, leaves = average_linkage(D)
    parts, cut, lowest_split, highest_kept = cut_tree(Z, n, k)
    pos = {leaf: p for p, leaf in enumerate(order)}
    ranked = sorted(enumerate(parts), key=lambda t: (-len(leaves[t[1]]), t[0]))
    groups = []
    member_all = set(range(n))
    for gid, (_, node) in enumerate(ranked, start=1):
        mem = leaves[node]
        plain, _ = medoid(mem, D)
        complete = [i for i in mem if not (gaps or {}).get(ids[i])] if gaps is not None else None
        med, med_mean = medoid(mem, D, complete or None)
        within = [D[a, b] for i, a in enumerate(mem) for b in mem[i + 1:]]
        outside = sorted(member_all - set(mem))
        across = [D[a, b] for a in mem for b in outside]
        wm = float(np.mean(within)) if within else None
        am = float(np.mean(across)) if across else None
        groups.append({"group_id": gid, "label": f"Structural group {gid}", "node": node, "members": [ids[i] for i in mem],
                       "n_instances": len(mem), "n_entries": len({entries[ids[i]] for i in mem}),
                       "representative": ids[med], "medoid": ids[plain],
                       "representative_rule": "medoid" if gaps is None or (med == plain and not gaps.get(ids[med])) else "medoid among members without internal coordinate gaps" if complete else "medoid (every member has internal coordinate gaps)",
                       "representative_internal_gaps": None if gaps is None else list(gaps.get(ids[med], [])),
                       "representative_mean_dissimilarity": round(med_mean, 6),
                       "leaf_range": [pos[mem[0]], pos[mem[-1]]],
                       "within": {"mean": None if wm is None else round(wm, 6), "max": round(float(max(within)), 6) if within else None, "n_pairs": len(within)},
                       "across": {"mean": None if am is None else round(am, 6), "n_pairs": len(across)},
                       "across_to_within_ratio": round(am / wm, 2) if wm and am else None,
                       "small": len(mem) < SMALL_GROUP, "instance_level_view": 2 <= len(mem) <= DETAIL_MAX})
    heights = sorted(z[2] for z in Z)
    gaps = [heights[i + 1] - heights[i] for i in range(len(heights) - 1)]
    top = heights[-10:]
    sil = {str(kk): round(silhouette(D, partition(Z, n, kk, leaves)), 4) for kk in range(2, 11)}
    return {"linkage": [[int(a), int(b), round(h, 6), int(s)] for a, b, h, s in Z], "leaf_order": order,
            "cophenetic_correlation": round(cophenetic_correlation(Z, leaves, D), 4), "max_height": round(heights[-1], 6),
            "cut": {"k": k, "height": round(cut, 6), "lowest_merge_undone": round(lowest_split, 6), "highest_merge_kept": round(highest_kept, 6),
                    "gap": round(lowest_split - highest_kept, 6), "rule": "undo the k - 1 highest merges; height midway between the lowest merge undone and the highest kept",
                    "is_largest_gap_of_tree": abs((lowest_split - highest_kept) - max(gaps)) < 1e-12,
                    "top_merge_heights": [round(h, 6) for h in top], "silhouette_by_k": sil},
            "groups": groups, "assignments": {ids[i]: g["group_id"] for g in groups for i in leaves[g["node"]]}}


# ----------------------------------------------------------------------------- build
def research_evidence(summary_path: Optional[Path], metric: str) -> dict[str, Any]:
    """Stability evidence computed by the research for this metric. Two summary layouts exist: per-threshold (RF00167:
    thresholds/<T>/metrics, stability_T<T>, threshold_sensitivity_ARI) and per-run (RF00162, RF00059: runs/<primary>/metrics,
    stability, run_sensitivity_ARI). Both store the duplicate floor, cophenetic correlation and per-k partitions per metric."""
    if not summary_path or not summary_path.exists():
        return {"status": "unavailable"}
    s = json.load(open(summary_path))
    key = {"discrepancy": "disc_base", "tm": "one_minus_tm_fixed"}[metric]
    def _filtered(d: dict[str, Any]) -> Optional[dict[str, Any]]:
        out = {k.split("|", 1)[1]: v for k, v in d.items() if k.split("|", 1)[0] == key}
        return out or None
    if "thresholds" in s:
        t = str(s.get("primary_T", 64))
        m = s["thresholds"][t]["metrics"].get(key, {})
        ev: dict[str, Any] = {"source": str(summary_path), "layout": "per_threshold", "primary": f"T{t}",
                              "noise_floor_p95": m.get("noise_floor", {}).get("p95"), "noise_floor_n_pairs": m.get("noise_floor", {}).get("n"),
                              "cophenetic_correlation": m.get("cophenetic_correlation"),
                              "silhouette_by_k": {k: v.get("silhouette") for k, v in m.get("partitions", {}).items()},
                              "threshold_sensitivity_ARI": _filtered(s.get("threshold_sensitivity_ARI", {}))}
        st = s.get(f"stability_T{t}", {})
        sep = s.get(f"apo_vs_bound_separability_T{t}", {})
    else:
        run = s.get("primary_run")
        m = s.get("runs", {}).get(run, {}).get("metrics", {}).get(key, {})
        ev = {"source": str(summary_path), "layout": "per_run", "primary": run,
              "noise_floor_p95": m.get("noise_floor", {}).get("p95"), "noise_floor_n_pairs": m.get("noise_floor", {}).get("n"),
              "cophenetic_correlation": m.get("cophenetic_correlation"),
              "silhouette_by_k": {k: v.get("silhouette") for k, v in m.get("partitions", {}).items()} or None,
              "run_sensitivity_ARI": _filtered(s.get("run_sensitivity_ARI", {}))}
        st = s.get("stability", {})
        sep = {}
    if key in st.get("position_bootstrap", {}):
        ev["position_bootstrap"] = st["position_bootstrap"][key]
    jk = {el: {kk: v for kk, v in d.items() if kk.startswith(key)} for el, d in st.get("element_jackknife", {}).items()}
    if any(jk.values()):
        ev["element_jackknife"] = jk
    if metric == "discrepancy" and sep:
        ev["apo_vs_bound_separability"] = sep.get("all", {}).get("disc_base")
    return ev


def verify_against_mappings(idx: InstanceIndex, ids: list[str], shared: list[int]) -> dict[str, Any]:
    """Every clustered instance must be a mapped instance of this checkout with a match-state residue with coordinates at every shared position."""
    bad = []
    for sn in ids:
        if sn not in idx.by_seq:
            raise ValueError(f"{sn} is not an instance of this family")
        for p in shared:
            r = idx.match_residue(sn, p)
            if r is None or not r.get("has_coordinates"):
                bad.append((sn, p))
    if bad:
        raise ValueError(f"shared positions without an observed match residue in the main-app mapping: {bad[:5]}")
    return {"n_instances_checked": len(ids), "n_positions_checked": len(shared), "result": "every clustered instance has a match-state residue with coordinates at every shared position"}


def build_structure_clustering(acc: str, out_dir: Path, research_root: Path, write: bool = True) -> dict[str, Any]:
    cfg = json.load(open(CONFIG_PATH))[acc]
    run = research_root / cfg["research_run"]
    summary = research_root / cfg["research_summary"] if cfg.get("research_summary") else None
    fam_dir = out_dir / "families" / acc
    idx = InstanceIndex(fam_dir)
    man = json.load(open(run / "manifest.json"))
    ids = [i if isinstance(i, str) else i["seq_name"] for i in man["instances"]]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate instance ids in the manifest")
    shared = [int(p) for p in man["shared_positions"]]
    check = verify_against_mappings(idx, ids, shared)
    entries = {sn: idx.by_seq[sn]["pdb_id"] for sn in ids}
    gaps = {sn: internal_gaps(json.load(open(fam_dir / "mappings" / f"{sn}.json"))) for sn in ids}
    excluded = [{"seq_name": e["seq_name"], "pdb_id": idx.by_seq[e["seq_name"]]["pdb_id"], "n_observed": e.get("n_observed"), "reason": e.get("reason"),
                 "reason_text": f"{e.get('n_observed')} of the {man['threshold_min_observed']} required core positions observed"} for e in man.get("excluded", [])]
    if set(ids) | {e["seq_name"] for e in excluded} != set(idx.by_seq):
        raise ValueError("clustered + excluded instances do not cover every mapped instance")
    metrics: dict[str, Any] = {}
    sources: dict[str, Any] = {}
    for key, meta in METRICS.items():
        path = run / "matrices" / meta["file"]
        M = read_matrix(path, ids)
        src = {"file": str(path), "sha256": _sha256(path)}
        check_matrix(M, 1.0 if key == "tm" else 0.0, meta["file"])
        dist = 1.0 - M if key == "tm" else M
        np.fill_diagonal(dist, 0.0)
        mc = cfg["metrics"][key]
        res = cluster_metric(ids, entries, dist, int(mc["k"]), gaps)
        iu = np.triu_indices(len(ids), 1)
        ev = research_evidence(summary, key)
        metrics[key] = {k: v for k, v in meta.items() if k != "file"} | {
            "values_upper": [round(float(x), 6) for x in M[iu]], "values_layout": "row-major upper triangle (i < j) in dataset.instances order; diagonal "
            + ("1" if key == "tm" else "0"),
            "value_range": {"min": round(float(M[iu].min()), 6), "max": round(float(M[iu].max()), 6)},
            "distance_range": {"min": round(float(dist[iu].min()), 6), "max": round(float(dist[iu].max()), 6)},
            "clustering": {"method": "agglomerative hierarchical clustering", "linkage_method": "average (UPGMA)", "distance": meta["distance"],
                           **res, "validation": mc["validation"], "rationale": mc["rationale"],
                           # a family can omit the floor (configured, with the research's reason); its research value stays in research_evidence
                           "floor": {"value": None if cfg.get("omit_floor") else ev.get("noise_floor_p95"),
                                     "definition": "95th percentile of the duplicate pairs (NCS copies, isomorphous re-determinations); smaller differences are not interpretable"
                                     + (f". Not shown: {cfg['floor_note']}" if cfg.get("omit_floor") and cfg.get("floor_note") else ""),
                                     "n_pairs": ev.get("noise_floor_n_pairs"), "source": "research clustering_summary.json"},
                           "research_evidence": ev}}
        if key == "discrepancy":   # Structural Dissimilarity mode: the discrepancy matrix in TPPL order (matrix only, independent of the clustering)
            metrics[key]["seriation"] = seriation_record(ids, dist)
        sources[key] = src
    mp = idx.mapping(ids[0])["provenance"]
    out = {"schema_version": SCHEMA_VERSION, "rfam_acc": acc, "visualisation": {"detail_max": DETAIL_MAX, "small_group": SMALL_GROUP},
           "dataset": {"n_mapped_instances": len(idx.by_seq), "n_mapped_entries": len({r["pdb_id"] for r in idx.instances}),
                       "n_clustered_instances": len(ids), "n_clustered_entries": len(set(entries.values())), "n_excluded_instances": len(excluded),
                       "inclusion_rule": f"observed (complete base, C4', C3') at >= {man['threshold_min_observed']} consensus positions of the aptamer core; "
                                         f"comparisons use the {len(shared)} positions observed in every included instance",
                       "coverage_threshold": man["threshold_min_observed"], "fixed_correspondence": man.get("fixed_correspondence"),
                       "shared_positions": shared, "n_shared_positions": len(shared), "instances": ids, "excluded": excluded,
                       "mapping_check": check},
           "default_metric": "discrepancy", "metrics": metrics,
           "provenance": {"research_run": str(run), "research_manifest_generated": man.get("generated"), "research_software": man.get("software"),
                          "research_reports": [str(research_root / "experiments/rna_structural_clustering" / r) for r in cfg.get("research_reports", [])],
                          "matrix_sources": sources, "config": str(CONFIG_PATH.name), "algorithm_reference": "antibody_antigen_interface/app-react/src/groups/clustering.js (ported)",
                          "rfam_release": mp.get("rfam_release"), "cm_checksum": mp.get("cm_checksum"), "mapping_pipeline_version": mp.get("pipeline_version"),
                          "mapping_generated": mp.get("generated"), "pipeline_version": __version__,
                          "generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                          "note": "pairwise values copied from the research run; no structural comparison was recomputed"}}
    if write:
        validate_json(out, "structure_clustering")
        with open(fam_dir / "structure_clustering.json", "w") as f:
            json.dump(out, f, indent=1)
    return out
