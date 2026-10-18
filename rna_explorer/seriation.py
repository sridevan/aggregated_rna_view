"""Tree-penalized path length (TPPL) seriation of a pairwise dissimilarity matrix, for the Structures tab's Structural Dissimilarity mode.

The algorithm is the reference implementation, vendored unmodified as `ordering_similarity.py` from
https://github.com/sridevan/ordering_analysis (commit 37f63c9, itself adapted from BGSU RNA-3D-correspondence; Aliyev & Zirbel 2023,
EJOR 305(2):617-629) and called as in its run.py: treePenalizedPathLength(dist, 10, 39873). This wrapper only fixes the input order,
refuses incomplete input and keeps the global random state untouched:
- rows/columns are put in sorted identifier order before the call, so the result never depends on how the caller listed the instances;
- the reference imputes missing (NaN or negative) values with the largest observed value; a missing comparison is not a maximal
  difference, so incomplete matrices are refused instead (the clustering pipeline already rejects non-finite matrices);
- scipy is required: without it the reference silently drops the tree penalty and returns a different ordering.
"""
from __future__ import annotations
import hashlib
import random
from pathlib import Path
from typing import Any
import numpy as np

REFERENCE = {"repository": "https://github.com/sridevan/ordering_analysis", "commit": "37f63c930e3b75e5cf450ad7be8ad6f76d227990",
             "file": "ordering_similarity.py", "function": "treePenalizedPathLength (tree penalty from average linkage, multiple greedy insertion, orientPath)",
             "citation": "Aliyev, D. A. & Zirbel, C. L. (2023) Seriation using tree-penalized path length. EJOR 305(2), 617-629. doi:10.1016/j.ejor.2022.06.026"}
REPETITIONS = 10
SEED = 39873
VENDORED = Path(__file__).with_name("ordering_similarity.py")


def vendored_sha256() -> str:
    return hashlib.sha256(VENDORED.read_bytes()).hexdigest()


def tppl_order(ids: list[str], D: np.ndarray, repetitions: int = REPETITIONS, seed: int = SEED) -> list[str]:
    """Identifiers in TPPL order (the same order applies to rows and columns)."""
    import scipy.cluster.hierarchy  # noqa: F401  (required; see module docstring)
    from . import ordering_similarity as ref
    D = np.asarray(D, dtype=float)
    n = len(ids)
    if D.shape != (n, n):
        raise ValueError(f"matrix shape {D.shape} does not match {n} identifiers")
    if len(set(ids)) != n:
        raise ValueError("duplicate identifiers")
    if not np.all(np.isfinite(D)) or (D < 0).any():
        raise ValueError("incomplete or negative dissimilarities: refusing to impute missing comparisons")
    if np.max(np.abs(D - D.T)) > 1e-9 or np.max(np.abs(np.diag(D))) > 1e-9:
        raise ValueError("matrix is not a symmetric dissimilarity with a zero diagonal")
    perm = sorted(range(n), key=lambda i: ids[i])
    sids = [ids[i] for i in perm]
    S = D[np.ix_(perm, perm)]
    state = random.getstate()
    try:
        order = [int(i) for i in ref.treePenalizedPathLength(S, repetitions, seed)]
    finally:
        random.setstate(state)
    if sorted(order) != list(range(n)):
        raise ValueError("seriation did not return a permutation")
    return [sids[i] for i in order]


def seriation_record(ids: list[str], D: np.ndarray) -> dict[str, Any]:
    order = tppl_order(ids, D)
    return {"method": "tree-penalized path length (TPPL) seriation", "repetitions": REPETITIONS, "seed": SEED,
            "input_order": "instance identifiers sorted lexicographically before the call",
            "order": order, "missing_pairs": 0,
            "missing_values": "none: the matrix must be complete (missing comparisons are never imputed as maximal differences)",
            "reference": REFERENCE | {"vendored_sha256": vendored_sha256()}}
