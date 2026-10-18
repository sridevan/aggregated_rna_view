"""TPPL seriation (Structural Dissimilarity mode): the vendored reference, the wrapper's input handling, and the stored family orders."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pytest
from rna_explorer.seriation import REFERENCE, SEED, tppl_order, vendored_sha256

pytest.importorskip("scipy.cluster.hierarchy")
from rna_explorer import ordering_similarity as ref  # noqa: E402

FAMILIES = Path(__file__).resolve().parents[1] / "data" / "families"
# sridevan/ordering_analysis @ 37f63c9, ordering_similarity.py, byte for byte
REFERENCE_SHA256 = "1215f5802987ecc118bf6bf0b45a5a07ca46c93946774e86f71bd164d1f8ebc2"
# the reference's own test dataset (testOrdering: 50 uniform points in 3D, seed 2276393), ordered by the unmodified reference
REF_ORDER_100 = [15, 37, 43, 21, 11, 6, 8, 25, 48, 23, 44, 24, 13, 35, 1, 45, 26, 5, 7, 29, 18, 19, 31, 39, 40, 0, 14, 12, 20, 47, 36, 3, 2, 10, 17, 9, 41, 34,
                 27, 38, 30, 4, 49, 33, 28, 46, 16, 42, 32, 22]
REF_ORDER_10 = [22, 32, 42, 1, 3, 36, 47, 0, 14, 12, 20, 45, 26, 7, 5, 18, 29, 35, 13, 24, 2, 17, 10, 38, 27, 34, 9, 41, 30, 4, 49, 28, 33, 16, 46, 40, 39, 31,
                19, 23, 44, 48, 25, 8, 6, 11, 21, 43, 37, 15]
IDS = [f"p{i:02d}" for i in range(50)]   # sorted order == index order


def _uniform():
    return ref.generateUniformDataset(50, 3, 2276393)[1]


def _square(upper, n):
    D = np.zeros((n, n)); k = 0
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = upper[k]; k += 1
    return D


def test_vendored_file_is_the_unmodified_reference():
    assert vendored_sha256() == REFERENCE_SHA256 == hashlib.sha256((Path(ref.__file__)).read_bytes()).hexdigest()
    assert REFERENCE["commit"].startswith("37f63c9") and SEED == 39873


def test_reference_ordering_on_its_test_dataset():
    assert [int(i) for i in ref.treePenalizedPathLength(_uniform(), 100, 39873)] == REF_ORDER_100
    assert tppl_order(IDS, _uniform(), repetitions=100) == [IDS[i] for i in REF_ORDER_100]
    assert tppl_order(IDS, _uniform()) == [IDS[i] for i in REF_ORDER_10]          # the app's call: 10 repetitions, seed 39873


def test_deterministic_permutation_and_independent_of_input_order():
    D = _uniform()
    a, b = tppl_order(IDS, D), tppl_order(IDS, D)
    assert a == b and sorted(a) == IDS and len(set(a)) == 50
    p = np.random.default_rng(7).permutation(50)
    assert tppl_order([IDS[i] for i in p], D[np.ix_(p, p)]) == a
    # rows and columns refer to the right instances: reordering by the result reproduces each pair's value
    pos = {s: i for i, s in enumerate(IDS)}
    R = D[np.ix_([pos[s] for s in a], [pos[s] for s in a])]
    for r in range(0, 50, 7):
        for c in range(0, 50, 5):
            assert R[r, c] == D[pos[a[r]], pos[a[c]]]


def test_global_random_state_untouched():
    import random
    random.seed(5); before = random.random(); random.seed(5)
    tppl_order(IDS, _uniform())
    assert random.random() == before


@pytest.mark.parametrize("bad", ["nan", "negative", "asymmetric", "diagonal", "shape", "duplicate"])
def test_incomplete_or_invalid_input_is_refused(bad):
    D = _uniform(); ids = list(IDS)
    if bad == "nan":
        D[3, 9] = D[9, 3] = np.nan          # the reference would impute the maximum; a missing comparison is not a maximal difference
    elif bad == "negative":
        D[3, 9] = D[9, 3] = -1.0
    elif bad == "asymmetric":
        D[3, 9] += 0.5
    elif bad == "diagonal":
        D[4, 4] = 0.1
    elif bad == "shape":
        ids = ids[:-1]
    else:
        ids[1] = ids[0]
    with pytest.raises(ValueError):
        tppl_order(ids, D)


def test_small_matrices():
    assert tppl_order(["b", "a"], np.array([[0, 1.0], [1.0, 0]])) == ["a", "b"]
    assert sorted(tppl_order(["c", "a", "b"], np.array([[0, 1, 2.0], [1, 0, 3], [2, 3, 0]]))) == ["a", "b", "c"]


@pytest.mark.parametrize("acc", ["RF00167", "RF00162", "RF00059", "RF00005", "RF00177"])
def test_stored_family_orders(acc):
    d = json.load(open(FAMILIES / acc / "structure_clustering.json"))
    ids = d["dataset"]["instances"]; m = d["metrics"]["discrepancy"]; s = m["seriation"]
    assert sorted(s["order"]) == sorted(ids) and len(s["order"]) == len(ids)
    assert s["seed"] == 39873 and s["repetitions"] == 10 and s["missing_pairs"] == 0 and s["reference"]["vendored_sha256"] == REFERENCE_SHA256
    D = _square(m["values_upper"], len(ids))
    assert np.allclose(D, D.T) and np.all(np.diag(D) == 0) and np.all(np.isfinite(D))
    assert tppl_order(ids, D) == s["order"]      # reproducible from the stored values
