import math
import numpy as np
import pytest
from rna_explorer.orientation import kabsch, select_pairs, orient_instance, fit_quality


def rot(axis, deg):
    a = np.array(axis, float); a /= np.linalg.norm(a); th = math.radians(deg)
    k = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(th) * k + (1 - math.cos(th)) * (k @ k)


def cloud(n=25, seed=1):
    return np.random.default_rng(seed).normal(size=(n, 3)) * 10


def test_identity_and_known_rotation():
    p = cloud()
    k = kabsch(p, p)
    assert np.allclose(k.rotation, np.eye(3), atol=1e-9) and np.allclose(k.translation, 0, atol=1e-9) and k.rmsd < 1e-9
    r, t = rot([1, 2, 3], 37), np.array([5.0, -2.0, 11.0])
    q = (r @ p.T).T + t                                       # reference = R p + t
    k2 = kabsch(p, q)
    assert np.allclose(k2.rotation, r, atol=1e-8) and np.allclose(k2.translation, t, atol=1e-8) and k2.rmsd < 1e-8
    assert abs(k2.determinant - 1) < 1e-9
    # convention check: applying R,t to target reproduces reference
    assert np.allclose((np.array(k2.rotation) @ p.T).T + np.array(k2.translation), q)


def test_reflection_rejected_and_rmsd_reproducible():
    p = cloud(seed=2)
    q = p * np.array([-1, 1, 1])                             # mirror image: not reachable by a proper rotation
    k = kabsch(p, q)
    assert abs(k.determinant - 1) < 1e-9 and k.rmsd > 1.0
    noisy = p + np.random.default_rng(3).normal(scale=0.5, size=p.shape)
    a, b = kabsch(p, noisy), kabsch(p, noisy)
    assert a.rmsd == b.rmsd and 0.3 < a.rmsd < 1.2
    with pytest.raises(ValueError):
        kabsch(p[:2], p[:2])


def _mapping(pdb, seq_ids, cons, pp="*", coords=True):
    return {"pdb_id": pdb, "entity_id": 1, "label_asym_id": "A", "auth_asym_id": "A",
            "residues": [{"polymer_seq_position": s, "alignment_state": "match", "rfam_consensus_pos": c, "posterior_prob_char": pp, "has_coordinates": coords} for s, c in zip(seq_ids, cons)]}


def test_pair_selection_missing_atoms_low_pp_and_identifiers():
    cons = list(range(1, 13))
    ref = _mapping("1ref", list(range(1, 13)), cons)
    tgt = _mapping("2tgt", list(range(5, 17)), cons)          # same consensus, different polymer numbering
    pts = cloud(12, seed=4)
    ref_xyz = {s: {"C4'": tuple(pts[i])} for i, s in enumerate(range(1, 13))}
    tgt_xyz = {s: {"C4'": tuple(pts[i])} for i, s in enumerate(range(5, 17))}
    del tgt_xyz[9]["C4'"]; tgt_xyz[9]["P"] = tuple(pts[4] + 0.1); ref_xyz[5]["P"] = tuple(pts[4] + 0.1)   # fallback atom (P) present in both residues
    tgt_xyz.pop(10)                                                         # residue without any atom
    sel = select_pairs(ref, tgt, ref_xyz, tgt_xyz)
    assert sel.n_shared_consensus == 12 and len(sel.pairs) == 11 and sel.excluded["missing_atom"] == 1
    assert sel.atoms_used == {"C4'": 10, "P": 1}
    assert all(t == r + 4 for _, r, t, _ in sel.pairs)                     # target label_seq_ids preserved (not re-numbered)
    rec = orient_instance(ref, tgt, ref_xyz, tgt_xyz, 12, 12)
    assert rec["fit_quality"] in ("good", "moderate") and rec["n_pairs"] == 11 and rec["check"]["label_seq_id"] == 5
    # low PP in the target is excluded under the strict filter; relaxed when too few remain
    tgt_low = _mapping("2tgt", list(range(5, 17)), cons, pp="9")
    sel2 = select_pairs(ref, tgt_low, ref_xyz, tgt_xyz)
    assert sel2.pp_filter == "9" and len(sel2.pairs) == 11
    tgt_bad = _mapping("2tgt", list(range(5, 17)), cons, pp="5")
    sel3 = select_pairs(ref, tgt_bad, ref_xyz, tgt_xyz)
    assert len(sel3.pairs) == 0 and orient_instance(ref, tgt_bad, ref_xyz, tgt_xyz, 12, 12)["fit_quality"] == "insufficient"
    assert fit_quality(1.0, 5) == "insufficient" and fit_quality(5.0, 30) == "poor"
