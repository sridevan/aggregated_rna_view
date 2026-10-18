"""Independent numerical check of the FreeSASA per-nucleotide values: a self-contained numpy Shrake-Rupley implementation
(golden-spiral sphere points, KD-tree-free brute force) on the same atoms/radii, compared per nucleotide and per chain."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
from rna_explorer.sasa import load_chain, RADII, PROBE_RADIUS
from rna_explorer.cache import Cache
from rna_explorer.pdbe import get_updated_cif

def sphere_points(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)

def shrake_rupley(xyz, radii, probe=PROBE_RADIUS, n=2000):
    pts = sphere_points(n); R = radii + probe; out = np.zeros(len(xyz))
    for i in range(len(xyz)):
        d = np.linalg.norm(xyz - xyz[i], axis=1); nb = np.where((d < R[i] + R) & (d > 0))[0]
        surf = xyz[i] + R[i] * pts
        buried = np.zeros(n, dtype=bool)
        for j in nb:
            buried |= np.linalg.norm(surf - xyz[j], axis=1) < R[j]
        out[i] = 4 * np.pi * R[i] ** 2 * (~buried).mean()
    return out

def main(seq_name="1y26_1_X_A"):
    pilot = json.load(open(Path(__file__).parent / "sasa.json"))
    m = next(x for x in pilot["instances"] if x["seq_name"] == seq_name)
    cache = Cache(ROOT / "cache"); chain = load_chain(get_updated_cif(cache, m["pdb_id"]), m["label_asym_id"])
    atoms = [(r["label_seq_id"], a) for r in chain["residues"] for a in r["atoms"] if a["element"] in RADII]
    xyz = np.array([[a["x"], a["y"], a["z"]] for _, a in atoms]); radii = np.array([RADII[a["element"]] for _, a in atoms])
    sasa = shrake_rupley(xyz, radii)
    per_nt = {}
    for (i, _), s in zip(atoms, sasa): per_nt[i] = per_nt.get(i, 0.0) + s
    ref = {n["label_seq_id"]: n["sasa_A2"] for n in m["nucleotides"] if n["sasa_A2"] is not None}
    diffs = [(i, ref[i], per_nt[i]) for i in ref if i in per_nt]
    rel = [abs(a - b) / max(a, 1e-9) for _, a, b in diffs]
    out = {"seq_name": seq_name, "n_atoms": len(atoms), "freesasa_total_A2": round(sum(ref.values()), 2), "numpy_total_A2": round(float(sum(per_nt.values())), 2),
           "total_relative_difference": round(abs(sum(ref.values()) - sum(per_nt.values())) / sum(ref.values()), 5),
           "per_nucleotide_max_abs_diff_A2": round(max(abs(a - b) for _, a, b in diffs), 3), "per_nucleotide_median_rel_diff": round(float(np.median(rel)), 5), "per_nucleotide_max_rel_diff": round(max(rel), 5),
           "n_points_numpy": 2000, "worst": sorted(((round(abs(a - b), 2), i, a, round(b, 2)) for i, a, b in diffs), reverse=True)[:3]}
    print(json.dumps(out, indent=1)); json.dump(out, open(Path(__file__).parent / f"independent_check_{seq_name}.json", "w"), indent=1)

if __name__ == "__main__":
    main(*sys.argv[1:])
