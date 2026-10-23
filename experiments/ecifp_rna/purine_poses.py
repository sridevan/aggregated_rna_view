"""Do the purine-riboswitch ligands sit in the same pose? Pocket-superposed comparison of ligand ring placement.

For every ligand copy in the main RF00167 site (out/RF00167.json, site of ADE):
1. coordinates of the ligand and of the contacted RNA chain (model 1) from the PDBe ModelServer (streamed, not stored);
2. superposition onto the reference pocket (1Y26 X, ADE 90) on C1' C4' P of the Rfam positions whose C4' lies within 12 A of the
   reference ligand and that the chain has (Kabsch);
3. ligand ring atoms identified on the CCD definition with SMARTS: purine core N1 C2 N3 C4 C5 C6 N7 C8 N9 (any C/N at the ring
   positions, so 7-deaza and 8-aza analogues match); otherwise a 6-membered aromatic ring (pyrimidines) or no core;
4. per copy: pocket RMSD, ligand centroid shift and, for purine cores, core RMSD to the reference core and RMSD of the
   best ring flip (core atoms permuted by each of the in-plane 180-degree flips), so a flipped pose is reported as such.

usage: purine_poses.py <repo> <cache_dir> <out_json>
"""
import json, sys, collections
from pathlib import Path
import numpy as np, requests, gemmi
from rdkit import Chem
from pdbeccdutils.core import ccd_reader

REPO, CACHE, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
MS = "https://www.ebi.ac.uk/pdbe/model-server/v1/{pdb}/atoms?{q}&encoding=cif"
REF = ("1y26", "X", 90)
CORE = ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"]
PURINE = Chem.MolFromSmarts("[#6,#7;R]1~[#6;R]~[#6,#7;R]~[#6;R]2~[#6;R](~[#6;R]~1)~[#6,#7;R]~[#6,#7;R]~[#6,#7;R]~2")
# in-plane flips of the purine core as index permutations (identity excluded): Watson-Crick/Hoogsteen-like swaps are
# not symmetric for a bicycle, so only the two mirror flips across the long and short axes are considered
FLIPS = {"flip_long_axis": [2, 1, 0, 3, 4, 5, 8, 7, 6], "flip_short_axis": None}


def structure(pdb, q):
    r = requests.get(MS.format(pdb=pdb, q=q), timeout=120)
    r.raise_for_status()
    st = gemmi.make_structure_from_block(gemmi.cif.read_string(r.text).sole_block())
    return st[0]


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return lambda X: (X - pc) @ R.T + qc


def core_names(ccd):
    p = CACHE / f"ccd_{ccd}.cif"
    mol = ccd_reader.read_pdb_cif_file(str(p)).component.mol
    m = mol.GetSubstructMatch(PURINE)
    name = lambda i: mol.GetAtomWithIdx(i).GetProp("name")
    if m:
        return "purine", [name(i) for i in m]
    ring6 = [r for r in mol.GetRingInfo().AtomRings() if len(r) == 6 and all(mol.GetAtomWithIdx(i).GetSymbol() in "CN" for i in r)]
    return ("ring6", [name(i) for i in ring6[0]]) if ring6 else ("none", [])


def main():
    fam = REPO / "data/families/RF00167"
    res = json.load(open(REPO / "experiments/ecifp_rna/out/RF00167.json"))
    ok = [r for r in res["rows"] if r["status"] == "ok"]
    site = collections.Counter(r["site"] for r in ok if r["ccd"] == "ADE").most_common(1)[0][0]
    rows = [r for r in ok if r["site"] == site]
    inst = json.load(open(fam / "instances.json"))
    sn = {(i["pdb_id"], i["auth_asym_id"]): i["seq_name"] for i in inst}

    def rna_chain(r):
        it = json.load(open(CACHE / f"bli_{r['pdb']}_{r['chain']}_{r['num']}.json"))[r["pdb"]][0]["interactions"]
        c = collections.Counter(x["end"]["chain_id"] for x in it if (r["pdb"], x["end"]["chain_id"]) in sn)
        return c.most_common(1)[0][0]

    def pocket(pdb, ch):
        m = json.load(open(fam / "mappings" / f"{sn[(pdb, ch)]}.json"))
        pos = {(str(x["author_residue_number"]), (x.get("author_insertion_code") or "").strip()): x["rfam_consensus_pos"] for x in m["residues"] if x["rfam_consensus_pos"]}
        model = structure(pdb, f"auth_asym_id={ch}&model_nums=1")
        out = {}
        for chn in model:
            for rr in chn:
                p = pos.get((str(rr.seqid.num), rr.seqid.icode.strip()))
                if p is None:
                    continue
                for a in rr:
                    if a.name in ("C1'", "C4'", "P") and (a.altloc in ("\0", "A")):
                        out[(p, a.name)] = np.array(a.pos.tolist())
        return out

    def ligand(r):
        model = structure(r["pdb"], f"auth_asym_id={r['chain']}&auth_seq_id={r['num']}&model_nums=1")
        return {a.name: np.array(a.pos.tolist()) for chn in model for rr in chn for a in rr if a.altloc in ("\0", "A") and a.element.name != "H"}

    ref = next(r for r in rows if (r["pdb"], r["chain"], r["num"]) == REF)
    ref_p, ref_l = pocket(REF[0], rna_chain(ref)), ligand(ref)
    ref_core = np.array([ref_l[n] for n in core_names("ADE")[1]])
    lc = np.array(list(ref_l.values())).mean(0)
    near = sorted({p for (p, a), x in ref_p.items() if a == "C4'" and np.linalg.norm(x - lc) <= 12})
    out = []
    for r in rows:
        rec = dict(pdb=r["pdb"], chain=r["chain"], num=r["num"], ccd=r["ccd"], mode=r["mode"])
        try:
            ch = rna_chain(r); P = pocket(r["pdb"], ch); L = ligand(r)
            keys = [k for k in ref_p if k[0] in near and k in P]
            f = kabsch(np.array([P[k] for k in keys]), np.array([ref_p[k] for k in keys]))
            rec["pocket_atoms"] = len(keys)
            rec["pocket_rmsd"] = round(float(np.sqrt(((f(np.array([P[k] for k in keys])) - np.array([ref_p[k] for k in keys])) ** 2).sum(1).mean())), 2)
            Lx = {n: f(x[None])[0] for n, x in L.items()}
            rec["centroid_shift"] = round(float(np.linalg.norm(np.array(list(Lx.values())).mean(0) - ref_core.mean(0))), 2)
            kind, names = core_names(r["ccd"])
            rec["core"] = kind
            if kind == "purine" and all(n in Lx for n in names):
                C = np.array([Lx[n] for n in names])
                rmsd = lambda A: round(float(np.sqrt(((A - ref_core) ** 2).sum(1).mean())), 2)
                rec["core_rmsd"] = rmsd(C)
                rec["flip_rmsd"] = rmsd(C[FLIPS["flip_long_axis"]])
            elif kind == "ring6" and all(n in Lx for n in names):
                C = np.array([Lx[n] for n in names])
                rec["ring6_centroid_to_ref_6ring"] = round(float(np.linalg.norm(C.mean(0) - ref_core[:6].mean(0))), 2)
        except Exception as e:   # noqa: BLE001 - reported per copy
            rec["error"] = str(e)[:200]
        out.append(rec)
        print(rec, flush=True)
    json.dump(dict(reference=REF, pocket_positions=near, rows=out), open(OUT, "w"), indent=1)


main()
