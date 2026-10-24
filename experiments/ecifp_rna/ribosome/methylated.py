"""Resistance-methylated ribosomes in the PDB: where does the added methyl sit relative to the drug, and what changes compared with the
unmethylated ribosome bound to the same drug?

Entries (found by inventory of modified rRNA nucleotides, 2026-10-24; PDBe API modified_AA_or_NA on rRNA entries of RF00177/RF02541):
- Erm   m6,6A2058 (CCD MA6 at 23S 2058): 6XHV (apo), 7RQ9, 8FC4, 8FC5, 8FC6, 8UD7, 8VTX
- Cfr   m2,8A2503 (CCD 8AH at 23S 2503): 7LVK, 8G29, 8G2A (apo), 7S1K, 8G2B, 8G2C, 8UD8
- NpmA  m1A1408   (CCD 1MA at 16S 1408): 9B4Z, 9B51
(7CPJ m7G1405 is a stalled E. coli ribosome with no resistance enzyme in the study: treated as a modelling artefact, not used.)

Coordinates: PDBe ModelServer residueSurroundings (radius 10 A, the server maximum) around every copy of the drug, model 1, so no full
mmCIF is stored. Per drug copy:
- methyl carbons = heavy atoms of the modified residue that standard A lacks; nearest drug heavy atom to each;
- parent atom (N6 / C8 / N1) to drug, in both the methylated and the control structure;
- virtual methyls placed on the parent atom in the control (ideal in-plane geometry), nearest drug atom = the clash a pure steric model
  would predict; checked against the real methyl positions in the methylated structure;
- the pair superposed on RNA atoms within 10 A of the drug present in both (same chain id, residue number, atom name): drug RMSD over
  same-name atoms, largest per-atom drug shift, base shift of the modified nucleotide and of the residues lining the drug;
- RNA residues within 4 A of the drug in each structure (lost / gained);
- drug occupancy and mean B relative to the mean B of RNA atoms within 5 A.

usage: methylated.py <cache_dir> <out_json>   (needs gemmi and numpy)
"""
import json, sys, time, urllib.request
from pathlib import Path
import numpy as np
import gemmi

CACHE, OUT = Path(sys.argv[1]), Path(sys.argv[2])
CACHE.mkdir(parents=True, exist_ok=True)

# (methylated entry, control entry or None, drug CCD, enzyme)
PAIRS = [
    ("8fc4", "8fc1", "ERY", "Erm"), ("8fc5", "8fc2", "ZIT", "Erm"), ("8fc6", "8fc3", "TEL", "Erm"),
    ("8vtx", "8vtw", "A1AE1", "Erm"), ("7rq9", "7rq8", "6IF", "Erm"), ("8ud7", "8ud6", "WC9", "Erm"),
    ("8g2b", "7rq8", "6IF", "Cfr"), ("8g2c", "8g2d", "TYK", "Cfr"), ("8ud8", "8ud6", "WC9", "Cfr"), ("7s1k", "7s1i", "RD8", "Cfr"),
    ("9b4z", "9b50", "84G", "NpmA"), ("9b51", None, "GET", "NpmA"),
]
SITE = {"Erm": (2058, "MA6", "N6"), "Cfr": (2503, "8AH", "C8"), "NpmA": (1408, "1MA", "N1")}
STD_A = {"P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
         "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"}
BASE = ["N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"]


def fragment(pdb, ccd):
    p = CACHE / f"{pdb}_{ccd}.cif"
    if not p.exists():
        url = f"https://www.ebi.ac.uk/pdbe/model-server/v1/{pdb}/residueSurroundings?label_comp_id={ccd}&radius=10&encoding=cif&model_nums=1"
        for k in range(4):
            try:
                with urllib.request.urlopen(url, timeout=300) as r:
                    p.write_bytes(r.read())
                break
            except Exception:
                time.sleep(5 + 5 * k)
    st = gemmi.read_structure(str(p))
    st.setup_entities()
    return st[0]


def heavy(res):
    return [a for a in res if a.element.name not in ("H", "D")]


def xyz(a):
    return np.array([a.pos.x, a.pos.y, a.pos.z])


def unit(v):
    return v / np.linalg.norm(v)


def virtual_methyls(res, parent):
    """Ideal in-plane methyl carbon(s) on the parent atom of an A (or modified A) residue."""
    A = {a.name: xyz(a) for a in res}
    if parent == "C8":                                     # C8-CH3 along the external bisector of N7-C8-N9
        return [A["C8"] + 1.50 * unit(unit(A["C8"] - A["N7"]) + unit(A["C8"] - A["N9"]))]
    if parent == "N1":                                     # N1-CH3 along the external bisector of C2-N1-C6
        return [A["N1"] + 1.47 * unit(unit(A["N1"] - A["C2"]) + unit(A["N1"] - A["C6"]))]
    # N6: two methyls, C6-N6-C 120 deg, in the base plane
    u = unit(A["N6"] - A["C6"]); n = unit(np.cross(A["C5"] - A["C6"], A["N1"] - A["C6"]))
    out = []
    for s in (1, -1):
        th = np.radians(60) * s
        v = u * np.cos(th) + np.cross(n, u) * np.sin(th)
        out.append(A["N6"] + 1.45 * v)
    return out


def nearest(points, atoms):
    best = (9e9, None)
    for p in points:
        for a in atoms:
            d = float(np.linalg.norm(p - xyz(a)))
            if d < best[0]:
                best = (d, a.name)
    return best


def copies(model, ccd, site_num):
    """Each drug copy with the site residue of its own ribosome (nearest residue numbered site_num)."""
    out = []
    for ch in model:
        for r in ch:
            if r.name != ccd:
                continue
            lig = heavy(r)
            c = np.mean([xyz(a) for a in lig], axis=0)
            site = min(((float(np.linalg.norm(np.mean([xyz(a) for a in heavy(s)], axis=0) - c)), ch2.name, s)
                        for ch2 in model for s in ch2 if s.seqid.num == site_num and s.name in ("A", "MA6", "8AH", "1MA", "2MA", "6MZ")),
                       default=(None, None, None), key=lambda t: t[0] if t[0] is not None else 9e9)
            out.append(dict(chain=ch.name, num=r.seqid.num, lig=r, site_chain=site[1], site=site[2], site_dist=site[0]))
    return out


def rna_atoms(model, skip):
    d = {}
    for ch in model:
        for r in ch:
            if r.name in skip or len(r) < 5 or gemmi.find_tabulated_residue(r.name) is None:
                continue
            if not gemmi.find_tabulated_residue(r.name).is_nucleic_acid() and r.name not in ("MA6", "8AH", "1MA", "2MA", "6MZ"):
                continue
            for a in heavy(r):
                d[(ch.name, r.seqid.num, r.seqid.icode, a.name)] = (xyz(a), a.b_iso, r.name)
    return d


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return lambda X: (X - pc) @ R.T + qc


def lining(model, lig, cut):
    L = np.array([xyz(a) for a in heavy(lig)])
    out = set()
    for ch in model:
        for r in ch:
            if r is lig or r.name == lig.name or r.name == "HOH" or len(r) < 5:
                continue
            X = np.array([xyz(a) for a in heavy(r)])
            if len(X) and np.min(np.linalg.norm(X[:, None] - L[None], axis=2)) <= cut:
                out.add((ch.name, r.seqid.num, r.name))
    return out


def describe(model, cp, enzyme, methylated):
    num, mod, parent = SITE[enzyme]
    lig, site = cp["lig"], cp["site"]
    L = heavy(lig)
    row = dict(drug_chain=cp["chain"], drug_num=cp["num"], site_chain=cp["site_chain"], site_name=site.name if site else None,
               occupancy=round(float(np.mean([a.occ for a in L])), 2))
    pocket = [a.b_iso for ch in model for r in ch if r.name != lig.name and len(r) >= 5 for a in heavy(r)
              if min(np.linalg.norm(xyz(a) - xyz(b)) for b in L) <= 5]
    row["b_ratio"] = round(float(np.mean([a.b_iso for a in L]) / np.mean(pocket)), 2) if pocket else None
    if site is None:
        return row
    d, an = nearest([xyz(site[parent][0])], L)
    row["parent_to_drug"] = (round(d, 2), an)
    row["virtual_methyls_to_drug"] = [(round(x, 2), a) for x, a in (nearest([v], L) for v in virtual_methyls(site, parent))]
    if methylated:
        extra = [a for a in heavy(site) if a.name not in STD_A]
        row["methyls"] = []
        for a in extra:
            bonded = min(((float(np.linalg.norm(xyz(a) - xyz(b))), b.name) for b in heavy(site) if b.name in STD_A))[1]
            dd, aa = nearest([xyz(a)], L)
            row["methyls"].append(dict(atom=a.name, on=bonded, to_drug=round(dd, 2), drug_atom=aa, occ=round(a.occ, 2)))
        on_parent = [xyz(a) for a in extra if min(((float(np.linalg.norm(xyz(a) - xyz(b))), b.name) for b in heavy(site) if b.name in STD_A))[1] == parent]
        vm = virtual_methyls(site, parent)
        row["virtual_vs_real_A"] = [round(min(float(np.linalg.norm(v - r)) for r in on_parent), 2) for v in vm] if on_parent else None
    return row


def main():
    rows, secondary = [], {}
    for meth, ctrl, ccd, enz in PAIRS:
        num, mod, parent = SITE[enz]
        M = fragment(meth, ccd)
        # copies more than 15 A from the site residue are secondary sites (e.g. aminoglycosides at H69): counted, not compared
        allm = copies(M, ccd, num)
        mc = [c for c in allm if c["site_dist"] is not None and c["site_dist"] <= 15]
        secondary[meth] = len(allm) - len(mc)
        C = fragment(ctrl, ccd) if ctrl else None
        cc = [c for c in copies(C, ccd, num) if c["site_dist"] is not None and c["site_dist"] <= 15] if C else []
        for cp in mc:
            row = dict(methylated=meth, control=ctrl, drug=ccd, enzyme=enz, site=f"{'16S' if enz == 'NpmA' else '23S'} {num}",
                       meth=describe(M, cp, enz, True))
            # partner copy in the control: the drug copy whose site residue sits on the same chain id, else the first
            partner = next((c for c in cc if c["site_chain"] == cp["site_chain"]), cc[0] if cc else None)
            if partner and cp["site"] is not None and partner["site"] is not None:
                row["ctrl"] = describe(C, partner, enz, False)
                am, ac = rna_atoms(M, {ccd}), rna_atoms(C, {ccd})
                common = [k for k in am if k in ac and k[3] in STD_A]
                P = np.array([am[k][0] for k in common]); Q = np.array([ac[k][0] for k in common])
                f = kabsch(P, Q)
                row["superposition"] = dict(atoms=len(common), rmsd=round(float(np.sqrt(np.mean(np.sum((f(P) - Q) ** 2, 1)))), 2))
                lm = {a.name: xyz(a) for a in heavy(cp["lig"])}; lc = {a.name: xyz(a) for a in heavy(partner["lig"])}
                same = sorted(set(lm) & set(lc))
                if same:
                    dv = np.linalg.norm(f(np.array([lm[n] for n in same])) - np.array([lc[n] for n in same]), axis=1)
                    row["drug_shift"] = dict(atoms=len(same), of=len(lm), rmsd=round(float(np.sqrt(np.mean(dv ** 2))), 2),
                                             max=round(float(dv.max()), 2), max_atom=same[int(dv.argmax())])
                # per-residue base shifts for residues lining the drug in either structure
                lin = lining(M, cp["lig"], 4.0) | {(c, n, r) for c, n, r in lining(C, partner["lig"], 4.0)}
                shifts = []
                for ch, n, rn in sorted(lin, key=lambda t: t[1]):
                    ks = [k for k in common if k[0] == ch and k[1] == n and k[3] in BASE]
                    if ks:
                        v = np.linalg.norm(f(np.array([am[k][0] for k in ks])) - np.array([ac[k][0] for k in ks]), axis=1)
                        shifts.append((ch, n, rn, round(float(np.sqrt(np.mean(v ** 2))), 2)))
                row["lining_base_shift"] = shifts
                sm = {(c, n) for c, n, _ in lining(M, cp["lig"], 4.0)}; sc = {(c, n) for c, n, _ in lining(C, partner["lig"], 4.0)}
                row["contacts_4A"] = dict(meth=len(sm), ctrl=len(sc), lost=sorted(n for c, n in sc - sm), gained=sorted(n for c, n in sm - sc))
            rows.append(row)
            m = row["meth"]
            print(meth, ccd, enz, "copy", cp["chain"], "| methyls", [(x["on"], x["to_drug"]) for x in m.get("methyls", [])],
                  "| virt/real", m.get("virtual_vs_real_A"), "| occ", m["occupancy"], "B", m["b_ratio"],
                  "| ctrl parent", row.get("ctrl", {}).get("parent_to_drug"), "virt", row.get("ctrl", {}).get("virtual_methyls_to_drug"),
                  "occ", row.get("ctrl", {}).get("occupancy"), "B", row.get("ctrl", {}).get("b_ratio"),
                  "| shift", row.get("drug_shift"), "| sup", row.get("superposition"))
    json.dump(dict(source="PDBe ModelServer residueSurroundings radius 10 A, model 1; retrieved 2026-10-24", pairs=PAIRS, secondary_copies=secondary, rows=rows),
              open(OUT, "w"), indent=1, default=str)


main()
