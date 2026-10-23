"""How do the SAM riboswitch classes recognise SAM? Contact chemistry per ligand moiety and RNA edge, across all classes.

Input: the family directories from sam_classes.py (out/sam/<ACC>/) and data/families/RF00162 (SAM-I).
1. Ligand copies and PDBe Arpeggio contacts per family (ecifp_rna.collect); sites on the family's own Rfam positions and binding
   modes (ecifp_rna.cluster), per family (positions are not comparable across classes with different folds).
2. Every contact split by ligand moiety (atom names of SAM/SAH/SFG: adenine, ribose, sulfonium centre, amino acid) and by RNA part
   (base_edges.edges: Watson-Crick/Hoogsteen/sugar edge, corners, base face, ribose, phosphate; same rules as base_edges.py).
3. Ligand conformation from coordinates (PDBe ModelServer, the copy's residue only): glycosidic torsion chi (O4'-C1'-N9-C4) and the
   distance from the adenine N6 to the amino-acid CA (short = folded/compact, long = extended).
4. Water- and ion-mediated contacts from coordinates (ModelServer residueSurroundings, 7 A, model 1): a ligand N/O atom within 3.5 A of
   a water O that is within 3.5 A of an N/O atom of a family RNA chain; a metal within the fingeRNAt cutoff (Mg 3.2, Na 3.6, K 3.9,
   other 3.5 A) of a ligand atom, and whether the same metal is within the cutoff of an RNA atom (inner-sphere bridge).
5. ECIFP between all primary-site copies of all classes (position-free, so comparable across folds).

usage: sam_recognition.py <repo> <cache_dir> <out_dir>
"""
import json, sys, collections, itertools
from pathlib import Path
import numpy as np, requests, gemmi

REPO, CACHE, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
sys.path[:0] = [str(Path(__file__).parent)]
import ecifp_rna as E
from base_edges import edges, EDGE_ALL, EDGE, PARENT, RIBOSE, PHOS, DROP

CLASSES = {"RF00162": "SAM-I", "RF00634": "SAM-IV", "RF01725": "SAM-I/IV", "RF00521": "SAM-II", "RF01767": "SAM-III",
           "RF01826": "SAM-V", "RF02885": "SAM-VI", "RF01727": "SAM/SAH", "RF01057": "SAH"}
SAM_LIKE = {"SAM", "SAH", "SFG"}
MOIETY = {**{a: "adenine" for a in ("N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9")},
          **{a: "ribose" for a in ("C1'", "C2'", "O2'", "C3'", "O3'", "C4'", "O4'")},
          **{a: "S+ centre" for a in ("C5'", "SD", "CE", "CD", "NE")},     # SFG: CD replaces SD, NE (amine) replaces the methyl
          **{a: "amino acid" for a in ("N", "CA", "C", "O", "OXT", "CB", "CG")}}
MS = "https://www.ebi.ac.uk/pdbe/model-server/v1/{pdb}/{q}&model_nums=1&encoding=cif"
METAL_CUT = {"MG": 3.2, "NA": 3.6, "K": 3.9}
WATER = {"HOH", "DOD"}


def model(pdb, q, name):
    p = CACHE / f"ms_{name}.cif"
    if not p.exists():
        r = requests.get(MS.format(pdb=pdb, q=q), timeout=120); r.raise_for_status(); p.write_text(r.text)
    st = gemmi.make_structure_from_block(gemmi.cif.read_string(p.read_text()).sole_block())
    return st[0]


def torsion(a, b, c, d):
    return round(float(np.degrees(gemmi.calculate_dihedral(a, b, c, d))), 0)


def conformation(r):
    m = model(r["pdb"], f"atoms?auth_asym_id={r['chain']}&auth_seq_id={r['num']}&label_comp_id={r['ccd']}", f"{r['pdb']}_{r['chain']}_{r['num']}")
    at = {a.name: a.pos for ch in m for res in ch for a in res if a.altloc in ("\0", "A")}
    out = {}
    if all(k in at for k in ("O4'", "C1'", "N9", "C4")):
        chi = torsion(at["O4'"], at["C1'"], at["N9"], at["C4"])
        out["chi"] = chi
        out["glycosidic"] = "syn" if -30 <= chi <= 120 else "anti"     # syn about +60; anti about -160..-60 (O4'-C1'-N9-C4)
    if "N6" in at and "CA" in at:
        out["n6_ca"] = round(at["N6"].dist(at["CA"]), 1)
    out["modelled"] = sorted({MOIETY.get(a, "other") for a in at if not a.startswith("H")})
    out["amino_acid_complete"] = all(k in at for k in ("N", "CA", "C", "O", "CB", "CG"))
    return out


def centre_environment(r, fam_chains):
    """RNA N/O atoms within 4.0 A of the sulfur (SD; CD in SFG) and of the methyl carbon (CE), any contact type: the electrostatic
    recognition of the sulfonium cation (SAM) vs the neutral thioether (SAH) is by carbonyl/ring oxygens and nitrogens nearby"""
    m = model(r["pdb"], f"residueSurroundings?auth_asym_id={r['chain']}&auth_seq_id={r['num']}&label_comp_id={r['ccd']}&radius=7",
              f"sur_{r['pdb']}_{r['chain']}_{r['num']}")
    lig = {a.name: a.pos for ch in m for res in ch for a in res if ch.name == r["chain"] and res.seqid.num == r["num"] and res.name == r["ccd"]
           and a.altloc in ("\0", "A")}
    out = {}
    for c in ("SD", "CD", "CE"):
        if c not in lig or (c == "CD" and r["ccd"] != "SFG"):
            continue
        near = sorted((round(a.pos.dist(lig[c]), 1), f"{res.name}{res.seqid.num}:{a.name}", part(res.name, a.name))
                      for ch in m if ch.name in fam_chains for res in ch if res.name != r["ccd"] and res.name not in WATER and len(res) > 1 for a in res
                      if a.element.name in ("O", "N") and a.altloc in ("\0", "A") and a.pos.dist(lig[c]) <= 4.0)
        out[c] = near
    return out


def bridges(r, fam_chains):
    """water- and metal-mediated ligand-RNA contacts around one copy"""
    m = model(r["pdb"], f"residueSurroundings?auth_asym_id={r['chain']}&auth_seq_id={r['num']}&label_comp_id={r['ccd']}&radius=7",
              f"sur_{r['pdb']}_{r['chain']}_{r['num']}")
    lig, wat, met, rna = [], [], [], []
    for ch in m:
        for res in ch:
            for a in res:
                if a.altloc not in ("\0", "A"):
                    continue
                if ch.name == r["chain"] and res.seqid.num == r["num"] and res.name == r["ccd"]:
                    lig.append((a.name, a.element.name, a.pos))
                elif res.name in WATER and a.element.name == "O":
                    wat.append((f"{res.name}{res.seqid.num}", a.pos))
                elif a.element.is_metal and len(res) == 1:
                    met.append((res.name, a.pos))
                elif ch.name in fam_chains:
                    rna.append((ch.name, res.name, res.seqid.num, (res.seqid.icode or " ").strip(), a.name, a.element.name, a.pos))
    pol_l = [x for x in lig if x[1] in ("N", "O")]
    pol_r = [x for x in rna if x[5] in ("N", "O")]
    w_out, m_out = [], []
    for wn, wp in wat:
        lw = [l[0] for l in pol_l if l[2].dist(wp) <= 3.5]
        rw = [x for x in pol_r if x[6].dist(wp) <= 3.5]
        if lw and rw:
            w_out.append(dict(water=wn, ligand_atoms=lw, rna=[f"{x[0]}:{x[1]}{x[2]}{x[3]}:{x[4]}" for x in rw],
                              rna_parts=sorted({part(x[1], x[4]) for x in rw})))
    for mn, mp in met:
        cut = METAL_CUT.get(mn, 3.5)
        lm = [l[0] for l in lig if l[2].dist(mp) <= cut]
        if lm:
            rm = [x for x in rna if x[6].dist(mp) <= cut]
            m_out.append(dict(metal=mn, ligand_atoms=lm, inner_sphere_rna=[f"{x[0]}:{x[1]}{x[2]}{x[3]}:{x[4]}" for x in rm]))
    return w_out, m_out


def part(resname, atom):
    """RNA part of one atom, geometry only (used for water bridges, where the contact is polar by construction)"""
    if atom in RIBOSE or atom == "O2'":
        return "ribose" if atom != "O2'" else "sugar edge (O2')"
    if atom in PHOS:
        return "phosphate"
    e = EDGE.get(PARENT.get(resname, resname), {}).get(atom, "")
    return {"W": "Watson-Crick edge", "H": "Hoogsteen edge", "S": "sugar edge", "WH": "W/H corner", "WS": "W/S corner"}.get(e, "base other")


EDGE_NAME = {"W": "Watson-Crick edge", "H": "Hoogsteen edge", "S": "sugar edge", "F": "base face", "R": "ribose", "P": "phosphate",
             "WH": "W/H corner", "WS": "W/S corner"}


def moiety_contacts(r):
    """(ligand moiety, RNA part) pairs and the Rfam position of each, for one copy (Arpeggio contacts, base_edges rules)"""
    it = json.load(open(CACHE / f"bli_{r['pdb']}_{r['chain']}_{r['num']}.json"))[r["pdb"]][0]["interactions"]
    pairs = set()
    for x in it:
        types = set(x["interaction_details"])
        if not types - DROP:
            continue
        e = x["end"]
        if e["chain_id"] not in r["_chains"]:
            continue
        ks = edges(e["chem_comp_id"], e["atom_names"], types)
        if len(e["atom_names"]) == 1 and e["atom_names"][0] in EDGE.get(PARENT.get(e["chem_comp_id"], e["chem_comp_id"]), {}):
            two = EDGE[PARENT.get(e["chem_comp_id"], e["chem_comp_id"])][e["atom_names"][0]]
            if len(two) == 2 and ks:
                ks = {two}                                   # corner atom: report as a corner (no per-base suppression here)
        mo = {MOIETY.get(a, "other") for a in x["ligand_atoms"]}
        for k in ks:
            for m in mo:
                pairs.add((m, EDGE_NAME.get(k, k)))
    return pairs


def main():
    E.configure(CACHE, OUT)
    import os; os.chdir(REPO)
    allrows, report = [], dict(classes={})
    for acc, name in CLASSES.items():
        fam = None if acc == "RF00162" else OUT / "sam" / acc
        rows = E.collect(acc, fam)
        famdir = fam or REPO / "data/families" / acc
        chains = collections.defaultdict(set)
        for i in json.load(open(Path(famdir) / "instances.json")):
            chains[i["pdb_id"]].add(i["auth_asym_id"])
        ok, ecifp, lig_only, ov = E.cluster(rows) if any(r["status"] == "ok" for r in rows) else ([], None, None, None)
        sam = [r for r in ok if r["ccd"] in SAM_LIKE]
        main_site = collections.Counter(r["site"] for r in sam).most_common(1)[0][0] if sam else None
        prim = [r for r in sam if r["site"] == main_site]
        for r in prim:
            r["_chains"] = chains[r["pdb"]]
            r["class"] = name
            r["moiety_parts"] = sorted(moiety_contacts(r))
            r["conformation"] = conformation(r)
            r["water_bridges"], r["metals"] = bridges(r, chains[r["pdb"]])
            r["centre_env"] = centre_environment(r, chains[r["pdb"]])
            del r["_chains"]
        n = len(prim)
        mp = collections.Counter(p for r in prim for p in map(tuple, r["moiety_parts"]))
        pos = collections.Counter(p for r in prim for p in r["rfam_positions"])
        wb = collections.Counter((m, prt) for r in prim for w in r["water_bridges"] for m in {MOIETY.get(a, "other") for a in w["ligand_atoms"]} for prt in w["rna_parts"])
        wb_copies = collections.Counter()
        for r in prim:
            for k in {(m, prt) for w in r["water_bridges"] for m in {MOIETY.get(a, "other") for a in w["ligand_atoms"]} for prt in w["rna_parts"]}:
                wb_copies[k] += 1
        mets = collections.Counter()
        for r in prim:
            for k in {(x["metal"], MOIETY.get(a, "other"), "inner-sphere RNA" if x["inner_sphere_rna"] else "ligand only") for x in r["metals"] for a in x["ligand_atoms"]}:
                mets[k] += 1
        conf = [r["conformation"] for r in prim]
        report["classes"][acc] = dict(
            name=name, statuses=dict(collections.Counter(r["status"] for r in rows)),
            ligands_all=dict(collections.Counter(r["ccd"] for r in rows)), sites={f"{s}:{c}": v for (s, c), v in collections.Counter((r["site"], r["ccd"]) for r in ok).most_common()} if ok else {},
            primary_site=main_site, n=n, ligands=dict(collections.Counter(r["ccd"] for r in prim)),
            entries=sorted({r["pdb"] for r in prim}),
            positions={str(p): k for p, k in sorted(pos.items())},
            moiety_parts={f"{m} | {p}": k for (m, p), k in sorted(mp.items(), key=lambda x: (x[0][0], -x[1]))},
            water_bridges={f"{m} | {p}": k for (m, p), k in sorted(wb_copies.items(), key=lambda x: (x[0][0], -x[1]))},
            metals={" | ".join(k): v for k, v in mets.items()},
            conformation=dict(chi=[c.get("chi") for c in conf], glycosidic=dict(collections.Counter(c.get("glycosidic") for c in conf)),
                              n6_ca=[c.get("n6_ca") for c in conf],
                              amino_acid_modelled=sum(1 for c in conf if c.get("amino_acid_complete"))),
            centre_env={f"{r['pdb']}_{r['chain']}_{r['num']}_{r['ccd']}": r["centre_env"] for r in prim},
            modes=dict(collections.Counter(r["mode"] for r in prim)))
        allrows += prim
        print(f"== {acc} {name}: statuses {report['classes'][acc]['statuses']}; primary site {n} copies {report['classes'][acc]['ligands']}")
    # ECIFP across classes on the primary-site copies
    L = [set(r["lig_bits"]) for r in allrows]; P = [set(r["rna_bits"]) for r in allrows]
    k = len(allrows)
    ec = np.array([[E.tani(L[i], L[j]) ** 0.4 * E.tani(P[i], P[j]) ** 0.6 for j in range(k)] for i in range(k)])
    po = np.array([[E.tani(P[i], P[j]) for j in range(k)] for i in range(k)])
    names = list(dict.fromkeys(r["class"] for r in allrows))
    idx = {c: [i for i, r in enumerate(allrows) if r["class"] == c] for c in names}
    def med(M, a, b):
        v = [M[i, j] for i in idx[a] for j in idx[b] if i != j]
        return round(float(np.median(v)), 2) if v else None
    report["class_ecifp_median"] = {a: {b: med(ec, a, b) for b in names} for a in names}
    report["class_polymer_tanimoto_median"] = {a: {b: med(po, a, b) for b in names} for a in names}
    report["rows"] = allrows
    report["labels"] = [f"{r['class']}_{r['pdb']}_{r['chain']}_{r['num']}_{r['ccd']}" for r in allrows]
    report["ecifp"] = ec.round(3).tolist()
    json.dump(report, open(OUT / "sam_recognition.json", "w"), indent=1)
    print("ECIFP class medians:")
    print("          " + " ".join(f"{c[:8]:>8}" for c in names))
    for a in names:
        print(f"{a[:9]:<9} " + " ".join(f"{str(report['class_ecifp_median'][a][b]):>8}" for b in names))


main()
