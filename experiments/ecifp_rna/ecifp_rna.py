"""Feasibility test: ECIFP (PDBeurope/ecifp, a21de95) applied to RNA-ligand complexes of our riboswitch families.

Replicates ecifp's fingerprint exactly (PDBe Arpeggio contacts from bound_ligand_interactions, contact types minus
proximal/clash/vdw_clash/vdw, water and other bound components dropped; Morgan radius 2, 1024 bits from the PDBe CCD
component; ligand and polymer bits = bits of the contacting atoms; polymer residues merged by chem_comp_id; similarity =
Tanimoto_lig^0.4 * Tanimoto_pol^0.6) and adds the RNA analogue of ecifp's UniProt-numbered binding sites: contacted
residues mapped to Rfam consensus positions through our validated mappings.

usage: ecifp_rna.py <repo> <cache_dir> <out_dir> RF00167 [RF...]
"""
import json, sys, os, collections, tempfile, itertools
from pathlib import Path
import numpy as np, requests
from concurrent.futures import ThreadPoolExecutor
from rdkit.Chem import rdFingerprintGenerator
from pdbeccdutils.core import ccd_reader
from scipy.cluster.hierarchy import linkage, fcluster

CACHE = OUT = None   # set by main() or configure()


def configure(cache, out):
    global CACHE, OUT
    CACHE, OUT = Path(cache), Path(out)
    CACHE.mkdir(parents=True, exist_ok=True); OUT.mkdir(parents=True, exist_ok=True)
API = "https://www.ebi.ac.uk/pdbe/api/v2/pdb"
DROP = {"proximal", "clash", "vdw_clash", "vdw"}
R, SIZE = 2, 1024


def get_json(url, name):
    p = CACHE / name
    if p.exists():
        return json.loads(p.read_text())
    r = requests.get(url, timeout=60)
    d = r.json() if r.status_code == 200 else None
    p.write_text(json.dumps(d))
    return d


_comp = {}
def component(ccd):
    if ccd not in _comp:
        d = "prd" if ccd.startswith("PRD") else "clc" if ccd.startswith("CLC") else "ccd"
        sub = ccd[-1] if d != "ccd" else ccd[0]
        p = CACHE / f"ccd_{ccd}.cif"
        if not p.exists():
            r = requests.get(f"https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/{d}/{sub}/{ccd}/{ccd}.cif", timeout=60)
            p.write_text(r.text if r.status_code == 200 else "")
        try:
            _comp[ccd] = ccd_reader.read_pdb_cif_file(str(p)).component.mol
        except Exception:
            _comp[ccd] = None
    return _comp[ccd]


_gen = rdFingerprintGenerator.GetMorganGenerator(radius=R, fpSize=SIZE)
def bits(ccd, atoms):
    mol = component(ccd)
    if mol is None:
        return None
    ao = rdFingerprintGenerator.AdditionalOutput(); ao.AllocateAtomToBits()
    _gen.GetFingerprint(mol, additionalOutput=ao)
    a2b = ao.GetAtomToBits()
    on = set()
    for at in mol.GetAtoms():
        if at.HasProp("name") and at.GetProp("name") in atoms:
            on |= set(a2b[at.GetIdx()])
    return on


def tani(a, b):
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def collect(acc, fam=None):
    """ligand-copy rows (fingerprints, contacted positions) for one family directory"""
    fam = Path(fam) if fam else Path("data/families") / acc
    inst = json.load(open(fam / "instances.json"))
    lig = json.load(open(fam / "ligands.json"))
    skip = set(itertools.chain.from_iterable(lig["classification"]["classes"].values()))
    chains = collections.defaultdict(dict)   # pdb -> auth chain -> seq_name
    for i in inst:
        chains[i["pdb_id"].lower()][i["auth_asym_id"]] = i["seq_name"]
    # author residue -> Rfam consensus position, per chain
    pos = {}
    for pdb, ch in chains.items():
        for a, sn in ch.items():
            m = json.load(open(fam / "mappings" / f"{sn}.json"))
            for r in m["residues"]:
                pos[(pdb, a, str(r["author_residue_number"]), r.get("author_insertion_code") or "")] = r["rfam_consensus_pos"]
    copies = []
    def entry(pdb):
        bm = get_json(f"{API}/bound_molecules/{pdb}", f"bm_{pdb}.json") or {}
        out = []
        for b in bm.get(pdb, []):
            comp = [(l["chain_id"], l["author_residue_number"], l["chem_comp_id"]) for l in b["composition"]["ligands"]]
            for (c, n, ccd) in comp:
                if ccd in skip:
                    continue
                it = get_json(f"{API}/bound_ligand_interactions/{pdb}/{c}/{n}", f"bli_{pdb}_{c}_{n}.json")
                out.append((pdb, b.get("bm_id"), c, n, ccd, comp, it))
        return out
    with ThreadPoolExecutor(6) as ex:
        for res in ex.map(entry, sorted(chains)):
            copies.extend(res)
    rows = []
    for pdb, bmid, c, n, ccd, comp, it in copies:
        if not it or pdb not in it:
            rows.append(dict(pdb=pdb, chain=c, num=n, ccd=ccd, status="no_interactions")); continue
        bound = {(x[0], x[1], x[2]) for x in comp}
        lat, pat, other, rpos, offmap = set(), collections.defaultdict(set), collections.Counter(), set(), 0
        for x in it[pdb][0]["interactions"]:
            if not [t for t in x["interaction_details"] if t not in DROP]:
                continue
            e = x["end"]
            key = (e["chain_id"], e["author_residue_number"], e["chem_comp_id"])
            if key in bound or e["chem_comp_id"] in ("HOH", "DOD"):
                continue
            if e["chain_id"] not in chains[pdb]:
                other[e["chain_id"]] += 1; continue          # non-family polymer (protein, other RNA)
            lat |= set(x["ligand_atoms"]); pat[e["chem_comp_id"]] |= set(e["atom_names"])
            p = pos.get((pdb, e["chain_id"], str(e["author_residue_number"]), (e.get("author_insertion_code") or "").strip()))
            if p: rpos.add(p)
            else: offmap += 1
        lb = bits(ccd, lat)
        pb = set().union(*[bits(k, v) or set() for k, v in pat.items()]) if pat else set()
        rows.append(dict(pdb=pdb, chain=c, num=n, ccd=ccd, status="ok" if lat else "no_family_contacts", n_lig_atoms=len(lat),
                         lig_bits=sorted(lb or []), rna_bits=sorted(pb), rfam_positions=sorted(rpos), contacts_unmapped=offmap,
                         other_chain_contacts=dict(other), family=acc))
    return rows


def cluster(rows):
    """site (Rfam-position overlap, single linkage 0.4) and mode (ECIFP, average linkage 0.5) labels on the 'ok' rows; matrices"""
    ok = [r for r in rows if r["status"] == "ok"]
    n = len(ok)
    L = [set(r["lig_bits"]) for r in ok]; P = [set(r["rna_bits"]) for r in ok]; S = [set(r["rfam_positions"]) for r in ok]
    ecifp = np.array([[tani(L[i], L[j]) ** 0.4 * tani(P[i], P[j]) ** 0.6 for j in range(n)] for i in range(n)])
    lig_only = np.array([[tani(L[i], L[j]) for j in range(n)] for i in range(n)])
    ov = np.array([[len(S[i] & S[j]) / min(len(S[i]), len(S[j])) if S[i] and S[j] else 0 for j in range(n)] for i in range(n)])
    site = fcluster(linkage((1 - ov)[np.triu_indices(n, 1)], "single"), 0.4, "distance") if n > 1 else np.ones(n, int)
    mode = fcluster(linkage((1 - ecifp)[np.triu_indices(n, 1)], "average"), 0.5, "distance") if n > 1 else np.ones(n, int)
    for r, s, m in zip(ok, site, mode):
        r["site"], r["mode"] = int(s), int(m)
    return ok, ecifp, lig_only, ov


def run(acc):
    rows = collect(acc)
    ok, ecifp, lig_only, ov = cluster(rows)
    n, site = len(ok), [r["site"] for r in ok]
    json.dump(dict(acc=acc, params=dict(radius=R, size=SIZE, weights=[0.4, 0.6], site_cut=0.4, mode_cut=0.5),
                   rows=rows, labels=[f"{r['pdb']}_{r['chain']}_{r['num']}_{r['ccd']}" for r in ok],
                   ecifp=ecifp.round(3).tolist(), lig_only=lig_only.round(3).tolist(), site_overlap=ov.round(3).tolist()),
              open(OUT / f"{acc}.json", "w"))
    print(f"== {acc}: {len(rows)} ligand copies, {n} with family-RNA contacts; statuses {dict(collections.Counter(r['status'] for r in rows))}")
    for s in sorted(set(site)):
        mem = [r for r in ok if r["site"] == s]
        core = collections.Counter(p for r in mem for p in r["rfam_positions"])
        print(f" site {s}: {len(mem)} copies, ligands {dict(collections.Counter(r['ccd'] for r in mem))}; positions in >=50%: {sorted(p for p, k in core.items() if k >= len(mem) / 2)}")
        for m in sorted({r['mode'] for r in mem}):
            mm = [r for r in mem if r["mode"] == m]
            print(f"   mode {m}: {len(mm)} {dict(collections.Counter(r['ccd'] for r in mm))}")


if __name__ == "__main__":
    REPO = Path(sys.argv[1]); configure(sys.argv[2], sys.argv[3])
    os.chdir(REPO)
    for a in sys.argv[4:]:
        run(a)
