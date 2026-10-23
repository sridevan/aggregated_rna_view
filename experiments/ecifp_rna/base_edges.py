"""Ligand contacts by nucleotide edge (Leontis-Westhof): Watson-Crick (W), Hoogsteen (H), sugar (S), base face (stacking/pi:
Arpeggio ring-centroid contacts), ribose, phosphate. Atoms on two edges (e.g. A N6: W and H) count for both.
Input: the cached PDBe bound_ligand_interactions used by ecifp_rna.py and its out/<ACC>.json rows (site labels, Rfam positions).
usage: base_edges.py <cache_dir> <out_dir> ACC [...]"""
import json, sys, collections
from pathlib import Path

EDGE = {  # base atom -> edges
    "A": {"N6": "WH", "N1": "W", "C2": "WS", "N3": "S", "N7": "H", "C8": "H"},
    "G": {"O6": "WH", "N1": "W", "N2": "WS", "N3": "S", "N7": "H", "C8": "H"},
    "C": {"N4": "WH", "N3": "W", "O2": "WS", "C5": "H", "C6": "H"},
    "U": {"O4": "WH", "N3": "W", "O2": "WS", "C5": "H", "C6": "H"},
}
PARENT = {"PSU": "U", "5MU": "U", "H2U": "U", "5MC": "C", "OMC": "C", "1MA": "A", "2MG": "G", "7MG": "G", "M2G": "G", "OMG": "G", "OMU": "U"}
RIBOSE = {"C1'", "C2'", "C3'", "C4'", "O4'", "C5'"}
PHOS = {"P", "OP1", "OP2", "OP3", "O5'", "O3'"}
DROP = {"proximal", "clash", "vdw_clash", "vdw"}
HB = {"hbond", "polar", "ionic", "metal_complex"}   # contacts that define an edge (hydrogen-bond-type); weak CH contacts excluded
EDGE_ALL = {a: e for d in EDGE.values() for a, e in d.items() if len(e) == 1}
EDGE_ALL.update({"O2'": "S"})
NAME = {"W/H": "W/H corner", "W/S": "W/S corner", "W/": "W", "W": "Watson-Crick edge", "H": "Hoogsteen edge", "S": "sugar edge", "F": "base face", "R": "ribose", "P": "phosphate"}


def edges(resname, atoms, types):
    if len(atoms) > 1:                        # Arpeggio ring centroid: pi/stacking contact with the base face
        return {"F"}
    a = atoms[0]
    if a in RIBOSE:                           # backbone: any contact that passes the ECIFP filter
        return {"R"}
    if a in PHOS:
        return {"P"}
    if not types & HB:                        # base edges (and O2', sugar edge): hydrogen-bond-type contacts only
        return set()
    if a == "O2'":
        return {"S"}
    return set(EDGE.get(PARENT.get(resname, resname), {}).get(a, ""))


if __name__ == "__main__":
  cache, out = Path(sys.argv[1]), Path(sys.argv[2])
  for acc in sys.argv[3:]:
        d = json.load(open(out / f"{acc}.json"))
        ok = [r for r in d["rows"] if r["status"] == "ok"]
        inst = json.load(open(f"data/families/{acc}/instances.json"))
        pos = {}
        for i in inst:
            for r in json.load(open(f"data/families/{acc}/mappings/{i['seq_name']}.json"))["residues"]:
                pos[(i["pdb_id"], i["auth_asym_id"], str(r["author_residue_number"]), (r.get("author_insertion_code") or "").strip())] = r["rfam_consensus_pos"]
        per = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))  # site -> position -> edge -> copies
        nsite = collections.Counter(r["site"] for r in ok)
        for r in ok:
            it = json.load(open(cache / f"bli_{r['pdb']}_{r['chain']}_{r['num']}.json"))[r["pdb"]][0]["interactions"]
            seen = collections.defaultdict(set)
            for x in it:
                if not [t for t in x["interaction_details"] if t not in DROP]:
                    continue
                e = x["end"]
                p = pos.get((r["pdb"], e["chain_id"], str(e["author_residue_number"]), (e.get("author_insertion_code") or "").strip()))
                if p:
                    seen[p] |= {(k, e["atom_names"][0] if len(e["atom_names"]) == 1 else "ring") for k in edges(e["chem_comp_id"], e["atom_names"], set(x["interaction_details"]))}
            for p, ka in list(seen.items()):
                # an atom on two edges (W/H or W/S corner) counts for an edge only when no single-edge atom of the base was contacted
                base = {k for k, a in ka if k in "WHS"}
                single = {k for k, a in ka if k in "WHS" and len(EDGE_ALL.get(a, "")) == 1}
                es = {k for k, a in ka if k not in "WHS"} | (single if single else {"W/" + "".join(sorted(base - {"W"}))} if base else set())
                seen[p] = es
            for p, es in seen.items():
                for k in es:
                    per[r["site"]][p][k] += 1
                r.setdefault("edges", {})[str(p)] = sorted(es)
        json.dump(dict(acc=acc, edge_names=NAME, rows=ok), open(out / f"{acc}_edges.json", "w"))
        print(f"== {acc}")
        for s, n in nsite.most_common():
            if n < 3:
                continue
            print(f" site {s} ({n} copies; ligands {dict(collections.Counter(r['ccd'] for r in ok if r['site'] == s).most_common(4))}):")
            for p in sorted(per[s], key=lambda p: -sum(per[s][p].values()))[:10]:
                c = per[s][p]
                print(f"   pos {p}: " + ", ".join(f"{NAME.get(k, k)} {v}/{n}" for k, v in c.most_common()))
