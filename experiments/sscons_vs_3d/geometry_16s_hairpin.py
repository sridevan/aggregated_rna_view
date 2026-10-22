"""Geometry of the 16S consensus hairpin pairs 75-86, 76-85, 77-84, 78-83 (E. coli G79-C90, A80-U89, A81-U88, G82-C87) and a reference
pair of the same helix region that FR3D annotates (73-88, A77-U92), in a sample of chains; same measure as geometry.py (PDBe updated
mmCIF streamed into memory, never stored). Output: followup/hairpin_16s_geometry.csv"""
import csv, json, sys, urllib.request
from pathlib import Path
import gemmi

ROOT = Path(__file__).resolve().parents[2]; OUT = Path(__file__).parent / "followup"
SAMPLE = ["7k00_1_A_A", "4ybb_1_AA_A", "8aye_7_A_G", "9tnl_2_A_B", "5afi_1_a_A", "4u27_1_AA_A", "4v9d_1_AA_A"]
PAIRS = [(73, 88), (75, 86), (76, 85), (77, 84), (78, 83)]


def at(r, n):
    x = r.find_atom(n, "*"); return x.pos if x else None


d = ROOT / "data/families/RF00177"; meta = {x["seq_name"]: x for x in json.load(open(d / "instances.json"))}
rows = []
for sn in SAMPLE:
    if sn not in meta: print("not in family", sn, file=sys.stderr); continue
    mt = meta[sn]; m = json.load(open(d / mt["mapping_file"])); bypos = {x["rfam_consensus_pos"]: x for x in m["residues"] if x["alignment_state"] == "match"}
    txt = urllib.request.urlopen(f"https://www.ebi.ac.uk/pdbe/entry-files/download/{mt['pdb_id'].lower()}_updated.cif", timeout=600).read().decode()
    st = gemmi.make_structure_from_block(gemmi.cif.read_string(txt).sole_block()); st.setup_entities(); del txt
    res = {r.label_seq: r for ch in st[0] for r in ch if r.subchain == mt["label_asym_id"]}
    obs = json.load(open(d / f"interactions_by_instance/{sn}.json"))["observations"]
    cls = {tuple(sorted((o[4], o[5]))): o[1] for o in obs if o[0] == "basepair" and isinstance(o[4], int) and isinstance(o[5], int)}
    for i, j in PAIRS:
        ri, rj = bypos.get(i), bypos.get(j); a = res.get(ri["polymer_seq_position"]) if ri else None; b = res.get(rj["polymer_seq_position"]) if rj else None
        if not a or not b: rows.append({"seq_name": sn, "pair": f"{i}-{j}", "res": "missing"}); continue
        pa, pb = ri["nucleotide_parent"], rj["nucleotide_parent"]
        p, q = at(a, "N1" if pa in "AG" else "N3"), at(b, "N1" if pb in "AG" else "N3"); c1, c2 = at(a, "C1'"), at(b, "C1'")
        rows.append({"seq_name": sn, "pdb": mt["pdb_id"].upper(), "method": mt["experimental_method"], "resolution": mt["resolution"], "pair": f"{i}-{j}",
                     "res": f"{a.name}{a.seqid.num}-{b.name}{b.seqid.num}", "fr3d": cls.get((i, j), "-"),
                     "wc_N_N": round(p.dist(q), 2) if p and q else None, "c1_c1": round(c1.dist(c2), 2) if c1 and c2 else None})
    print(sn, file=sys.stderr)
with open(OUT / "hairpin_16s_geometry.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["seq_name", "pdb", "method", "resolution", "pair", "res", "fr3d", "wc_N_N", "c1_c1"]); w.writeheader(); w.writerows(rows)
