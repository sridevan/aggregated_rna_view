"""Geometry of SS_cons pair observations that carry no FR3D base pair (followup/unpaired_causes.csv): Watson-Crick-edge N-N distance
(purine N1 / pyrimidine N3) and C1'-C1' distance from the PDBe updated mmCIF, streamed into memory (model 1; files are never stored).

Categories: wc_like (N-N <= 3.6 A and C1'-C1' 8.5-12.5 A: bases positioned as a Watson-Crick-edge pair that FR3D did not annotate),
distorted (N-N <= 5.0 A otherwise), apart (N-N > 5.0 A), no_atoms. Also measured for a reference sample: SS_cons pairs FR3D annotates cWW
in the same chains (calibration). Usage: python3 experiments/sscons_vs_3d/geometry.py [ACC ...]   (default: the five families except 16S)
"""
import csv, json, sys, urllib.request
from collections import defaultdict
from pathlib import Path
import gemmi

ROOT = Path(__file__).resolve().parents[2]; OUT = Path(__file__).parent / "followup"
FAMS = sys.argv[1:] or ["RF00167", "RF00162", "RF00059", "RF00005", "RF00011"]
sys.path.insert(0, str(Path(__file__).parent)); from analyse import ss_pairs  # noqa: E402


def load_structure(pdb):
    txt = urllib.request.urlopen(f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb}_updated.cif", timeout=300).read().decode()
    st = gemmi.make_structure_from_block(gemmi.cif.read_string(txt).sole_block()); st.setup_entities(); return st


def measure(a, b, pa, pb):
    def at(r, n):
        x = r.find_atom(n, "*"); return x.pos if x else None
    p, q = at(a, "N1" if pa in "AG" else "N3"), at(b, "N1" if pb in "AG" else "N3"); c1, c2 = at(a, "C1'"), at(b, "C1'")
    nn = round(p.dist(q), 2) if p and q else None; cc = round(c1.dist(c2), 2) if c1 and c2 else None
    cat = ("no_atoms" if nn is None or cc is None else "wc_like" if nn <= 3.6 and 8.5 <= cc <= 12.5 else "distorted" if nn <= 5.0 else "apart")
    return nn, cc, cat


rows = []
unp = [r for r in csv.DictReader(open(OUT / "unpaired_causes.csv")) if r["acc"] in FAMS]
by_chain = defaultdict(list)
for r in unp: by_chain[(r["acc"], r["seq_name"])].append(r)
for (acc, sn), rr in sorted(by_chain.items()):
    d = ROOT / f"data/families/{acc}"; meta = {x["seq_name"]: x for x in json.load(open(d / "instances.json"))}[sn]
    m = json.load(open(d / meta["mapping_file"])); bypos = {x["rfam_consensus_pos"]: x for x in m["residues"] if x["alignment_state"] == "match"}
    try:
        st = load_structure(meta["pdb_id"].lower())
    except Exception as e:  # noqa: BLE001
        print(sn, "failed", e, file=sys.stderr); continue
    res = {r.label_seq: r for ch in st[0] for r in ch if r.subchain == meta["label_asym_id"]}
    obs = json.load(open(d / f"interactions_by_instance/{sn}.json"))["observations"]
    cww = {tuple(sorted((o[4], o[5]))) for o in obs if o[0] == "basepair" and o[1] == "cWW" and isinstance(o[4], int) and isinstance(o[5], int)}
    fam = json.load(open(d / "family.json")); ref = [p for p in ss_pairs(fam["ss_cons"]) if p in cww][:6]
    for kind, items in (("unpaired", [(int(r["pos_1"]), int(r["pos_2"]), r) for r in rr]), ("reference_cWW", [(i, j, None) for i, j in ref])):
        for i, j, r in items:
            ri, rj = bypos.get(i), bypos.get(j)
            a, b = res.get(ri["polymer_seq_position"]) if ri else None, res.get(rj["polymer_seq_position"]) if rj else None
            if not a or not b: continue
            nn, cc, cat = measure(a, b, ri["nucleotide_parent"], rj["nucleotide_parent"])
            rows.append({"acc": acc, "seq_name": sn, "pdb": meta["pdb_id"].upper(), "method": meta["experimental_method"], "resolution": meta["resolution"],
                         "kind": kind, "pos_1": i, "pos_2": j, "res": f"{a.name}{a.seqid.num}-{b.name}{b.seqid.num}", "nts": ri["nucleotide_parent"] + rj["nucleotide_parent"],
                         "complementary": int((ri["nucleotide_parent"] + rj["nucleotide_parent"]) in {"AU", "UA", "GC", "CG", "GU", "UG"}),
                         "wc_N_N": nn, "c1_c1": cc, "geometry": cat, "cause": r["cause"] if r else ""})
    print(sn, len(rr), file=sys.stderr)
with open(OUT / "unpaired_geometry.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
