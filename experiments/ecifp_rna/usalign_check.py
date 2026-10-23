"""Structure-based correspondence (US-align) vs CM-based correspondence for clan siblings, per main-family Rfam position.

For each (reference chain of the main family, query chain) pair:
1. both chains (model 1, polymer residues only) from the PDBe ModelServer, written as temporary PDB files (not kept);
2. US-align -mol RNA -outfmt 1 (sequence-independent structural alignment, C3' atoms) with its superposition matrix (-m);
3. alignment strings mapped back to residues (author number + insertion code), residue by residue;
4. per main-family position (reference residue -> Rfam position from our mapping): the query residue paired by US-align, its C3'
   distance after the US-align superposition, and the query residue the CM cross mapping puts at that position; agree / differ /
   one-sided.
Also runs US-align -mm 4 (multiple structure alignment) over each set and reports whether the multiple alignment pairs the same
residues with the reference as the pairwise runs.

usage: usalign_check.py <repo> <usalign_binary> <tmp_dir> <out_json>
"""
import json, sys, subprocess, tempfile, collections
from pathlib import Path
import numpy as np, requests, gemmi

REPO, USALIGN, TMP, OUT = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
TMP.mkdir(parents=True, exist_ok=True)
MS = "https://www.ebi.ac.uk/pdbe/model-server/v1/{pdb}/atoms?auth_asym_id={ch}&model_nums=1&encoding=cif"
SETS = {
    "RF00167": dict(ref=("1y26", "X"), pocket=[25, 26, 27, 28, 29, 51, 53, 54, 55, 56, 57, 79, 80, 81, 82],
                    queries=[("RF00167", "4fe5", "B"), ("RF01510", "3ski", "A"), ("RF01510", "3skl", "A")]),
    "RF00162": dict(ref=("2gis", "A"), pocket=[7, 8, 46, 47, 48, 49, 69, 70, 71, 72, 102, 103, 104, 105],
                    queries=[("RF00162", "4kqy", "A"), ("RF00634", "6uet", "A"), ("RF00634", "6ues", "A"),
                             ("RF01725", "4l81", "A"), ("RF01725", "4oqu", "A")]),
}


def chain_pdb(pdb, ch):
    p = TMP / f"{pdb}_{ch}.pdb"
    st = gemmi.make_structure_from_block(gemmi.cif.read_string(requests.get(MS.format(pdb=pdb, ch=ch), timeout=120).text).sole_block())
    st.remove_ligands_and_waters(); st.remove_hydrogens()
    for m in list(st)[1:]:
        del st[1]
    st.setup_entities()
    p.write_text(st.make_pdb_string())
    res = [(r.seqid.num, r.seqid.icode.strip(), r.name, np.array(r.find_atom("C3'", "*").pos.tolist()))
           for c in st[0] for r in c if r.find_atom("C3'", "*") and r.het_flag != "H"]   # US-align reads ATOM records only
    return p, res


def usalign(q, r, extra=()):
    mat = TMP / "matrix.txt"
    out = subprocess.run([USALIGN, str(q), str(r), "-mol", "RNA", "-outfmt", "1", "-m", str(mat), *extra], capture_output=True, text=True, check=True).stdout
    seqs = [l.strip() for l in out.splitlines() if l and not l.startswith(">") and not l.startswith("#") and not l.startswith("$")]
    head = [l for l in out.splitlines() if l.startswith(">")]
    a_q, a_r = seqs[0], seqs[1]
    m = [l.split() for l in mat.read_text().splitlines() if l.strip()[:1] in "012"]
    t = np.array([float(x[1]) for x in m[:3]]); U = np.array([[float(v) for v in x[2:5]] for x in m[:3]])
    return a_q, a_r, (t, U), head


def mapping(acc, sn_path):
    return {(str(x["author_residue_number"]), (x.get("author_insertion_code") or "").strip()): x["rfam_consensus_pos"]
            for x in json.load(open(sn_path))["residues"]}


def seq_name(acc, pdb, ch, main):
    d = REPO / "data/families" / acc if acc == main else REPO / "experiments/ecifp_rna/out/clan" / f"{acc}_on_{main}"
    inst = json.load(open(d / "instances.json"))
    sn = next(i["seq_name"] for i in inst if i["pdb_id"] == pdb and i["auth_asym_id"] == ch)
    return d / "mappings" / f"{sn}.json"


def pairs(a_q, a_r, rq, rr):
    """residue-index pairs from the two aligned strings; checks the residue count"""
    if len(a_q.replace("-", "")) != len(rq) or len(a_r.replace("-", "")) != len(rr):
        raise ValueError(f"alignment length mismatch: {len(a_q.replace('-', ''))}/{len(rq)}, {len(a_r.replace('-', ''))}/{len(rr)}")
    out, i, j = {}, 0, 0
    for x, y in zip(a_q, a_r):
        if x != "-" and y != "-":
            out[j] = i
        i += x != "-"; j += y != "-"
    return out   # reference residue index -> query residue index


report = {}
for main, s in SETS.items():
    rp, rr = chain_pdb(*s["ref"])
    ref_map = mapping(main, seq_name(main, *s["ref"], main))
    ref_pos = {k: ref_map.get((str(n), ic)) for k, (n, ic, _, _) in enumerate(rr)}
    rows, files = [], [rp]
    for acc, pdb, ch in s["queries"]:
        qp, rq = chain_pdb(pdb, ch); files.append(qp)
        a_q, a_r, (t, U), head = usalign(qp, rp)
        pr = pairs(a_q, a_r, rq, rr)
        cm = mapping(acc, seq_name(acc, pdb, ch, main))           # query residue -> main-family position (CM cross map)
        cm_by_pos = {p: k for k, p in cm.items() if p}
        per = []
        for k, p in ref_pos.items():
            if not p:
                continue
            j = pr.get(k)
            us = f"{rq[j][2]}{rq[j][0]}{rq[j][1]}" if j is not None else None
            d = float(np.linalg.norm(U @ rq[j][3] + t - rr[k][3])) if j is not None else None
            c = cm_by_pos.get(p)
            cmr = next((f"{r[2]}{r[0]}{r[1]}" for r in rq if (str(r[0]), r[1]) == c), None) if c else None
            per.append(dict(pos=p, ref=f"{rr[k][2]}{rr[k][0]}{rr[k][1]}", usalign=us, dist=round(d, 1) if d is not None else None, cm=cmr,
                            status="agree" if us and us == cmr else "differ" if us and cmr else "usalign_only" if us else "cm_only" if cmr else "neither"))
        tm = [h for h in head if "TM-score" in h]
        st = collections.Counter(x["status"] for x in per)
        pk = [x for x in per if x["pos"] in s["pocket"]]
        rows.append(dict(family=acc, pdb=pdb, chain=ch, header=head, counts=dict(st), pocket=pk,
                         differ=[x for x in per if x["status"] != "agree"], close_agree=sum(1 for x in per if x["status"] == "agree" and x["dist"] <= 3)))
        print(f"{main} <- {acc} {pdb} {ch}: {dict(st)} | {' '.join(h.split('TM-score=')[-1][:6] for h in tm)}")
        for x in pk:
            print(f"    pocket {x['pos']:>4} {x['ref']:>6}: US-align {x['usalign']} ({x['dist']} A) | CM {x['cm']} -> {x['status']}")
    # multiple structure alignment over the set: same pairs with the reference?
    lst = TMP / f"list_{main}.txt"; lst.write_text("\n".join(f.name for f in files) + "\n")
    msta = subprocess.run([USALIGN, "-dir", str(TMP) + "/", str(lst), "-mm", "4", "-mol", "RNA", "-outfmt", "1"], capture_output=True, text=True).stdout
    report[main] = dict(ref=s["ref"], pocket_positions=s["pocket"], pairwise=rows, msta_raw=msta)
    print(f"  MSTA ({len(files)} chains) output lines: {len(msta.splitlines())}")
json.dump(report, open(OUT, "w"), indent=1)
