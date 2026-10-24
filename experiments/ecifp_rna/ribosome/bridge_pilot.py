"""Pilot: can SSU rRNA positions be bridged between the bacterial (RF00177) and eukaryotic (RF01960) Rfam families, so that
aminoglycoside contacts from bacterial, mitochondrial and eukaryotic cytosolic ribosomes land on one axis?

Small set (one SSU rRNA chain per entry, aminoglycoside complexes; 2026-10-24):
  bacterial  7K00 A (E. coli, paromomycin, 1.98 A), 2VQE A (T. thermophilus, paromomycin, 2.5 A)
  mito       9VXM AA (human 12S, neomycin, 2.9 A)
  cytosolic  9VXN S2 (human 18S, neomycin), 5NDV 2 (S. cerevisiae, paromomycin), 9FXO S1 (Leishmania major, paromomycin),
             9G6J A (Candida albicans, G418), 9SUL AA (Trypanosoma brucei, G418)

Steps
1. Own-family mapping with the pipeline's functions (whole-chain segment, cmalign to the family CM, build_chain_mapping):
   bacterial and mitochondrial chains to RF00177, cytosolic chains to RF01960.
2. Bridge A (sequence/CM): cytosolic chains also aligned to RF00177 (cross-family cmalign).
   Bridge B (structure): US-align -mol RNA (C3') of every chain onto 7K00 A; each residue takes the RF00177 position of the 7K00
   residue it is paired with (pairs further than 4 A after superposition are not used).
   Agreement of A and B per residue; mitochondrial and bacterial chains also checked (own CM vs US-align route).
3. Aminoglycoside copies (PDBe bound_ligand_interactions, non-vdW contacts; inventory cache) on the SSU chain: contacted residues
   placed on RF00177 positions through each route; the copy at the decoding site is the one contacting 1408/1492/1493 (E. coli
   numbering via 7K00). Base identity per species at the known determinant positions.

usage: bridge_pilot.py <repo> <scratch> <usalign> <inventory_cache> <out_dir> <cpu>
"""
import json, sys, subprocess, collections
from pathlib import Path
import numpy as np, gemmi

REPO, SCRATCH, USALIGN, INVCACHE, OUT, CPU = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4]), Path(sys.argv[5]), int(sys.argv[6])
sys.path[:0] = [str(REPO)]
from rna_explorer.cache import Cache
from rna_explorer.config import find_cmalign
from rna_explorer.pdbe import get_molecules, get_residue_listing, get_chain_cif
from rna_explorer.rfam import get_cm
from rna_explorer.sequences import chain_instance_from_pdbe, unique_sequences, write_fasta
from rna_explorer.mmcif import atom_site_observed
from rna_explorer.mapping import build_chain_mapping
from rna_explorer.pipeline import align_uniques

SET = [  # pdb, auth chain, group, species, own family
    ("7k00", "A", "bacterial", "E. coli", "RF00177"), ("2vqe", "A", "bacterial", "T. thermophilus", "RF00177"),
    ("9vxm", "AA", "mitochondrial", "H. sapiens 12S", "RF00177"),
    ("9vxn", "S2", "cytosolic", "H. sapiens 18S", "RF01960"), ("5ndv", "2", "cytosolic", "S. cerevisiae", "RF01960"),
    ("9fxo", "S1", "cytosolic", "L. major", "RF01960"), ("9g6j", "A", "cytosolic", "C. albicans", "RF01960"),
    ("9sul", "AA", "cytosolic", "T. brucei", "RF01960"),
]
REF = ("7k00", "A")
DETERMINANTS = [1405, 1406, 1407, 1408, 1409, 1410, 1490, 1491, 1492, 1493, 1494, 1495, 1496]   # E. coli numbering, h44 A site
AMINOGLYCOSIDES = {"PAR", "NMY", "GET", "84G", "LLL", "KAN", "TOY", "SIS", "AKN", "AM2", "LUJ"}
TMP = SCRATCH / "bridge_tmp"; TMP.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)


def instances(cache, entries):
    mols = get_molecules(cache, [p for p, *_ in entries])
    out = {}
    for pdb, ch, *_ in entries:
        listing = get_residue_listing(cache, pdb)
        ents = {int(e["entity_id"]): e for e in mols[pdb]}
        for m in listing["molecules"]:
            if "polyribonucleotide" not in str(ents[int(m["entity_id"])].get("molecule_type", "")):   # ligands share the author chain id
                continue
            for c in m["chains"]:
                if c["chain_id"] != ch:
                    continue
                rs = sorted(c["residues"], key=lambda r: r["residue_number"])
                seg = dict(entity_id=m["entity_id"], chain_id=ch, struct_asym_id=c["struct_asym_id"],
                           start=dict(residue_number=rs[0]["residue_number"], author_residue_number=rs[0]["author_residue_number"]),
                           end=dict(residue_number=rs[-1]["residue_number"], author_residue_number=rs[-1]["author_residue_number"]))
                cif = get_chain_cif(cache, pdb, seg["struct_asym_id"])
                out[(pdb, ch)] = (chain_instance_from_pdbe(pdb, seg, ents[int(m["entity_id"])], listing,
                                                           atom_site_observed(cif, seg["struct_asym_id"], seg["entity_id"]), None), cif)
    return out


def cm_map(cache, cmalign, insts, acc, tag):
    """residue (author number, icode) -> consensus position of acc, per chain; plus fit"""
    d = TMP / f"{tag}_{acc}"; d.mkdir(exist_ok=True)
    fasta = d / "u.fasta"; write_fasta(unique_sequences([i for i, _ in insts.values()]), fasta)
    aligned, cm, scores = align_uniques(cmalign, Path(get_cm(cache, acc).path), fasta, d, CPU)
    out = {}
    for key, (inst, _) in insts.items():
        m = build_chain_mapping(inst, aligned[inst.sequence_key], cm, acc)
        res = {(r.author_residue_number, (r.author_insertion_code or "").strip()): (r.rfam_consensus_pos, r.residue_name, r.alignment_state)
               for r in m.rows}
        sc = scores[inst.sequence_key]
        out[key] = dict(res=res, fit=dict(bits=sc["bit_score"], avg_pp=sc["avg_pp"], cm=[sc["cm_from"], sc["cm_to"]]),
                        states=dict(collections.Counter(r.alignment_state for r in m.rows)))
    return out


def chain_pdb(cif, pdb, ch):
    st = gemmi.read_structure(str(cif)); st.remove_ligands_and_waters(); st.remove_hydrogens()
    for m in list(st)[1:]:
        del st[1]
    st.setup_entities()
    for c in st[0]:
        c.name = "A"          # US-align / PDB format: one-character chain id
    p = TMP / f"{pdb}_{ch}.pdb"; p.write_text(st.make_pdb_string())
    res = [(r.seqid.num, r.seqid.icode.strip(), r.name, np.array(r.find_atom("C3'", "*").pos.tolist()))
           for c in st[0] for r in c if r.find_atom("C3'", "*") and r.het_flag != "H"]
    return p, res


def usalign(q, r):
    mat = TMP / "matrix.txt"
    out = subprocess.run([USALIGN, str(q), str(r), "-mol", "RNA", "-outfmt", "1", "-m", str(mat)], capture_output=True, text=True, check=True).stdout
    seqs = [l.strip() for l in out.splitlines() if l and l[0] not in ">#$"]
    head = [l for l in out.splitlines() if l.startswith(">")]
    m = [l.split() for l in mat.read_text().splitlines() if l.strip()[:1] in "012"]
    t = np.array([float(x[1]) for x in m[:3]]); U = np.array([[float(v) for v in x[2:5]] for x in m[:3]])
    return seqs[0], seqs[1], t, U, head


def pairs(a_q, a_r, nq, nr):
    if len(a_q.replace("-", "")) != nq or len(a_r.replace("-", "")) != nr:
        raise ValueError("alignment length mismatch")
    out, i, j = {}, 0, 0
    for x, y in zip(a_q, a_r):
        if x != "-" and y != "-":
            out[i] = j
        i += x != "-"; j += y != "-"
    return out   # query index -> reference index


def copies(pdb, ch):
    """aminoglycoside copies contacting this chain: list of (ccd, copy id, contacted author residues on ch)"""
    out = []
    for f in sorted(INVCACHE.glob(f"bli_{pdb}_*.json")):
        d = json.loads(f.read_text())
        if not isinstance(d, dict) or pdb not in d:
            continue
        x = d[pdb][0]
        ccd = x.get("ligand", {}).get("chem_comp_id") if isinstance(x.get("ligand"), dict) else None
        good = [i for i in x["interactions"] if set(i["interaction_details"]) - {"proximal", "vdw", "vdw_clash", "clash"}]
        ccd = ccd or (good[0]["begin"]["chem_comp_id"] if good else None)
        if ccd not in AMINOGLYCOSIDES:
            continue
        res = sorted({(i["end"]["author_residue_number"], (i["end"].get("author_insertion_code") or "").strip()) for i in good if i["end"]["chain_id"] == ch})
        if res:
            out.append(dict(ccd=ccd, copy=f.stem[4:], residues=res))
    return out


def main():
    cache = Cache(SCRATCH / "http"); cmalign = find_cmalign()
    insts = instances(cache, SET)
    own = {}
    for acc in ("RF00177", "RF01960"):
        sub = {k: v for k, v in insts.items() if next(s[4] for s in SET if s[:2] == k) == acc}
        own.update(cm_map(cache, cmalign, sub, acc, "own"))
    cross = cm_map(cache, cmalign, {k: v for k, v in insts.items() if next(s[4] for s in SET if s[:2] == k) == "RF01960"}, "RF00177", "cross")
    # reference: 7K00 residue -> RF00177 position, and E. coli number -> position
    ref_pos = {k: v[0] for k, v in own[REF]["res"].items()}
    ecoli_to_pos = {n: p for (n, ic), p in ref_pos.items() if not ic and p}
    pos_to_ecoli = {p: n for n, p in ecoli_to_pos.items()}
    rp, rr = chain_pdb(insts[REF][1], *REF)
    rows = []
    for pdb, ch, grp, sp, acc in SET:
        key = (pdb, ch)
        qp, rq = chain_pdb(insts[key][1], pdb, ch)
        a_q, a_r, t, U, head = usalign(qp, rp)
        pr = pairs(a_q, a_r, len(rq), len(rr))
        us = {}
        for i, j in pr.items():
            d = float(np.linalg.norm(U @ rq[i][3] + t - rr[j][3]))
            if d <= 4.0:
                us[(rq[i][0], rq[i][1])] = (ref_pos.get((rr[j][0], rr[j][1])), round(d, 1))
        cm = own[key]["res"] if acc == "RF00177" else cross[key]["res"]
        both = [(k, cm[k][0], us[k][0]) for k in cm if cm[k][0] and k in us and us[k][0]]
        agree = sum(1 for _, a, b in both if a == b)
        near = sum(1 for _, a, b in both if abs(a - b) <= 1)
        # determinant positions: residue and base on each route
        det = {}
        for n in DETERMINANTS:
            p = ecoli_to_pos.get(n)
            r_cm = next((f"{cm[k][1]}{k[0]}{k[1]}" for k in cm if cm[k][0] == p), None)
            r_us = next((f"{next(x[2] for x in rq if (x[0], x[1]) == k)}{k[0]}{k[1]}" for k in us if us[k][0] == p), None)
            det[n] = dict(pos=p, cm=r_cm, usalign=r_us)
        # aminoglycoside copies on this chain -> E. coli numbering through each route
        cps = []
        for c in copies(pdb, ch):
            via_cm = sorted({pos_to_ecoli.get(cm[k][0]) for k in c["residues"] if k in cm and cm[k][0] and pos_to_ecoli.get(cm[k][0])})
            via_us = sorted({pos_to_ecoli.get(us[k][0]) for k in c["residues"] if k in us and us[k][0] and pos_to_ecoli.get(us[k][0])})
            unmapped = [k for k in c["residues"] if not (k in cm and cm[k][0])]
            cps.append(dict(ccd=c["ccd"], copy=c["copy"], own=[f"{n}{ic}" for n, ic in c["residues"]], ecoli_via_cm=via_cm, ecoli_via_usalign=via_us,
                            unmapped_cm=[f"{n}{ic}" for n, ic in unmapped], decoding_site=bool({1408, 1492, 1493} & (set(via_cm) | set(via_us)))))
        tm = [h.split("TM-score=")[-1].split()[0] for h in head if "TM-score=" in h]
        row = dict(pdb=pdb, chain=ch, group=grp, species=sp, own_family=acc, own_fit=own[key]["fit"], own_states=own[key]["states"],
                   cross_fit=cross[key]["fit"] if key in cross else None, cross_states=cross[key]["states"] if key in cross else None,
                   usalign=dict(tm=tm, paired_within_4A=len(us), chain_residues=len(rq)),
                   route_agreement=dict(both=len(both), same=agree, within_1=near), determinants=det, copies=cps)
        rows.append(row)
        print(f"{pdb} {ch:<3} {sp:<16} own {acc} bits {own[key]['fit']['bits']} pp {own[key]['fit']['avg_pp']} | "
              f"cross {row['cross_fit']} | US TM {tm} paired {len(us)}/{len(rq)} | CM vs US same {agree}/{len(both)} (+-1 {near})")
        print("    A-site:", " ".join(f"{n}:{v['cm']}/{v['usalign']}" for n, v in det.items()))
        for c in cps:
            if c["decoding_site"]:
                print(f"    {c['ccd']} {c['copy']} decoding site -> E. coli via CM {c['ecoli_via_cm']} | via US {c['ecoli_via_usalign']} | unmapped {c['unmapped_cm']}")
        print(f"    copies on chain: {len(cps)}, at decoding site: {sum(c['decoding_site'] for c in cps)}")
    json.dump(dict(set=SET, reference=REF, determinants_ecoli=DETERMINANTS,
                   sources="PDBe residue_listing/molecules/ModelServer chains; Rfam CMs RF00177, RF01960 (cmalign via the pipeline); US-align; PDBe bound_ligand_interactions (inventory cache); retrieved 2026-10-24",
                   rows=rows), open(OUT / "bridge_pilot.json", "w"), indent=1, default=str)


main()
