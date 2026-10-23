"""ECIFP within an Rfam clan: sibling-family chains placed on the main family's consensus positions.

For each sibling family S of a main family M (same clan):
1. S's PDBe-mapped chains are built exactly as in the pipeline (PDBe nucleic_mappings segment, molecules sequence, residue listing,
   coordinates from the PDBe ModelServer chain file), in a scratch cache;
2. the unique segment sequences are aligned with cmalign twice: to M's CM (the cross mapping, giving M positions) and to S's own CM
   (for comparison of fit: bit score and posterior probability);
3. the cross mapping is written as a family-like directory (instances.json, mappings/, ligands.json with M's ligand classes) and the
   ligand copies of S and M go through ecifp_rna.collect + cluster together, so sites are defined on M's positions;
4. 3D check of the correspondence: each S chain is superposed on M's reference pocket (C1' C4' P of M positions within 12 A of the
   reference ligand) through the cross-mapped positions; a low RMSD supports the correspondence, a high one rejects it.

usage: clan_cross.py <repo> <scratch> <out_dir> <cpu> MAIN:SIB[,SIB] [...]
"""
import json, sys, collections
from pathlib import Path
import numpy as np, requests, gemmi

REPO, SCRATCH, OUT, CPU = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), int(sys.argv[4])
sys.path[:0] = [str(REPO), str(Path(__file__).parent)]
from rna_explorer.cache import Cache
from rna_explorer.config import find_cmalign
from rna_explorer.pdbe import get_nucleic_mappings, get_molecules, get_residue_listing, get_chain_cif
from rna_explorer.rfam import get_cm
from rna_explorer.sequences import chain_instance_from_pdbe, unique_sequences, write_fasta
from rna_explorer.mmcif import atom_site_observed
from rna_explorer.mapping import build_chain_mapping
from rna_explorer.pipeline import align_uniques
import ecifp_rna as E

REFS = {"RF00167": ("1y26", "X", 90), "RF00162": ("2gis", "A", 301)}
MS = "https://www.ebi.ac.uk/pdbe/model-server/v1/{pdb}/atoms?{q}&encoding=cif"


def cross_map(main, sib, cache, cmalign):
    d = OUT / f"{sib}_on_{main}"
    (d / "mappings").mkdir(parents=True, exist_ok=True)
    maps = get_nucleic_mappings(cache, sib)
    mols = get_molecules(cache, sorted(maps))
    insts = []
    for pdb, segs in sorted(maps.items()):
        listing = get_residue_listing(cache, pdb)
        ents = {int(e["entity_id"]): e for e in mols.get(pdb, [])}
        for seg in segs:
            cif = get_chain_cif(cache, pdb, seg["struct_asym_id"])
            insts.append(chain_instance_from_pdbe(pdb, seg, ents[int(seg["entity_id"])], listing, atom_site_observed(cif, seg["struct_asym_id"], seg["entity_id"]), None))
    fasta = d / "unique_sequences.fasta"
    write_fasta(unique_sequences(insts), fasta)
    fits = {}
    for acc in (main, sib):
        work = d / f"work_{acc}"; work.mkdir(exist_ok=True)
        fits[acc] = align_uniques(cmalign, Path(get_cm(cache, acc).path), fasta, work, CPU)
    aligned, cm, scores = fits[main]
    own = fits[sib][2]
    chains = []
    for inst in insts:
        m = build_chain_mapping(inst, aligned[inst.sequence_key], cm, main)
        res = [dict(polymer_seq_position=r.polymer_seq_position, author_residue_number=r.author_residue_number,
                    author_insertion_code=r.author_insertion_code, residue_name=r.residue_name, alignment_state=r.alignment_state,
                    rfam_consensus_pos=r.rfam_consensus_pos, posterior_prob=r.posterior_prob_value) for r in m.rows]
        json.dump(dict(seq_name=inst.seq_name, aligned_to=main, residues=res), open(d / "mappings" / f"{inst.seq_name}.json", "w"))
        sc, so = scores[inst.sequence_key], own[inst.sequence_key]
        chains.append(dict(pdb_id=inst.pdb_id, auth_asym_id=inst.auth_asym_id, label_asym_id=inst.label_asym_id, seq_name=inst.seq_name,
                           segment=[inst.map_start, inst.map_end], fit_main=dict(bits=sc["bit_score"], avg_pp=sc["avg_pp"], cm=[sc["cm_from"], sc["cm_to"]]),
                           fit_own=dict(bits=so["bit_score"], avg_pp=so["avg_pp"], cm=[so["cm_from"], so["cm_to"]]),
                           states=dict(collections.Counter(r["alignment_state"] for r in res))))
    json.dump(chains, open(d / "instances.json", "w"), indent=1)
    lig = json.load(open(REPO / "data/families" / main / "ligands.json"))
    json.dump(dict(classification=lig["classification"]), open(d / "ligands.json", "w"))
    return d, chains


def model(pdb, q):
    r = requests.get(MS.format(pdb=pdb, q=q), timeout=120); r.raise_for_status()
    return gemmi.make_structure_from_block(gemmi.cif.read_string(r.text).sole_block())[0]


def backbone(pdb, chain, mapping):
    pos = {(str(x["author_residue_number"]), (x.get("author_insertion_code") or "").strip()): (x["rfam_consensus_pos"], x.get("posterior_prob"))
           for x in mapping["residues"] if x["rfam_consensus_pos"]}
    out = {}
    for ch in model(pdb, f"auth_asym_id={chain}&model_nums=1"):
        for rr in ch:
            p = pos.get((str(rr.seqid.num), rr.seqid.icode.strip()))
            if p:
                for a in rr:
                    if a.name in ("C1'", "C4'", "P") and a.altloc in ("\0", "A"):
                        out[(p[0], a.name)] = (np.array(a.pos.tolist()), p[1])
    return out


def rmsd_fit(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    X = (P - pc) @ R.T + qc
    return float(np.sqrt(((X - Q) ** 2).sum(1).mean())), np.sqrt(((X - Q) ** 2).sum(1))


def pocket_check(main, sibdir, chains):
    pdb, ch, num = REFS[main]
    inst = {(i["pdb_id"], i["auth_asym_id"]): i["seq_name"] for i in json.load(open(REPO / "data/families" / main / "instances.json"))}
    ref = backbone(pdb, ch, json.load(open(REPO / "data/families" / main / "mappings" / f"{inst[(pdb, ch)]}.json")))
    lig = np.array([a.pos.tolist() for c in model(pdb, f"auth_asym_id={ch}&auth_seq_id={num}&model_nums=1") for r in c for a in r]).mean(0)
    near = sorted({p for (p, a), (x, _) in ref.items() if a == "C4'" and np.linalg.norm(x - lig) <= 12})
    out = []
    for c in chains:
        b = backbone(c["pdb_id"], c["auth_asym_id"], json.load(open(sibdir / "mappings" / f"{c['seq_name']}.json")))
        keys = [k for k in ref if k[0] in near and k in b]
        pp = [b[k][1] for k in keys if k[1] == "C4'" and b[k][1] is not None]
        fit = rmsd_fit(np.array([b[k][0] for k in keys]), np.array([ref[k][0] for k in keys])) if len(keys) >= 9 else (None, [])
        per = {f"{k[0]}{'' if k[1] == 'C4' + chr(39) else ':' + k[1]}": round(float(x), 1) for k, x in zip(keys, fit[1]) if k[1] == "C4'"}
        out.append(dict(seq_name=c["seq_name"], pocket_positions=len(near), positions_mapped=len({k[0] for k in keys}), atoms=len(keys),
                        rmsd=round(fit[0], 2) if fit[0] is not None else None, c4_deviation=per,
                        median_pp=round(float(np.median(pp)), 3) if pp else None))
    return near, out


def main():
    cache = Cache(SCRATCH / "http")
    cmalign = find_cmalign()
    E.configure(SCRATCH / "ecifp_cache", OUT)
    import os; os.chdir(REPO)
    for spec in sys.argv[5:]:
        main_acc, sibs = spec.split(":")
        rows = E.collect(main_acc)
        report = dict(main=main_acc, siblings={})
        for sib in sibs.split(","):
            d, chains = cross_map(main_acc, sib, cache, cmalign)
            near, pocket = pocket_check(main_acc, d, chains)
            srows = E.collect(sib, d)
            rows += srows
            report["siblings"][sib] = dict(chains=chains, pocket_positions=near, pocket=pocket,
                                           ligand_statuses=dict(collections.Counter(r["status"] for r in srows)))
        ok, ecifp, lig_only, ov = E.cluster(rows)
        report.update(rows=rows, labels=[f"{r['family']}_{r['pdb']}_{r['chain']}_{r['num']}_{r['ccd']}" for r in ok],
                      ecifp=ecifp.round(3).tolist(), lig_only=lig_only.round(3).tolist(), site_overlap=ov.round(3).tolist())
        json.dump(report, open(OUT / f"clan_{main_acc}.json", "w"))
        print(f"=== {main_acc} + {sibs}")
        for sib, s in report["siblings"].items():
            print(f" {sib}: ligand copies {s['ligand_statuses']}; pocket positions {s['pocket_positions']}")
            for c, p in zip(s["chains"], s["pocket"]):
                print(f"   {c['seq_name']}: on {main_acc} {c['fit_main']} | own {c['fit_own']} | states {c['states']} | pocket {p}")
        for st in sorted({r['site'] for r in ok}):
            mem = [r for r in ok if r["site"] == st]
            if not any(r["family"] != main_acc for r in mem):
                continue
            print(f" site {st}: {dict(collections.Counter((r['family'], r['ccd']) for r in mem))}")
            for m in sorted({r['mode'] for r in mem}):
                print(f"   mode {m}: {dict(collections.Counter((r['family'], r['ccd']) for r in mem if r['mode'] == m))}")
        idx = {i: r for i, r in enumerate(ok)}
        for sib in report["siblings"]:
            S = [i for i, r in idx.items() if r["family"] == sib]
            M = [i for i, r in idx.items() if r["family"] == main_acc]
            for i in S:
                best = max(M, key=lambda j: ecifp[i, j])
                print(f"  {idx[i]['pdb']} {idx[i]['ccd']} site {idx[i]['site']}: best {main_acc} match {idx[best]['pdb']} {idx[best]['ccd']} "
                      f"ECIFP {ecifp[i, best]:.2f} (lig {lig_only[i, best]:.2f}, site overlap {ov[i, best]:.2f}); positions {idx[i]['rfam_positions']}")


main()
