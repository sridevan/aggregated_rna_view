"""SAM riboswitch classes: every PDB chain of each SAM/SAH riboswitch family placed on its own family's Rfam consensus positions.

Families and entries (2026-10-23): PDBe nucleic_mappings for RF00162, RF00634, RF01725, RF01727 (6HAG), RF01826, RF01057 (3NPQ);
the other entries have no PDBe Rfam mapping and were found by RCSB full-text search and confirmed by cmscan against the nine
SAM/SAH CMs (each chain's best hit; see REPORT.md, SAM classes). For those, the segment is the whole chain (cmalign local mode
still leaves flanks as local ends) and the instance records segment_source = "whole_chain_cmscan".

Per family the pipeline's own functions build the chain instances (molecules sequence, residue listing, ModelServer chain file for
coordinate availability), align the unique sequences to the family CM (cmalign) and write a family-like directory
out/sam/<ACC>/ (instances.json, mappings/, ligands.json with the RF00162 ligand classification) that ecifp_rna.collect reads.
RF00162 (SAM-I) is not rebuilt: data/families/RF00162 is used as is.

usage: sam_classes.py <repo> <scratch> <out_dir> <cpu>
"""
import json, sys, collections
from pathlib import Path

REPO, SCRATCH, OUT, CPU = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), int(sys.argv[4])
sys.path[:0] = [str(REPO)]
from rna_explorer.cache import Cache
from rna_explorer.config import find_cmalign
from rna_explorer.pdbe import get_nucleic_mappings, get_molecules, get_residue_listing, get_chain_cif
from rna_explorer.rfam import get_cm
from rna_explorer.sequences import chain_instance_from_pdbe, unique_sequences, write_fasta
from rna_explorer.mmcif import atom_site_observed
from rna_explorer.mapping import build_chain_mapping
from rna_explorer.pipeline import align_uniques

FAMILIES = {   # accession -> (class name, entries)
    "RF00634": ("SAM-IV", ["6ues", "6uet", "6wlq", "6wlr"]),
    "RF01725": ("SAM-I/IV", ["4l81", "4oqu"]),
    "RF00521": ("SAM-II", ["2qwy"]),
    "RF01767": ("SAM-III", ["3e5c", "3e5e", "3e5f", "6c27"]),
    "RF01826": ("SAM-V", ["6fz0"]),
    "RF02885": ("SAM-VI", ["6las", "6lau", "6lax", "6laz"]),
    "RF01727": ("SAM/SAH", ["6hag", "6yl5", "6ylb", "6ymi", "6ymj", "6ymk", "6yml", "6ymm"]),
    "RF01057": ("SAH", ["3npn", "3npq"]),
}


def whole_chain_segments(pdb, mols, listing):
    """one segment per RNA chain (whole chain), in the nucleic_mappings shape"""
    rna = {int(e["entity_id"]) for e in mols if "ribonucleotide" in e.get("molecule_type", "").lower() and "deoxy" not in e.get("molecule_type", "").lower()}
    segs = []
    for m in listing.get("molecules", []):
        if int(m["entity_id"]) not in rna:
            continue
        for c in m["chains"]:
            rs = sorted(c["residues"], key=lambda r: r["residue_number"])
            segs.append(dict(entity_id=m["entity_id"], chain_id=c["chain_id"], struct_asym_id=c["struct_asym_id"],
                             start=dict(residue_number=rs[0]["residue_number"], author_residue_number=rs[0]["author_residue_number"]),
                             end=dict(residue_number=rs[-1]["residue_number"], author_residue_number=rs[-1]["author_residue_number"])))
    return segs


def build(acc, name, entries, cache, cmalign):
    d = OUT / acc
    (d / "mappings").mkdir(parents=True, exist_ok=True)
    pdbe = get_nucleic_mappings(cache, acc)
    mols = get_molecules(cache, entries)
    insts, source = [], {}
    for pdb in entries:
        listing = get_residue_listing(cache, pdb)
        ents = {int(e["entity_id"]): e for e in mols.get(pdb, [])}
        segs, src = (pdbe[pdb], "pdbe_nucleic_mappings") if pdb in pdbe else (whole_chain_segments(pdb, mols.get(pdb, []), listing), "whole_chain_cmscan")
        for seg in segs:
            cif = get_chain_cif(cache, pdb, seg["struct_asym_id"])
            inst = chain_instance_from_pdbe(pdb, seg, ents[int(seg["entity_id"])], listing, atom_site_observed(cif, seg["struct_asym_id"], seg["entity_id"]), None)
            insts.append(inst); source[inst.seq_name] = src
    fasta = d / "unique_sequences.fasta"
    write_fasta(unique_sequences(insts), fasta)
    work = d / "work"; work.mkdir(exist_ok=True)
    aligned, cm, scores = align_uniques(cmalign, Path(get_cm(cache, acc).path), fasta, work, CPU)
    chains = []
    for inst in insts:
        m = build_chain_mapping(inst, aligned[inst.sequence_key], cm, acc)
        res = [dict(polymer_seq_position=r.polymer_seq_position, author_residue_number=r.author_residue_number,
                    author_insertion_code=r.author_insertion_code, residue_name=r.residue_name, alignment_state=r.alignment_state,
                    rfam_consensus_pos=r.rfam_consensus_pos, posterior_prob=r.posterior_prob_value) for r in m.rows]
        json.dump(dict(seq_name=inst.seq_name, aligned_to=acc, residues=res), open(d / "mappings" / f"{inst.seq_name}.json", "w"))
        sc = scores[inst.sequence_key]
        chains.append(dict(pdb_id=inst.pdb_id, auth_asym_id=inst.auth_asym_id, label_asym_id=inst.label_asym_id, seq_name=inst.seq_name,
                           segment=[inst.map_start, inst.map_end], segment_source=source[inst.seq_name], molecule_name=inst.molecule_name,
                           fit=dict(bits=sc["bit_score"], avg_pp=sc["avg_pp"], cm=[sc["cm_from"], sc["cm_to"]]),
                           states=dict(collections.Counter(r["alignment_state"] for r in res))))
    json.dump(chains, open(d / "instances.json", "w"), indent=1)
    lig = json.load(open(REPO / "data/families/RF00162/ligands.json"))
    json.dump(dict(classification=lig["classification"]), open(d / "ligands.json", "w"))
    json.dump(dict(acc=acc, name=name, entries=entries, clen=cm.clen if hasattr(cm, "clen") else None), open(d / "family.json", "w"))
    print(f"== {acc} {name}: {len(chains)} chains")
    for c in chains:
        print(f"   {c['seq_name']:<14} {c['segment_source']:<22} bits {c['fit']['bits']:>6} pp {c['fit']['avg_pp']} cm {c['fit']['cm']} {c['states']}")


if __name__ == "__main__":
    cache = Cache(SCRATCH / "http")
    cmalign = find_cmalign()
    for acc, (name, entries) in FAMILIES.items():
        build(acc, name, entries, cache, cmalign)
