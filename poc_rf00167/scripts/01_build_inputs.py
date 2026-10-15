"""Build cmalign input FASTA + chain metadata from PDBe API JSON (POC, RF00167)."""
import json, sys, re
from pathlib import Path
DATA = Path(sys.argv[1]); OUT = Path(sys.argv[2]); OUT.mkdir(exist_ok=True)
SELECTION = [("6p2h","A"),("1y26","X"),("3ds7","B"),("4fe5","B"),("9nbh","A")]  # (pdb_id, author chain_id)

maps = json.load(open(DATA/"rf00167_nucleic_mappings.json"))["RF00167"]["PDB"]

def cif_source(cif_path, entity_id):
    """Minimal mmCIF scan for organism of an entity (non-loop, key-value form only)."""
    txt = open(cif_path).read()
    out = {}
    for cat, orgkey, taxkey in [("_entity_src_gen","pdbx_gene_src_scientific_name","pdbx_gene_src_ncbi_taxonomy_id"),
                                ("_entity_src_nat","pdbx_organism_scientific","pdbx_ncbi_taxonomy_id"),
                                ("_pdbx_entity_src_syn","organism_scientific","ncbi_taxonomy_id")]:
        ent = re.search(rf"^{cat}\.entity_id\s+(\S+)", txt, re.M)
        if ent and ent.group(1) == str(entity_id):
            org = re.search(rf"^{cat}\.{orgkey}\s+(.+)$", txt, re.M)
            tax = re.search(rf"^{cat}\.{taxkey}\s+(\S+)", txt, re.M)
            out = {"cif_src_category": cat,
                   "cif_organism": org.group(1).strip().strip("'\"") if org else None,
                   "cif_tax_id": tax.group(1) if tax else None}
            break
    return out

chains = []
fa_mapped, fa_full = [], []
for pdb, chain in SELECTION:
    seg = [s for s in maps[pdb] if s["chain_id"] == chain]
    assert len(seg) == 1, (pdb, chain, seg)
    seg = seg[0]
    ent = seg["entity_id"]
    mol = [e for e in json.load(open(DATA/f"molecules_{pdb}.json"))[pdb] if e["entity_id"] == ent][0]
    res = json.load(open(DATA/f"residues_{pdb}.json"))[pdb]["molecules"]
    rl = [c for m in res if m["entity_id"] == ent for c in m["chains"] if c["chain_id"] == chain and c["struct_asym_id"] == seg["struct_asym_id"]]
    assert len(rl) == 1
    residues = sorted(rl[0]["residues"], key=lambda r: r["residue_number"])
    seq = mol["sequence"]
    assert len(seq) == len(residues) == mol["length"]
    s, e = seg["start"]["residue_number"], seg["end"]["residue_number"]
    mapped = seq[s-1:e]
    src = mol["source"][0] if mol["source"] else {}
    name = f"{pdb}_{ent}_{chain}_{seg['struct_asym_id']}"
    chains.append({"pdb_id": pdb, "entity_id": ent, "auth_chain_id": chain, "label_asym_id": seg["struct_asym_id"],
                   "seq_name": name, "sequence": seq, "pdb_sequence": mol["pdb_sequence"],
                   "modified_positions": mol["pdb_sequence_indices_with_multiple_residues"],
                   "map_start": s, "map_end": e, "map_auth_start": seg["start"]["author_residue_number"],
                   "map_auth_end": seg["end"]["author_residue_number"], "mapped_sequence": mapped,
                   "molecule_name": mol["molecule_name"], "sample_preparation": mol.get("sample_preparation"),
                   "api_organism": src.get("organism_scientific_name"), "api_tax_id": src.get("tax_id"),
                   "api_expression_host": src.get("expression_host_scientific_name"),
                   **cif_source(DATA/f"{pdb}_updated.cif", ent),
                   "residues": residues})
    fa_mapped.append(f">{name}/{s}-{e}\n{mapped}\n")
    fa_full.append(f">{name}/1-{len(seq)}\n{seq}\n")
(OUT/"mapped_regions.fa").write_text("".join(fa_mapped))
(OUT/"full_sequences.fa").write_text("".join(fa_full))
json.dump(chains, open(OUT/"chains.json","w"), indent=1)
for c in chains:
    print(f"{c['seq_name']:16s} len={len(c['sequence']):3d} mapped={c['map_start']}-{c['map_end']} ({len(c['mapped_sequence'])} nt) auth={c['map_auth_start']}-{c['map_auth_end']} | API: {c['api_organism']} ({c['api_tax_id']}) | CIF: {c.get('cif_organism')} ({c.get('cif_tax_id')}) [{c.get('cif_src_category')}]")
