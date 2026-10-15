"""Step 3: build cmalign input (PDBe-mapped regions, parent-base normalised) + chain metadata. Reuses the RF00167 design."""
import json, sys, hashlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent)); from common import entity_source
DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(exist_ok=True)
FAM = "RF00177"
# (pdb_id, auth chain) and rationale
SELECTION = [
 ("7k00","A",  "Conventional bacterial 16S: E. coli K-12 70S, 2.0 A cryo-EM; E. coli numbering is the community reference; carries modelled modified nucleotides"),
 ("4y4o","1a", "Conventional bacterial 16S: Thermus thermophilus HB8 70S X-ray 2.3 A; second well-characterised bacterium, 2 copies in ASU (only 1a taken)"),
 ("6eri","BA", "Chloroplast 16S: spinach chloroplast ribosome 3.0 A; RF00177 mapping confirmed (Rfam hit 1-1523, bits 1517)"),
 ("7p2e","A",  "Mitochondrial 12S: human mitoribosome SSU 2.x A; explicitly mapped to RF00177 (Rfam hit model 9-1522, bits 251); shortest full-length instance class (946 nt) -> many deletions"),
 ("5mrc","aa", "Mitochondrial 15S: S. cerevisiae mitoribosome 3.25 A; 1627 nt mapped, low bit score (260) -> divergent with large insertions"),
 ("6yw5","aa", "Mitochondrial 16S: Neurospora crassa mitoribosome SSU 2.85 A; longest mapped region in RF00177 (1808 nt) -> largest insertions"),
 ("9gyt","2",  "Mitochondrial 18S: Brassica oleracea mitoribosome 2.11 A; 1760 nt mapped yet high bit score (1232) -> insertions in a confidently aligned background; 6 modified nucleotides"),
 ("8apn","B3", "Fragmented mitochondrial SSU rRNA: Polytomella magna mtSSU-3 (379 nt) maps to model 488-980 only -> partial model coverage, terminal insertions (anchor 0 test)"),
]
if len(sys.argv) > 3:   # optional: JSON list of [pdb_id, auth_chain, rationale] overrides the default selection
    SELECTION = [tuple(x) for x in json.load(open(sys.argv[3]))]
maps = json.load(open(DATA/f"api/{FAM.lower()}_nucleic_mappings.json"))[FAM]["PDB"]
chains, fa_mapped, fa_full = [], [], []
for pdb, chain, why in SELECTION:
    seg = [s for s in maps[pdb] if s["chain_id"] == chain]; assert len(seg) == 1, (pdb, chain); seg = seg[0]
    ent = seg["entity_id"]
    mol = [e for e in json.load(open(DATA/f"api/molecules_{pdb}.json"))[pdb] if e["entity_id"] == ent][0]
    rl = [c for m in json.load(open(DATA/f"api/residues_{pdb}.json"))[pdb]["molecules"] if m["entity_id"] == ent
          for c in m["chains"] if c["chain_id"] == chain and c["struct_asym_id"] == seg["struct_asym_id"]]
    assert len(rl) == 1, (pdb, chain, len(rl))
    residues = sorted(rl[0]["residues"], key=lambda r: r["residue_number"])
    seq = mol["sequence"]  # PDBe one-letter sequence, modified residues already given as parent base
    assert len(seq) == len(residues) == mol["length"], (pdb, len(seq), len(residues), mol["length"])
    s, e = seg["start"]["residue_number"], seg["end"]["residue_number"]
    modpos = mol.get("pdb_sequence_indices_with_multiple_residues") or {}
    mod_in_region = {int(k): v for k, v in modpos.items() if s <= int(k) <= e}
    nonstd = sorted({seq[i] for i in range(s-1, e)} - set("ACGU"))
    src = (mol.get("source") or [{}])[0]
    name = f"{pdb}_{ent}_{chain}_{seg['struct_asym_id']}"
    mapped = seq[s-1:e]
    chains.append({"pdb_id": pdb, "entity_id": ent, "auth_chain_id": chain, "label_asym_id": seg["struct_asym_id"], "seq_name": name, "rationale": why,
                   "sequence": seq, "pdb_sequence": mol["pdb_sequence"], "modified_positions": modpos, "modified_positions_in_mapped_region": mod_in_region,
                   "non_ACGU_letters_in_mapped_region": nonstd, "map_start": s, "map_end": e,
                   "map_auth_start": seg["start"]["author_residue_number"], "map_auth_end": seg["end"]["author_residue_number"],
                   "mapped_sequence": mapped, "mapped_sequence_md5": hashlib.md5(mapped.encode()).hexdigest(),
                   "molecule_name": mol["molecule_name"], "api_organism": src.get("organism_scientific_name"), "api_tax_id": src.get("tax_id"),
                   **entity_source(DATA/f"cif/{pdb}_updated.cif", ent), "residues": residues})
    fa_mapped.append(f">{name}/{s}-{e}\n{mapped}\n"); fa_full.append(f">{name}/1-{len(seq)}\n{seq}\n")
(OUT/"rf00177_sequences.fasta").write_text("".join(fa_mapped)); (OUT/"rf00177_full_sequences.fasta").write_text("".join(fa_full))
json.dump(chains, open(OUT/"chains.json", "w"), indent=1)
for c in chains:
    print(f"{c['seq_name']:18s} len={len(c['sequence']):5d} mapped={c['map_start']}-{c['map_end']} ({len(c['mapped_sequence'])}) auth={c['map_auth_start']}-{c['map_auth_end']} mod_in_region={len(c['modified_positions_in_mapped_region'])} nonACGU={c['non_ACGU_letters_in_mapped_region']} | API: {c['api_organism']} | CIF: {c['cif_organism']} ({c['cif_tax_id']}) [{c['cif_src_category']}]")
    if c["modified_positions_in_mapped_region"]:
        print("    modified:", {k: v["three_letter_code"] + "->" + v["one_letter_code"] for k, v in list(c["modified_positions_in_mapped_region"].items())[:12]})
