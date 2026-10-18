# Arpeggio pilot on deposited archive coordinates (RF00167 + 70S ribosome 16S)

Date: 2026-10-15. Status: pilot only; `rna_explorer/`, `frontend/` and `data/families/` unchanged. Scripts: `scripts/pilot.py`, `scripts/symmetry_check.py` (run with the Arpeggio venv python; see `../ARPEGGIO_FEASIBILITY.md` for the environment: pdbe-arpeggio 1.4.4, openbabel-wheel 3.1.1.23, gemmi 0.7.5).

## Coordinate context (recorded in every output as `coordinate_context`)
- Source: wwPDB archive mmCIF, deposited coordinates, `https://files.rcsb.org/download/{PDB}.cif` (not the PDBe "updated" files), stored in `cif/`.
- Model 1, asymmetric unit exactly as deposited. **No biological assembly generated, no symmetry expansion.** Hydrogens are added by Arpeggio/OpenBabel at run time.
- Contact definitions: `any` = Arpeggio INTER contact between a ligand atom and an RNA atom (includes `proximal` pairs within the 5 Å interacting cutoff); `specific` = at least one non-proximal class (hbond, polar, weak_hbond, weak_polar, vdw, vdw_clash, ionic, metal_complex, aromatic, carbonyl, …).

## Workflow (per PDB entry)
1. Target RNA chains = Rfam-mapped instances (`instances.json` + `mappings/` through `InstanceIndex`, read-only). For 7K00 the RF00177 POC numbering CSV (`poc_rf00177/results/ecoli_7k00_numbering_to_consensus.csv`) is the mapping and label_seq_id comes from the deposited mmCIF.
2. Nonpolymer ligand instances read from the mmCIF with gemmi: residues of entities of type non-polymer, water excluded; modified nucleotides (polymer entities) and branched entities are not ligand instances. Classified with `rna_explorer/annotation_config.json` (ligand / metal_ion / additive).
3. One Arpeggio run per entry with every ligand instance in a selection file (`work/<pdb>/selection.txt`, author ids).
4. INTER contacts kept only when one side is a ligand instance and the other a residue of a target chain (water and ligand–ligand pairs dropped).
5. Atom level retained (`atom_contacts`: ligand instance, ligand atom, RNA chain/residue/atom, distance, classes); residue level (`residue_contacts`: n atom pairs, min distance, class union, `specific`).
6. RNA residues resolved to (pdb_id, label_asym_id, label_seq_id), consensus position, alignment state, PP; insertions / scaffold residues keep `rfam_consensus_pos: null` with a status.
7. Aggregation per RNA chain instance (`per_chain[seq_name].positions`) and separately per entry (`entry.positions`, distinct chains per position); family level `RF00167_chain_level_tracks.json` with counting unit *RNA chain instance* and explicit denominators (observable, evaluated, positive, not contacted, not evaluated).
8. Ligand classes kept separate everywhere (`by_class`, three tracks per position).

Outputs: `results/<pdb>.json` (66 + 7k00), `results/RF00167_residue_contacts.csv`, `results/RF00167_chain_level_tracks.json`, `results/symmetry_check_RF00167.{json,csv}`, `results/7k00.json`, `results/7k00_residue_contacts.csv`, logs in `results/*.log`, raw Arpeggio JSON in `work/<pdb>/`.

## RF00167 results (66 entries, 84 chain instances, all evaluated)

| Quantity | Value |
|---|---|
| Runtime (Arpeggio) | 67.7 s total, mean 1.0 s, max 12.7 s per entry |
| Ligand instances in deposited files | 622 (metal_ion 505, ligand 82, additive 35) |
| Atom-level ligand–RNA contacts | 42,283 |
| Residue-level contacts | 3,431 (mapped 3,313; outside mapped region 106; insert 12; unresolved 0) |
| Ligand-class residue contacts | 728 (specific 510), on 75 chain instances in 63 entries |

Chain-level ligand track vs the Milestone 3B entity-level track (pocket positions; `n_inst any / specific`, local entries, 3B entry-entities):

| Consensus | any | specific | entries (local) | 3B entry-entities | f chain (any) | f entity (3B) |
|---|---|---|---|---|---|---|
| 25 | 74 | 73 | 62 | 62 | 0.937 | 0.969 |
| 26 | 75 | 73 | 63 | 62 | 0.938 | 0.954 |
| 51 | 67 | 55 | 59 | 54 | 0.893 | 0.857 |
| 54 | 74 | 20 | 62 | 7 | 0.925 | 0.108 |
| 55 | 74 | 74 | 62 | 62 | 0.925 | 0.954 |
| 56 | 74 | 64 | 62 | 60 | 0.925 | 0.923 |
| 79 | 74 | 22 | 62 | 24 | 0.925 | 0.369 |
| 80 | 73 | 71 | 61 | 61 | 0.913 | 0.939 |
| 81 | 74 | 61 | 62 | 60 | 0.937 | 0.938 |

Reading: the `any` definition (5 Å proximity) systematically adds the stacking neighbours at consensus 54 and 79 to nearly every instance, while `specific` reproduces the PDBe/entity-level pocket. The chain-level track must therefore expose both counts and default to `specific`, or state the 5 Å definition explicitly.

Validation against PDBe graph-api (preferred assembly, Arpeggio with PDBe preprocessing): for the 544 ligand instances PDBe processed, **0 PDBe RNA partners are missing locally**; 1,502 local partners are absent from PDBe, of which 1,454 are proximal-only (metal 1,304, ligand 126, additive 24) and 48 specific (additive 21, metal 20, ligand 7). 78 ligand instances (12.5 %) are not processed by PDBe at all (outside the preferred assembly) and are covered only by this pilot.

## Contacts potentially omitted (no assembly, no symmetry expansion)
`scripts/symmetry_check.py` places every ligand atom in the crystal lattice (gemmi NeighborSearch with cell + space group) and reports target-chain RNA atoms of *symmetry images* within 5 Å; cryo-EM entries (9NBC, 9NBH, 9NDD) are not applicable.

| Quantity | Value |
|---|---|
| X-ray entries checked | 63 |
| Entries with ≥ 1 symmetry-image ligand–RNA contact | 48 |
| Symmetry-image residue contacts | 409 (metal_ion 396: NCO 347, MG 36, K 9, MN 4; additive 6: ACT, BR; ligand 7) |
| Ligand-class cases | 9IWF GTP A/101 ↔ B/21, B/49 (3.2 Å) and B/101 ↔ B/49; 9LKE GTP A/102 ↔ A/35 (3.1 Å), A/71; 9V4Y ANG A/102 ↔ A/1 (2.9 Å), A/35 (2.6 Å) |

Interpretation: the omissions are crystal-packing contacts of ions and additives, plus the terminal-nucleotide GTP artefacts and a second 8-aminoguanine copy (9V4Y A/102) bridging to a symmetry mate. No pocket ligand of RF00167 loses a contact: the PDBe comparison (which includes symmetry mates) found no RNA partner absent locally. For families with ligands at crystal interfaces the check must be run and the flagged contacts shown as "symmetry-mate contact, not in deposited asymmetric unit".

## 70S ribosome pilot: 7K00 (E. coli 70S, cryo-EM 2.0 Å), 16S rRNA chain A (RF00177)
Mapping: `poc_rf00177/results/ecoli_7k00_numbering_to_consensus.csv` (1,542 author positions, contiguous; label_seq_id from the deposited mmCIF; CLEN 1533). Deposited file `cif/7K00.cif` (17 MB, 8,587 HETATM records, one model).

| Quantity | Value |
|---|---|
| Nonpolymer ligand instances (whole 70S, water excluded) | 329: MG 309, SPD 16, ZN 2, SPM 1, PAR 1 (paromomycin) → metal_ion 311, additive 17, ligand 1 |
| Arpeggio runtime | **24,134 s (6 h 42 min)** for one run with all 329 instances selected; 81,104 raw contacts |
| 16S chain A residue-level contacts | 472 (all mapped to consensus; 0 insert/unmapped): MG 439, SPD 19, PAR 14 |
| Paromomycin (PAR) contacts | h44 decoding site: author 1405–1409 and 1488–1496 → consensus 1396–1400 and 1479–1487; 11 of 14 residues with specific interactions (min distances 2.5–3.1 Å) |

The biology is recovered exactly (A1408, G1491/G1494 and the 1405–1409 / 1488–1496 strands of helix 44), the mapping to RF00177 consensus positions is complete, and the per-chain / per-entry aggregation works unchanged for a 1,542-nt chain. The runtime, however, rules out the whole-structure, all-ligands protocol for ribosomes: pdbe-arpeggio's binding-site expansion and neighbour search over ~150k atoms is repeated for every selected residue (309 Mg²⁺ dominate). Options for production, in order of preference: (1) pre-filter ligand instances with a gemmi neighbour search (keep only instances with any atom within 5 Å of a target chain) and pass only those; (2) run Arpeggio on a sub-structure containing the target chain, its neighbouring chains and the ligands (gemmi can write it; identifiers unchanged); (3) run per ligand class, deferring Mg²⁺ (which are mostly captured adequately as "metal_ion" proximity by the RCSB neighbour records). The MG/SPD contacts are retained in `results/7k00.json` with class labels and can be aggregated once a runtime strategy is chosen.

## Recommendation
Deposited asymmetric-unit coordinates are a sound default for the family pipeline: small RNAs run in ~1 s per entry, every ligand copy is covered, identifiers map without loss, and the symmetry check quantifies what the choice omits (crystal-contact ions and additives; no pocket ligand in RF00167). For ribosomes a ligand pre-filter is mandatory before Arpeggio. Before the production pipeline changes: decide on `specific` vs `any` as the chain-level contact definition (the tables above show that `any` inflates stacking neighbours), keep the entity-level PDBe track as a separate layer, and surface symmetry-image contacts as a flagged category rather than silently dropping them.
