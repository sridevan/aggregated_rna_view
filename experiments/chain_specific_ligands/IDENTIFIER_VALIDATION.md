# Identifier validation

Canonical id: `(pdb_id, label_asym_id, label_seq_id)`; consensus position from `data/families/RF00167/mappings/<seq_name>.json` via `rna_explorer.annotation_mapping.InstanceIndex` (read-only; no new correspondence pipeline). Script: `scripts/compare_contacts.py`; table: `CONTACT_COMPARISON.csv` (499 rows).

| Source | Identifier form | Resolution route | Rows | mapped | outside_mapped_region | unresolved |
|---|---|---|---|---|---|---|
| PDBe v2 binding_sites | polymer index (= label_seq_id), no chain | `resolve_seq_index` on the entity's chains; chain left **blank** (entity level) | 70 | 62 | 8 (9NBH/9NDD scaffold) | 0 |
| PDBe graph-api bound_molecule_interactions | author chain + author number + insertion code | `resolve_unit` (author chain → instance by `auth_asym_id`; author number + ins code → residue) | 70 | 62 | 8 | 0 |
| PDBe graph-api bound_ligand_interactions | author ids + atom names + distance | same | 30 | 30 | 0 | 0 |
| RCSB ligand_neighbors | `ligand_asym_id`, RNA `seq_id` (label), `auth_seq_id`, atom ids, model, alt ids | `resolve_seq_index` on the instance (label asym from the URL) | 281 | 223 | 58 (scaffold / MG on 9NBH, 9NDD) | 0 |
| local Arpeggio | author chain/number/ins code/atom, `label_comp_id` | `resolve_unit` | 48 | 48 | 0 | 0 |

Preserved per row: author chain, author number, insertion code (blank in all test entries), ligand instance (author chain + number, `struct_asym_id`, `bm_id`), interaction types, minimum distance, ligand/RNA atom names where the source gives them. Model id: RCSB gives `ligand_model_id` (always 1 here); PDBe interaction records carry none (PDBe processes model 1 of the preferred assembly); local Arpeggio ran on model 1 of the asymmetric unit. Alternate conformations: `ligand_monomers.alternate_conformers` = 0 for every instance; RCSB `alt_id` fields absent in these records. Symmetry/assembly: not exposed by any interaction record; PDBe's assembly choice is known only through `api/pdb/entry/summary` (`preferred: true`).

Ligand instance ↔ label asym: `ligand_monomers.struct_asym_id` equals the RCSB `asym_id` for all 59 ligand instances (e.g. 1Y26 ADE X/90 = G, 3RKF DX4 B/91 = O, 9NDD OXG A/201 = E), so PDBe (author) and RCSB (label) records join without ambiguity.

Consistency checks
- 3RKF v2 label_seq_id 7, 8, 33, 37, 60, 61 → chain B author 21, 22, 47, 51, 74, 75 = the Arpeggio bm1 partners exactly (v2 is the chain B set, despite `bestChainId = A`).
- 1Y26: v2 (7 residues) = graph-api bm = graph-api atom-level = consensus 25, 26, 51, 55, 56, 80, 81. RCSB adds consensus 79 (A73, 3.5 Å); local Arpeggio adds 79 and 54 (C50, 3.9 Å, proximal only).
- 9NDD chain A: PDBe reports A65 (consensus 52, polar + weak H-bond 3.49 Å) while RCSB reports A66 (consensus 53, 3.53 Å); both are within 5 Å but each source lists only one. Local Arpeggio was not run on 9NDD (cryo-EM scaffold; not needed for feasibility).
- No author/label numbering mismatch and no unmatched identifier in 499 rows. Residues on the 9NBH/9NDD scaffolds resolve to the instance but have no consensus position (`outside_mapped_region`), as in Milestone 3B.

Ambiguities to carry into the schema: (a) entity-level v2 rows have no chain by construction; (b) for multi-copy entities only the PDBe-processed copies are chain-resolved by PDBe, the rest only by RCSB proximity or local Arpeggio; (c) PDBe interactions may include symmetry-related partners that the asymmetric-unit sources cannot see (none observed in the four entries).
