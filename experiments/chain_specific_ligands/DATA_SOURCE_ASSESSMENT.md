# Data source assessment: chain-specific RNA–ligand interactions

Retrieval date: 2026-10-15. All responses in `cache/` (names = URL with `/`, `?`, `=` replaced by `_`); `scripts/fetch_sources.sh` reproduces them.

## 1. PDBe

| Endpoint | Status | Content | Identifiers | Chain-specific? |
|---|---|---|---|---|
| `api/v2/pdb/entry/binding_sites/{pdb}/{entity}` (Milestone 3B) | 200 | per CCD id: residue index ranges with `boundMoleculeId`, `bestChainId`, scaffold | polymer index = label_seq_id; **no RNA chain**; bm ids | **No** (entity level). The residues are the union of the bound molecules processed by PDBe (see §1.3). `bestChainId` is unreliable: 3RKF reports `A` while every underlying bm (bm1–bm6) sits on chain B |
| `graph-api/pdb/bound_molecules/{pdb}` | 200 | bound molecules (`bm_id`) with composition: CCD, **author chain**, author residue number, insertion code, entity | author ids only | lists ligand instances, one bm per instance (connections for multi-component ligands) |
| `api/pdb/entry/ligand_monomers/{pdb}` | 200 | every ligand instance in the entry: author chain + number, `struct_asym_id` (= label_asym_id), `entity_id`, `residue_number`, `alternate_conformers`, `bm_id` (**null** when PDBe did not process the instance), `chem_comp_name`, weight | author **and** label ids | the key for joining bm ids to label asym ids |
| `graph-api/pdb/bound_molecule_interactions/{pdb}/{bm_id}` | 200 | Arpeggio interactions per bound molecule, **residue level**: `begin` (ligand residue), `end` (partner residue incl. waters), classes `atom_atom` (hbond, polar, weak_hbond, weak_polar, vdw, vdw_clash, aromatic, metal_complex, …), `atom_plane` (CARBONPI, CATIONPI, DONORPI), `plane_plane` (FF, OF, EE), `group_*`; `atom_names: null`; **no distances** | author chain + number + insertion code | **Yes** |
| `graph-api/pdb/bound_ligand_interactions/{pdb}/{chain}/{resnum}` | 200 for processed ligands, 404 otherwise | Arpeggio interactions per ligand residue, **atom level**: `ligand_atoms`, partner `atom_names`, `interaction_type` (atom-atom / atom-plane / plane-plane), `interaction_details`, **`distance`** | author ids | **Yes** |
| `graph-api/pdb/bound_ligand_interactions/{pdb}/{bm_id}` | 404 | – | – | the documented path takes chain/resnum, not bm id |
| `api/pdb/entry/binding_sites/{pdb}` (v1), `api/v2/pdb/entry/bound_molecules/{pdb}`, `graph-api/pdb/ligand_sites/{pdb}`, `graph-api/residue_mapping/{pdb}`, `graph-api/pdb/bound_excluded_molecules/{pdb}`, `graph-api/pdb/ligand_monomers/{pdb}` | 404 | not operational / do not exist | | |
| `api/pdb/entry/summary/{pdb}`, `api/pdb/entry/assembly/{pdb}` | 200 | assemblies with `preferred` flag, chains per assembly | author chain ids | explains the coverage limit below |

Documentation: `cache/pdbe_graph_api_openapi.json` (OpenAPI, 568 KB) lists `/pdb/bound_molecules/{pdb_id}`, `/pdb/bound_molecule_interactions/{pdb_id}/{bm_id}`, `/pdb/bound_ligand_interactions/{pdb_id}/{chain_id}/{seq_id}`, `/compound/interaction/{hetcode}`, `/pdb/entry/binding_sites/{pdb_id}/{entity_id}`, `/pdb/entry/ligand_monomers/{pdb_id}`. No endpoint exposes model, symmetry operator or assembly id in the interaction records; alternate conformers appear only as a count in `ligand_monomers`.

### 1.3 Coverage limitation: PDBe processes bound molecules of the preferred assembly only
3RKF has four author-defined monomer assemblies (one chain each); PDBe marks assembly 2 (chain B) as preferred. `bound_molecules` returns only bm1–bm6 = DX4 B/91 and five NCO on chain B; the 24 ligand instances on chains A, C, D have `bm_id: null` in `ligand_monomers`, no `bound_molecule_interactions` record, and `bound_ligand_interactions/3rkf/A/91` returns 404. The v2 `binding_sites` residues for 3RKF DX4 are exactly the chain B Arpeggio partners (author 21, 22, 47, 51, 74, 75 = label_seq_id 7, 8, 33, 37, 60, 61). 9NDD (one assembly with four chains) has all four OXG instances processed; 9NBH all 18 MG. Consequently the 3B record `chain_attribution: unresolved` for 3RKF is resolvable (chain B) from `ligand_monomers`/`bound_molecules`, but the other copies are simply absent from PDBe.

## 2. RCSB

| Resource | Status | Content | Identifiers | Nature |
|---|---|---|---|---|
| Data API `core/polymer_entity_instance/{PDB}/{asym}` → `rcsb_ligand_neighbors[]` | 200 | per RNA instance: ligand atom ↔ RNA atom pairs with `distance`, `ligand_asym_id`, `ligand_comp_id`, `ligand_entity_id`, `ligand_model_id`, `ligand_alt_id`, `ligand_is_bound` (covalent/metal flag), RNA `seq_id` (= label_seq_id), `auth_seq_id`, `comp_id`, `atom_id`, `alt_id` | **label ids** (asym_id, seq_id) + author number | **distance-based proximity**; empirically one record per (ligand instance, RNA residue) giving the closest atom pair; 281 records across 12 RNA instances, max distance 4.849 Å → 5 Å cutoff (schema has no description text) |
| Data API `core/nonpolymer_entity_instance/{PDB}/{asym}` → `rcsb_target_neighbors[]` | 200 | the same pairs from the ligand side (`target_asym_id`, `target_seq_id`, `target_atom_id`, `distance`, `target_is_bound`, `target_model_id`) | label ids | same data, ligand-centric; covers **every** ligand instance (all four 3RKF DX4 copies) |
| Data API `core/entry/{PDB}` | 200 | entity/instance inventories (`non_polymer_entity_ids`, `assembly_ids`, `model_ids`) | | inventory only |
| ModelServer `models.rcsb.org/v1/{pdb}/ligand?label_comp_id=ADE&radius=5&encoding=cif` | 200 | mmCIF of the ligand plus atoms within the radius (1Y26: 10 atom records) | label ids in the CIF | raw coordinates for a residue surrounding query; no interaction annotation |

No RCSB endpoint provides interaction types; `ligand_is_bound` is the only classification (5 rows here, all MG/NCO coordination).

## 3. Three kinds of evidence
1. **Explicit precomputed interactions** (PDBe graph-api, Arpeggio): geometric criteria per interaction class (hydrogen bonds, π-stacking, cation-π, vdW, clashes), hydrogen placement by PDBe's pipeline, computed on the preferred assembly (symmetry mates included). Chain- and atom-resolved, but only for ligand instances in the preferred assembly.
2. **Distance-based proximity** (RCSB neighbours, local Arpeggio "proximal"): every residue with an atom within a cutoff (5 Å); no interaction semantics; includes stacking neighbours and packing contacts that explicit classes may omit (e.g. 1Y26 A73 N6···N6 3.5 Å, consensus 79, reported by RCSB and local Arpeggio but not by PDBe). Computed on the asymmetric unit; covers all ligand instances.
3. **Entity-level binding-site annotation** (PDBe v2 binding_sites): union of the PDBe-processed bound molecules' residues, reported on the polymer sequence without chain. Correct for single-copy entities, resolvable through `ligand_monomers` bm ids for multi-copy entities, but blind to copies outside the preferred assembly.

## 4. Test entries

| Entry | RNA chains (label=auth) | Ligand instances (`ligand_monomers`) | PDBe processed (bm) | RCSB neighbours | Notes |
|---|---|---|---|---|---|
| 1Y26 | A=X | ADE X/90, MG ×5 | all 6 | 18 pairs on chain A | single copy; all sources agree on the 7 Arpeggio partners; RCSB and local Arpeggio add A73 (consensus 79) |
| 3RKF | A, B, C, D | DX4 ×4 (A/B/C/D 91), NCO ×26 | 6, chain B only | all 30 instances | preferred assembly = chain B; chains A, C, D absent from PDBe |
| 9NDD | A, B, C, D | OXG ×4 (201) | all 4 | all 4 | one assembly; chain-specific contacts available; scaffold residues map to `outside_mapped_region` |
| 9NBH | A, B, C, D | MG ×18 | all 18 | all 18 | no ligand-class molecule; aptamer core unmodelled; MG sites on scaffold (consensus null) or P1 end (consensus 1) |

Example chain-specific records: `cache/www.ebi.ac.uk_pdbe_graph-api_pdb_bound_ligand_interactions_9ndd_A_201.json` (OXG A/201 ↔ A/65 polar, weak_hbond, 3.49 Å), `cache/data.rcsb.org_rest_v1_core_nonpolymer_entity_instance_3RKF_E.json` (DX4 instance E = chain A copy, 8 RNA residues on target_asym_id A).
