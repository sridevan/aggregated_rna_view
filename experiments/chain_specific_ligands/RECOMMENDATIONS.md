# Recommendations: chain-specific RNA–ligand contacts

## Checkpoint answers
1. **Can PDBe provide chain-specific contacts directly?** Yes, through `graph-api/pdb/bound_molecule_interactions/{pdb}/{bm_id}` (residue level, interaction classes) and `graph-api/pdb/bound_ligand_interactions/{pdb}/{chain}/{resnum}` (atom level with distances), joined to label asym ids via `api/pdb/entry/ligand_monomers/{pdb}`. Coverage is limited to ligand instances in PDBe's **preferred assembly**: 3RKF chains A, C, D (24 of 30 ligand instances) have no PDBe interactions at all. Author identifiers only.
2. **Can RCSB?** Yes, `rcsb_ligand_neighbors` / `rcsb_target_neighbors` give every ligand instance ↔ RNA residue pair within 5 Å with label ids, atom ids and distances, for all copies. They are **proximity contacts**, not typed interactions; the only classification is `ligand_is_bound` (covalent / metal coordination).
3. **Most reliable identifiers:** RCSB (native `asym_id` + `seq_id` = our canonical ids, plus model and alt ids). PDBe author ids resolved without error through the existing mapping (0 unresolved in 499 rows), and `ligand_monomers.struct_asym_id` equals the RCSB asym id for all 59 ligand instances, so the two can be joined.
4. **Comparison with Milestone 3B (v2 binding_sites):** for single-copy entities v2 equals the PDBe Arpeggio residue set (1Y26: 7/7). For multi-copy entities v2 is the union over the PDBe-processed copies only (3RKF: chain B; 9NDD: all four) and `bestChainId` can be wrong (3RKF says A, data are chain B). RCSB proximity is a superset of the Arpeggio set by 1–2 residues per ligand (stacked or packing neighbours such as consensus 79 in 1Y26) and occasionally picks a different residue at the 5 Å margin (9NDD A65 vs A66). No source disagrees on the core pocket (consensus 25, 26, 51, 55, 56, 80, 81).
5. **Is local Arpeggio necessary?** No for residue-level chain attribution. Optional for typed, atom-level interactions on copies PDBe skipped; cheap (≈1 s per entry) but the default CLI run does not reproduce PDBe's interaction classes without further preprocessing.
6. **What changes in 3B:** see below. Entity-level counts stay; a chain-level evidence layer is added; `bestChainId` is demoted to metadata with a known-unreliable note.

## Preferred source
Primary: **RCSB ligand neighbours** for chain-specific residue sets (complete, label ids, distances). Secondary: **PDBe graph-api interactions** for interaction classes where available (preferred assembly), joined by ligand instance. Keep **PDBe v2 binding_sites** as the entity-level layer already in production. Local Arpeggio: deferred; revisit if typed interactions are required for all copies.

## Proposed schema changes (Milestone 3B → 3C/3D)
`ligands.json` per bound-molecule record (one per entry/entity/CCD) gains:
```
"ligand_instances": [ { "ligand_label_asym_id": "O", "ligand_auth_asym_id": "B", "ligand_auth_resnum": 91, "pdbe_bm_id": "bm1" | null,
                        "rna_contacts": [ { "seq_name": "3rkf_1_B_B", "label_asym_id": "B", "label_seq_id": 7, "rfam_consensus_pos": 25,
                                            "evidence": { "rcsb_proximity": {"min_distance": 3.3, "ligand_atom": "N9", "rna_atom": "O2", "cutoff": 5.0, "model": 1},
                                                          "pdbe_arpeggio": {"classes": ["aromatic","hbond"], "min_distance": 2.9} | null } } ],
                        "sources": { "rcsb_ligand_neighbors": {"status": "evaluated", "url", "retrieved"}, "pdbe_bound_ligand_interactions": {"status": "evaluated|not_in_preferred_assembly", ...} } } ],
"chain_attribution": "resolved_single_chain" | "resolved_by_instance_evidence" | "unresolved"
```
Per-instance `annotations/<seq_name>.json` residues get `bound_molecules[].attribution ∈ {chain_evidence, entity_unresolved}` with the evidence source, so the Mol* and 2D viewers highlight only residues with chain-level evidence for the selected instance.

`tracks.json` `ligand_binding_sites` keeps the entry-entity counts as primary and adds a second block `chain_level` with its own counting unit (instances with evidence), its own denominators (observable ∧ instance evaluated by RCSB/PDBe) and `evidence_definition` (proximity 5 Å vs Arpeggio classes). The two blocks are never summed; the frontend labels which one is shown. Switch the default family track to chain-level only after the per-position frequencies are shown to be stable against the entity-level ones across all 66 entries (not only the 4 tested).

## Pipeline integration
1. `annotation_sources.py`: add fetchers for PDBe `ligand_monomers`, `bound_molecules`, `bound_ligand_interactions/{pdb}/{chain}/{resnum}` and RCSB `polymer_entity_instance` (one call per RNA instance; 84 calls) with cache seeding as today.
2. `annotations.py`: build `ligand_instances` by joining RCSB `ligand_asym_id` with `ligand_monomers.struct_asym_id`; resolve RNA residues with `resolve_seq_index` (RCSB) and `resolve_unit` (PDBe); apply the existing class/artefact config per CCD; record `not_in_preferred_assembly` when `bm_id` is null.
3. Reconcile with the v2 record: every v2 residue should appear in at least one instance's contacts; report any that do not (traceability), and keep v2 as the entity-level layer.
4. Tests: synthetic multi-copy entity with one processed and one unprocessed ligand copy; identifier join; denominators for both counting units; schema validation.
5. Open items for review: 5 Å proximity as the chain-level contact definition (vs Arpeggio classes); treatment of symmetry-mate contacts (not visible in asymmetric-unit sources); whether to run local Arpeggio for unprocessed copies after calibrating its output against PDBe.
