# Annotation schema (Milestone 3B, implemented)

This replaces the Milestone 3A proposal. Decisions taken for 3B are marked **[3B decision]**. Machine-readable schemas live in `rna_explorer/schemas/*.schema.json`; every file is validated on write.

Principles: keep every source observation with its original identifier and provenance; aggregate only through the validated residue mapping; store numerators and denominators explicitly; never merge categories that overlap (J3/J4) or sources with different observation levels (entity vs chain); never force insertions, local ends or unmapped nucleotides onto consensus positions.

## Files (`data/families/<ACC>/`)

| File | Content | Schema |
|---|---|---|
| `reference_structural_elements.json` | one categorical family-level track transferred from the BGSU reference chain | `reference_structural_elements` |
| `loops.json` | complete Motif Atlas loop catalogue with nucleotide membership and motif groups | `loops` |
| `ligands.json` | bound-molecule records (entity level), per-CCD summary, classification config used | `ligands` |
| `tracks.json` | per-consensus-position aggregation for all tracks with denominators | `tracks` |
| `annotations/<seq_name>.json` | per-instance nucleotide annotations and source status | `annotation_instance` |
| `annotation_summary.json` | build counts, unresolved identifiers, file sizes | (none) |

## A. Reference structural elements
- `reference`: PDB, chain (author and label), entity, seq_name, number of annotated nucleotides, citation lines from the family file.
- `positions[]`: `rfam_consensus_pos`, `element` (first label), `all_labels`, reference unit id, reference label_seq_id, nucleotide, PP of the reference residue.
- `track[]`: length CLEN, label or null. `elements[]`: label with consensus spans.
- `reference_nucleotides_not_on_consensus[]`: reference residues that are insertions/local ends (retained, never forced). `conflicts[]`: two labels for one position (none in RF00167).
- `instance_transfer`: counts of instances for which the BGSU per-chain transfer was evaluated/unavailable and the agreement of those labels with the family track. **[3B decision]** The family track is derived from the reference file through *our* mapping of the reference chain; BGSU's per-chain labels are kept per instance for highlighting and checked for agreement (5,283/5,283 agree), never counted as independent observations.

## B. Loops
Per loop: `loop_id`, `loop_type` (FR3D category preserved), `atlas_status` (`valid` | `problematic` | `not_listed`), `motif_group`, **[3B decision]** `motif_group_assignment` (`member` = instance of the group in the current release; `homologous` / `geometric` = mapped to the group through `motif_group_via_loop_id`), `annotation_text`, `problem` (Atlas reason for problematic loops), `location`, `chains[]` (author + label chain, seq_name, in_family), `multi_chain`, `instances[]`, `strands[]` (contiguous segments; discontinuous loops keep every strand), `nucleotides[]` (unit id, canonical id, alignment state, consensus position or null, mapping_status, PP, coordinates, `resolution` = resolved | chain_not_in_family | residue_not_found), `consensus_positions[]`, status/PP counts, source URLs, retrieval date and Motif Atlas release.
Counting rule: one loop id = one observation whatever the number of chains.

## C. Bound molecules
Per record (one per entry, entity and CCD id): class (`ligand` | `metal_ion` | `additive` | `water`, from `annotation_config.json`), `class_rule`, PDBe `scaffoldId`/`cofactorId`, `observation_level: entity`, `best_chain_label_asym_id` (metadata only), `n_entity_chains`, `entity_chains[]`, `chain_attribution` (`resolved` only for single-chain entities, else `unresolved`), `bound_molecule_ids[]`, `contacts[]` deduplicated per label_seq_id (bound-molecule ids merged per contact; consensus position, mapping status, PP, which chains have coordinates), `flags[]`, `in_default_track`, `exclusion_reason`.
**[3B decision]** Artefact rules in config: `ccd_is_polymer_residue` (excluded from default track) and `contacts_only_polymer_termini` (window 2; excluded from default track; flagged, retained). `by_ccd` summary gives entry/entity counts and consensus positions with entry counts.

## D. Tracks (`positions[]`, length CLEN)
Each position carries `observable` (instances with a modelled match residue; entries; entry-entities) and:
- `reference_structural_elements`: `element`, reference id, `instance_transfer` availability counts.
- `hairpin_loops`, `internal_loops`, `junction_loops_J3`, `junction_loops_J4`: `n_loop_observations`, `n_instances`, `n_entries`, `n_distinct_motif_groups`, `by_motif_group{group: loops/instances/entries}` (+ `unassigned`), `n_loops_without_motif_group`, `loop_ids[]` (traceability), PP `confidence`, `denominators`, `frequency`.
- `junction_loops_union`: set union of J3 and J4 (instances, entries, loop ids), never the sum.
- `ligand_binding_sites`: counting unit `(pdb_id, entity_id)`; `n_entry_entities`, `n_entries`, `entry_entities[]`, `by_ccd_n_entry_entities` (default classes only), `by_class_n_entry_entities` (all classes), `n_entry_entities_excluded_from_default`, PP `confidence`, `denominators`, `frequency`.

Denominators (chain-level tracks): `n_observable_instances`, `n_source_evaluated_instances`, `n_observable_and_evaluated`, `n_positive`, `n_positive_observable_and_evaluated`, `n_not_annotated_by_source`, `n_source_unavailable`, `n_positive_not_observable`. `frequency = n_positive_observable_and_evaluated / n_observable_and_evaluated`, null when the denominator is 0. Entity-level ligand denominators use the same fields over entry-entities. Terminology: "not annotated by the source", never "confirmed absence".

## E. Per-instance records
`sources` (status `evaluated` | `unavailable` per source with URL, retrieval date, attribution scope), `loops[]` involving the chain, `bound_molecules[]` of the entity (attribution `chain` | `entity_unresolved`, `is_best_chain`), `residues[]` (every polymer residue: canonical ids, alignment state, consensus position or null, PP, coordinates, `reference_element`, `loops[]`, `bound_molecules[]`). Residues outside the mapped region and insertions keep their annotations with `rfam_consensus_pos: null`.

## F. Versioning
Every file: `schema_version`, Rfam release and CM checksum, mapping pipeline version and date, annotation pipeline version, config version, generation time, source URLs and retrieval dates, Motif Atlas release (4.15 at build time).
