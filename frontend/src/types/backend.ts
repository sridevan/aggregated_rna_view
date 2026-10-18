/** Types mirroring rna_explorer/schemas/*.schema.json (Milestone 1 data contracts). */
export type AlignmentState = 'match' | 'insert' | 'local_end' | 'not_submitted_flank'
export type ConsensusState = 'match' | 'deletion' | 'uncovered'

export interface Provenance {
  rfam_acc: string; rfam_release: string; cm_name?: string; cm_checksum: string; cm_clen: number; cm_build_date?: string
  infernal_version: string; cmalign_options: string; cmalign_mode?: string; pipeline_version?: string; generated: string
  pdbe_mapping_source?: string; rfam_hit_source?: string
}
export interface CmalignScore { name?: string; length?: number; cm_from: number; cm_to: number; trunc: string; bit_score: number; avg_pp: number }

export interface InstanceRow {
  pdb_id: string; entity_id: number; label_asym_id: string; auth_asym_id: string; seq_name: string; sequence_key: string; mapping_file: string
  organism: string | null; tax_id: number | string | null; organism_source: string | null; molecule_name: string
  experimental_method: string | null; resolution: number | null; release_year: number | null; title: string | null
  chain_length: number; map_start: number; map_end: number; map_auth_start: number | string | null; map_auth_end: number | string | null
  mapped_length: number; n_modified_in_mapped: number
  rfam_bit_score: number | null; rfam_evalue: string | number | null; rfam_cm_start: number | null; rfam_cm_end: number | null
  clan_families_hit: { acc: string; bits: number }[]; clan_best_family: string | null; is_clan_best: boolean | null
  cmalign: CmalignScore | null
  n_match: number; n_insert: number; n_local_end: number; n_deleted_in_span: number; n_uncovered: number; consensus_span: [number, number]
  fraction_clen_matched: number; n_mapped_with_coordinates: number; n_match_with_coordinates: number; pp_star_fraction: number | null
  pp_category_counts: Record<string, number>; validation_ok: boolean | null
  organism_species?: string | null; species_tax_id?: string | null; organism_rank?: string | null; pdbe_2d_available?: boolean | null
}
export interface DisplayRepresentative extends InstanceRow { selection_rule: string; eligible_count: number; min_clen_fraction_used: number }
export interface FamilyOverview {
  rfam_acc: string; rfam_id: string; description: string; comment?: string | null; rfam_url?: string; display_representative?: DisplayRepresentative | null; clen: number; consensus_sequence: string; ss_cons: string; clan_members: string[]
  curation?: Record<string, unknown> | null
  counts: { pdb_entries: number; instances: number; unique_mapped_sequences: number; instances_clan_best: number; instances_validation_ok: number; skipped_segments: number }
  experimental_methods: Record<string, number>; representative: InstanceRow | null
  coverage_summary: { sequence_coverage_mean: number | null; structural_coverage_mean: number | null; positions_with_sequence_coverage: number; positions_with_structural_coverage: number; positions_clen: number }
  provenance: Provenance
}
export interface LevelStats {
  n_match: number[]; n_deletion: number[]; n_uncovered: number[]; coverage: (number | null)[]
  pp_n: number[]; pp_sum: number[]; pp_mean: (number | null)[]; pp_min: (number | null)[]; pp_high: number[]; pp_medium: number[]; pp_low: number[]
}
export interface SequenceLevelStats extends LevelStats { n_unique_sequences: number }
export interface StructureLevelStats extends LevelStats { n_instances: number; n_match_with_coordinates: number[]; structural_coverage: (number | null)[] }
export interface ConsensusStats {
  rfam_acc: string; clen: number; rf: string[]; ss_cons: string[]
  sequence_level: SequenceLevelStats; structure_level: StructureLevelStats
  pp_bins: Record<string, number>; pp_display_thresholds: { high: number; medium: number }; notes?: string[]; provenance: Provenance
}
export interface MappedResidue {
  polymer_seq_position: number; author_residue_number: number | string | null; author_insertion_code: string; residue_name: string; nucleotide_parent: string
  is_modified: boolean; alignment_state: AlignmentState; rfam_consensus_pos: number | null; insertion_anchor: number | null; insertion_offset: number | null
  posterior_prob_char: string; posterior_prob: number | null; has_coordinates: boolean
}
export interface ConsensusColumn { rfam_consensus_pos: number; state: ConsensusState; polymer_seq_position: number | null }
export interface ChainMapping {
  pdb_id: string; entity_id: number; label_asym_id: string; auth_asym_id: string; rfam_acc: string; sequence_key: string
  mapped_region: { start: number; end: number }; alignment: CmalignScore | null; summary?: Record<string, unknown>
  residues: MappedResidue[]; consensus_columns: ConsensusColumn[]; provenance: Provenance
}

/** orientation.json (rna_explorer/orientation.py): rigid-body transforms of each instance onto the family's display representative.
 *  Convention: x_ref = R · x_target + t (column vectors), R row-major, det(R) = +1. */
export interface OrientationTransform {
  target: { pdb_id: string; entity_id: number; label_asym_id: string; auth_asym_id: string }
  rotation: number[][] | null; translation: number[] | null; rmsd: number | null; determinant: number | null
  n_pairs: number; n_shared_consensus: number; n_shared_consensus_with_coords: number
  reference_coverage: number | null; target_coverage: number | null; pp_filter_used: string
  atoms_used: Record<string, number>; excluded: Record<string, number>; paired_consensus_positions: number[]
  fit_quality: 'good' | 'moderate' | 'poor' | 'insufficient'; is_identity: boolean
  check: { consensus_pos: number; label_seq_id: number; atom: string; original: number[]; transformed: number[] } | null
}
export interface OrientationDataset {
  schema_version: string; rfam_acc: string
  reference: { seq_name: string; pdb_id: string; entity_id: number; label_asym_id: string; auth_asym_id: string; selection_rule: string }
  method: Record<string, unknown>; transforms: Record<string, OrientationTransform>; provenance: Record<string, unknown>
}

/** tracks.json (rna_explorer/annotations.py, Milestone 3B): per-consensus-position aggregation with explicit denominators. */
export interface LoopDenominators { n_observable_instances: number; n_source_evaluated_instances: number; n_observable_and_evaluated: number; n_positive: number; n_positive_observable_and_evaluated: number; n_not_annotated_by_source: number; n_source_unavailable: number; n_positive_not_observable: number }
export interface LoopTrackPosition {
  n_loop_observations: number; n_instances: number; n_entries: number; n_distinct_motif_groups: number
  by_motif_group: Record<string, { n_loops: number; n_instances: number; n_entries: number }>; n_loops_without_motif_group: number; loop_ids: string[]
  confidence: { high: number; medium: number; low: number }; denominators: LoopDenominators; frequency: number | null
}
export interface JunctionUnionPosition { n_loop_observations: number; n_instances: number; n_entries: number; denominators: LoopDenominators; frequency: number | null; rule: string
  /** junction types in the union and instances per type (absent in files built before J5 support) */
  loop_types?: string[]; by_type_n_instances?: Record<string, number> }
export interface LigandTrackPosition {
  n_entry_entities: number; n_entries: number; entry_entities: string[]; by_ccd_n_entry_entities: Record<string, number>
  /** per bound ccd in contact here, counted per PDB entry; denominator = entries with the bound ccd and an observed nucleotide here */
  by_ccd_entries: Record<string, { n_positive: number; n_positive_observable_and_evaluated: number; n_observable_and_evaluated: number; n_not_observable: number; n_source_unavailable: number; frequency: number | null }>
  by_class_n_entry_entities: Record<string, number>
  n_entry_entities_excluded_from_default: number; confidence: { high: number; medium: number; low: number }
  denominators: { counting_unit: string; n_observable_entry_entities: number; n_source_evaluated_entry_entities: number; n_observable_and_evaluated: number; n_positive: number; n_positive_observable_and_evaluated: number; n_not_annotated_by_source: number; n_source_unavailable: number; n_positive_not_observable: number }
  frequency: number | null
}
export interface TrackPosition {
  rfam_consensus_pos: number; observable: { n_instances: number; n_entries: number; n_entry_entities: number }
  reference_structural_elements: { element: string | null; reference: string }
  hairpin_loops: LoopTrackPosition; internal_loops: LoopTrackPosition; junction_loops_J3: LoopTrackPosition; junction_loops_J4: LoopTrackPosition
  junction_loops_union: JunctionUnionPosition; ligand_binding_sites: LigandTrackPosition
}
export interface TracksDataset {
  schema_version: string; rfam_acc: string; clen: number; n_instances: number; n_entries: number; n_entry_entities: number
  source_availability: Record<string, Record<string, number>>; terminology: Record<string, string>; tracks: Record<string, unknown>; positions: TrackPosition[]; provenance: Record<string, unknown>
}
/** ligands.json (Milestone 3B): bound-molecule catalogue (entity level) and per-CCD summary. */
export interface LigandSummary { name: string | null; ligand_class: string; n_entry_entities: number; n_entries: number; entries: string[]; n_bound_molecule_instances: number; n_entries_in_default_track: number; consensus_positions_n_entries: Record<string, number>; flag_counts: Record<string, number> }
export interface BoundMoleculeRecord { pdb_id: string; entity_id: number; ccd_id: string; name: string | null; ligand_class: string; chain_attribution: string; n_entity_chains: number; consensus_positions: number[]; flags: string[]; in_default_track: boolean; exclusion_reason: string | null }
export interface LigandsDataset { schema_version: string; rfam_acc: string; clen: number; source: Record<string, unknown>; classification: { default_track_classes: string[] }; counts: Record<string, unknown>; by_ccd: Record<string, LigandSummary>; bound_molecules: BoundMoleculeRecord[]; provenance: Record<string, unknown> }

/** sasa_tracks.json (rna_explorer/sasa.py): median SASA of the isolated RNA chain per consensus position. */
export interface SasaTrackPosition {
  rfam_consensus_pos: number; n_contributing_instances: number; n_mapped_instances_with_coordinates: number; n_excluded_incomplete: number; n_excluded_unsupported: number
  nucleotide_identity: Record<string, number>; contributing_instances?: string[]
  median_sasa_A2: number | null; mean_sasa_A2: number | null; min_sasa_A2: number | null; max_sasa_A2: number | null; sd_sasa_A2: number | null; q1_sasa_A2: number | null; q3_sasa_A2: number | null
}
export interface SasaTracksDataset {
  schema_version: string; rfam_acc: string; clen: number; track_id: string; category: string
  measurement: { name: string; unit: string; statistic?: string; context?: string }; aggregation_rules?: Record<string, string>
  colour_domain_suggestion_A2?: { min: number | null; max: number | null; rule?: string }
  summary: { n_positions_with_value: number; n_instances_measured?: number; median_of_medians_A2?: number | null; max_median_A2?: number | null }
  positions: SasaTrackPosition[]; method: Record<string, unknown>; provenance: Record<string, unknown>
}

/** pocket_tracks.json (rna_explorer/pockets.py): fpocketR pocket occurrence per consensus position, three filter variants. */
export interface PocketTrackPosition {
  rfam_consensus_pos: number; n_supporting_instances: number; n_evaluable_instances: number; frequency: number | null; n_unevaluable_instances: number
  unevaluable_reasons?: Record<string, number>; n_supporting_entries: number; n_distinct_pockets: number; pocket_ids?: string[]; score_summary?: { min: number; median: number; max: number } | null
}
export interface PocketTracksDataset {
  schema_version: string; rfam_acc: string; clen: number; track_id: string; category: string; measurement: { name: string; definition: string; evaluable?: string }
  filters: Record<string, string>; n_instances: number; status_counts: Record<string, number>
  tracks: Record<string, PocketTrackPosition[]>; method: Record<string, unknown>; provenance: Record<string, unknown>
}

/** modification_tracks.json (rna_explorer/modifications.py): nonstandard deposited nucleotide components per consensus position. */
export interface ModificationTrackPosition {
  rfam_consensus_pos: number; n_positive: number; n_negative: number; n_missing_identity: number; missing_reasons: Record<string, number>; n_no_residue: number; n_source_unavailable: number
  denominator: number; frequency: number | null; n_supporting_entries: number; n_denominator_without_coordinates: number
  ccd: { ccd_id: string; n_instances: number; n_entries: number }[]; supporting_observations: string[]
}
export interface ModificationComponent {
  chem_comp_name: string | null; chem_comp_type: string | null; parent: string | null; parent_source: 'CCD' | 'NAKB' | null; ccd_parent: string | null; nakb_parent: string[] | null
  nakb_change_locations: string[] | null; n_observations: number; n_instances: number; n_entries: number; n_on_consensus_axis: number; consensus_positions: number[]; pdbe_reported: number | null
}
export interface ModificationTracksDataset {
  schema_version: string; rfam_acc: string; clen: number; track_id: string; category: string; measurement: { name: string; definition: string; unit?: string }
  aggregation_rules: Record<string, string>; summary: { n_instances: number; n_instances_ok?: number; n_nonstandard: number; n_on_consensus_axis: number; n_off_axis: number; off_axis_observations: unknown[] }
  components: Record<string, ModificationComponent>; positions: ModificationTrackPosition[]; method: Record<string, unknown>; provenance: Record<string, unknown>
}

/** structure_clustering.json (rna_explorer/structure_clustering.py): structural groups from precomputed pairwise matrices. */
export type StructureMetricId = 'discrepancy' | 'tm'
export interface StructureGroupRecord {
  group_id: number; label: string; node: number; members: string[]; n_instances: number; n_entries: number; representative: string; representative_rule: string
  /** plain medoid (may differ when it has internal coordinate gaps) and the representative's unmodelled internal nucleotides; absent in older files */
  medoid?: string; representative_internal_gaps?: string[] | null
  representative_mean_dissimilarity: number; leaf_range: [number, number]; within: { mean: number | null; max: number | null; n_pairs: number }
  across: { mean: number | null; n_pairs: number }; across_to_within_ratio: number | null; small: boolean; instance_level_view: boolean
}
export interface StructureMetricClustering {
  method: string; linkage_method: string; distance: string; linkage: [number, number, number, number][]; leaf_order: number[]; cophenetic_correlation: number; max_height: number
  cut: { k: number; height: number; lowest_merge_undone: number; highest_merge_kept: number; gap: number; rule: string; is_largest_gap_of_tree: boolean; top_merge_heights: number[]; silhouette_by_k: Record<string, number> }
  groups: StructureGroupRecord[]; assignments: Record<string, number>; validation: 'stable' | 'provisional'; rationale: string
  floor: { value: number | null; definition: string; n_pairs?: number | null; source?: string }; research_evidence: Record<string, unknown>
}
export interface StructureMetric {
  label: string; short: string; kind: 'dissimilarity' | 'similarity'; stored_quantity: string; distance: string; higher_is_more_similar: boolean; unit: string; definition: string
  values_upper: number[]; values_layout: string; value_range: { min: number; max: number }; distance_range: { min: number; max: number }; clustering: StructureMetricClustering
}
export interface StructureClusteringDataset {
  schema_version: string; rfam_acc: string; visualisation: { detail_max: number; small_group: number }
  dataset: { n_mapped_instances: number; n_mapped_entries: number; n_clustered_instances: number; n_clustered_entries: number; n_excluded_instances: number; inclusion_rule: string
    coverage_threshold: number; shared_positions: number[]; n_shared_positions: number; instances: string[]; excluded: { seq_name: string; pdb_id: string; n_observed: number; reason: string; reason_text: string }[] }
  default_metric: StructureMetricId; metrics: Record<StructureMetricId, StructureMetric>; provenance: Record<string, unknown>
}
/** modifications.json observations (subset used by the Structures tab's associated features). */
export interface ModificationObservation { observation_id: string; seq_name: string; pdb_id: string; ccd_id: string; chem_comp_name: string | null; on_consensus_axis: boolean; outside_mapped_region: boolean }
export interface ModificationsDataset { observations: ModificationObservation[]; instances?: { seq_name: string; status: string }[] }

/** loops.json (rna_explorer/annotations.py): BGSU RNA 3D Motif Atlas loops with their complete nucleotide membership. */
export interface LoopNucleotide { seq_name: string | null; pdb_id: string; label_asym_id: string; auth_asym_id: string; label_seq_id: number; author_residue_number: number | string; insertion_code: string; residue_name: string
  alignment_state: string; rfam_consensus_pos: number | null; has_coordinates: boolean; mapping_status: string; unit_id: string }
export interface LoopRecord { loop_id: string; pdb_id: string; loop_type: 'HL' | 'IL' | 'J3' | 'J4' | string; motif_group: string | null
  /** BGSU assignment of the motif group (member of the release | homologous | geometric, i.e. mapped through another loop); BGSU annotation of the group */
  motif_group_assignment?: string | null; motif_group_via_loop_id?: string | null; annotation_text?: string | null; atlas_status?: string; source?: { motif_atlas_release?: string }
  n_nucleotides: number; multi_chain: boolean; instances: string[]
  strands: { auth_asym_id: string; start_unit_id: string; end_unit_id: string; nucleotide_indices: [number, number]; n_nucleotides: number }[]; nucleotides: LoopNucleotide[] }
export interface LoopsDataset { loops: LoopRecord[] }
