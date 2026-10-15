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
