import type { ChainMapping, MappedResidue } from '../types/backend'

/** A structural residue addressed for Mol* (label identifiers) with author identifiers carried separately for display. */
export interface ResidueTarget {
  rfam_consensus_pos: number | null
  label_asym_id: string
  label_seq_id: number
  auth_asym_id: string
  author_residue_number: number | string | null
  author_insertion_code: string
  residue_name: string
  nucleotide_parent: string
  alignment_state: MappedResidue['alignment_state']
  insertion_anchor: number | null
  insertion_offset: number | null
  posterior_prob_char: string
  has_coordinates: boolean
}
export interface Translation {
  /** match residues with experimental coordinates: the only consensus-derived residues sent to the 3D viewer */
  highlightable: ResidueTarget[]
  /** match residues without coordinates (sequence-level only) */
  unobserved: ResidueTarget[]
  deleted: number[]
  uncovered: number[]
  missing: number[]
}
export interface MappingIndex {
  byConsensus: Map<number, { state: 'match' | 'deletion' | 'uncovered'; polymer: number | null }>
  byPolymer: Map<number, MappedResidue>
}
export function indexMapping(m: ChainMapping): MappingIndex {
  return {
    byConsensus: new Map(m.consensus_columns.map((c) => [c.rfam_consensus_pos, { state: c.state, polymer: c.polymer_seq_position }])),
    byPolymer: new Map(m.residues.map((r) => [r.polymer_seq_position, r])),
  }
}
export function toTarget(m: ChainMapping, r: MappedResidue): ResidueTarget {
  return {
    rfam_consensus_pos: r.rfam_consensus_pos, label_asym_id: m.label_asym_id, label_seq_id: r.polymer_seq_position, auth_asym_id: m.auth_asym_id,
    author_residue_number: r.author_residue_number, author_insertion_code: r.author_insertion_code, residue_name: r.residue_name,
    nucleotide_parent: r.nucleotide_parent, alignment_state: r.alignment_state, insertion_anchor: r.insertion_anchor, insertion_offset: r.insertion_offset,
    posterior_prob_char: r.posterior_prob_char, has_coordinates: r.has_coordinates,
  }
}
/** Family (consensus) coordinates -> structural residues of one chain. Pure; independent of any viewer. */
export function translateConsensus(m: ChainMapping, positions: Iterable<number>, idx: MappingIndex = indexMapping(m)): Translation {
  const t: Translation = { highlightable: [], unobserved: [], deleted: [], uncovered: [], missing: [] }
  for (const pos of [...positions].sort((a, b) => a - b)) {
    const c = idx.byConsensus.get(pos)
    if (!c) { t.missing.push(pos); continue }
    if (c.state === 'deletion') { t.deleted.push(pos); continue }
    if (c.state === 'uncovered') { t.uncovered.push(pos); continue }
    const r = c.polymer == null ? undefined : idx.byPolymer.get(c.polymer)
    if (!r || r.alignment_state !== 'match' || r.rfam_consensus_pos !== pos) { t.missing.push(pos); continue }
    ;(r.has_coordinates ? t.highlightable : t.unobserved).push(toTarget(m, r))
  }
  return t
}
/** Inverse: a chain residue (label_seq_id) -> its consensus position, or its instance-level identity when it has none. */
export function residueToConsensus(m: ChainMapping, labelSeqId: number, idx: MappingIndex = indexMapping(m)): ResidueTarget | null {
  const r = idx.byPolymer.get(labelSeqId)
  return r ? toTarget(m, r) : null
}
export type PositionState = 'match_observed' | 'match_unobserved' | 'deletion' | 'uncovered' | 'missing'
/** Per-consensus-position state of a chain for drawing the instance track. */
export function positionStates(m: ChainMapping, clen: number, idx: MappingIndex = indexMapping(m)): PositionState[] {
  const out: PositionState[] = []
  for (let pos = 1; pos <= clen; pos++) {
    const c = idx.byConsensus.get(pos)
    if (!c) { out.push('missing'); continue }
    if (c.state !== 'match') { out.push(c.state); continue }
    const r = c.polymer == null ? undefined : idx.byPolymer.get(c.polymer)
    out.push(r?.has_coordinates ? 'match_observed' : 'match_unobserved')
  }
  return out
}
/** Insert / local-end residues grouped by anchor (markers between consensus positions). */
export function insertionsByAnchor(m: ChainMapping): Map<number, { insert: number; local_end: number }> {
  const out = new Map<number, { insert: number; local_end: number }>()
  for (const r of m.residues) {
    if ((r.alignment_state === 'insert' || r.alignment_state === 'local_end') && r.insertion_anchor != null) {
      const e = out.get(r.insertion_anchor) ?? { insert: 0, local_end: 0 }
      e[r.alignment_state] += 1
      out.set(r.insertion_anchor, e)
    }
  }
  return out
}
