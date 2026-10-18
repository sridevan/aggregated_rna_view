import type { ChainMapping, MappedResidue } from '../types/backend'
import { indexMapping, type MappingIndex, type PositionState } from './consensus'

/*
 * Position numbering (display only). Rfam consensus match positions stay the internal coordinate system everywhere; a numbering reference
 * only relabels positions with the author residue identifiers of one RNA chain, looked up from that chain's validated mapping
 * (consensus_columns -> residues). No offsets, no re-alignment: a position without a match residue in the reference has no label.
 */

/** Author residue identifier of a deposited residue (residue name + author number + insertion code, e.g. "G949", "A306B", "PSU516");
 *  null when the mapping has no author number (never invented). */
export function authorResidueId(r: Pick<MappedResidue, 'residue_name' | 'author_residue_number' | 'author_insertion_code'>): string | null {
  if (r.author_residue_number == null || r.author_residue_number === '') return null
  return `${r.residue_name}${r.author_residue_number}${r.author_insertion_code ?? ''}`
}
/** Author number + insertion code only (axis labels: "949", "306B") */
export function authorNumber(r: Pick<MappedResidue, 'author_residue_number' | 'author_insertion_code'>): string | null {
  return r.author_residue_number == null || r.author_residue_number === '' ? null : `${r.author_residue_number}${r.author_insertion_code ?? ''}`
}

/** What one chain has at one Rfam match position. `residue` only for a match (the deposited residue aligned to the position). */
export interface ChainAtPosition { state: PositionState; residue: MappedResidue | null }
export function chainAt(idx: MappingIndex, pos: number): ChainAtPosition {
  const c = idx.byConsensus.get(pos)
  if (!c) return { state: 'missing', residue: null }
  if (c.state !== 'match') return { state: c.state, residue: null }
  const r = c.polymer == null ? undefined : idx.byPolymer.get(c.polymer)
  if (!r || r.alignment_state !== 'match' || r.rfam_consensus_pos !== pos) return { state: 'missing', residue: null }
  return { state: r.has_coordinates ? 'match_observed' : 'match_unobserved', residue: r }
}

/** Display value of a chain's residue at a position: author identifier, or the reason there is none (never a number that is not in the data). */
export function residueValue(a: ChainAtPosition): string {
  if (a.residue) return authorResidueId(a.residue) ?? 'Unavailable'
  if (a.state === 'deletion') return 'Not present'
  return 'Unmapped'
}
export const STATUS_TEXT: Record<PositionState, string> = {
  match_observed: 'Coordinates available', match_unobserved: 'Coordinates unavailable', deletion: 'Deletion', uncovered: 'Outside aligned region', missing: 'Not available in the mapping',
}

/** A numbering reference: one chain, its label for display and, per Rfam position (index pos - 1), its author number or null. */
export interface Numbering {
  seqName: string
  /** "5J7L A" (PDB id + author chain) */
  label: string
  idx: MappingIndex
  /** author number (+ insertion code) of the reference residue at each Rfam position; null = no match residue / no author number */
  numbers: (string | null)[]
}
export function buildNumbering(seqName: string, pdbId: string, authAsymId: string, m: ChainMapping, clen: number): Numbering {
  const idx = indexMapping(m)
  const numbers = Array.from({ length: clen }, (_, k) => { const a = chainAt(idx, k + 1); return a.residue ? authorNumber(a.residue) : null })
  return { seqName, label: `${pdbId.toUpperCase()} ${authAsymId}`, idx, numbers }
}

const STEPS = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
/**
 * Axis ticks for a numbering reference at a zoom level: Rfam positions whose reference author number is a round number (step chosen like
 * the Rfam axis: labels about `minPx` apart), thinned so labels never overlap; long stretches without a round number (irregular numbering,
 * insertion codes, gaps) are filled with any mapped residue. Author numbers are never assumed continuous or increasing.
 */
export function referenceTicks(numbers: (string | null)[], cellPx: number, minPx = 60): { pos: number; text: string }[] {
  const step = STEPS.find((t) => t * cellPx >= minPx) ?? STEPS[STEPS.length - 1]
  const width = (s: string) => s.length * 5.2 + 10                     // label width (px) incl. a margin
  const centre = (pos: number) => (pos - 0.5) * cellPx
  const picked: { pos: number; text: string }[] = []
  const fits = (pos: number, text: string) => picked.every((p) => Math.abs(centre(p.pos) - centre(pos)) >= Math.max(minPx * 0.8, (width(p.text) + width(text)) / 2))
  numbers.forEach((n, k) => { if (n != null && /^-?\d+$/.test(n) && +n % step === 0 && fits(k + 1, n)) picked.push({ pos: k + 1, text: n }) })
  // fill stretches without a round number (irregular numbering, insertion codes, gaps): any mapped residue at least 1.5 spacings from every label
  numbers.forEach((n, k) => { if (n != null && picked.every((p) => Math.abs(centre(p.pos) - centre(k + 1)) >= 1.5 * minPx) && fits(k + 1, n)) picked.push({ pos: k + 1, text: n }) })
  return picked.sort((a, b) => a.pos - b.pos)
}
