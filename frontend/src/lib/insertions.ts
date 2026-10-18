import type { ChainMapping } from '../types/backend'
import type { GroupSelection } from '../state/selection'

/*
 * Rfam insertion events of one chain, read from the validated per-chain mapping (cmalign alignment states). Only residues whose
 * alignment_state is 'insert' (CM insert state) or 'local_end' (skipped by a local-end transition) are insertion residues;
 * 'not_submitted_flank' residues lie outside the Rfam-mapped region and are never insertions; gaps in author numbering and author
 * insertion codes play no role. insertion_anchor = the consensus position preceding the insert column (0 = before position 1,
 * CLEN = after the last); insertion_offset = 1..n within the run. Inserted residues never receive a consensus position.
 */
export interface InsertionResidue {
  label_seq_id: number; author_residue_number: number | string | null; author_insertion_code: string
  residue_name: string; nucleotide_parent: string; offset: number; has_coordinates: boolean; posterior_prob_char: string
}
export interface InsertionRun {
  id: string                                  // `${anchor}:${first label_seq_id}` (unique within a chain)
  kind: 'insert' | 'local_end'
  anchor: number
  left: number | null                         // flanking consensus match positions (null at the model ends)
  right: number | null
  residues: InsertionResidue[]                // in sequence order
  sequence: string                            // parent nucleotides in sequence order
}
export interface InsertionMarker { id: string; boundaries: number[]; runs: InsertionRun[]; length: number }

/** Ordered insertion runs: consecutive residues (label_seq_id +1) with the same state and anchor and offsets 1, 2, ... form one run. */
export function insertionRuns(m: ChainMapping, clen: number): InsertionRun[] {
  const runs: InsertionRun[] = []
  let cur: InsertionRun | null = null
  for (const r of [...m.residues].sort((a, b) => a.polymer_seq_position - b.polymer_seq_position)) {
    const kind = r.alignment_state
    if ((kind !== 'insert' && kind !== 'local_end') || r.insertion_anchor == null) { cur = null; continue }
    const res: InsertionResidue = { label_seq_id: r.polymer_seq_position, author_residue_number: r.author_residue_number, author_insertion_code: r.author_insertion_code,
      residue_name: r.residue_name, nucleotide_parent: r.nucleotide_parent, offset: r.insertion_offset ?? 0, has_coordinates: r.has_coordinates, posterior_prob_char: r.posterior_prob_char }
    const last = cur?.residues[cur.residues.length - 1]
    if (cur && last && cur.kind === kind && cur.anchor === r.insertion_anchor && res.label_seq_id === last.label_seq_id + 1 && res.offset === last.offset + 1) {
      cur.residues.push(res); cur.sequence += res.nucleotide_parent
    } else {
      const a: number = r.insertion_anchor
      cur = { id: `${a}:${res.label_seq_id}`, kind, anchor: a, left: a >= 1 ? a : null, right: a < clen ? a + 1 : null, residues: [res], sequence: res.nucleotide_parent }
      runs.push(cur)
    }
  }
  return runs
}

/** Visible width of a "+N" label in track units (bold font 7.5: ~4.2 units per character, plus padding). */
export const markerWidth = (n: number) => 4 + 4.2 * (`+${n}`).length
/**
 * One marker per insertion boundary; runs at the same boundary stay separate inside the marker. Markers whose labels would overlap
 * at the shared grid scale (cell width `cell`) are aggregated into one marker; the consensus grid itself never moves.
 */
export function insertionMarkers(runs: InsertionRun[], cell: number): InsertionMarker[] {
  const byAnchor = new Map<number, InsertionRun[]>()
  for (const r of runs) byAnchor.set(r.anchor, [...(byAnchor.get(r.anchor) ?? []), r])
  const out: InsertionMarker[] = []
  for (const a of [...byAnchor.keys()].sort((x, y) => x - y)) {
    const rs = byAnchor.get(a)!, n = rs.reduce((s, r) => s + r.residues.length, 0)
    const prev = out[out.length - 1]
    const px = prev ? prev.boundaries.reduce((s, b) => s + b, 0) / prev.boundaries.length : 0
    if (prev && (a - px) * cell < (markerWidth(prev.length) + markerWidth(n)) / 2 + 1) {
      prev.boundaries.push(a); prev.runs.push(...rs); prev.length += n; prev.id = `${prev.boundaries[0]}-${a}`
    } else out.push({ id: String(a), boundaries: [a], runs: [...rs], length: n })
  }
  return out
}

/** Boundary wording without inventing a consensus position for the inserted residues. */
export function boundaryText(r: Pick<InsertionRun, 'left' | 'right'>): string {
  if (r.left != null && r.right != null) return `Between Rfam positions ${r.left} and ${r.right}`
  if (r.right != null) return `Before Rfam position ${r.right} (start of the model)`
  if (r.left != null) return `After Rfam position ${r.left} (end of the model)`
  return 'Position in the model unknown'
}
export const residueLabel = (x: InsertionResidue) => `${x.nucleotide_parent}${x.label_seq_id}`
export const authorLabel = (x: InsertionResidue) => `${x.author_residue_number ?? '?'}${x.author_insertion_code ?? ''}`

/** The run's residues as an instance-bound group (no consensus positions), for highlighting in 2D and 3D. */
export function insertionSelection(run: InsertionRun, pdbId: string, seqName: string): GroupSelection {
  return { kind: 'group', id: `insertion:${seqName}:${run.id}`, label: `Insertion +${run.residues.length}`, type: 'INS', source: { pdb_id: pdbId, seq_name: seqName },
    members: run.residues.map((x) => ({ unit_id: `${seqName}:${x.label_seq_id}`, rfam_consensus_pos: null, seq_name: seqName, label_seq_id: x.label_seq_id })),
    segments: [[0, run.residues.length - 1]] }
}
