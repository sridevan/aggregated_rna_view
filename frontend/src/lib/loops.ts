import type { LoopNucleotide, LoopRecord } from '../types/backend'
import type { GroupSelection } from '../state/selection'

/*
 * Per-instance BGSU loop tracks for Family Overview. Identity = original BGSU loop ID; membership = the loop's own nucleotide list (loops.json),
 * never a residue-number range. A block is drawn per original strand (segment) over the consensus match positions of that strand's members in
 * the selected chain; members without a match position (inserts, flanks) or in another chain stay in the membership but have no column.
 */
export type LoopTrack = 'HL' | 'IL' | 'JL'
export interface LoopBlock { from: number; to: number; strand: number }
export interface InstanceLoop {
  id: string; type: string; track: LoopTrack; typeName: string
  /** explicit BGSU motif-group assignment only (never inferred): group id, BGSU annotation text, assignment type, loop it was mapped through */
  motif: { group: string; name: string | null; assignment: string | null; via: string | null; release: string | null } | null
  nucleotides: LoopNucleotide[]          // complete original membership (all chains)
  inChain: LoopNucleotide[]              // members in the selected chain
  nOtherChain: number                    // members in another chain (inter-chain loops)
  nWithoutCoordinates: number            // members of the selected chain without coordinates
  nStrands: number; blocks: LoopBlock[]; lane: number
}
export const TRACK_OF = (t: string): LoopTrack => (t === 'HL' ? 'HL' : t === 'IL' ? 'IL' : 'JL')
export const TYPE_NAME: Record<string, string> = { HL: 'Hairpin loop', IL: 'Internal loop', J3: 'Three-way junction loop', J4: 'Four-way junction loop' }

export function instanceLoops(loops: LoopRecord[], seqName: string): InstanceLoop[] {
  const out: InstanceLoop[] = []
  for (const l of loops) {
    if (!l.nucleotides.some((n) => n.seq_name === seqName)) continue
    const blocks: LoopBlock[] = []
    l.strands.forEach((s, k) => {
      const ps = l.nucleotides.slice(s.nucleotide_indices[0], s.nucleotide_indices[1] + 1).filter((n) => n.seq_name === seqName && n.alignment_state === 'match' && n.rfam_consensus_pos != null).map((n) => n.rfam_consensus_pos!)
      if (ps.length) blocks.push({ from: Math.min(...ps), to: Math.max(...ps), strand: k })
    })
    const inChain = l.nucleotides.filter((n) => n.seq_name === seqName)
    out.push({ id: l.loop_id, type: l.loop_type, track: TRACK_OF(l.loop_type), typeName: TYPE_NAME[l.loop_type] ?? l.loop_type,
      motif: l.motif_group ? { group: l.motif_group, name: l.annotation_text ?? null, assignment: l.motif_group_assignment ?? null, via: l.motif_group_via_loop_id ?? null, release: l.source?.motif_atlas_release ?? null } : null, nucleotides: l.nucleotides, inChain,
      nOtherChain: l.nucleotides.length - inChain.length, nWithoutCoordinates: inChain.filter((n) => !n.has_coordinates).length, nStrands: l.strands.length, blocks, lane: 0 })
  }
  // lanes: within a track, a loop goes to the first lane where none of its blocks overlaps a block already there (overlaps are kept, never merged)
  for (const t of ['HL', 'IL', 'JL'] as const) {
    const lanes: LoopBlock[][] = []
    for (const lp of out.filter((x) => x.track === t && x.blocks.length).sort((a, b) => a.blocks[0].from - b.blocks[0].from || a.id.localeCompare(b.id))) {
      let k = lanes.findIndex((lane) => !lane.some((b) => lp.blocks.some((c) => c.from <= b.to && b.from <= c.to)))
      if (k < 0) { lanes.push([]); k = lanes.length - 1 }
      lanes[k].push(...lp.blocks); lp.lane = k
    }
  }
  return out.sort((a, b) => a.id.localeCompare(b.id))
}

/** Whole-loop selection: every original member (all segments), in the source order; segments = the original strands. */
export function loopSelection(lp: InstanceLoop, pdbId: string, seqName: string, strands: LoopRecord['strands']): GroupSelection {
  return { kind: 'group', id: lp.id, label: `${lp.typeName} · ${lp.id}`, type: lp.type, source: { pdb_id: pdbId, seq_name: seqName },
    members: lp.nucleotides.map((n) => ({ unit_id: n.unit_id, rfam_consensus_pos: n.alignment_state === 'match' ? n.rfam_consensus_pos : null, seq_name: n.seq_name ?? '', label_seq_id: n.label_seq_id })),
    segments: strands.map((s) => [s.nucleotide_indices[0], s.nucleotide_indices[1]] as [number, number]) }
}

/** Structured tooltip content of a loop: category, original BGSU loop id, counts, and the motif group only when BGSU assigns one explicitly. */
export function loopTooltip(lp: InstanceLoop): { title: string; id: string; counts: string; motif: string | null; motifNote: string | null; otherChain: string | null } {
  const n = lp.nucleotides.length
  const m = lp.motif
  return {
    title: `${lp.typeName} (${lp.type})`,
    id: lp.id,
    counts: `${n} nucleotide${n === 1 ? '' : 's'} · ${lp.nStrands} segment${lp.nStrands === 1 ? '' : 's'}`,
    motif: m ? `Motif: ${m.group}${m.name ? ` (${m.name})` : ''}` : null,
    motifNote: m ? [m.release ? `Motif Atlas ${m.release}` : '', m.assignment === 'member' ? 'loop is a group member' : m.assignment ? `${m.assignment} assignment${m.via ? ` via ${m.via}` : ''}` : ''].filter(Boolean).join(' · ') || null : null,
    otherChain: lp.nOtherChain ? `${lp.nOtherChain} in another chain` : null,
  }
}
