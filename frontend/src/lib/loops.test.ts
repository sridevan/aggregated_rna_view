import { describe, expect, it } from 'vitest'
import { instanceLoops, loopSelection, loopTooltip, TRACK_OF, typeName } from './loops'
import type { LoopNucleotide, LoopRecord } from '../types/backend'

const nt = (id: number, pos: number | null, state = 'match', sn = 'x_1_A_A', coords = true): LoopNucleotide => ({ seq_name: sn, pdb_id: 'x', label_asym_id: sn.split('_')[3], auth_asym_id: sn.split('_')[2], label_seq_id: id,
  author_residue_number: id + 10, insertion_code: '', residue_name: 'A', alignment_state: state, rfam_consensus_pos: pos, has_coordinates: coords, mapping_status: state === 'match' ? 'mapped' : state, unit_id: `X|1|A|A|${id + 10}` })
const loop = (id: string, type: string, strands: LoopNucleotide[][]): LoopRecord => {
  let k = 0
  return { loop_id: id, pdb_id: 'x', loop_type: type, motif_group: null, n_nucleotides: strands.flat().length, multi_chain: false, instances: ['x_1_A_A'],
    strands: strands.map((s) => { const r = { auth_asym_id: 'A', start_unit_id: s[0].unit_id, end_unit_id: s[s.length - 1].unit_id, nucleotide_indices: [k, k + s.length - 1] as [number, number], n_nucleotides: s.length }; k += s.length; return r }),
    nucleotides: strands.flat() }
}

describe('per-instance loop tracks', () => {
  const recs = [
    loop('HL_X_001', 'HL', [[nt(1, 10), nt(2, 11), nt(3, 12)]]),
    loop('HL_X_002', 'HL', [[nt(4, 13), nt(5, 14)]]),                                                   // adjacent to HL_X_001, different ID
    loop('IL_X_001', 'IL', [[nt(6, 20), nt(7, null, 'insert'), nt(8, 21)], [nt(30, 40), nt(31, 41, 'match', 'x_1_A_A', false)]]),   // two segments; an insert; one member without coordinates
    loop('J3_X_001', 'J3', [[nt(1, 10), nt(2, 11)], [nt(30, 40)], [nt(50, 60)]]),
    loop('J4_X_001', 'J4', [[nt(2, 11), nt(3, 12)], [nt(50, 60)]]),                                    // overlaps J3_X_001 in the same track
    loop('IL_X_002', 'IL', [[nt(90, null, 'not_submitted_flank')], [nt(91, null, 'not_submitted_flank')]]),   // no consensus position at all
    loop('IL_X_003', 'IL', [[nt(70, 30), nt(71, 31)], [nt(5, 50, 'match', 'x_1_C_C'), nt(6, 51, 'match', 'x_1_C_C')]]),   // inter-chain
    loop('HL_Y_001', 'HL', [[nt(1, 10, 'match', 'y_1_A_A')]]),                                           // another instance
  ]
  const L = instanceLoops(recs, 'x_1_A_A'), by = Object.fromEntries(L.map((l) => [l.id, l]))
  it('keeps the BGSU IDs of this instance only', () => { expect(L.map((l) => l.id)).toEqual(['HL_X_001', 'HL_X_002', 'IL_X_001', 'IL_X_002', 'IL_X_003', 'J3_X_001', 'J4_X_001']) })
  it('one block per original segment; adjacent loops stay distinct', () => {
    expect(by.HL_X_001.blocks).toEqual([{ from: 10, to: 12, strand: 0 }]); expect(by.HL_X_002.blocks).toEqual([{ from: 13, to: 14, strand: 0 }])
    expect(by.HL_X_001.lane).toBe(0); expect(by.HL_X_002.lane).toBe(0)                                  // adjacent, not overlapping: same lane, two blocks
    expect(by.IL_X_001.blocks).toEqual([{ from: 20, to: 21, strand: 0 }, { from: 40, to: 41, strand: 1 }])
  })
  it('overlapping loops of one track go to separate lanes; JL groups J3 and J4', () => {
    expect([by.J3_X_001.track, by.J4_X_001.track]).toEqual(['JL', 'JL']); expect(new Set([by.J3_X_001.lane, by.J4_X_001.lane])).toEqual(new Set([0, 1]))
  })
  it('membership is the original list: inserts, flanks, other-chain members and members without coordinates are kept', () => {
    expect(by.IL_X_001.nucleotides.map((n) => n.label_seq_id)).toEqual([6, 7, 8, 30, 31]); expect(by.IL_X_001.nWithoutCoordinates).toBe(1)
    expect(by.IL_X_002.blocks).toEqual([]); expect(by.IL_X_002.nucleotides.length).toBe(2)
    expect(by.IL_X_003.nOtherChain).toBe(2); expect(by.IL_X_003.blocks).toEqual([{ from: 30, to: 31, strand: 0 }])     // only this chain's segment is drawn
  })
  it('whole-loop selection: every member, original segments, no residue range', () => {
    const rec = recs[2], g = loopSelection(by.IL_X_001, 'x', 'x_1_A_A', rec.strands)
    expect(g.members.map((m) => m.label_seq_id)).toEqual([6, 7, 8, 30, 31]); expect(g.segments).toEqual([[0, 2], [3, 4]])
    expect(g.members.find((m) => m.label_seq_id === 7)!.rfam_consensus_pos).toBeNull()                  // insert: no consensus coordinate
    expect(g.members.some((m) => m.label_seq_id > 8 && m.label_seq_id < 30)).toBe(false)               // nothing from the intervening stem
    expect(g.source).toEqual({ pdb_id: 'x', seq_name: 'x_1_A_A' }); expect(g.label).toBe('Internal loop · IL_X_001')
  })
})

describe('structured loop tooltip', () => {
  const base = loop('IL_X_001', 'IL', [[nt(6, 20), nt(7, null, 'insert'), nt(8, 21)], [nt(30, 40)]])
  it('category, original loop id, nucleotide and segment counts; no motif line without an explicit assignment', () => {
    const t = loopTooltip(instanceLoops([base], 'x_1_A_A')[0])
    expect(t).toMatchObject({ title: 'Internal loop (IL)', id: 'IL_X_001', counts: '4 nucleotides · 2 segments', motif: null, motifNote: null, otherChain: null })
  })
  it('singular counts and junction category names', () => {
    const t = loopTooltip(instanceLoops([loop('J3_X_001', 'J3', [[nt(1, 10)]])], 'x_1_A_A')[0])
    expect(t.title).toBe('Three-way junction loop (J3)'); expect(t.counts).toBe('1 nucleotide · 1 segment')
  })
  it('explicit BGSU motif group, with its BGSU annotation and assignment basis; never the loop id', () => {
    const rec: LoopRecord = { ...base, motif_group: 'IL_19025.1', annotation_text: 'Isolated non-canonical cWW pair', motif_group_assignment: 'homologous', motif_group_via_loop_id: 'IL_Y_002', source: { motif_atlas_release: '4.15' } }
    const t = loopTooltip(instanceLoops([rec], 'x_1_A_A')[0])
    expect(t.motif).toBe('Motif: IL_19025.1 (Isolated non-canonical cWW pair)'); expect(t.motifNote).toBe('Motif Atlas 4.15 · homologous assignment via IL_Y_002')
    const bare = loopTooltip(instanceLoops([{ ...rec, annotation_text: null, motif_group_assignment: 'member', motif_group_via_loop_id: null }], 'x_1_A_A')[0])
    expect(bare.motif).toBe('Motif: IL_19025.1'); expect(bare.motifNote).toBe('Motif Atlas 4.15 · loop is a group member')
  })
  it('reports members in another chain', () => {
    const t = loopTooltip(instanceLoops([loop('IL_X_003', 'IL', [[nt(70, 30)], [nt(5, 50, 'match', 'x_1_C_C')]])], 'x_1_A_A')[0])
    expect(t.counts).toBe('2 nucleotides · 2 segments'); expect(t.otherChain).toBe('1 in another chain')
  })
})

describe('loop type names', () => {
  it('names every junction order, including higher-order junctions', () => {
    expect(typeName('J3')).toBe('Three-way junction loop'); expect(typeName('J5')).toBe('Five-way junction loop')
    expect(typeName('J7')).toBe('7-way junction loop'); expect(typeName('HL')).toBe('Hairpin loop'); expect(typeName('XX')).toBe('XX')
    expect(TRACK_OF('J5')).toBe('JL')
  })
})
