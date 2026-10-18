import { describe, expect, it } from 'vitest'
import { readFileSync, readdirSync } from 'node:fs'
import { boundaryText, insertionMarkers, insertionRuns, insertionSelection, markerWidth } from './insertions'
import type { ChainMapping, MappedResidue } from '../types/backend'

const DIR = new URL('../../../data/families/RF00167/mappings/', import.meta.url)
const load = (name: string) => JSON.parse(readFileSync(new URL(`${name}.json`, DIR), 'utf8')) as ChainMapping
const CLEN = 102

describe('RF00167 insertion events from the validated mappings', () => {
  it('9IWF A: one run of 2 insert-state nucleotides (CA) between consensus 77 and 78, with deposited identifiers', () => {
    const runs = insertionRuns(load('9iwf_1_A_A'), CLEN)
    expect(runs).toHaveLength(1)
    expect(runs[0]).toMatchObject({ kind: 'insert', anchor: 77, left: 77, right: 78, sequence: 'CA' })
    expect(runs[0].residues.map((r) => [r.label_seq_id, r.author_residue_number, r.author_insertion_code, r.offset, r.has_coordinates])).toEqual([[59, 60, '', 1, true], [60, 61, '', 2, true]])
  })
  it('9LJN A: same boundary, its own residue numbering (label 60-61, author 60-61); the leading flank G1 is not an insertion', () => {
    const m = load('9ljn_1_A_A'), runs = insertionRuns(m, CLEN)
    expect(runs.map((r) => [r.anchor, r.sequence, r.residues.map((x) => x.label_seq_id)])).toEqual([[77, 'CA', [60, 61]]])
    expect(m.residues.find((r) => r.polymer_seq_position === 1)!.alignment_state).toBe('not_submitted_flank')
  })
  it('construct extensions and terminal flanks (not_submitted_flank) never become insertions: 1Y26 X, 6P2H A, 9NBH A; no markers without insertions (4FE5 B)', () => {
    for (const sn of ['1y26_1_X_A', '6p2h_1_A_A', '9nbh_1_A_A', '4fe5_1_B_A']) expect(insertionRuns(load(sn), CLEN)).toEqual([])
  })
  it('every insert-state residue of RF00167 belongs to exactly one run, and the runs agree with the mapping summaries', () => {
    let nIns = 0, nRes = 0
    for (const f of readdirSync(DIR).filter((x) => x.endsWith('.json'))) {
      const m = load(f.slice(0, -5)), runs = insertionRuns(m, CLEN)
      nIns += m.residues.filter((r) => r.alignment_state === 'insert' || r.alignment_state === 'local_end').length
      nRes += runs.reduce((s, r) => s + r.residues.length, 0)
      const sm = m.summary as { n_insertion_sites: number; longest_insertion: number }
      expect(new Set(runs.map((r) => r.anchor)).size).toBe(sm.n_insertion_sites)
      expect(Math.max(0, ...runs.map((r) => r.residues.length))).toBe(sm.longest_insertion)
      for (const r of runs) expect(r.residues.every((x, i) => x.offset === i + 1)).toBe(true)
    }
    expect(nRes).toBe(nIns); expect(nIns).toBe(34)
  })
})

const res = (id: number, state: MappedResidue['alignment_state'], anchor: number | null, offset: number | null, nt = 'A', pos: number | null = null): MappedResidue => ({
  polymer_seq_position: id, author_residue_number: id + 100, author_insertion_code: id === 9 ? 'A' : '', residue_name: nt, nucleotide_parent: nt, is_modified: false,
  alignment_state: state, rfam_consensus_pos: pos, insertion_anchor: anchor, insertion_offset: offset, posterior_prob_char: '9', posterior_prob: 0.95, has_coordinates: id !== 7 })
const synth = (residues: MappedResidue[]) => ({ residues } as unknown as ChainMapping)

describe('insertion runs: synthetic edge cases', () => {
  const m = synth([
    res(1, 'insert', 0, 1, 'G'), res(2, 'insert', 0, 2, 'G'),                    // terminal insert before consensus 1
    res(3, 'match', null, null, 'C', 1), res(4, 'match', null, null, 'C', 2),
    res(5, 'insert', 2, 1, 'U'),                                                 // run 1 after 2
    res(6, 'local_end', 2, 1, 'A'), res(7, 'local_end', 2, 2, 'A'),              // different state at the same anchor: its own run
    res(8, 'match', null, null, 'G', 3),
    res(9, 'insert', 10, 1, 'C'),                                                // after a deletion of 4..10: anchor 10, not the last match
    res(11, 'insert', 10, 1, 'U'),                                               // numbering gap + offset restart: separate run, same anchor
    res(12, 'match', null, null, 'G', 11),
    res(13, 'insert', 11, 1, 'A'),                                               // adjacent boundary
    res(14, 'insert', 12, 1, 'A'), res(15, 'insert', 12, 2, 'A'), res(16, 'insert', 12, 3, 'A'), res(17, 'insert', 12, 4, 'A'), res(18, 'insert', 12, 5, 'A'), res(19, 'insert', 12, 6, 'A'), res(20, 'insert', 12, 7, 'A'), res(21, 'insert', 12, 8, 'A'), res(22, 'insert', 12, 9, 'A'), res(23, 'insert', 12, 10, 'A'),
    res(30, 'not_submitted_flank', null, null, 'G'),                             // flank: ignored
    res(31, 'insert', 102, 1, 'U'),                                              // terminal insert after the last position
  ])
  const runs = insertionRuns(m, CLEN)
  it('terminal runs get one-sided boundaries; no invented flanking position', () => {
    expect(runs[0]).toMatchObject({ anchor: 0, left: null, right: 1, sequence: 'GG' }); expect(boundaryText(runs[0])).toBe('Before Rfam position 1 (start of the model)')
    const last = runs[runs.length - 1]; expect(last).toMatchObject({ anchor: 102, left: 102, right: null }); expect(boundaryText(last)).toBe('After Rfam position 102 (end of the model)')
  })
  it('runs split on state, numbering gaps and offset restarts; never concatenated across them', () => {
    expect(runs.map((r) => [r.kind, r.anchor, r.sequence])).toEqual([['insert', 0, 'GG'], ['insert', 2, 'U'], ['local_end', 2, 'AA'], ['insert', 10, 'C'], ['insert', 10, 'U'], ['insert', 11, 'A'], ['insert', 12, 'AAAAAAAAAA'], ['insert', 102, 'U']])
    expect(runs[3].residues[0]).toMatchObject({ label_seq_id: 9, author_residue_number: 109, author_insertion_code: 'A' })
    expect(runs[2].residues.map((x) => x.has_coordinates)).toEqual([true, false])                                   // residue 7 has no coordinates
    expect(boundaryText(runs[3])).toBe('Between Rfam positions 10 and 11')
  })
  it('markers: one per boundary with its runs kept separate; overlapping labels at adjacent boundaries are aggregated; the run order is kept', () => {
    const mk = insertionMarkers(runs, 9)
    expect(mk.map((x) => [x.id, x.boundaries, x.length, x.runs.length])).toEqual([['0', [0], 2, 1], ['2', [2], 3, 2], ['10-12', [10, 11, 12], 13, 4], ['102', [102], 1, 1]])
    expect(insertionMarkers(runs, 40).map((x) => x.id)).toEqual(['0', '2', '10', '11', '12', '102'])   // wide cells: no collisions, nothing aggregated
    expect(markerWidth(2)).toBeLessThan(markerWidth(10))
  })
  it('selection group: instance-bound members without consensus positions', () => {
    const g = insertionSelection(runs[1], 'x', 'x_1_A_A')
    expect(g).toMatchObject({ kind: 'group', type: 'INS', id: 'insertion:x_1_A_A:2:5', source: { pdb_id: 'x', seq_name: 'x_1_A_A' } })
    expect(g.members).toEqual([{ unit_id: 'x_1_A_A:5', rfam_consensus_pos: null, seq_name: 'x_1_A_A', label_seq_id: 5 }])
  })
})
