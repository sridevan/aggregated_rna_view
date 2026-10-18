import { describe, expect, it } from 'vitest'
import type { ChainMapping, MappedResidue } from '../types/backend'
import { indexMapping } from './consensus'
import { authorResidueId, buildNumbering, chainAt, referenceTicks, residueValue, STATUS_TEXT } from './numbering'

/* SYNTHETIC TEST FIXTURE (not real data): a 10-position model and one chain with
 *   pos 1-2 outside the aligned region, 3 = G100, 4 = A101 (no coordinates), 5 deleted, 6 = U102A (insertion code), an Rfam insertion
 *   residue C103 between 6 and 7, 7 = G110 (author numbering gap), 8 = PSU111 (modified), 9 = A (no author number), 10 outside. */
const res = (p: number, num: number | null, icode: string, name: string, state: MappedResidue['alignment_state'], pos: number | null, coords = true): MappedResidue => ({
  polymer_seq_position: p, author_residue_number: num, author_insertion_code: icode, residue_name: name, nucleotide_parent: name === 'PSU' ? 'U' : name, is_modified: name === 'PSU',
  alignment_state: state, rfam_consensus_pos: pos, insertion_anchor: state === 'insert' ? 6 : null, insertion_offset: state === 'insert' ? 1 : null, posterior_prob_char: '*', posterior_prob: 1, has_coordinates: coords,
})
const residues = [res(1, 100, '', 'G', 'match', 3), res(2, 101, '', 'A', 'match', 4, false), res(3, 102, 'A', 'U', 'match', 6), res(4, 103, '', 'C', 'insert', null),
  res(5, 110, '', 'G', 'match', 7), res(6, 111, '', 'PSU', 'match', 8), res(7, null, '', 'A', 'match', 9)]
const cols = [[1, 'uncovered', null], [2, 'uncovered', null], [3, 'match', 1], [4, 'match', 2], [5, 'deletion', null], [6, 'match', 3], [7, 'match', 5], [8, 'match', 6], [9, 'match', 7], [10, 'uncovered', null]] as const
const fixture = { pdb_id: 'test', entity_id: 1, label_asym_id: 'A', auth_asym_id: 'X', rfam_acc: 'RF99999', sequence_key: 'k', mapped_region: { start: 1, end: 7 }, alignment: null,
  residues, consensus_columns: cols.map(([p, s, q]) => ({ rfam_consensus_pos: p, state: s, polymer_seq_position: q })), provenance: {} } as unknown as ChainMapping

describe('position numbering from explicit mappings (synthetic fixture)', () => {
  const idx = indexMapping(fixture)
  it('author identifiers: residue name + author number + insertion code; never invented', () => {
    expect(authorResidueId(residues[0])).toBe('G100')
    expect(authorResidueId(residues[2])).toBe('U102A')
    expect(authorResidueId(residues[5])).toBe('PSU111')
    expect(authorResidueId(residues[6])).toBeNull()
  })
  it('states A-D are distinct: coordinates, no coordinates, deletion, outside the aligned region', () => {
    expect([3, 4, 5, 1].map((p) => chainAt(idx, p).state)).toEqual(['match_observed', 'match_unobserved', 'deletion', 'uncovered'])
    expect([3, 4, 5, 1, 9].map((p) => residueValue(chainAt(idx, p)))).toEqual(['G100', 'A101', 'Not present', 'Unmapped', 'Unavailable'])
    expect(STATUS_TEXT.deletion).toBe('Deletion'); expect(STATUS_TEXT.uncovered).toBe('Outside aligned region'); expect(STATUS_TEXT.match_unobserved).toBe('Coordinates unavailable')
  })
  it('numbering keeps every Rfam column, labels only mapped positions, ignores Rfam insertion residues and author-number gaps', () => {
    const n = buildNumbering('test_1_A_X', 'test', 'X', fixture, 10)
    expect(n.label).toBe('TEST X')
    expect(n.numbers).toEqual([null, null, '100', '101', null, '102A', '110', '111', null, null])
    expect(n.numbers).toHaveLength(10)
    expect(n.numbers).not.toContain('103')          // the inserted residue has no Rfam position and gets no column
  })
})

describe('reference axis ticks', () => {
  it('round author numbers at the zoom step, real residues only, no overlap', () => {
    const numbers = Array.from({ length: 200 }, (_, k) => (k < 10 ? null : String(k + 31)))   // positions 11.. numbered 42..
    const t = referenceTicks(numbers, 10)
    expect(t.length).toBeGreaterThan(5)
    for (const x of t) expect(x.text).toBe(numbers[x.pos - 1])
    for (let i = 1; i < t.length; i++) expect((t[i].pos - t[i - 1].pos) * 10).toBeGreaterThanOrEqual(48)
    expect(t.filter((x) => +x.text % 10 === 0).length).toBeGreaterThanOrEqual(t.length - 2)
  })
  it('non-monotonic numbering and insertion codes are labels, not arithmetic', () => {
    const numbers = ['306', '306A', '306B', '306C', '307', null, '5', '6', '7', '8', '9', '10', ...Array(30).fill(null), '1000']
    const t = referenceTicks(numbers, 30)
    for (const x of t) expect(x.text).toBe(numbers[x.pos - 1])
    expect(t.map((x) => x.text)).toContain('1000')
  })
  it('fit-all of a long chain stays sparse', () => {
    const numbers = Array.from({ length: 1533 }, (_, k) => String(k + 1))
    const t = referenceTicks(numbers, 0.7)
    expect(t.length).toBeLessThan(20); expect(t.every((x) => +x.text % 100 === 0)).toBe(true)
  })
})
