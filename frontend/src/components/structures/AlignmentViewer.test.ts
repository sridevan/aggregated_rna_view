import { describe, expect, it } from 'vitest'
import { BACKBONE, residueLabel, tubeScale } from './AlignmentViewer'

describe('3D alignment view shows the overall chain geometry only', () => {
  it('passes only sugar-phosphate backbone atoms to 3Dmol (no base atoms, so no base ladders)', () => {
    for (const a of ['P', "C4'", "C1'", "O3'"]) expect(BACKBONE.has(a)).toBe(true)
    for (const a of ['N1', 'N3', 'N9', 'C2', 'C4', 'C5', 'C6', 'C8', 'O6', 'N6', 'N2', 'O2', 'O4', 'N4']) expect(BACKBONE.has(a)).toBe(false)
  })
  it('thicker tubes for large chains', () => {
    expect(tubeScale(70)).toBe(1); expect(tubeScale(500)).toBeGreaterThan(1); expect(tubeScale(1500)).toBeGreaterThan(tubeScale(500))
  })
})

describe('residue hover label', () => {
  it('PDB id upper case, author chain, residue name and author number, insertion code appended', () => {
    expect(residueLabel('5j7l', { chain: 'AA', resn: 'G', resi: '230', icode: '' })).toBe('5J7L.AA G230 (author)')
    expect(residueLabel('4v9d', { chain: 'BA', resn: 'PSU', resi: '516', icode: 'A' })).toBe('4V9D.BA PSU516A (author)')
  })
})
