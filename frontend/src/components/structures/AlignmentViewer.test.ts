import { describe, expect, it } from 'vitest'
import { residueLabel } from './AlignmentViewer'

describe('3D view helpers', () => {
  it('residue hover label: PDB id, full author chain, residue name, author number and insertion code', () => {
    expect(residueLabel('5j7l', { chain: 'AA', resn: 'G', resi: '230', icode: '' })).toBe('5J7L.AA G230 (author)')
    expect(residueLabel('4v9d', { chain: 'BA', resn: 'PSU', resi: '516', icode: 'A' })).toBe('4V9D.BA PSU516A (author)')
  })
})
