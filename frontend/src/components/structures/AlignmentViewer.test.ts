import { describe, expect, it } from 'vitest'
import { BACKBONE, drawableNucleotides, residueLabel, tubeScale } from './AlignmentViewer'

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

describe('modified nucleotides are drawn as their parent base (display copy only)', () => {
  const at = (resn: string, resi: string, name: string, het = false) => ({ het, name, resn, chain: 'A', resi, icode: '', x: 0, y: 0, z: 0, occ: 1, b: 0, el: name[0] })
  it('infers the parent from base atoms; standard residues, ligands and order untouched', () => {
    const atoms = [at('G', '1', 'P'), at('PSU', '2', 'P'), at('PSU', '2', 'N1'), at('PSU', '2', 'O4'), at('2MG', '3', 'P'), at('2MG', '3', 'N9'), at('2MG', '3', 'O6'),
      at('1MA', '4', 'P'), at('1MA', '4', 'N9'), at('1MA', '4', 'N6'), at('5MC', '5', 'P'), at('5MC', '5', 'N4'), at('4SU', '6', 'S4'), at('MG', '90', 'MG', true)]
    const d = drawableNucleotides(atoms)
    expect(d.map((a) => a.resn)).toEqual(['G', 'U', 'U', 'U', 'G', 'G', 'G', 'A', 'A', 'A', 'C', 'C', 'U', 'MG'])
    expect(atoms[1].resn).toBe('PSU')                          // the original atoms (used for hover labels) keep the deposited name
    expect(d.map((a) => a.name)).toEqual(atoms.map((a) => a.name))
  })
})
