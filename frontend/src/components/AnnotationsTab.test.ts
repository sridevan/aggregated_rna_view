import { describe, expect, it } from 'vitest'
import { pct, fitCellPx, tickStepFor, visibleRange } from './AnnotationsTab'

describe('pct (tooltip percentage formatting, display only)', () => {
  it('uses at most one decimal and drops a trailing .0', () => {
    expect(pct(66, 80)).toBe('82.5%')
    expect(pct(61, 80)).toBe('76.3%')
    expect(pct(12, 12)).toBe('100%')
    expect(pct(76, 80)).toBe('95%')
  })
  it('never divides by zero', () => { expect(pct(0, 0)).toBe('0%') })
})

import { scaleColor, type ScaleDef } from './AnnotationsTab'
describe('presence scale (sparse tracks)', () => {
  const s: ScaleDef = { id: 'p', kind: 'presence', label: 'x', light: '#c2185b', dark: '#c2185b', bins: 1, binLabels: ['observed'] }
  it('uses one fixed colour whatever the frequency', () => {
    expect(scaleColor(s, 0.012)).toBe('#c2185b'); expect(scaleColor(s, 0.026)).toBe('#c2185b'); expect(scaleColor(s, 1)).toBe('#c2185b')
  })
  it('leaves the frequency scales unchanged', () => {
    const f: ScaleDef = { id: 'f', kind: 'frequency', label: 'f', light: '#000000', dark: '#ffffff', bins: 5, binLabels: [] }
    expect(scaleColor(f, 0.012)).toBe('rgb(0,0,0)'); expect(scaleColor(f, 0.95)).toBe('rgb(255,255,255)')
  })
})

describe('annotation zoom window', () => {
  it('Fit all spreads the model over the plot width', () => {
    expect(fitCellPx(1020, 102)).toBe(10)
    expect(fitCellPx(0, 102)).toBe(9)          // before measurement
  })
  it('axis label step keeps labels >= 60 px apart', () => {
    expect(tickStepFor(10)).toBe(10)           // RF00167 at Fit all
    expect(tickStepFor(0.7)).toBe(100)         // 16S at Fit all
    expect(tickStepFor(32)).toBe(2)
    expect(tickStepFor(0.01)).toBe(1000)
  })
  it('visible range is clamped to the model and covers partial cells', () => {
    expect(visibleRange(0, 1000, 10, 102)).toEqual([1, 100])
    expect(visibleRange(15, 100, 10, 102)).toEqual([2, 12])
    expect(visibleRange(0, 2000, 10, 102)).toEqual([1, 102])
  })
})

import { defaultLigandCcds } from './AnnotationsTab'
import { depictionUrls, galleryItems, ligandPageUrl } from './LigandsTab'
import type { LigandsDataset } from '../types/backend'
describe('Ligands tab: default ligand set, gallery search and CCD depiction URLs', () => {
  const s = (name: string, n: number) => ({ name, ligand_class: 'ligand', n_entry_entities: n, n_entries: n, entries: [], n_bound_molecule_instances: n, n_entries_in_default_track: n, consensus_positions_n_entries: {}, flag_counts: {} })
  const bm = (pdb_id: string, ccd_id: string, in_default_track = true) => ({ pdb_id, entity_id: 1, ccd_id, name: null, ligand_class: 'ligand', chain_attribution: 'resolved', n_entity_chains: 1, consensus_positions: [25], flags: [], in_default_track, exclusion_reason: null })
  const lig = { by_ccd: { ADE: s('ADENINE', 2), HPA: s('HYPOXANTHINE', 2), GUN: s('GUANINE', 5), GTP: { ...s("GUANOSINE-5'-TRIPHOSPHATE", 0), n_entries: 4 }, MG: { ...s('MAGNESIUM ION', 0), ligand_class: 'metal_ion' } },
    bound_molecules: [bm('1abc', 'ADE'), bm('2abc', 'ADE'), bm('2abc', 'ADE'), bm('3abc', 'HPA'), bm('4abc', 'HPA'), bm('5abc', 'GUN'), bm('6abc', 'GTP', false)] } as unknown as LigandsDataset
  it('the gallery is the default track set, in track order (entries desc, then CCD id), counting distinct PDB entries only', () => {
    expect(defaultLigandCcds(lig).map(([c]) => c)).toEqual(['GUN', 'ADE', 'HPA'])
    const g = galleryItems(lig)
    expect(g.map((x) => x.ccd)).toEqual(['GUN', 'ADE', 'HPA'])
    expect(g.find((x) => x.ccd === 'ADE')).toEqual({ ccd: 'ADE', name: 'Adenine', nEntries: 2, entries: ['1abc', '2abc'] })
  })
  it('search matches the CCD id or the name, case-insensitively; blank shows all', () => {
    expect(galleryItems(lig, 'hpa').map((x) => x.ccd)).toEqual(['HPA'])
    expect(galleryItems(lig, 'GUAN').map((x) => x.ccd)).toEqual(['GUN'])
    expect(galleryItems(lig, '  ').length).toBe(3); expect(galleryItems(lig, 'xyz')).toEqual([])
  })
  it('PDBe chemistry image first, RCSB CCD image as fallback; PDBe ligand page', () => {
    expect(depictionUrls('HPA')).toEqual(['https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/HPA_200.svg', 'https://cdn.rcsb.org/images/ccd/unlabeled/H/HPA.svg'])
    expect(depictionUrls('A1JAI')[1]).toBe('https://cdn.rcsb.org/images/ccd/unlabeled/A/A1JAI.svg')
    expect(ligandPageUrl('6AP')).toBe('https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/6AP')
  })
})
