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

import { pageList, PAGE_SIZE } from './LigandsTab'
import { checkIdealSdf, idealSdfUrls, principalAxes } from '../lib/idealLigand'
describe('Ligands tab: gallery pages and ideal 3D coordinates', () => {
  it('9 per page; page numbers with gaps only when there are more than 7 pages', () => {
    expect(PAGE_SIZE).toBe(9)
    expect(pageList(1, 4)).toEqual([1, 2, 3, 4])
    expect(pageList(5, 10)).toEqual([1, '…', 4, 5, 6, '…', 10]); expect(pageList(1, 10)).toEqual([1, 2, '…', 10]); expect(pageList(10, 10)).toEqual([1, '…', 9, 10])
  })
  const mol = (atoms: number[][]) => ['X', '', '', `${String(atoms.length).padStart(3)}  0  0  0  0  0  0  0  0  0999 V2000`, ...atoms.map(([x, y, z]) => `${x.toFixed(4).padStart(10)}${y.toFixed(4).padStart(10)}${z.toFixed(4).padStart(10)} C   0  0`), 'M  END'].join('\n')
  it('accepts an ideal molfile, reports all-zero coordinates as missing and anything else as invalid', () => {
    expect(checkIdealSdf(mol([[0, 0, 0], [1.4, 0, 0], [0, 1.4, 0]]))).toMatchObject({ ok: true, nAtoms: 3 })
    expect(checkIdealSdf(mol([[0, 0, 0], [0, 0, 0]]))).toEqual({ ok: false, reason: 'no_ideal_coordinates' })
    expect(checkIdealSdf('<html>')).toEqual({ ok: false, reason: 'invalid' }); expect(checkIdealSdf(mol([[0, 0, 0]]).replace('M  END', ''))).toEqual({ ok: false, reason: 'invalid' })
    const v3000 = ['X', '', '', '  0  0  0     0  0            999 V3000', 'M  V30 BEGIN CTAB', 'M  V30 COUNTS 2 1 0 0 0', 'M  V30 BEGIN ATOM', 'M  V30 1 C 0 0 0 0', 'M  V30 2 O 1.2 0 0 0', 'M  V30 END ATOM', 'M  V30 END CTAB', 'M  END'].join('\n')
    expect(checkIdealSdf(v3000)).toMatchObject({ ok: true, nAtoms: 2 })
  })
  it('principal axes: a planar molecule is viewed along its normal, its longest extent first', () => {
    const { centre, axes } = principalAxes([[0, 0, 5], [4, 0, 5], [0, 1, 5], [4, 1, 5], [2, 0.5, 5]])
    expect(centre.map((v) => +v.toFixed(6))).toEqual([2, 0.5, 5])
    expect(Math.abs(axes[0][0])).toBeCloseTo(1, 5); expect(Math.abs(axes[1][1])).toBeCloseTo(1, 5); expect(Math.abs(axes[2][2])).toBeCloseTo(1, 5)
  })
  it('PDBe ideal SDF first, the identical RCSB file only as fallback', () => {
    expect(idealSdfUrls('ADE')).toEqual(['https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/ADE_ideal.sdf', 'https://files.rcsb.org/ligands/download/ADE_ideal.sdf'])
  })
})

import { seqDiffLabel, seqDiffTracks } from './AnnotationsTab'
import type { SequenceDifferenceTracksDataset } from '../types/backend'
describe('reported sequence differences (mmCIF _struct_ref_seq_dif)', () => {
  const o = (pdb: string, ch: string, ref: string | null, dep: string | null, auth: string, kind: 'substitution' | 'added_residue' | 'absent_residue', cls: string) =>
    ({ pdb_id: pdb, auth_asym_id: ch, label_asym_id: ch, seq_name: `${pdb}_1_${ch}_${ch}`, label_seq_id: 1, author_residue_number: auth, author_insertion_code: '', reference: ref, deposited: dep, change_kind: kind, classification: cls })
  const sd = {
    rfam_acc: 'RFX', clen: 3, category: 'Reported sequence differences',
    classifications: { 'Engineered mutation': { n_observations: 3, n_on_consensus_axis: 3, n_instances: 3, n_entries: 2, change_kinds: {}, consensus_positions: [2] },
                       Conflict: { n_observations: 1, n_on_consensus_axis: 1, n_instances: 1, n_entries: 1, change_kinds: {}, consensus_positions: [3] },
                       'SEE REMARK 999': { n_observations: 2, n_on_consensus_axis: 0, n_instances: 1, n_entries: 1, change_kinds: {}, consensus_positions: [] } },
    positions: [{ rfam_consensus_pos: 1, classifications: {} },
                { rfam_consensus_pos: 2, classifications: { 'Engineered mutation': { n_observations: 3, n_instances: 3, n_entries: 2, changes: { 'C>U': 2, 'U>C': 1 }, change_kinds: { substitution: 3 }, observation_ids: ['a', 'b', 'c'] } } },
                { rfam_consensus_pos: 3, classifications: { Conflict: { n_observations: 1, n_instances: 1, n_entries: 1, changes: { 'A>G': 1 }, change_kinds: { substitution: 1 }, observation_ids: ['d'] } } }],
    off_axis: [], counts: { n_entries: 3, n_entries_evaluated: 3, n_instances_evaluated: 4, n_observations: 6, n_on_consensus_axis: 4, n_off_axis: 2, n_unresolved_rows: 0 },
    observations: { a: o('3fo4', 'A', 'C', 'U', '74', 'substitution', 'Engineered mutation'), b: o('2b57', 'A', 'C', 'U', '74', 'substitution', 'Engineered mutation'),
                    c: o('2b57', 'B', 'U', 'C', '74', 'substitution', 'Engineered mutation'), d: o('5c7u', 'B', 'A', 'G', '15', 'substitution', 'Conflict') },
  } as unknown as SequenceDifferenceTracksDataset
  it('one track per classification on the axis (verbatim labels), most instances first; fixed colour per label', () => {
    const t = seqDiffTracks(sd)
    expect(t.map((x) => x.label)).toEqual(['Engineered mutation', 'Conflict'])                    // SEE REMARK 999 is off-axis only: no track
    expect(t[0].scale.dark).toBe('#6d28d9'); expect(t[1].scale.dark).toBe('#b45309'); expect(t[0].scale.kind).toBe('presence')
  })
  it('marks only where a difference is recorded (absence is not a negative); tooltip, details and per-chain observations from the data', () => {
    const [em] = seqDiffTracks(sd)
    expect(em.at(1)).toBeNull(); expect(em.at(3)).toBeNull()
    const m = em.at(2) as Exclude<ReturnType<typeof em.at>, null | 'unevaluable'>
    expect(m.tip.lines).toEqual(['Engineered mutation: C→U ×2, U→C', 'Reported in 3 RNA instances (2 PDB entries)', '3FO4 A C74U', '2B57 A C74U', '2B57 B U74C'])
    expect(m.observations).toEqual([{ label: '3FO4 A C74U', seqName: '3fo4_1_A_A' }, { label: '2B57 A C74U', seqName: '2b57_1_A_A' }, { label: '2B57 B U74C', seqName: '2b57_1_B_B' }])
    expect(m.summary.find(([k]) => k === 'Reference')![1]).toContain('not the Rfam consensus')
  })
  it('labels substitutions, added and absent residues', () => {
    expect(seqDiffLabel(o('2b57', 'A', null, 'G', '16', 'added_residue', 'Engineered mutation'))).toBe('2B57 A G16 added')
    expect(seqDiffLabel(o('4kqy', 'A', 'U', null, '', 'absent_residue', 'deletion'))).toBe('4KQY A U absent')
    expect(seqDiffLabel({ ...o('9xyz', 'A', 'C', 'U', '74', 'substitution', 'x'), author_insertion_code: 'A' })).toBe('9XYZ A C74A→U')
    expect(seqDiffLabel(o('1ehz', 'A', 'U', 'H2U', '16', 'substitution', 'modified residue'))).toBe('1EHZ A U16→H2U')
  })
})
