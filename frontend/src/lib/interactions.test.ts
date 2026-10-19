import { readFileSync } from 'node:fs'
import { describe, expect, it } from 'vitest'
import { arcPath, baseCorrespondence, categoryAvailability, CATEGORY_STYLE, categoryColor, CIRCULAR_MAX_CLEN, defaultComparePair, DEFAULT_PAIR_MAX_RESOLUTION, helixCoverage, matchesPositions, medianHelixCoverage, multiColor, parseSearch, defaultFamilyOrder, groupByEndpoints, groupOccurrence, showCircular, circleLayout, compareInstances, indexInstance, kabsch, pairKey, pairLabel, passes, positionLabel, resolveLookup, statusIn, supportingInstances, DEFAULT_FILTER } from './interactions'
import { buildNumbering } from './numbering'
import type { ChainMapping, InteractionsDataset, StructureClusteringDataset } from '../types/backend'

const FAM = new URL('../../../data/families/RF00167/', import.meta.url)
const ds = JSON.parse(readFileSync(new URL('interactions.json', FAM), 'utf8')) as InteractionsDataset
const mapping = (sn: string) => JSON.parse(readFileSync(new URL(`mappings/${sn}.json`, FAM), 'utf8')) as ChainMapping
const numberingOf = (sn: string) => { const i = ds.instances.find((x) => x.seq_name === sn)!; return buildNumbering(sn, i.pdb_id, i.auth_asym_id, mapping(sn), ds.clen) }

// a synthetic 6-position family: instance x (all eligible), y (position 5 not modelled), z (annotation unavailable)
const toy = {
  clen: 6, instances: [
    { seq_name: 'x', pdb_id: 'X', auth_asym_id: 'A', label_asym_id: 'A', annotation: 'available', annotation_note: null, eligibility: 'EEEEEE', nucleotides: 'GGACUC', n_interactions: 3 },
    { seq_name: 'y', pdb_id: 'Y', auth_asym_id: 'A', label_asym_id: 'A', annotation: 'available', annotation_note: null, eligibility: 'EEEEME', nucleotides: 'GGACUC', n_interactions: 2 },
    { seq_name: 'z', pdb_id: 'Z', auth_asym_id: 'A', label_asym_id: 'A', annotation: 'unavailable', annotation_note: 'no BGSU FR3D interaction records for this chain', eligibility: 'EEEEEE', nucleotides: 'GGACUC', n_interactions: 0 },
  ],
  observations: {
    x: [['basepair', 'cWW', 1, 6, 1, 6, 'u', 'v', 'verified', ''], ['basepair', 'cWW', 2, 5, 2, 5, 'u', 'v', 'verified', ''], ['stacking', 's35', 1, 2, 1, 2, 'u', 'v', 'verified', ''], ['basepair', 'tWW', 3, 9, 3, '3^1', 'u', 'v', 'verified', '']],
    y: [['basepair', 'tWW', 1, 6, 1, 6, 'u', 'v', 'verified', ''], ['stacking', 'ns35', 1, 2, 1, 2, 'u', 'v', 'verified', '']],
    z: [],
  },
  pairs: [],
} as unknown as InteractionsDataset

describe('instance status', () => {
  const x = indexInstance(toy, 'x')!, y = indexInstance(toy, 'y')!, z = indexInstance(toy, 'z')!
  it('present / absent / not comparable, with reasons and alternative annotations', () => {
    expect(statusIn(x, 'basepair', 1, 6)).toMatchObject({ state: 'present', classes: ['cWW'] })
    expect(statusIn(y, 'basepair', 2, 5)).toMatchObject({ state: 'not_comparable' })
    expect(statusIn(y, 'basepair', 2, 5).reason).toContain('Rfam 5: nucleotide not modelled')
    expect(statusIn(z, 'basepair', 1, 6).reason).toContain('FR3D annotation unavailable')
    // near stacking is never "present"; it is listed as another annotation of an absent interaction
    expect(statusIn(y, 'stacking', 1, 2)).toMatchObject({ state: 'absent', other: ['Bs ns35'] })
    expect(statusIn(y, 'basepair', 1, 2)).toMatchObject({ state: 'absent', other: ['Bs ns35'] })
    expect(x.nOffConsensus).toBe(1)
  })
  it('comparison: reclassified (cWW vs tWW), only A, not comparable (missing nucleotide is never a loss)', () => {
    const rows = Object.fromEntries(compareInstances(x, y).map((r) => [r.id, r.status]))
    expect(rows).toEqual({ [pairKey('basepair', 1, 6)]: 'reclassified', [pairKey('basepair', 2, 5)]: 'not_comparable', [pairKey('stacking', 1, 2)]: 'only_a' })
    expect(compareInstances(x, x).every((r) => r.status === 'shared')).toBe(true)
    expect(compareInstances(x, z).every((r) => r.status === 'not_comparable')).toBe(true)
  })
  it('filters: base pairs by default; canonical and structure classes apply to base pairs only', () => {
    const bp = { category: 'basepair', structure_class: 'tertiary', region_1: 'L2', region_2: 'L3' } as const
    expect(passes(DEFAULT_FILTER, bp, false)).toBe(true)
    expect(passes(DEFAULT_FILTER, { ...bp, category: 'stacking', structure_class: null }, false)).toBe(false)
    expect(passes({ ...DEFAULT_FILTER, canonical: 'canonical' }, bp, false)).toBe(false)
    expect(passes({ ...DEFAULT_FILTER, structure: 'secondary' }, bp, true)).toBe(false)
    expect(passes({ ...DEFAULT_FILTER, region: 'L3' }, bp, true)).toBe(true)
    expect(passes({ ...DEFAULT_FILTER, region: 'P1' }, bp, true)).toBe(false)
  })
})

describe('RF00167 data', () => {
  const n4fe5 = numberingOf('4fe5_1_B_A')
  it('family-wide frequencies are the precomputed ones and do not depend on the numbering', () => {
    const p = ds.pairs.find((q) => q.id === 'basepair:41-66')!
    expect(pairLabel(41, 66, n4fe5)).toBe('G37–C61')
    expect(pairLabel(41, 66, null)).toBe('Rfam 41–Rfam 66')
    expect([p.n_observed, p.n_eligible]).toEqual([80, 80])
    // recount from the per-instance observations
    expect(ds.instances.filter((i) => statusIn(indexInstance(ds, i.seq_name)!, 'basepair', 41, 66).state === 'present').length).toBe(80)
  })
  it('lookup through the reference numbering and Rfam positions; explicit errors', () => {
    expect(resolveLookup('G37', n4fe5, ds.clen)).toEqual({ kind: 'position', positions: [41] })
    expect(resolveLookup('37', n4fe5, ds.clen)).toEqual({ kind: 'position', positions: [41] })
    expect(resolveLookup('G37-C61', n4fe5, ds.clen)).toEqual({ kind: 'pair', positions: [41, 66] })
    expect(resolveLookup('C61 G37', n4fe5, ds.clen)).toEqual({ kind: 'pair', positions: [41, 66] })
    expect(resolveLookup('Rfam 41', n4fe5, ds.clen)).toEqual({ kind: 'position', positions: [41] })
    expect(resolveLookup('41-66', null, ds.clen)).toEqual({ kind: 'pair', positions: [41, 66] })
    expect(resolveLookup('A37', n4fe5, ds.clen)).toMatchObject({ kind: 'error', message: expect.stringContaining('is G37') })
    expect(resolveLookup('9999', n4fe5, ds.clen)).toMatchObject({ kind: 'error' })
    expect(resolveLookup('Rfam 500', n4fe5, ds.clen)).toMatchObject({ kind: 'error' })
    expect(positionLabel(18, n4fe5)).toMatchObject({ fromReference: false, text: 'Rfam 18' })   // no 4FE5 B residue at Rfam 18
  })
  it('5E54 B (apo) vs 5SWE (adenine-bound): the Stagno et al. 2017 junction differences are recovered from the data', () => {
    const nb = numberingOf('5e54_1_B_B')
    const pos = (q: string) => (resolveLookup(q, nb, ds.clen) as { positions: [number, number] }).positions
    const rows = new Map(compareInstances(indexInstance(ds, '5e54_1_B_B')!, indexInstance(ds, '5swe_1_X_A')!).map((r) => [r.id, r]))
    const st = (q: string) => { const [i, j] = pos(q); return rows.get(pairKey('basepair', i, j))?.status }
    for (const q of ['U22-A52', 'U49-A76', 'C50-U75', 'A23-G46']) expect(st(q)).toBe('only_b')        // bound junction pairs
    for (const q of ['U47-A52', 'U49-A73', 'A23-U74', 'U48-U75']) expect(st(q)).toBe('only_a')        // apo pairs
    expect(st('U20-A76')).toBe('shared')
  })
  it('supporting instances are actual observations, preferred instance first', () => {
    const p = ds.pairs.find((q) => q.id === 'basepair:41-66')!
    const s = supportingInstances(ds, p, () => 2, ['4fe5_1_B_A'])
    expect(s).toHaveLength(80); expect(s[0].seqName).toBe('4fe5_1_B_A')
  })
})

describe('geometry', () => {
  it('deterministic clockwise layout; arcs end on their positions', () => {
    const L = circleLayout(1, 102, 500, 200)
    const [x, y] = L.point(1); expect(x).toBeGreaterThan(250); expect(y).toBeLessThan(60)       // 5' end just right of 12 o'clock
    expect(L.angle(50)).toBeGreaterThan(L.angle(10))
    const d = arcPath(L, 10, 60), [x1, y1] = L.point(10), [x2, y2] = L.point(60)
    expect(d.startsWith(`M${x1.toFixed(1)},${y1.toFixed(1)}A`)).toBe(true); expect(d.endsWith(`${x2.toFixed(1)},${y2.toFixed(1)}`)).toBe(true)
    expect(arcPath(circleLayout(1, 102, 500, 200), 3, 4)).toContain('A')
    // per category: its own colour family and line pattern; intensity darkens monotonically with frequency; rare interactions stay visible
    const rel = (h: string) => { const c = [1, 3, 5].map((o) => parseInt(h.slice(o, o + 2), 16) / 255).map((x) => (x <= 0.03928 ? x / 12.92 : ((x + 0.055) / 1.055) ** 2.4)); return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2] }
    const contrast = (h: string) => 1.05 / (rel(h) + 0.05)
    expect(Object.fromEntries(Object.entries(CATEGORY_STYLE).map(([k, v]) => [k, [v.family, v.pattern]]))).toEqual({ basepair: ['blue', 'solid'], stacking: ['green', 'dashed'], baseribose: ['orange', 'dotted'], basephosphate: ['purple', 'dash-dot'] })
    expect(new Set(Object.values(CATEGORY_STYLE).map((v) => v.dash ?? '')).size).toBe(4)
    for (const c of ['basepair', 'stacking', 'baseribose', 'basephosphate'] as const) {
      expect(categoryColor(c, 0)).toBe(CATEGORY_STYLE[c].ramp[0]); expect(categoryColor(c, 1)).toBe(CATEGORY_STYLE[c].ramp[4]); expect(categoryColor(c, null)).toBe(CATEGORY_STYLE[c].ramp[0])
      const ls = [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1].map((f) => rel(categoryColor(c, f)))
      expect(ls.every((v, k) => k === 0 || v < ls[k - 1])).toBe(true)
      expect(contrast(categoryColor(c, 0))).toBeGreaterThanOrEqual(2)          // lowest frequency still visible on white
    }
    expect(contrast(multiColor(0))).toBeGreaterThanOrEqual(2)
  })
  it('kabsch recovers a rigid transform; base correspondence by nucleotide identity', () => {
    const P: [number, number, number][] = [[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3], [1, 1, 1]]
    const c = Math.cos(0.7), s = Math.sin(0.7), Q = P.map(([x, y, z]) => [c * x - s * y + 5, s * x + c * y - 2, z + 1] as [number, number, number])
    const f = kabsch(P, Q)!; expect(f.rmsd).toBeLessThan(1e-6); expect(f.R[0][0]).toBeCloseTo(c, 6); expect(f.t[0]).toBeCloseTo(5, 6)
    const names = (s: string) => new Set(s.split(' '))
    const G = { names: names("N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4 C1'"), parent: 'G' }, A = { names: names("N9 C8 N7 C5 C6 N6 N1 C2 N3 C4 C1'"), parent: 'A' }
    const U = { names: names("N1 C2 O2 N3 C4 O4 C5 C6 C1'"), parent: 'U' }
    expect(baseCorrespondence(G, G)).toHaveLength(11); expect(baseCorrespondence(G, A)).toHaveLength(9)
    expect(baseCorrespondence(G, U)).toEqual([['N9', 'N1'], ['C4', 'C2'], ['C8', 'C6'], ["C1'", "C1'"]])
    expect(baseCorrespondence(U, G)).toEqual([['N1', 'N9'], ['C2', 'C4'], ['C6', 'C8'], ["C1'", "C1'"]])
  })
})

describe('unified search', () => {
  const n4 = numberingOf('4fe5_1_B_A'), terms = ['Helix P1', 'J1/2', 'J2/3', 'L2', 'L3', 'Base pairing']
  it('residue queries resolve through the selected numbering or Rfam positions; codes and regions are text', () => {
    expect(parseSearch('G37', n4, ds.clen, terms)).toEqual({ kind: 'residue', positions: [41] })
    expect(parseSearch('37', n4, ds.clen, terms)).toEqual({ kind: 'residue', positions: [41] })
    expect(parseSearch('G37-C61', n4, ds.clen, terms)).toEqual({ kind: 'residue', positions: [41, 66] })
    expect(parseSearch('Rfam 41', n4, ds.clen, terms)).toEqual({ kind: 'residue', positions: [41] })
    expect(parseSearch('Rfam 41-66', n4, ds.clen, terms)).toEqual({ kind: 'residue', positions: [41, 66] })
    expect(parseSearch('41-66', null, ds.clen, terms)).toEqual({ kind: 'residue', positions: [41, 66] })
    for (const t of ['cWW', 's35', '9BR', '0BPh', 'tHS']) expect(parseSearch(t, n4, ds.clen, terms)).toEqual({ kind: 'text', text: t.toLowerCase() })
    for (const t of ['J2/3', 'L2', 'P1', 'helix', 'base pairing']) expect(parseSearch(t, n4, ds.clen, terms).kind).toBe('text')
    expect(parseSearch('', n4, ds.clen, terms)).toEqual({ kind: 'none' })
  })
  it('invalid or ambiguous residues give a message, never a guess', () => {
    expect(parseSearch('A37', n4, ds.clen, terms)).toMatchObject({ kind: 'error', message: expect.stringContaining('is G37') })
    expect(parseSearch('9999', n4, ds.clen, terms)).toMatchObject({ kind: 'error', message: expect.stringContaining('No residue') })
    expect(parseSearch('Rfam 500', n4, ds.clen, terms)).toMatchObject({ kind: 'error' })
  })
  it('position matching: one position = involves it; two = exactly that pair in either direction', () => {
    expect(matchesPositions({ pos_1: 41, pos_2: 66 }, [41])).toBe(true); expect(matchesPositions({ pos_1: 66, pos_2: 41 }, [41, 66])).toBe(true)
    expect(matchesPositions({ pos_1: 41, pos_2: 67 }, [41, 66])).toBe(false)
  })
})

describe('four BGSU categories, grouping, ordering, cutoff, default pair', () => {
  // instance w: base of 3 complete, ribose of 1 complete, phosphate of 1 absent (5' end); v: ribose of 1 incomplete
  const dir = {
    clen: 4, summary: { by_category: { basepair: { n_entries_with_file: 1 }, stacking: { n_entries_with_file: 1 }, baseribose: { n_entries_with_file: 1 }, basephosphate: { n_entries_with_file: 0 } } },
    instances: [
      { seq_name: 'w', pdb_id: 'W', auth_asym_id: 'A', label_asym_id: 'A', annotation: 'available', annotation_note: null, eligibility: 'EEEE', sugar_eligibility: 'EEEE', phosphate_eligibility: 'IEEE', unavailable_categories: [], nucleotides: 'GACU', n_interactions: 2 },
      { seq_name: 'v', pdb_id: 'V', auth_asym_id: 'A', label_asym_id: 'A', annotation: 'available', annotation_note: null, eligibility: 'EEEE', sugar_eligibility: 'IEEE', phosphate_eligibility: 'EEEE', unavailable_categories: ['basephosphate'], nucleotides: 'GACU', n_interactions: 0 },
    ],
    observations: { w: [['baseribose', '0BR', 3, 1, 3, 1, 'u', 'v', 'verified', ''], ['basepair', 'tSS', 1, 3, 1, 3, 'u', 'v', 'verified', '']], v: [] },
    pairs: [
      { id: 'baseribose:3-1', category: 'baseribose', pos_1: 3, pos_2: 1, directional: true, frequency: 1, n_observed: 1, n_eligible: 1 },
      { id: 'basepair:1-3', category: 'basepair', pos_1: 1, pos_2: 3, directional: false, frequency: 0.5, n_observed: 1, n_eligible: 2 },
    ],
  } as unknown as InteractionsDataset
  it('category-specific eligibility: ribose / phosphate of the contacted nucleotide; unavailable categories are not comparable', () => {
    const w = indexInstance(dir, 'w')!, v = indexInstance(dir, 'v')!
    expect(statusIn(w, 'baseribose', 3, 1)).toMatchObject({ state: 'present', classes: ['0BR'] })
    expect(statusIn(w, 'baseribose', 1, 3)).toMatchObject({ state: 'absent' })                 // the other direction is a different interaction
    expect(statusIn(v, 'baseribose', 3, 1)).toMatchObject({ state: 'not_comparable', reason: 'Rfam 1: incomplete ribose' })
    expect(statusIn(w, 'basephosphate', 3, 1).reason).toBe('Rfam 1: incomplete or absent phosphate')
    expect(statusIn(v, 'basephosphate', 3, 1).reason).toContain('unavailable')
    // the base pair between the same nucleotides lists the base-ribose contact as another annotation, and vice versa
    expect(statusIn(w, 'basepair', 1, 3).other).toEqual(['Br 0BR'])
    expect(categoryAvailability(dir, 'basephosphate')).toBe('none'); expect(categoryAvailability(dir, 'stacking')).toBe('empty'); expect(categoryAvailability(dir, 'baseribose')).toBe('ok')
  })
  it('multiple interactions: grouped by unordered endpoints, members kept; explicit occurrence (not a sum or average)', () => {
    const g = groupByEndpoints(dir.pairs); expect(g).toHaveLength(1); expect(g[0].members.map((m) => m.id)).toEqual(['baseribose:3-1', 'basepair:1-3'])
    // w: eligible for both, present -> counted; v: eligible for the base pair only, absent -> in the denominator only
    expect(groupOccurrence(dir, dir.pairs)).toEqual({ k: 1, n: 2, frequency: 0.5 })
  })
  it('default family order: entries desc, then frequency desc, eligible desc, positions', () => {
    // entries first; then frequency, eligible chains, positions
    const rows = [{ entries: 2, frequency: 1, n: 4, i: 1, j: 9 }, { entries: 13, frequency: 0.16, n: 80, i: 2, j: 8 }, { entries: 65, frequency: 1, n: 80, i: 5, j: 6 }, { entries: 60, frequency: 0.91, n: 80, i: 3, j: 7 }, { entries: 65, frequency: 1, n: 80, i: 4, j: 6 }, { entries: 65, frequency: 0.9, n: 80, i: 1, j: 2 }]
    expect([...rows].sort(defaultFamilyOrder).map((r) => `${r.entries}/${r.frequency}/${r.i}`)).toEqual(['65/1/4', '65/1/5', '65/0.9/1', '60/0.91/3', '13/0.16/2', '2/1/1'])
  })
  it('circular diagrams up to 300 consensus positions', () => {
    expect(CIRCULAR_MAX_CLEN).toBe(300); expect(showCircular(299)).toBe(true); expect(showCircular(300)).toBe(true); expect(showCircular(301)).toBe(false)
  })
  it('default pair: highest valid discrepancy; null values skipped; ties by stable ids; unannotated or non-comparable pairs excluded', () => {
    const sc = (ids: string[], upper: (number | null)[]) => ({ dataset: { instances: ids }, metrics: { discrepancy: { values_upper: upper } } }) as unknown as StructureClusteringDataset
    const res2 = () => 2
    expect(defaultComparePair(toy, sc(['x', 'y', 'z'], [0.2, null, 0.3]), res2)).toEqual({ a: 'x', b: 'y', value: 0.2 })   // y-z has no annotation in z, x-z null
    expect(defaultComparePair(toy, sc(['y', 'x'], [0.5]), res2)).toEqual({ a: 'x', b: 'y', value: 0.5 })                  // A = smaller id
    // resolution: both chains must be <= 3.0 Å; a chain without a resolution is excluded
    expect(DEFAULT_PAIR_MAX_RESOLUTION).toBe(3)
    expect(defaultComparePair(toy, sc(['x', 'y'], [0.5]), (sn) => (sn === 'y' ? 3.25 : 2))).toBeNull()
    expect(defaultComparePair(toy, sc(['x', 'y'], [0.5]), (sn) => (sn === 'y' ? 3.0 : 2))).toEqual({ a: 'x', b: 'y', value: 0.5 })
    expect(defaultComparePair(toy, sc(['x', 'y'], [0.5]))).toBeNull()
    const four = { ...toy, instances: [...toy.instances, { ...toy.instances[0], seq_name: 'w' }], observations: { ...toy.observations, w: toy.observations!.x } } as unknown as InteractionsDataset
    expect(defaultComparePair(four, sc(['w', 'x', 'y'], [0.4, 0.4, 0.1]), res2)).toEqual({ a: 'w', b: 'x', value: 0.4 })   // tie: (w,x) before (w,y)
    expect(defaultComparePair(toy, null)).toBeNull()
  })
  it('RF00167: the default pair is the largest discrepancy among chains <= 3.0 Å (not the 3.25 Å 9V50 A of the unrestricted maximum)', () => {
    const sc = JSON.parse(readFileSync(new URL('structure_clustering.json', FAM), 'utf8')) as StructureClusteringDataset
    const rows = JSON.parse(readFileSync(new URL('instances.json', FAM), 'utf8')) as { seq_name: string; resolution: number | null }[]
    const res = (sn: string) => rows.find((r) => r.seq_name === sn)?.resolution ?? null
    const ids = (sc.dataset as unknown as { instances: string[] }).instances, up = (sc.metrics.discrepancy as unknown as { values_upper: (number | null)[] }).values_upper
    let best = -1, k = 0
    for (let r = 0; r < ids.length; r++) for (let c = r + 1; c < ids.length; c++, k++) { const v = up[k]; if (v != null && (res(ids[r]) ?? 99) <= 3 && (res(ids[c]) ?? 99) <= 3) best = Math.max(best, v) }
    const p = defaultComparePair(ds, sc, res)!
    expect(p.value).toBe(best); expect(p).toEqual({ a: '5e54_1_B_B', b: '9ljn_1_A_A', value: best })
    expect(defaultComparePair(ds, sc, res, null)!.b).toBe('9v50_1_A_A')                 // without the resolution rule
  })
  it('helix coverage: consensus helix base pairs annotated / eligible per chain; the family median', () => {
    expect(helixCoverage(ds, '5e54_1_B_B')).toEqual({ k: 19, n: 19, fraction: 1 })
    expect(helixCoverage(ds, '9v50_1_A_A')).toEqual({ k: 13, n: 20, fraction: 0.65 })
    expect(medianHelixCoverage(ds)).toBe(0.95)
  })
  it('RF00167: base-pair and stacking pairs are unchanged; four categories present', () => {
    expect(ds.pairs.filter((q) => q.category === 'basepair' || q.category === 'stacking')).toHaveLength(199)
    for (const c of ['basepair', 'stacking', 'baseribose', 'basephosphate'] as const) expect(categoryAvailability(ds, c)).toBe('ok')
    expect(ds.pairs.filter((q) => q.directional).every((q) => q.category === 'baseribose' || q.category === 'basephosphate')).toBe(true)
  })
})
