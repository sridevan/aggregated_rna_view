import { readFileSync, existsSync } from 'node:fs'
import { describe, expect, it } from 'vitest'
import { buildTree, clusteringFor, GROUP_COLORS, squareMatrix, UNPLACED, viridis, instanceLabel } from './structureGroups'
import { categoricalRows, medianOf } from '../components/structures/AssociatedFeatures'
import type { StructureClusteringDataset } from '../types/backend'

describe('tree and matrix helpers', () => {
  it('rebuilds the drawing order from a scipy-style linkage', () => {
    const t = buildTree([[0, 1, 1, 2], [3, 2, 4, 3]], 3)
    expect(t.leaves).toEqual([0, 1, 2]); expect(t.height).toBe(4); expect(t.left!.leaves).toEqual([0, 1])
  })
  it('expands the upper triangle with the right diagonal', () => {
    expect(squareMatrix([0.5, 0.7, 0.9], 3, 1)).toEqual([[1, 0.5, 0.7], [0.5, 1, 0.9], [0.7, 0.9, 1]])
  })
  it('viridis ends and labels', () => {
    expect(viridis(0)).toBe('rgb(68,1,84)'); expect(viridis(1)).toBe('rgb(253,231,37)'); expect(instanceLabel('4fe5_1_B_A')).toBe('4FE5 B')
  })
})

describe('associated features counting', () => {
  const cols = [{ id: 'g1', label: 'G1', members: ['a', 'b', 'c'], color: '#000', selected: true }, { id: 'g2', label: 'G2', members: ['d', 'e'], color: '#111', selected: false }]
  it('counts instances once per feature, with known-status denominators, in a stable order', () => {
    const f: Record<string, string[]> = { a: ['ADE', 'ADE'], b: ['ADE', 'GUN'], c: [], d: ['GUN'], e: ['HPA'] }
    const known = (sn: string) => sn !== 'c'                                     // c: status unknown -> not a negative, not in n
    const rows = categoricalRows(cols, known, (sn) => new Set(f[sn]))
    expect(rows.map((r) => [r.key, r.cells.map((x) => `${x.k}/${x.n}`)])).toEqual([['ADE', ['2/2', '0/2']], ['GUN', ['1/2', '1/2']], ['HPA', ['0/2', '1/2']]])
  })
  it('medians skip missing values and report how many instances contributed', () => {
    expect(medianOf(['a', 'b', 'c', 'd'], (s) => ({ a: 2, b: null, c: 4, d: 3 } as Record<string, number | null>)[s])).toEqual({ median: 3, n: 3 })
    expect(medianOf(['a'], () => null)).toEqual({ median: null, n: 0 })
  })
})

const FILE = new URL('../../public/families/RF00167/structure_clustering.json', import.meta.url)
describe.skipIf(!existsSync(FILE))('RF00167 clustering file', () => {
  const d = JSON.parse(readFileSync(FILE, 'utf8')) as StructureClusteringDataset
  it('discrepancy: 72 + 4, group colours, 1 - TM only for the TM metric', () => {
    const c = clusteringFor(d, 'discrepancy')
    expect(c.groups.map((g) => g.n_instances)).toEqual([72, 4]); expect(c.groups.map((g) => g.color)).toEqual([GROUP_COLORS[0], GROUP_COLORS[1]])
    expect(c.order.length).toBe(76); expect(new Set(c.order).size).toBe(76)
    expect(c.dist[0][1]).toBe(c.value[0][1])
    const t = clusteringFor(d, 'tm')
    expect(t.groups.map((g) => g.n_instances)).toEqual([51, 17, 4, 4]); expect(t.dist[0][1]).toBeCloseTo(1 - t.value[0][1], 12); expect(t.value[3][3]).toBe(1)
    expect(t.groups.every((g) => g.color !== UNPLACED)).toBe(true)
    for (const g of t.groups) expect(t.order.slice(g.from, g.to + 1)).toEqual(g.members)   // groups are contiguous in leaf order
  })
})

import { groupShade, pathBetween, structureColours, viewerColours, SUPERPOSE_COLOURS } from './structureGroups'
describe('comparison route and structure colours', () => {
  it('route goes up to the lowest common ancestor and down again', () => {
    const t = buildTree([[0, 1, 1, 2], [2, 3, 1.5, 2], [4, 5, 3, 4]], 4)          // ((0,1),(2,3))
    expect([...pathBetween(t, 4, 5)!].sort()).toEqual([4, 5, 6]); expect([...pathBetween(t, 0, 1)!].sort()).toEqual([0, 1, 4])
    expect(pathBetween(t, 0, 99)).toBeNull()
  })
  it('representatives take the exact group colour; members distinct shades of their own group', () => {
    const c = { groups: [{ group_id: 1, representative: 'a', color: '#e8685a' }, { group_id: 2, representative: 'x', color: '#b0327a' }], groupOf: { a: 1, b: 1, c: 1, x: 2, y: 2 } } as never
    const col = structureColours(['a', 'b', 'c', 'x', 'y'], c)
    expect(col.a).toBe('#e8685a'); expect(col.x).toBe('#b0327a')
    expect(new Set(Object.values(col)).size).toBe(5); expect(col.b).toBe(groupShade('#e8685a', 0)); expect(col.y).toBe(groupShade('#b0327a', 0))
  })
  it('viewer: one structure keeps its group colour; a superposed pair always takes the fixed contrasting pair, also within one group', () => {
    const c = { groups: [{ group_id: 1, representative: 'a', color: '#e8685a' }, { group_id: 2, representative: 'x', color: '#b0327a' }], groupOf: { a: 1, b: 1, x: 2 } } as never
    expect(viewerColours(['b'], c)).toEqual({ b: groupShade('#e8685a', 0) }); expect(viewerColours(['a'], c)).toEqual({ a: '#e8685a' })
    expect(viewerColours(['x', 'a'], c)).toEqual({ x: '#a569ff', a: '#bf8400' }); expect(viewerColours(['a', 'b'], c)).toEqual({ a: '#a569ff', b: '#bf8400' })
    expect(SUPERPOSE_COLOURS).toEqual(['#a569ff', '#bf8400'])
  })
})

