import { describe, expect, it } from 'vitest'
import { dissimilarityMatrix } from './structureGroups'
import type { StructureClusteringDataset } from '../types/backend'

// 4 instances in dataset order a,b,c,d; upper triangle row-major: ab ac ad bc bd cd
const make = (upper: (number | null)[], order: string[] | null) => ({
  dataset: { instances: ['a', 'b', 'c', 'd'] },
  metrics: { discrepancy: { values_upper: upper, distance_range: { min: 0.1, max: 0.6 }, ...(order ? { seriation: { method: 'tppl', repetitions: 10, seed: 39873, order, missing_pairs: 0, reference: {} } } : {}) } },
}) as unknown as StructureClusteringDataset
const U = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

describe('dissimilarityMatrix', () => {
  it('reorders rows and columns by the stored seriation and keeps every pair on its instances', () => {
    const m = dissimilarityMatrix(make(U, ['c', 'a', 'd', 'b']))!
    expect(m.ids).toEqual(['c', 'a', 'd', 'b'])
    const v = (x: string, y: string) => m.values[m.ids.indexOf(x)][m.ids.indexOf(y)]
    expect([v('a', 'b'), v('a', 'c'), v('a', 'd'), v('b', 'c'), v('b', 'd'), v('c', 'd')]).toEqual(U)
    for (let r = 0; r < 4; r++) { expect(m.values[r][r]).toBe(0); for (let c = 0; c < 4; c++) expect(m.values[r][c]).toBe(m.values[c][r]) }   // symmetric, zero diagonal
    expect(m.max).toBe(0.6); expect(m.nMissing).toBe(0)
  })
  it('the same input gives the same matrix; the order does not depend on the input instance order beyond the stored seriation', () => {
    expect(dissimilarityMatrix(make(U, ['b', 'd', 'a', 'c']))).toEqual(dissimilarityMatrix(make(U, ['b', 'd', 'a', 'c'])))
  })
  it('missing values stay null (never the maximum) and are counted once per pair', () => {
    const m = dissimilarityMatrix(make([0.1, null, 0.3, Number.NaN, 0.5, 0.6], ['a', 'b', 'c', 'd']))!
    expect(m.values[0][2]).toBeNull(); expect(m.values[2][0]).toBeNull(); expect(m.values[1][2]).toBeNull()
    expect(m.nMissing).toBe(2); expect(m.max).toBe(0.6)
  })
  it('no seriation: no matrix; an order that does not list every instance once is rejected', () => {
    expect(dissimilarityMatrix(make(U, null))).toBeNull()
    expect(() => dissimilarityMatrix(make(U, ['a', 'b', 'c']))).toThrow()
    expect(() => dissimilarityMatrix(make(U, ['a', 'b', 'c', 'c']))).toThrow()
    expect(() => dissimilarityMatrix(make(U, ['a', 'b', 'c', 'x']))).toThrow()
  })
})
