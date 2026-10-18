import { describe, expect, it } from 'vitest'
import { alignmentInfo, comparisonSummary } from './StructuresTab'

const ids = ['1abc_1_A_A', '2xyz_1_B_B']
const fit = (sn: string) => (sn.startsWith('1abc') ? { quality: 'good', rmsd: 1.234 } : { quality: 'poor', rmsd: 4.63 })

describe('comparisonSummary', () => {
  it('reads the stored pairwise discrepancy and the shared positions, nothing else', () => {
    const c = { ids, value: [[0, 0.5937], [0.5937, 0]], metric: { kind: 'dissimilarity', label: 'Geometric discrepancy' } } as unknown as Parameters<typeof comparisonSummary>[0]
    const s = comparisonSummary(c, ids[0], ids[1], 57)
    expect(s.value).toBe(0.5937)
    expect(s.metric).toBe('Geometric discrepancy 0.594 · 57 shared Rfam positions')
  })
  it('TM-score: shows the stored similarity and 1 − TM-score; an unknown chain gives no value', () => {
    const c = { ids, value: [[1, 0.71], [0.71, 1]], metric: { kind: 'similarity', label: 'TM-score' } } as unknown as Parameters<typeof comparisonSummary>[0]
    expect(comparisonSummary(c, ids[0], ids[1], 43).metric).toBe('TM-score 0.710 (1 − TM-score 0.290) · 43 shared Rfam positions')
    expect(comparisonSummary(c, ids[0], '9zzz_1_A_A', 43)).toEqual({ metric: 'pairwise value not available · 43 shared Rfam positions', value: null })
  })
})

describe('alignmentInfo', () => {
  it('lists the reference, each displayed chain fit and the distinction from the metric alignment', () => {
    const t = alignmentInfo(ids, fit, '4FE5 B', 'dissimilarity')
    expect(t).toContain('fitted onto 4FE5 B using C4′ atoms')
    expect(t).toContain('• 1ABC A: good fit, RMSD 1.23 Å\n• 2XYZ B: poor fit, RMSD 4.63 Å')
    expect(t).toContain('may differ from the alignment used to calculate the geometric discrepancy')
  })
  it('TM-score wording; a chain without a transform and an unknown reference are named', () => {
    const t = alignmentInfo(ids, (sn) => (sn.startsWith('1abc') ? fit(sn) : null), null, 'similarity')
    expect(t).toContain('• 2XYZ B: no alignment'); expect(t).toContain('the family reference'); expect(t).toContain('calculate the TM-score')
  })
})
