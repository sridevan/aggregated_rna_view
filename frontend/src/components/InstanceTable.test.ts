import { describe, expect, it } from 'vitest'
import { groupRank, type TableGroups } from './InstanceTable'

describe('default Structural group order', () => {
  // dendrogram leaf order: b(1) a(1, rep) c(1) | d(2, rep) e(2); x not clustered
  const info: Record<string, { id: number; rep: boolean; order: number }> = { b: { id: 1, rep: false, order: 0 }, a: { id: 1, rep: true, order: 1 }, c: { id: 1, rep: false, order: 2 }, d: { id: 2, rep: true, order: 3 }, e: { id: 2, rep: false, order: 4 } }
  const groups = (selectedId: number): Pick<TableGroups, 'of' | 'selectedId'> => ({ selectedId,
    of: (sn) => (info[sn] ? { kind: 'group', id: info[sn].id, label: '', color: '', representative: info[sn].rep, order: info[sn].order } : { kind: 'excluded', reason: 'n' }) })
  const order = (sel: number) => ['a', 'b', 'c', 'd', 'e', 'x'].sort((p, q) => { const rp = groupRank(groups(sel), p), rq = groupRank(groups(sel), q); return rp == null ? 1 : rq == null ? -1 : rp - rq })
  it('active group first with its representative on top, then the other groups in dendrogram order; not clustered last', () => {
    expect(order(1)).toEqual(['a', 'b', 'c', 'd', 'e', 'x'])
    expect(order(2)).toEqual(['d', 'e', 'b', 'a', 'c', 'x'])                 // other groups keep dendrogram order (their representatives are not lifted)
    expect(groupRank(groups(1), 'x')).toBeNull()
  })
  it('ranks are numeric (a string sort would put 10 before 9)', () => {
    expect(typeof groupRank(groups(1), 'c')).toBe('number'); expect(groupRank(groups(2), 'c')! > groupRank(groups(2), 'e')!).toBe(true)
  })
})
