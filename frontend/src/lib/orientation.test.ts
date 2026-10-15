import { describe, expect, it } from 'vitest'
import { readFileSync } from 'node:fs'
import { Mat4 } from 'molstar/lib/mol-math/linear-algebra'
import { applyMat4, isApplicable, toMat4 } from './orientation'
import type { OrientationDataset } from '../types/backend'

const data = JSON.parse(readFileSync(new URL('../../../data/families/RF00167/orientation.json', import.meta.url), 'utf8')) as OrientationDataset

describe('orientation matrix convention (backend numpy vs Mol* Mat4)', () => {
  it('reference transform is the identity', () => {
    const ref = data.transforms[data.reference.seq_name]
    expect(ref.is_identity).toBe(true)
    expect(Mat4.isIdentity(toMat4(ref)!)).toBe(true)
  })
  it('Mat4 built from R,t reproduces the backend check point for every applicable transform', () => {
    let n = 0
    for (const t of Object.values(data.transforms)) {
      if (!isApplicable(t) || !t.check) continue
      const m = toMat4(t)!
      const got = applyMat4(m, t.check.original)
      for (let i = 0; i < 3; i++) expect(Math.abs(got[i] - t.check.transformed[i])).toBeLessThan(1e-6)
      expect(Math.abs(t.determinant! - 1)).toBeLessThan(1e-6)
      n++
    }
    expect(n).toBeGreaterThanOrEqual(80)
  })
  it('insufficient fits are not applied', () => {
    expect(isApplicable(data.transforms['9nbh_1_A_A'])).toBe(false)
    expect(isApplicable(data.transforms['1y26_1_X_A'])).toBe(true)
  })
})
