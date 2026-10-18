import { describe, expect, it } from 'vitest'
import { contains, focusBox } from './PdbRnaViewer2D'

const full = { x: 0, y: 0, w: 1000, h: 800 }
describe('2D focus window', () => {
  it('a single nucleotide gets the minimum context window, centred on it', () => {
    expect(focusBox({ x: 500, y: 400, w: 10, h: 10 }, 260, full)).toEqual({ x: 375, y: 275, w: 260, h: 260 })
  })
  it('a large selection is padded by 25%, never beyond the whole drawing', () => {
    const b = focusBox({ x: 100, y: 100, w: 400, h: 100 }, 260, full)
    expect(b.w).toBe(500); expect(b.h).toBe(260)
    expect(focusBox({ x: 0, y: 0, w: 1000, h: 800 }, 260, full)).toMatchObject({ w: 1000, h: 800 })
  })
  it('containment', () => {
    expect(contains(full, { x: 10, y: 10, w: 5, h: 5 })).toBe(true)
    expect(contains({ x: 0, y: 0, w: 100, h: 100 }, { x: 95, y: 10, w: 10, h: 5 })).toBe(false)
  })
})
