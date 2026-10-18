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
