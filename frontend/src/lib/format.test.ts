import { describe, expect, it } from 'vitest'
import { shortMethod } from './format'

describe('shortMethod (display only)', () => {
  it('abbreviates the three common methods, case-insensitively, and leaves others unchanged', () => {
    expect(shortMethod('X-ray diffraction')).toBe('X-ray')
    expect(shortMethod('Electron Microscopy')).toBe('EM'); expect(shortMethod('electron microscopy')).toBe('EM')
    expect(shortMethod('Solution NMR')).toBe('NMR')
    expect(shortMethod('Neutron diffraction')).toBe('Neutron diffraction'); expect(shortMethod(null)).toBeNull()
  })
})
