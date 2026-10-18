import { describe, expect, it } from 'vitest'
import { accFromPath, AVAILABLE_ACCS, canonicalPath, FAMILY_BLURBS } from './route'

describe('family accession in the URL path', () => {
  it('reads /{rfam_acc}, case-insensitive, optional trailing slash', () => {
    expect(accFromPath('/RF00167')).toBe('RF00167'); expect(accFromPath('/rf00167/')).toBe('RF00167'); expect(accFromPath('/RF00001')).toBe('RF00001')
  })
  it('anything else is not an accession', () => {
    for (const p of ['/', '', '/RF167', '/RF00167/structures', '/families/RF00167', '/PF00001']) expect(accFromPath(p)).toBeNull()
  })
  it('canonical path: upper-case accession; the home page "/" for the root and anything else', () => {
    expect(canonicalPath('/rf00167')).toBe('/RF00167'); expect(canonicalPath('/')).toBe('/'); expect(canonicalPath('/whatever')).toBe('/'); expect(canonicalPath('')).toBe('/')
  })
  it('every available family has a home-page description', () => {
    for (const a of AVAILABLE_ACCS) expect(FAMILY_BLURBS[a]?.length).toBeGreaterThan(20)
  })
})
