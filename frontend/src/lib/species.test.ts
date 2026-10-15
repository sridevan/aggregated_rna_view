import { describe, expect, it } from 'vitest'
import { readFileSync } from 'node:fs'
import { classify, summariseSpecies } from './species'
import type { InstanceRow } from '../types/backend'

const row = (pdb: string, chain: string, organism: string | null, species: string | null, speciesTax: string | null, tax: string | null): InstanceRow =>
  ({ pdb_id: pdb, auth_asym_id: chain, seq_name: `${pdb}_1_${chain}_${chain}`, organism, organism_species: species, species_tax_id: speciesTax, tax_id: tax } as unknown as InstanceRow)

describe('species summary', () => {
  const rows = [
    row('1aaa', 'A', 'Bacillus subtilis subsp. subtilis str. 168', 'Bacillus subtilis', '1423', '224308'),
    row('1aaa', 'B', 'Bacillus subtilis subsp. subtilis str. 168', 'Bacillus subtilis', '1423', '224308'),   // same entry, same species -> 1 entry, 2 instances
    row('1bbb', 'A', 'Bacillus subtilis', 'Bacillus subtilis', '1423', '1423'),
    row('1ccc', 'X', 'Vibrio vulnificus', 'Vibrio vulnificus', '672', '672'),
    row('1ccc', 'Y', 'Bacillus subtilis', 'Bacillus subtilis', '1423', '1423'),                               // entry with two species counts once for each
    row('1zzz', 'A', 'Vibrio vulnificus', 'Vibrio vulnificus', '672', '672'),
    row('1ddd', 'A', 'synthetic construct', 'synthetic construct', '32630', '32630'),
    row('1eee', 'A', 'unclassified Paenibacillus', null, null, '77965'),                                      // tax id but no species-rank name
    row('1fff', 'A', null, null, null, null),                                                                  // no taxonomy at all
    ...['g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q'].map((s, i) => row(`2${s}${s}${s}`, 'A', `Species ${s}`, `Species ${s}`, `${1000 + i}`, `${1000 + i}`)),
  ]
  it('classifies synthetic, unresolved and unknown separately', () => {
    const by = (pdb: string) => rows.find((r) => r.pdb_id === pdb)!
    expect(classify(by('1ddd')).category).toBe('synthetic'); expect(classify(by('1eee')).category).toBe('unresolved'); expect(classify(by('1fff')).category).toBe('unknown'); expect(classify(by('1aaa')).category).toBe('species')
  })
  it('deduplicates PDB entries per species and counts instances', () => {
    const s = summariseSpecies(rows, 'entries', 10)
    const bs = s.rows.find((r) => r.label === 'Bacillus subtilis')!
    expect(bs.entries).toBe(3); expect(bs.instances).toBe(4); expect(bs.italic).toBe(true); expect(bs.sourceNames).toContain('Bacillus subtilis subsp. subtilis str. 168')
    expect(s.rows.find((r) => r.label === 'Vibrio vulnificus')!.entries).toBe(2)
    expect(s.totals).toEqual({ entries: 18, instances: 20 })
    expect(s.identifiedSpecies).toBe(13)
    expect(s.rows.map((r) => r.category)).toEqual([...Array(11).fill('species'), 'synthetic', 'unresolved', 'unknown'])   // top 10 + Other + 3 special rows
    const other = s.rows.find((r) => r.key === 'other')!
    expect(other.groupedSpecies).toBe(3); expect(other.entries).toBe(3); expect(other.label).toBe('Other identified species (3)'); expect(other.italic).toBe(false)
    expect(s.rows.find((r) => r.category === 'unresolved')!.sourceNames).toEqual(['unclassified Paenibacillus'])
  })
  it('orders by the selected unit', () => {
    expect(summariseSpecies(rows, 'instances', 10).rows[0].label).toBe('Bacillus subtilis')
    expect(summariseSpecies(rows, 'entries', 2).rows.find((r) => r.key === 'other')!.groupedSpecies).toBe(11)
  })
  it('works on the real RF00167 data', () => {
    const inst = JSON.parse(readFileSync(new URL('../../../data/families/RF00167/instances.json', import.meta.url), 'utf8')) as InstanceRow[]
    const s = summariseSpecies(inst)
    expect(s.totals.entries).toBe(66)
    expect(s.rows[0].label).toBe('Bacillus subtilis')
    expect(s.rows.find((r) => r.category === 'synthetic')).toBeTruthy()
    expect(s.rows.reduce((n, r) => n + (r.category === 'species' && r.key !== 'other' ? 1 : 0), 0) + (s.rows.find((r) => r.key === 'other')?.groupedSpecies ?? 0)).toBe(s.identifiedSpecies)
  })
})
