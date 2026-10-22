import { describe, expect, it } from 'vitest'
import { loadDrugbank, parsePdbeDrugbank, parseRcsbDrugbank, shortName } from './drugbank'

const PDBE_ADE = `data_ADE
_pdbe_chem_comp_drugbank_details.comp_id ADE
_pdbe_chem_comp_drugbank_details.drugbank_id DB00173
_pdbe_chem_comp_drugbank_details.type 'small molecule'
_pdbe_chem_comp_drugbank_details.name Adenine
_pdbe_chem_comp_drugbank_details.description 'A purine base.'
#
`
const res = (body: string, ok = true) => Promise.resolve({ ok, text: () => Promise.resolve(body), json: () => Promise.resolve(JSON.parse(body)) } as Response)

describe('DrugBank cross-references (PDBe first, RCSB only when PDBe has none)', () => {
  it('parses the PDBe key-value record, quoted names and a looped category', () => {
    expect(parsePdbeDrugbank(PDBE_ADE)).toEqual({ id: 'DB00173', name: 'Adenine' })
    expect(parsePdbeDrugbank(PDBE_ADE.replace('name Adenine', "name 'S-adenosyl-L-homocysteine'"))).toEqual({ id: 'DB00173', name: 'S-adenosyl-L-homocysteine' })
    expect(parsePdbeDrugbank('data_X\nloop_\n_pdbe_chem_comp_drugbank_details.comp_id\n_pdbe_chem_comp_drugbank_details.drugbank_id\nX DB01234\n#\n')).toEqual({ id: 'DB01234', name: null })
    expect(parsePdbeDrugbank('data_X\n_chem_comp.id X\n')).toBeNull()
  })
  it('parses the RCSB related-resource record', () => {
    expect(parseRcsbDrugbank({ rcsb_chem_comp_related: [{ resource_name: 'ChEBI', resource_accession_code: '1' }, { resource_name: 'DrugBank', resource_accession_code: 'DB00118' }] })).toEqual({ id: 'DB00118', name: null })
    expect(parseRcsbDrugbank({})).toBeNull()
  })
  it('uses PDBe when it has a mapping, RCSB only when PDBe has none, and reports what was checked', async () => {
    const calls: string[] = []
    const f = (body: Record<string, [string, boolean]>) => ((u: string) => { calls.push(u); const k = Object.keys(body).find((x) => u.includes(x))!; return res(...body[k]) }) as unknown as typeof fetch
    expect(await loadDrugbank('T1', f({ 'ebi.ac.uk': [PDBE_ADE, true], 'rcsb.org': ['{}', true] }))).toMatchObject({ id: 'DB00173', source: 'PDBe' })
    expect(calls.filter((u) => u.includes('rcsb'))).toEqual([])
    expect(await loadDrugbank('T2', f({ 'ebi.ac.uk': ['data_T2\n', true], 'rcsb.org': [JSON.stringify({ rcsb_chem_comp_related: [{ resource_name: 'DrugBank', resource_accession_code: 'DB00118' }] }), true] }))).toMatchObject({ id: 'DB00118', source: 'RCSB' })
    expect(await loadDrugbank('T3', f({ 'ebi.ac.uk': ['data_T3\n', true], 'rcsb.org': ['{}', true] }))).toEqual({ id: null, checked: ['PDBe', 'RCSB'] })
    expect(await loadDrugbank('T4', (() => Promise.reject(new Error('offline'))) as unknown as typeof fetch)).toEqual({ id: null, checked: [] })
  })
  it('shortens names to 20 characters with an ellipsis', () => {
    expect(shortName('Hypoxanthine')).toBe('Hypoxanthine'); expect(shortName('9h-purine-2,6-diamine')).toBe('9h-purine-2,6-diamin…'); expect(shortName('x'.repeat(20))).toBe('x'.repeat(20))
  })
})
